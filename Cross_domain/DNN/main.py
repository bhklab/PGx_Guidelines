import argparse
import os
import sys
import torch
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import torch.optim as optim
sys.setrecursionlimit(1000000)
import warnings
from sklearn.model_selection import KFold, GridSearchCV
import pandas as pd
from math import sqrt
import sklearn.preprocessing as sk
from sklearn import metrics
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import random
from random import randint
from sklearn.model_selection import StratifiedKFold
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics import mean_squared_error
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import torch.utils.data
import itertools 

from Prep_data import *
from Net import *
from Func import *

warnings.filterwarnings("ignore")
torch.set_num_threads(64)

def main():
    train_arg_parser = argparse.ArgumentParser()
    train_arg_parser.add_argument("--drug", type=str, default='Gemcitabine', help='input drug to train a model') 
    train_arg_parser.add_argument("--data_root", type=str, default='Data/', help="path to molecular and pharmacological data")        
    train_arg_parser.add_argument("--save_logs", type=str, default='./PGxResults/logs/', help='path of folder to write log')
    train_arg_parser.add_argument("--save_models", type=str, default='./PGxResults/models/', help='folder for saving model')
    train_arg_parser.add_argument("--save_results", type=str, default='./PGxResults/results/', help='folder for saving model')
    train_arg_parser.add_argument("--mode", type=str, default='cross', help='within or cross-domain evaluation mode')
    #train_arg_parser.add_argument("--mbs", type=int, default=32, help='mini-batch size')
    train_arg_parser.add_argument("--hd", type=int, default=2, help='strcuture of the network')
    train_arg_parser.add_argument("--ldr", type=float, default=0.5, help='dropout')
    train_arg_parser.add_argument("--idr", type=float, default=0.5, help='input dropout')
    train_arg_parser.add_argument("--wd", type=float, default=0.5, help='weight decay')
    train_arg_parser.add_argument("--lr", type=float, default=0.001, help='learning rate')
    train_arg_parser.add_argument("--epoch", type=int, default=30, help='number of epochs')
    train_arg_parser.add_argument("--fold", type=int, default=7, help='k-fold cross validation')              
    train_arg_parser.add_argument("--seed", type=int, default=42, help='set the random seed')          
    train_arg_parser.add_argument('--gpu', type=int, default=0, help='set using GPU or not')
        
    args = train_arg_parser.parse_args()

    paths = [args.save_logs, args.save_models, args.save_results]

    for path in paths:
        if not os.path.isdir(path):
            os.makedirs(path)
        with open(path +'/args.txt', 'w') as f:
            f.write(str(args))    
    
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    
    device = torch.device("cuda:" + str(args.gpu) if torch.cuda.is_available() else "cpu")
    
    
    if args.mode == 'cross':
        X_tr, Y_tr, X_ts_1, Y_ts_1, X_ts_2, Y_ts_2 = prep_data(args)
        
    if args.mode == 'within':
        X_tr, Y_tr = prep_data(args) 
          
            
    kfold = KFold(n_splits=args.fold, random_state=args.seed)
    k = 0
    best_pr = 0
    
    loss_fun = torch.nn.MSELoss()
    total_val = []
    total_aac = []
    
    for train_index, val_index in kfold.split(X_tr, Y_tr):
        train_loss = []
        train_pr = []
        val_loss = []
        val_pr = []
        
        k = k + 1

        X_train = X_tr[train_index,:]
        X_val =  X_tr[val_index,:]

        y_train = Y_tr[train_index]
        y_val = Y_tr[val_index]

        scaler = sk.StandardScaler()
        scaler.fit(X_train)
        X_train_N = scaler.transform(X_train)
        X_val_N = scaler.transform(X_val)

        TX_val_N = torch.FloatTensor(X_val_N)
        Ty_val = torch.FloatTensor(y_val)

        trainDataset = torch.utils.data.TensorDataset(torch.FloatTensor(X_train_N), torch.FloatTensor(y_train))

        trainLoader = torch.utils.data.DataLoader(dataset = trainDataset, batch_size=32, shuffle=True, num_workers=1)

        Model = Network(args, X_train_N)

        optimizer = optim.SGD(Model.parameters(), lr=args.lr, weight_decay = args.wd)
    
        f = open(os.path.join(args.save_logs, 'args.txt'), mode='a')
        f.write('fold:{}\n'.format(k))
        f.close()  
    
        for ite in range(args.epoch):

            epoch_loss, epoch_pr = train(args, Model, loss_fun, optimizer, trainLoader)

            train_loss.append(epoch_loss)      
            train_pr.append(epoch_pr)

            epoch_val_loss, epoch_Val_pr,_ = validate_workflow(args, Model, loss_fun, TX_val_N, Ty_val)
            val_loss.append(epoch_val_loss)
            val_pr.append(epoch_Val_pr)

            f = open(os.path.join(args.save_logs, 'args.txt'), mode='a')
            f.write('iteration:{}, train epoch loss:{}\n'.format(ite, epoch_loss))
            f.write('iteration:{}, validation epoch loss:{}\n'.format(ite, epoch_val_loss))
            f.close()                          
            
            if epoch_Val_pr > best_pr: 
                best_pr = epoch_Val_pr
                f = open(os.path.join(args.save_results, 'Best_val.txt'), mode='a')
                f.write('iteration:{}, best validation correlation:{}\n'.format(ite, best_pr))
                f.close()
                torch.save(Model.state_dict(), os.path.join(args.save_models, 'Best_Model.pt'))
                scaler_b = scaler
                
        plots(args, k, train_loss, train_pr, val_loss, val_pr)
        Model.load_state_dict(torch.load(os.path.join(args.save_models, 'Best_Model.pt')))
        Model.eval()
        _,_, preds= validate_workflow(args, Model, loss_fun, TX_val_N, Ty_val)
        total_val.append(preds.detach().numpy().flatten())
        total_aac.append(Ty_val.detach().numpy())
        
    Model.load_state_dict(torch.load(os.path.join(args.save_models, 'Best_Model.pt')))
    Model.eval()

    test_1, test_2, test_s1, test_s2, test_k1, test_k2 = heldout_test(args, Model, X_ts_1, Y_ts_1, X_ts_2, Y_ts_2, scaler_b)   
    
    final_pred = list(itertools.chain.from_iterable(total_val))
    final_labels = list(itertools.chain.from_iterable(total_aac))    
    f = open(os.path.join(args.save_results, 'Total_val.txt'), mode='a')
    f.write('Total validation Pearson:{}\n'.format(pearsonr(final_pred, final_labels)))
    f.write('Total validation Spearman:{}\n'.format(spearmanr(final_pred, final_labels)))
    f.write('Total validation Kendall:{}\n'.format(kendalltau(final_pred, final_labels)))
    f.write('---------------------------------\n')    
    f.close()

    f = open(os.path.join(args.save_results, 'Target.txt'), mode='a')
    f.write('Test Pearson GDSC:{}\n'.format(test_1))
    f.write('Test Pearson gCSI:{}\n'.format(test_2))
    f.write('---------------------------------\n')
    f.write('Test Spearman GDSC:{}\n'.format(test_s1))
    f.write('Test Spearman gCSI:{}\n'.format(test_s2))
    f.write('---------------------------------\n')    
    f.write('Test Kendall GDSC:{}\n'.format(test_k1))
    f.write('Test Kendall gCSI:{}\n'.format(test_k2))    
    f.close()

        
if __name__ == "__main__":
    main()    
        
        
        
   
    
    
    
    