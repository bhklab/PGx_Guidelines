import numpy as np
from scipy.stats import pearsonr, spearmanr, kendalltau
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import torch
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data
import os 


def train(args, Model, loss, optimizer, DataLoader_Tr):
    epoch_cost = []
    epoch_PR = []
    Model.train()
    
    for i, (exprs, AAC) in enumerate(DataLoader_Tr):
        Model.train()

        AAC_pred = Model(exprs)

        ls = loss(AAC_pred, AAC.view(-1,1))

        r_train,_ = pearsonr(AAC_pred.detach().numpy().flatten(), AAC.detach().numpy())

        optimizer.zero_grad()

        ls.backward()

        optimizer.step()

        epoch_cost.append(ls.item()) 
        epoch_PR.append(r_train)
        
    return np.mean(epoch_cost), np.mean(epoch_PR)

def validate_workflow(args, Model, loss, TX_val_N, Ty_val):
    Model.eval()

    AAC_pred_test = Model(TX_val_N)
    ls = loss(AAC_pred_test, Ty_val.view(-1,1))

    r_val,_ = pearsonr(AAC_pred_test.detach().numpy().flatten(), Ty_val.detach().numpy())
    sr_val,_ = spearmanr(AAC_pred_test.detach().numpy().flatten(), Ty_val.detach().numpy())
    kr_val,_ = kendalltau(AAC_pred_test.detach().numpy().flatten(), Ty_val.detach().numpy())

    return ls.item(), r_val, AAC_pred_test

def heldout_test(args, Best_M, X_ts_1, Y_ts_1, X_ts_2, Y_ts_2, scaler):
    
    X_ts1_N = scaler.transform(X_ts_1)
    X_ts2_N = scaler.transform(X_ts_2)
    
    TX_ts1_N = torch.FloatTensor(X_ts1_N)
    TX_ts2_N = torch.FloatTensor(X_ts2_N)
    
    AAC_ts1 = Best_M(TX_ts1_N)
    AAC_ts2 = Best_M(TX_ts2_N)
    
    r_ts1,_ = pearsonr(AAC_ts1.detach().numpy().flatten(), Y_ts_1) 
    r_ts2,_ = pearsonr(AAC_ts2.detach().numpy().flatten(), Y_ts_2) 
    sr_ts1,_ = spearmanr(AAC_ts1.detach().numpy().flatten(), Y_ts_1) 
    sr_ts2,_ = spearmanr(AAC_ts2.detach().numpy().flatten(), Y_ts_2) 
    kr_ts1,_ = kendalltau(AAC_ts1.detach().numpy().flatten(), Y_ts_1) 
    kr_ts2,_ = kendalltau(AAC_ts2.detach().numpy().flatten(), Y_ts_2)     
    
    return r_ts1, r_ts2, sr_ts1, sr_ts2, kr_ts1, kr_ts2 
    
def plots(args, k, train_loss, train_pr, val_loss, val_pr):
    
    if not os.path.isdir(args.save_results +"/"+ "plots/"):
        os.makedirs(args.save_results +"/"+ "plots/")
    new_path = args.save_results +"/"+ "plots/"
    
    title_loss = 'Loss Train fold = {}, drug ={}, lr = {}, epoch = {}, hd = {}. wd = {}, ldr = {}, idr = {}'.\
                      format(k, args.drug, args.lr, args.epoch, args.hd, args.wd, args.ldr, args.idr)
    plt.plot(np.squeeze(train_loss), '-r', np.squeeze(val_loss), '-b')
    plt.title(title_loss)
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['train', 'validation'], loc='upper right')
    plt.savefig(os.path.join(new_path, title_loss + '.png'), dpi = 300)
    plt.close()

    title_PR = 'Correlation Train fold = {}, drug ={}, lr = {}, epoch = {}, hd = {}. wd = {}, ldr = {}, idr = {}'.\
                      format(k, args.drug, args.lr, args.epoch, args.hd, args.wd, args.ldr, args.idr)
    plt.plot(np.squeeze(train_pr), '-r', np.squeeze(val_pr), '-b')
    plt.title(title_PR)
    plt.ylabel('Pearson')
    plt.xlabel('Epoch')
    plt.legend(['train', 'validation'], loc='upper right')
    plt.savefig(os.path.join(new_path, title_PR + '.png'), dpi = 300)
    plt.close()    
    
    