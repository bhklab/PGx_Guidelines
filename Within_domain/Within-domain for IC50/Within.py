import argparse
import os
import sys
import numpy as np
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import math
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
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import itertools

from withinfunc import *

def main():
    train_arg_parser = argparse.ArgumentParser()
    train_arg_parser.add_argument("--drug", type=str, default='Gemcitabine', help='input drug to train a model') 
    train_arg_parser.add_argument("--method", type=str, default='ridge', help="method to train")
    train_arg_parser.add_argument("--metric", type=str, default='rawic50', help="method to train")
        
    args = train_arg_parser.parse_args()

    alph = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    l1 = [0.1, 0.3, 0.5, 0.7, 0.9]
    alph_en = [1e-3, 1e-2, 1e-1, 1, 1e2, 1e3]
    n_estimators = [100, 500, 1000]
    depth = [10, 50]
    mtry = [1/4]
    folds = 10
    seeds = 42
    outer_folds = 10
    
    drug = args.drug
    method = args.method
    metric = args.metric
    
#     print(drug)
#     print(method)
#     print(metric)

    

    CTRP_exprs = pd.read_csv("Data_All/CTRP.exprsALL.tsv", sep = "\t", index_col=0)
    
    if metric == "logtrncic50":
        CTRP_ic50 = pd.read_csv("Data_All/CTRP.logIC50.tsv", sep = "\t", index_col=0)
        
    if metric == "logic50":
        CTRP_ic50 = pd.read_csv("Data_All/CTRP.onlylogIC50.tsv", sep = "\t", index_col=0)
        
    if metric == "rawic50":
        CTRP_ic50 = pd.read_csv("Data_All/CTRP.rawIC50.tsv", sep = "\t", index_col=0)
            
    if metric == "trncic50":        
        CTRP_ic50 = pd.read_csv("Data_All/CTRP.trcIC50.tsv", sep = "\t", index_col=0)

    CTRP_info = pd.read_csv("Data_All/CTRP.infoALL.tsv", sep = "\t", index_col=0)
    idx_other_ctrp = CTRP_info.index[CTRP_info["Tumor"] == 1]
    
#     print(CTRP_ic50.shape)
    
    save_pr = []
    save_mse = []
    save_ken = []
    save_sp = []    
    for it in range(10):
        kf = KFold(n_splits=outer_folds, shuffle=True)
        pred_drug = []
        final_pred = []
        labels = []
        final_labels = []  
        CTRP_ic50_drug = CTRP_ic50.loc[drug].dropna()
        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_ic50_drug.index)
        idx_ctrp = [x for x in idx_ctrp if x not in idx_other_ctrp]    
        CTRP_exprs_ic50_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_ic50_drug = CTRP_ic50_drug.loc[idx_ctrp]
        k = 0
        for train_index, test_index in kf.split(CTRP_exprs_ic50_drug):
            k+=1
            X_train_N = CTRP_exprs_ic50_drug.values[train_index,:]
            y_train = CTRP_ic50_drug.values[train_index]
            
            if method == "en":
                pred = EN(X_train_N, y_train, CTRP_exprs_ic50_drug.values[test_index,:], alph_en, l1, folds, seeds)
            if method =="rf":
                pred = RF(X_train_N, y_train, CTRP_exprs_ic50_drug.values[test_index,:], n_estimators, depth, mtry, folds, seeds)
            if method == "ridge":
                pred = RG(X_train_N, y_train, CTRP_exprs_ic50_drug.values[test_index,:], alph, folds, seeds)
            
            pred_drug.append(pred)
            labels.append(CTRP_ic50_drug.values[test_index])
        final_pred = list(itertools.chain.from_iterable(pred_drug))
        final_labels = list(itertools.chain.from_iterable(labels))

        save_pr.append(pearsonr(final_pred, final_labels)[0])
        save_mse.append(sqrt(mean_squared_error(final_pred, final_labels)))
        save_ken.append(kendalltau(final_pred, final_labels)[0])
        save_sp.append(spearmanr(final_pred, final_labels)[0])
        
        
    if method == "en":
        with open(metric+'en-res.txt', mode='a') as f:
            f.write('drug:{}\n'.format(drug))
            f.write('pearson mean:{}, pearson std:{}\n'.format(np.mean(save_pr), np.std(save_pr)))
            f.write('kendall mean:{}, kendall std:{}\n'.format(np.mean(save_ken), np.std(save_pr)))
            f.write('spearman mean:{}, spearman std:{}\n'.format(np.mean(save_sp), np.std(save_pr)))
            f.write('RMSE mean:{}, RMSE std:{}\n'.format(np.mean(save_mse), np.std(save_pr)))
            f.write('*********************************************\n')
            f.close()
        with open(metric+'en-vec.txt', mode='a') as f2:
            f2.write('drug:{}\n'.format(drug))
            f2.write('pearson:{}\n'.format(save_pr))
            f2.write('kendall:{}\n'.format(save_ken))
            f2.write('spearman:{}\n'.format(save_sp))
            f2.write('rmse:{}\n'.format(save_mse))
            f2.close()  
    if method =="rf":
        with open(metric+'rf-res.txt', mode='a') as f:
            f.write('drug:{}\n'.format(drug))
            f.write('pearson mean:{}, pearson std:{}\n'.format(np.mean(save_pr), np.std(save_pr)))
            f.write('kendall mean:{}, kendall std:{}\n'.format(np.mean(save_ken), np.std(save_pr)))
            f.write('spearman mean:{}, spearman std:{}\n'.format(np.mean(save_sp), np.std(save_pr)))
            f.write('RMSE mean:{}, RMSE std:{}\n'.format(np.mean(save_mse), np.std(save_pr)))
            f.write('*********************************************\n')
            f.close()
        with open(metric+'rf-vec.txt', mode='a') as f2:
            f2.write('drug:{}\n'.format(drug))
            f2.write('pearson:{}\n'.format(save_pr))
            f2.write('kendall:{}\n'.format(save_ken))
            f2.write('spearman:{}\n'.format(save_sp))
            f2.write('rmse:{}\n'.format(save_mse))
            f2.close() 
    if method =="ridge":
        with open(metric+'ridge-res.txt', mode='a') as f:
            f.write('drug:{}\n'.format(drug))
            f.write('pearson mean:{}, pearson std:{}\n'.format(np.mean(save_pr), np.std(save_pr)))
            f.write('kendall mean:{}, kendall std:{}\n'.format(np.mean(save_ken), np.std(save_pr)))
            f.write('spearman mean:{}, spearman std:{}\n'.format(np.mean(save_sp), np.std(save_pr)))
            f.write('RMSE mean:{}, RMSE std:{}\n'.format(np.mean(save_mse), np.std(save_pr)))
            f.write('*********************************************\n')
            f.close()
        with open(metric+'ridge-vec.txt', mode='a') as f2:
            f2.write('drug:{}\n'.format(drug))
            f2.write('pearson:{}\n'.format(save_pr))
            f2.write('kendall:{}\n'.format(save_ken))
            f2.write('spearman:{}\n'.format(save_sp))
            f2.write('rmse:{}\n'.format(save_mse))
            f2.close() 
    
if __name__ == "__main__":
    main()     