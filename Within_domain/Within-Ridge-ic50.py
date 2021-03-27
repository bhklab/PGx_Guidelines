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
from sklearn.model_selection import KFold
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics import mean_squared_error
import itertools
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

def RG(X_train, Y_train, X_test, alpha, cv, seed):
    Y_pred_0 = np.zeros([X_test.shape[0], 1])
    RG_pipe = Pipeline([('scaler', StandardScaler()),('RdG', Ridge())])
    model = GridSearchCV(RG_pipe, param_grid={"RdG__alpha": alpha}, scoring='neg_mean_squared_error', cv=KFold(n_splits=cv, shuffle=True, random_state=seed))
    y_train = Y_train
    x_train = X_train
    model.fit(x_train, y_train)
    Y_pred_0 = model.predict(X_test)
    return Y_pred_0


alph = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
folds = 10
outer_folds = 10
seeds = 42

drugs = ["Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib", 
         "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib"]

drugs = ["Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib", 
         "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib"]

CTRP_exprs = pd.read_csv("Data_All/CTRP.exprsALL.tsv", sep = "\t", index_col=0)

CTRP_aac = pd.read_csv("Data_All/CTRP.aacALL.tsv", sep = "\t", index_col=0)

CTRP_ic50 = pd.read_csv("Data_All/CTRP.logIC50.tsv", sep = "\t", index_col=0)

CTRP_info = pd.read_csv("Data_All/CTRP.infoALL.tsv", sep = "\t", index_col=0)
idx_other_ctrp = CTRP_info.index[CTRP_info["Tumor"] == 1]

with open('ridge-res-ic50.txt', mode='a') as f:
    f.write('*********************************************\n')
    f.write('IC50\n')
    f.write('*********************************************\n')
    f.close()
    
with open('ridge-vec-ic50.txt', mode='a') as f2:
    f2.write('*********************************************\n')
    f2.write('IC50\n')
    f2.write('*********************************************\n')
    f2.close()  
    

for drug in drugs:
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

            pred = RG(X_train_N, y_train, CTRP_exprs_ic50_drug.values[test_index,:],
                           alph, folds, seeds)
            pred_drug.append(pred)
            labels.append(CTRP_ic50_drug.values[test_index])
        final_pred = list(itertools.chain.from_iterable(pred_drug))
        final_labels = list(itertools.chain.from_iterable(labels))

        save_pr.append(pearsonr(final_pred, final_labels)[0])
        save_mse.append(sqrt(mean_squared_error(final_pred, final_labels)))
        save_ken.append(kendalltau(final_pred, final_labels)[0])
        save_sp.append(spearmanr(final_pred, final_labels)[0])
        
    with open('ridge-res-ic50.txt', mode='a') as f:
        f.write('drug:{}\n'.format(drug))
        f.write('pearson mean:{}, pearson std:{}\n'.format(np.mean(save_pr), np.std(save_pr)))
        f.write('kendall mean:{}, kendall std:{}\n'.format(np.mean(save_ken), np.std(save_pr)))
        f.write('spearman mean:{}, spearman std:{}\n'.format(np.mean(save_sp), np.std(save_pr)))
        f.write('RMSE mean:{}, RMSE std:{}\n'.format(np.mean(save_mse), np.std(save_pr)))
        f.write('*********************************************\n')
        f.close()
    with open('ridge-vec-ic50.txt', mode='a') as f2:
        f2.write('drug:{}\n'.format(drug))
        f2.write('pearson:{}\n'.format(save_pr))
        f2.write('kendall:{}\n'.format(save_ken))
        f2.write('spearman:{}\n'.format(save_sp))
        f2.write('rmse:{}\n'.format(save_mse))
        f2.close()        