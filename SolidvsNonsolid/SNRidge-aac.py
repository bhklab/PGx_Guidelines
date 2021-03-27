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
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score

def RG(X_train, Y_train, X_test, alpha, cv, seed):
    Y_pred_0 = np.zeros([X_test[0].shape[0], 1])
    Y_pred_1 = np.zeros([X_test[1].shape[0], 1])
    RG_pipe = Pipeline([('scaler', StandardScaler()),('RdG', Ridge())])
    model = GridSearchCV(RG_pipe, param_grid={"RdG__alpha": alpha}, scoring='neg_mean_squared_error', cv=KFold(n_splits=cv, shuffle=True, random_state=seed))
    y_train = Y_train
    x_train = X_train
    model.fit(x_train, y_train)
    Y_pred_0 = model.predict(X_test[0])
    Y_pred_1 = model.predict(X_test[1])
    return Y_pred_0, Y_pred_1 


drugs = ["Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib", 
         "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib"]

alph = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
folds = 10
outer_folds = 10
seeds = 42

CTRP_exprs = pd.read_csv("Data_All/CTRP.exprsALL.tsv", sep = "\t", index_col=0)
GDSC_exprs = pd.read_csv("Data_All/GDSCv2.exprsALL.tsv", sep = "\t", index_col=0)
gCSI_exprs = pd.read_csv("Data_All/gCSI.exprsALL.tsv", sep = "\t", index_col=0)

CTRP_exprs = CTRP_exprs.iloc[:-25,:]
GDSC_exprs = GDSC_exprs.iloc[:-25,:]
gCSI_exprs = gCSI_exprs.iloc[:-25,:]

CTRP_aac = pd.read_csv("Data_All/CTRP.aacALL.tsv", sep = "\t", index_col=0)
GDSC_aac = pd.read_csv("Data_All/GDSCv2.aacALL.tsv", sep = "\t", index_col=0)
gCSI_aac = pd.read_csv("Data_All/gCSI.aacALL.tsv", sep = "\t", index_col=0)

CTRP_ic50 = pd.read_csv("Data_All/CTRP.logIC50.tsv", sep = "\t", index_col=0)
GDSC_ic50 = pd.read_csv("Data_All/GDSC.logIC50.tsv", sep = "\t", index_col=0)
gCSI_ic50 = pd.read_csv("Data_All/gCSI.logIC50.tsv", sep = "\t", index_col=0)

CTRP_info = pd.read_csv("Data_All/CTRP.infoALL.tsv", sep = "\t", index_col=0)
idx_other_ctrp = CTRP_info.index[CTRP_info["Tumor"] == 1]
GDSC_info = pd.read_csv("Data_All/GDSCv2.infoALL.tsv", sep = "\t", index_col=0)
idx_other_gdsc = GDSC_info.index[GDSC_info["Tumor"] == 1]
gCSI_info = pd.read_csv("Data_All/gCSI.infoALL.tsv", sep = "\t", index_col=0)
idx_other_gcsi = gCSI_info.index[gCSI_info["Tumor"] == 1]



with open('ridge-sn-aac.txt', mode='a') as f:
    f.write('*********************************************\n')
    f.write('AAC\n')
    f.write('*********************************************\n')
    f.close()
    
with open('ridge-sn-aac.txt', mode='a') as f2:
    f2.write('*********************************************\n')
    f2.write('AAC\n')
    f2.write('*********************************************\n')
    f2.close()    

for drug in drugs:
    save_roc_gdsc = []
    save_aupr_gdsc = []
    save_roc_gcsi = []
    save_aupr_gcsi = []
    for it in range(10):
        N_nonsolid = len(CTRP_info.index[CTRP_info["Tumor"] == 1])
        solid_idx = CTRP_info.index[CTRP_info["Tumor"] == 0]
        idx_remove = np.random.choice(solid_idx, N_nonsolid, replace=False)

        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        GDSC_aac_drug = GDSC_aac.loc[drug].dropna()
        gCSI_aac_drug = gCSI_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_ctrp = [x for x in idx_ctrp if x not in idx_remove]
        idx_gdsc = GDSC_exprs.columns.intersection(GDSC_aac_drug.index)
        idx_gcsi = gCSI_exprs.columns.intersection(gCSI_aac_drug.index)

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        GDSC_exprs_drug = pd.DataFrame.transpose(GDSC_exprs.loc[:,idx_gdsc])
        GDSC_aac_drug = GDSC_aac_drug.loc[idx_gdsc]
        gCSI_exprs_drug = pd.DataFrame.transpose(gCSI_exprs.loc[:,idx_gcsi])
        gCSI_aac_drug = gCSI_aac_drug.loc[idx_gcsi]

        CTRP_info_drug = CTRP_info.loc[idx_ctrp,:]
        GDSC_info_drug = GDSC_info.loc[idx_gdsc,:]
        gCSI_info_drug = gCSI_info.loc[idx_gcsi,:]

        X_train_N = CTRP_exprs_drug.values
        y_train = CTRP_aac_drug.values

        pred_gdsc, pred_gcsi = RG(X_train_N, y_train, [GDSC_exprs_drug, gCSI_exprs_drug],
                       alph, folds, seeds)

        save_roc_gdsc.append(roc_auc_score(GDSC_info_drug["Tumor"].values, pred_gdsc, average='micro'))
        save_roc_gcsi.append(roc_auc_score(gCSI_info_drug["Tumor"].values, pred_gcsi, average='micro'))
        save_aupr_gdsc.append(average_precision_score(GDSC_info_drug["Tumor"].values, pred_gdsc, average='micro'))
        save_aupr_gcsi.append(average_precision_score(gCSI_info_drug["Tumor"].values, pred_gcsi, average='micro'))
        
    with open('ridge-sn-aac.txt', mode='a') as f:
        f.write('drug:{}\n'.format(drug))
        f.write('auroc gdsc mean:{}, auroc gdsc std:{}\n'.format(np.mean(save_roc_gdsc), np.std(save_roc_gdsc)))
        f.write('aupr gdsc mean:{}, aupr gdsc std:{}\n'.format(np.mean(save_aupr_gdsc), np.std(save_aupr_gdsc)))        
        f.write('auroc gcsi mean:{}, auroc gcsi std:{}\n'.format(np.mean(save_roc_gcsi), np.std(save_roc_gcsi)))
        f.write('aupr gcsi mean:{}, aupr gcsi std:{}\n'.format(np.mean(save_aupr_gcsi), np.std(save_aupr_gcsi)))
        f.write('*********************************************\n')
        f.close()
    with open('ridge-sn-aac.txt', mode='a') as f2:
        f2.write('drug:{}\n'.format(drug))
        f2.write('auroc gdsc:{}\n'.format(save_roc_gdsc))
        f2.write('aupr gdsc:{}\n'.format(save_aupr_gdsc))
        f2.write('auroc gcsi:{}\n'.format(save_roc_gcsi))
        f2.write('aupr gcsi:{}\n'.format(save_aupr_gcsi))
        f2.close()