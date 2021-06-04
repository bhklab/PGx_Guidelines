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


def RG(X_train, Y_train, X_test, alpha, cv, seed):
    Y_pred_0 = np.zeros([X_test.shape[0], 1])
    RG_pipe = Pipeline([('scaler', StandardScaler()),('RdG', Ridge())])
    model = GridSearchCV(RG_pipe, param_grid={"RdG__alpha": alpha}, scoring='neg_mean_squared_error', cv=KFold(n_splits=cv, shuffle=True, random_state=seed))
    y_train = Y_train
    x_train = X_train
    model.fit(x_train, y_train)
    Y_pred_0 = model.predict(X_test)
    return Y_pred_0

l1 = [0.1, 0.3, 0.5, 0.7, 0.9]
alph = [1e-2, 1e-1, 1, 1e2]
folds = 3
seeds = 42
alph_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

drugs = ["Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib", 
         "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib"]

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


col_All = ["AllAll_GDSC_Pearson", "AllAll_gCSI_Pearson", "AllAll_GDSC_Spearman","AllAll_gCSI_Spearman", "AllAll_GDSC_RMSE", "AllAll_gCSI_RMSE", "AllSolid_GDSC_Pearson", "AllSolid_gCSI_Pearson", "AllSolid_GDSC_Spearman", "AllSolid_gCSI_Spearman", "AllSolid_GDSC_RMSE", "AllSolid_gCSI_RMSE", "AllNonSolid_GDSC_Pearson", "AllNonSolid_gCSI_Pearson", "AllNonSolid_GDSC_Spearman", "AllNonSolid_gCSI_Spearman", "AllNonSolid_GDSC_RMSE", "AllNonSolid_gCSI_RMSE"]
df_All = pd.DataFrame(index=drugs, columns=col_All)

col_Solid = ["SolidAll_GDSC_Pearson", "SolidAll_gCSI_Pearson", "SolidAll_GDSC_Spearman", "SolidAll_gCSI_Spearman", "SolidAll_GDSC_RMSE", "SolidAll_gCSI_RMSE", "SolidSolid_GDSC_Pearson", "SolidSolid_gCSI_Pearson", "SolidSolid_GDSC_Spearman", "SolidSolid_gCSI_Spearman", "SolidSolid_GDSC_RMSE", "SolidSolid_gCSI_RMSE", "SolidNonSolid_GDSC_Pearson", "SolidNonSolid_gCSI_Pearson", "SolidNonSolid_GDSC_Spearman", "SolidNonSolid_gCSI_Spearman", "SolidNonSolid_GDSC_RMSE", "SolidNonSolid_gCSI_RMSE"]
df_Solid = pd.DataFrame(index=drugs, columns=col_Solid)

col_NonSolid = ["NonSolidAll_GDSC_Pearson", "NonSolidAll_gCSI_Pearson", "NonSolidAll_GDSC_Spearman", "NonSolidAll_gCSI_Spearman", "NonSolidAll_GDSC_RMSE", "NonSolidAll_gCSI_RMSE", "NonSolidSolid_GDSC_Pearson", "NonSolidSolid_gCSI_Pearson", "NonSolidSolid_GDSC_Spearman", "NonSolidSolid_gCSI_Spearman", "NonSolidSolid_GDSC_RMSE", "NonSolidSolid_gCSI_RMSE", "NonSolidNonSolid_GDSC_Pearson", "NonSolidNonSolid_gCSI_Pearson", "NonSolidNonSolid_GDSC_Spearman", "NonSolidNonSolid_gCSI_Spearman", "NonSolidNonSolid_GDSC_RMSE", "NonSolidNonSolid_gCSI_RMSE"]
df_NonSolid = pd.DataFrame(index=drugs, columns=col_NonSolid)

col_pr = "AllAll_GDSC_Pearson"
col_sp = "AllAll_GDSC_Spearman"
col_rmse = "AllAll_GDSC_RMSE"

for drug in drugs:
    save_pr = []
    save_rmse = []
    save_ken = []
    save_sp = [] 
    for it in range(10):
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        GDSC_aac_drug = GDSC_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_gdsc = GDSC_exprs.columns.intersection(GDSC_aac_drug.index)

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        GDSC_exprs_drug = pd.DataFrame.transpose(GDSC_exprs.loc[:,idx_gdsc])
        GDSC_aac_drug = GDSC_aac_drug.loc[idx_gdsc]
        
        idx = min(len(idx_ctrp), len(idx_gdsc))
        rand_ctrp = np.random.choice(CTRP_exprs_drug.index, idx, replace=False)
        rand_gdsc = np.random.choice(GDSC_exprs_drug.index, idx, replace=False)        

        X_train_N = CTRP_exprs_drug.loc[rand_ctrp].values
        y_train = CTRP_aac_drug[rand_ctrp].values

        pred_gdsc = RG(X_train_N, y_train, GDSC_exprs_drug.loc[rand_gdsc].values, alph_r, folds, seeds)
        
        save_pr.append(pearsonr(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)[0])
        save_sp.append(spearmanr(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)[0])
        save_rmse.append(sqrt(mean_squared_error(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)))
    
    df_All.loc[drug,col_pr] = np.mean(save_pr)
    df_All.loc[drug,col_sp] = np.mean(save_sp)
    df_All.loc[drug,col_rmse] = np.mean(save_rmse)
    
    
col_pr = "AllAll_gCSI_Pearson"
col_sp = "AllAll_gCSI_Spearman"
col_rmse = "AllAll_gCSI_RMSE"

for drug in drugs:
    save_pr = []
    save_rmse = []
    save_ken = []
    save_sp = [] 
    for it in range(10):
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        gCSI_aac_drug = gCSI_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_gcsi = gCSI_exprs.columns.intersection(gCSI_aac_drug.index)

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        gCSI_exprs_drug = pd.DataFrame.transpose(gCSI_exprs.loc[:,idx_gcsi])
        gCSI_aac_drug = gCSI_aac_drug.loc[idx_gcsi]
        
        idx = min(len(idx_ctrp), len(idx_gcsi))
        rand_ctrp = np.random.choice(CTRP_exprs_drug.index, idx, replace=False)
        rand_gcsi = np.random.choice(gCSI_exprs_drug.index, idx, replace=False)         

        X_train_N = CTRP_exprs_drug.loc[rand_ctrp].values
        y_train = CTRP_aac_drug[rand_ctrp].values

        pred_gcsi = RG(X_train_N, y_train, gCSI_exprs_drug.loc[rand_gcsi].values, alph_r, folds, seeds)
        
        save_pr.append(pearsonr(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)[0])
        save_sp.append(spearmanr(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)[0])
        save_rmse.append(sqrt(mean_squared_error(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)))
    
    df_All.loc[drug,col_pr] = np.mean(save_pr)
    df_All.loc[drug,col_sp] = np.mean(save_sp)
    df_All.loc[drug,col_rmse] = np.mean(save_rmse)    
     
col_pr = "AllSolid_GDSC_Pearson"
col_sp = "AllSolid_GDSC_Spearman"
col_rmse = "AllSolid_GDSC_RMSE"

for drug in drugs:
    save_pr = []
    save_rmse = []
    save_ken = []
    save_sp = [] 
    for it in range(10):
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        GDSC_aac_drug = GDSC_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_gdsc = GDSC_exprs.columns.intersection(GDSC_aac_drug.index)
        idx_gdsc = [x for x in idx_gdsc if x not in idx_other_gdsc]
        

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        GDSC_exprs_drug = pd.DataFrame.transpose(GDSC_exprs.loc[:,idx_gdsc])
        GDSC_aac_drug = GDSC_aac_drug.loc[idx_gdsc]
        
        idx = min(len(idx_ctrp), len(idx_gdsc))
        rand_ctrp = np.random.choice(CTRP_exprs_drug.index, idx, replace=False)
        rand_gdsc = np.random.choice(GDSC_exprs_drug.index, idx, replace=False)        

        X_train_N = CTRP_exprs_drug.loc[rand_ctrp].values
        y_train = CTRP_aac_drug[rand_ctrp].values

        pred_gdsc = RG(X_train_N, y_train, GDSC_exprs_drug.loc[rand_gdsc].values, alph_r, folds, seeds)
        
        save_pr.append(pearsonr(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)[0])
        save_sp.append(spearmanr(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)[0])
        save_rmse.append(sqrt(mean_squared_error(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)))
    
    df_All.loc[drug,col_pr] = np.mean(save_pr)
    df_All.loc[drug,col_sp] = np.mean(save_sp)
    df_All.loc[drug,col_rmse] = np.mean(save_rmse)
    
    
col_pr = "AllSolid_gCSI_Pearson"
col_sp = "AllSolid_gCSI_Spearman"
col_rmse = "AllSolid_gCSI_RMSE"

for drug in drugs:
    save_pr = []
    save_rmse = []
    save_ken = []
    save_sp = [] 
    for it in range(10):
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        gCSI_aac_drug = gCSI_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_gcsi = gCSI_exprs.columns.intersection(gCSI_aac_drug.index)
        idx_gcsi = [x for x in idx_gcsi if x not in idx_other_gcsi]
        

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        gCSI_exprs_drug = pd.DataFrame.transpose(gCSI_exprs.loc[:,idx_gcsi])
        gCSI_aac_drug = gCSI_aac_drug.loc[idx_gcsi]
        
        idx = min(len(idx_ctrp), len(idx_gcsi))
        rand_ctrp = np.random.choice(CTRP_exprs_drug.index, idx, replace=False)
        rand_gcsi = np.random.choice(gCSI_exprs_drug.index, idx, replace=False)         

        X_train_N = CTRP_exprs_drug.loc[rand_ctrp].values
        y_train = CTRP_aac_drug[rand_ctrp].values

        pred_gcsi = RG(X_train_N, y_train, gCSI_exprs_drug.loc[rand_gcsi].values, alph_r, folds, seeds)
        
        save_pr.append(pearsonr(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)[0])
        save_sp.append(spearmanr(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)[0])
        save_rmse.append(sqrt(mean_squared_error(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)))
    
    df_All.loc[drug,col_pr] = np.mean(save_pr)
    df_All.loc[drug,col_sp] = np.mean(save_sp)
    df_All.loc[drug,col_rmse] = np.mean(save_rmse)      

    
col_pr = "AllNonSolid_GDSC_Pearson"
col_sp = "AllNonSolid_GDSC_Spearman"
col_rmse = "AllNonSolid_GDSC_RMSE"

for drug in drugs:
    save_pr = []
    save_rmse = []
    save_ken = []
    save_sp = [] 
    for it in range(10):
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        GDSC_aac_drug = GDSC_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_gdsc = GDSC_exprs.columns.intersection(GDSC_aac_drug.index)
        idx_gdsc = [x for x in idx_gdsc if x in idx_other_gdsc]
        

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        GDSC_exprs_drug = pd.DataFrame.transpose(GDSC_exprs.loc[:,idx_gdsc])
        GDSC_aac_drug = GDSC_aac_drug.loc[idx_gdsc]
        
        idx = min(len(idx_ctrp), len(idx_gdsc))
        rand_ctrp = np.random.choice(CTRP_exprs_drug.index, idx, replace=False)
        rand_gdsc = np.random.choice(GDSC_exprs_drug.index, idx, replace=False)        

        X_train_N = CTRP_exprs_drug.loc[rand_ctrp].values
        y_train = CTRP_aac_drug[rand_ctrp].values

        pred_gdsc = RG(X_train_N, y_train, GDSC_exprs_drug.loc[rand_gdsc].values, alph_r, folds, seeds)
        
        save_pr.append(pearsonr(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)[0])
        save_sp.append(spearmanr(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)[0])
        save_rmse.append(sqrt(mean_squared_error(pred_gdsc, GDSC_aac_drug.loc[rand_gdsc].values)))
    
    df_All.loc[drug,col_pr] = np.mean(save_pr)
    df_All.loc[drug,col_sp] = np.mean(save_sp)
    df_All.loc[drug,col_rmse] = np.mean(save_rmse)
    
    
col_pr = "AllNonSolid_gCSI_Pearson"
col_sp = "AllNonSolid_gCSI_Spearman"
col_rmse = "AllNonSolid_gCSI_RMSE"

for drug in drugs:
    save_pr = []
    save_rmse = []
    save_ken = []
    save_sp = [] 
    for it in range(10):
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        gCSI_aac_drug = gCSI_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_gcsi = gCSI_exprs.columns.intersection(gCSI_aac_drug.index)
        idx_gcsi = [x for x in idx_gcsi if x in idx_other_gcsi]
        

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        gCSI_exprs_drug = pd.DataFrame.transpose(gCSI_exprs.loc[:,idx_gcsi])
        gCSI_aac_drug = gCSI_aac_drug.loc[idx_gcsi]
        
        idx = min(len(idx_ctrp), len(idx_gcsi))
        rand_ctrp = np.random.choice(CTRP_exprs_drug.index, idx, replace=False)
        rand_gcsi = np.random.choice(gCSI_exprs_drug.index, idx, replace=False)         

        X_train_N = CTRP_exprs_drug.loc[rand_ctrp].values
        y_train = CTRP_aac_drug[rand_ctrp].values

        pred_gcsi = RG(X_train_N, y_train, gCSI_exprs_drug.loc[rand_gcsi].values, alph_r, folds, seeds)
        
        save_pr.append(pearsonr(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)[0])
        save_sp.append(spearmanr(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)[0])
        save_rmse.append(sqrt(mean_squared_error(pred_gcsi, gCSI_aac_drug.loc[rand_gcsi].values)))
    
    df_All.loc[drug,col_pr] = np.mean(save_pr)
    df_All.loc[drug,col_sp] = np.mean(save_sp)
    df_All.loc[drug,col_rmse] = np.mean(save_rmse)     

df_All.to_csv('All_Results.csv')            
            
                        