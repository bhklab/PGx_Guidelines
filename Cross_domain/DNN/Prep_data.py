
import numpy as np
import pandas as pd

def prep_data(args):
    
    drug = args.drug
    
    if args.mode == 'cross':
        
        CTRP_exprs = pd.read_csv("Data_All/CTRP.exprsALL.tsv", sep = "\t", index_col=0)
        GDSC_exprs = pd.read_csv("Data_All/GDSCv2.exprsALL.tsv", sep = "\t", index_col=0)
        gCSI_exprs = pd.read_csv("Data_All/gCSI.exprsALL.tsv", sep = "\t", index_col=0)

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
        
        # change lines 32 to 34 for IC50
        
        CTRP_aac_drug = CTRP_aac.loc[drug].dropna()
        GDSC_aac_drug = GDSC_aac.loc[drug].dropna()
        gCSI_aac_drug = gCSI_aac.loc[drug].dropna()

        idx_ctrp = CTRP_exprs.columns.intersection(CTRP_aac_drug.index)
        idx_ctrp = [x for x in idx_ctrp if x not in idx_other_ctrp]    
        idx_gdsc = GDSC_exprs.columns.intersection(GDSC_aac_drug.index)
        idx_gdsc = [x for x in idx_gdsc if x not in idx_other_gdsc]    
        idx_gcsi = gCSI_exprs.columns.intersection(gCSI_aac_drug.index)
        idx_gcsi = [x for x in idx_gcsi if x not in idx_other_gcsi] 

        CTRP_exprs_drug = pd.DataFrame.transpose(CTRP_exprs.loc[:,idx_ctrp])
        CTRP_aac_drug = CTRP_aac_drug.loc[idx_ctrp]
        GDSC_exprs_drug = pd.DataFrame.transpose(GDSC_exprs.loc[:,idx_gdsc])
        GDSC_aac_drug = GDSC_aac_drug.loc[idx_gdsc]
        gCSI_exprs_drug = pd.DataFrame.transpose(gCSI_exprs.loc[:,idx_gcsi])
        gCSI_aac_drug = gCSI_aac_drug.loc[idx_gcsi]
        
        return CTRP_exprs_drug.values, CTRP_aac_drug.values, GDSC_exprs_drug.values, GDSC_aac_drug.values, gCSI_exprs_drug.values, gCSI_aac_drug.values
   