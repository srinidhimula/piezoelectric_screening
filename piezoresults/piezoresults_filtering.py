#!/usr/bin/env python3
"""Script to filter results of Piezoelectric constant 'e' obtained from highthroughput calculations.
Filters on magnetic moment 1) being equal with QMOF database and optimisation, 2) being equal in Optimisation and piezoelectric steps
Filters on Clamped and dynamic 'e' values wrt point group, ie 1) Check on some terms in 'e' should be zero as per symmetry
2) Check on some terms in 'e' should be nearly equal with a tolerance specified
"""
from piezofilters_pointgroup import *
import numpy as np
import pandas as pd

##Change the input results files and output result files here based on file you want to use##

df_qmofv13=pd.read_excel('../../qmof_sginfo1_v13.xlsx')
df_dynamic=pd.read_csv('../dynamic_set123_normdirec.csv')
df_clamped=pd.read_csv('../clamped_set123_normdirec.csv')
df_totale=pd.read_csv('../totale_set123_normdirec.csv')
csv_names=['../clamped_editedset123.csv','../dynamic_editedset123.csv']
df_checklist=[df_clamped,df_dynamic]
#print (df_clamped)

df_dynamic['info_optnetmagmom']=df_clamped['info_optnetmagmom'].values
df_dynamic['info_piezonetmagmom']=df_clamped['info_piezonetmagmom'].values
df_totale['info_optnetmagmom']=df_clamped['info_optnetmagmom'].values
df_totale['info_piezonetmagmom']=df_clamped['info_piezonetmagmom'].values


list_noncentro_pgs=['1','2','m','222','mm2','3','32','3m','6','4','6mm','4mm','-6','-6m2','-4','-42m','-43m','23','422','622']
list_centro_pgs=['-1','2/m','mmm','4/m','4/mmm','-3','-3m','6/m','6/mmm','m-3','m-3m']

def magmom_check(df_k1):
    df_qmofv13['outputs.pbe.net_magmom']=df_qmofv13['outputs.pbe.net_magmom'].apply(np.around)
    df_k1['info_optnetmagmom']=abs(df_k1['info_optnetmagmom']).apply(np.around)
    df_k1['info_piezonetmagmom']=abs(df_k1['info_piezonetmagmom']).apply(np.around)
    df_k1['info_magmom_matchopt&piezo']=np.where((df_k1['info_optnetmagmom']==df_k1['info_piezonetmagmom']),"equal_bothsteps","notequal_bothsteps")
    for column in df_k1["info_MOFCode"]:
        qmof_magmom=df_qmofv13.loc[df_qmofv13["name"]==column,"outputs.pbe.net_magmom"].iloc[0]
        obtained_optmagmom=df_k1.loc[df_k1["info_MOFCode"]==column,'info_optnetmagmom'].iloc[0]
        if obtained_optmagmom==qmof_magmom:
            df_k1.loc[df_k1["info_MOFCode"]==column,"info_magmom_matchwithqmof"]="equal"
            df_k1.loc[df_k1["info_MOFCode"]==column,"info_qmofmagmom"]=qmof_magmom
        if obtained_optmagmom!=qmof_magmom:
            #print ("matches",column)
            df_k1.loc[df_k1["info_MOFCode"]==column,
                      "info_magmom_matchwithqmof"]="notequal"
            df_k1.loc[df_k1["info_MOFCode"]==column,"info_qmofmagmom"]=qmof_magmom
    return df_k1

def check_filters():
    for csv_k,df_k in zip(csv_names,df_checklist):
        print (csv_k)
        df_k=magmom_check(df_k)
        for i in list_noncentro_pgs:
            pg_loc=df_k.index[df_k['info_PointGroup']==i].tolist()
            #print (pg_loc)
            df_updated1= echecks_pointgroup_zeroconv(df_k,i,pg_loc,pow(10,-3))
            df_updated2= echecks_pointgroup_depconv(df_updated1,i,pg_loc,pow(10,-3))
        for j in list_centro_pgs:
            pg_loc1=df_updated2.index[df_updated2['info_PointGroup']==j].tolist()
            for m in pg_loc1:
                df_updated2.at[m,"Status_dependentconv"]="Centro as per csd"
                df_updated2.at[m,"Status_zeroconvergence"]="Centro as per csd"
        df_updated2.to_csv(csv_k)

def mofs_passedfilters():
    df_clamped_passed=pd.read_csv(csv_names[0])
    df_dynamic_passed=pd.read_csv(csv_names[1])
    df_clamped_converged= df_clamped_passed.query('`info_magmom_matchopt&piezo`=="equal_bothsteps" & `info_magmom_matchwithqmof`=="equal" & (`Status_zeroconvergence`=="Success at all eik" or `Status_zeroconvergence`=="This sg doesnt have any")  & (`Status_dependentconv`=="Success at all eik" or `Status_dependentconv`=="This sg doesnt have any")')
    df_dynamic_converged= df_dynamic_passed.query('`info_magmom_matchopt&piezo`=="equal_bothsteps" & `info_magmom_matchwithqmof`=="equal" & (`Status_zeroconvergence`=="Success at all eik" or `Status_zeroconvergence`=="This sg doesnt have any")  & (`Status_dependentconv`=="Success at all eik" or `Status_dependentconv`=="This sg doesnt have any")')
    #print (df_clamped_converged)
    for column1 in df_totale["info_MOFCode"]:
        #print (df_qmofv13.loc[df_qmofv13["name"]==column1,"csd.Structure.Type"].iloc[0])
        df_totale.loc[df_totale["info_MOFCode"] == column1,
                      "info_ClampedNorm"] = df_clamped_passed.loc[df_clamped_passed["info_MOFCode"] == column1, "info_ClampedNorm"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"] == column1,
                      "info_DynamicNorm"] = df_dynamic_passed.loc[df_dynamic_passed["info_MOFCode"] == column1, "info_DynamicNorm"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"] == column1,
                      "Clampednorm_direction"] = df_clamped_passed.loc[df_clamped_passed["info_MOFCode"] == column1, "Clampednorm_direction"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"] == column1,
                      "Dynamicnorm_direction"] = df_dynamic_passed.loc[df_dynamic_passed["info_MOFCode"] == column1, "Dynamicnorm_direction"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"]==column1,
                    "info_smilesnodes"]=df_qmofv13.loc[df_qmofv13["name"]==column1,"info.mofid.smiles_nodes"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"]==column1,
                    "info_smileslinkers"]=df_qmofv13.loc[df_qmofv13["name"]==column1,"info.mofid.smiles_linkers"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"] == column1,
                      "info_natoms"] = df_qmofv13.loc[df_qmofv13["name"] == column1, "info.natoms"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"] == column1,
                      "info.density"] = df_qmofv13.loc[df_qmofv13["name"] == column1, "info.density"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"]==column1,
                    "info_Structuretype"]=df_qmofv13.loc[df_qmofv13["name"]==column1,"csd.Structure.Type"].iloc[0]
        df_totale.loc[df_totale["info_MOFCode"]==column1,
                    "info_qmofcsdequal?"]=df_qmofv13.loc[df_qmofv13["name"]==column1,"csd.qmof.equal?"].iloc[0]
        if (column1 in df_clamped_converged.values) and (column1 in df_dynamic_converged.values):
            df_totale.loc[df_totale["info_MOFCode"]==column1,
                        "info_ClampedDynamic_convergence"]="Converged"
        else:
            df_totale.loc[df_totale["info_MOFCode"]==column1,
                        "info_ClampedDynamic_convergence"]="NotConverged"
    df_totale.to_csv('../totale_editedset123.csv')
    return df_totale

"""check_filter function should be done for convergence checks on magmom equality in opt and piezo steps, magmom equality in qmof database and this work, 
convergence of 'eik' to zero and convergence of dependent 'eik' based on the point group. This is checked in Clamped and Dynamic eiks separately."""
check_filters()
"""The MOFs which satisfy the convergence criteria in dynamic and clamped eik components are then indicated as converged or not converged 
in the totale csv file"""
df_totale_temp=mofs_passedfilters()
