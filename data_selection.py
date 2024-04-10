"""Script to select data for highthroughput calculations from the QMOF database based on exclusion of some metals,
non-centro symmetric on CSD database, number of atoms"""
import numpy as np
import pandas as pd
qmof_v13file='../../qmof_sginfo1_v13.xlsx'
sheet='Sheet1'
calcsdone='polarmofs_done'
df_database=pd.read_excel(qmof_v13file,sheet_name=sheet,usecols='A:AR',engine='openpyxl')

"""Exclude some metals based on psuedopotential text """
metals_excluded=["Ag 02Apr2005","Au 04Oct2007","Pr_3 07Sep2000","Dy_3 06Sep2000","Er_3 06Sep2000","Eu_3 20Oct2008","Gd_3 06Sep2000","Ho_3 06Sep2000","Lu_3 06Sep2000","Nd_3 06Sep2000","U 06Sep2000","Sm_3 07Sep2000","Tb_3 06Sep2000","Th 07Sep2000","Tm_3 20Jan2003","U 06Sep2000"]
mofs_calcsdone=pd.read_excel(qmof_v13file,sheet_name=calcsdone,engine='openpyxl')
#print (mofs_calcsdone)
#print (df_database['info.source'])
"""Filtering by only selecting CSD deposited MOFs, noncentrosymmteric MOFs as per CSD and QMOF databases and in noncentro include only if spacegroup(sg) sg of csd=sg in qmof """
df_csd_noncentro_sgequal= df_database.query('`info.source`=="CSD" & `Structure.Type`!="Centro" & `csd.Structure.Type`!="Centro" & `csd.Structure.Type`!="#N/A" & `csd.qmof.equal?`=="equal"')
df_csd_noncentro_finalmetal=df_csd_noncentro_sgequal[df_csd_noncentro_sgequal.inputs_pbe_pseudopotentials.str.contains('|'.join(metals_excluded))==False]

df_calcnoncentro_atomless150=df_csd_noncentro_finalmetal.query('`info.natoms`<=150')
df_calcnoncentro_atomgrt150=df_csd_noncentro_finalmetal.query('`info.natoms`>150')

df_calcnoncentro_atomless150.to_csv("../InputMofsList/Candidates_atomless150.csv")
df_calcnoncentro_atomgrt150.to_csv("../InputMofsList/Candidates_atomgrt150.csv")
"""Filtering already performed calculations via mofs_calcdone"""
df_todolesthn150=df_calcnoncentro_atomless150[df_calcnoncentro_atomless150.name.str.contains('|'.join(list(mofs_calcsdone['name'])))==False]
df_todolesthn150.to_csv("../InputMofsList/Candidates_atomless150_excldone.csv")
df_todogrtthn150=df_calcnoncentro_atomgrt150[df_calcnoncentro_atomgrt150.name.str.contains('|'.join(list(mofs_calcsdone['name'])))==False]
df_todogrtthn150.to_csv("../InputMofsList/Candidates_atomgrt150_excldone.csv")
somemetals=["Cu 22Jun2005","Cd 06Sep2000","Zn 06Sep2000"]
df_todolesthn150_Cu=df_todolesthn150[df_todolesthn150.inputs_pbe_pseudopotentials.str.contains(somemetals[0])]
df_todolesthn150_Cd=df_todolesthn150[df_todolesthn150.inputs_pbe_pseudopotentials.str.contains(somemetals[1])]
df_todolesthn150_Zn=df_todolesthn150[df_todolesthn150.inputs_pbe_pseudopotentials.str.contains(somemetals[2])]
df_todolesthn150_rem=df_todolesthn150[df_todolesthn150.inputs_pbe_pseudopotentials.str.contains('|'.join(somemetals))==False]

#print (df_todolesthn150_Zn.merge(df_todolesthn150_Cu, how = 'inner' ,indicator=False))
#print (df_todolesthn150_rem)
#print (len(df_todolesthn150),len(df_todolesthn150_Zn),len(df_todolesthn150_Cu),len(df_todolesthn150_Cd),len(df_todolesthn150_rem))
#print (len(df_todolesthn150_Cu))
#df_todolesthn150_rem.to_csv('../InputMofsList/allexceptCuZnCd_454.txt',header=None,index=None, columns=("name","csd.spacegroup.num","csd.pointgroup"),sep='\t',mode='w')
#df_todolesthn150_Cu.to_csv('../InputMofsList/Culessatm150_279.txt',header=None,index=None, columns=("name","csd.spacegroup.num","csd.pointgroup"),sep='\t',mode='w')

