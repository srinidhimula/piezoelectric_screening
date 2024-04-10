import numpy as np
import sys
import pandas as pd
pd.set_option("display.max_columns", None)
from pymatgen.io.vasp.outputs import *
from pymatgen.core import *
#from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.local_env import *
from pymatgen.core.structure import *
from oximachinerunner import OximachineRunner

group1=[3,11,19,37,55]#["Li","Na","K","Rb","Cs"]
group2=[4,12,20,38,56]#["Be","Mg","Ca","Sr","Ba"]
group3=[21,39]#["Sc","Y"]
group4=[22,40,72]#["Ti","Zr","Hf"]
group5=[23,41,73]#["V","Nb","Ta"]
group6=[24,42,74]#["Cr","Mo","W"]
group7=[25,43,75]#["Mn","Tc","Re"]
group8=[26,44,76]#["Fe","Ru","Os"]
group9=[27,45,77]#["Co","Rh","Ir"]
group10=[28,46,78]#["Ni","Pd","Pt"]
group11=[29,47,79]#["Cu","Ag","Au"]
group12=[30,48,80]#["Zn","Cd","Hg"]
group13=[13,31,49,81]#["Al","Ga","In","Tl"]
group14=[32,50,82]#["Ge","Sn","Pb"]
group15=[33,51,83]#["As","Sb","Bi"]
group16=[52,84]#["Te","Po"]
groupfblock=[57,58,70,92]#["La","Ce","Yb","U"]
diff_groups=[group1,group2,group3,group4,group5,group6,group7,group8,group9,group10,group11,group12,group13,group14,group15,group16,groupfblock]
group_labels=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","fblock"]
df_structinfo=pd.DataFrame()
df_structinfo1=pd.DataFrame()

def struct_nodeinfo(contcarfile,where):
    """Determines metal atomic number,group number, coordination environment around it, coord type and if water is present.
    coordination number and type may not be reliable for all MOFs. Hence check results before using them."""
    struct=Structure.from_file(contcarfile) #read CONTCAR
    at_nos=np.asarray(struct.atomic_numbers) #get all atomic numbers
    groups_all=[]
    metal_index=[]
    for i,grp in enumerate(diff_groups):
        check_grp=np.where(np.isin(at_nos,grp))
        while i!=len(diff_groups):
            if check_grp[0].size!=0: #use predefined list for metals, metalloids etc and loop over to find its in a specific group
                groups_all.append(group_labels[i])
                metal_index.append(check_grp[0])
            break
    metal_index=np.hstack(metal_index)

    nn_object = CrystalNN()
    coord_num=[]
    coord_type=[]
    metal_atnos=[]
    dist_all=[]
    water_bonded=[]
    metal_atnos_uniq=[]
    coord_type_uniq=[]
    coord_num_uniq=[]
    distance_uniq=[]
    final_metalidx_uniq=[]
    #water_bonded_uniq=[]
    for l in range(len(metal_index)):
        atom_dist=[]
        metal_atnos.append(at_nos[metal_index[l]])
        coord_num.append(nn_object.get_cn(struct,metal_index[l]))
        coord_type.append(nn_object.get_cn_dict(struct,metal_index[l]))
        neighbors = nn_object.get_nn_info(struct,metal_index[l])
        # if "O" in nn_object.get_cn_dict(struct,metal_index[l]): #check if O is bonded to two H, ie is there a water molecule around
        #     for num1 in range(len(neighbors)):
        #         temp=list(neighbors[num1]['site'])
        #         print (temp)
        #         if str(temp[0]).split(" ")[-1]=="O":
        #             oxygen_index=neighbors[num1]['site_index']
        #             if str(nn_object.get_cn_dict(struct,oxygen_index))=="{'H': 2}":
        #                 water_bonded.append("Present")
        #             else:
        #                 water_bonded.append("Not Present")
        # else:
        #     water_bonded.append("Not Present")
        for num in neighbors:
            dist=num["site"].distance(struct.sites[metal_index[l]])
            atom_dist=np.append(atom_dist,dist)
            atom_dist=np.round(atom_dist,decimals=3)
        dist_all.append(atom_dist)
    metal_atnos=np.asarray(metal_atnos)
    coordination_data=list(zip(metal_atnos,coord_num,coord_type))
    coordination_data_uniq=[]
    indx_uniq=[]
    for i in range(len(coordination_data)):
        if coordination_data[i] not in coordination_data_uniq:
            coordination_data_uniq.append(coordination_data[i])
            indx_uniq.append(i)
    #print (indx_uniq,metal_atnos)
    for l in indx_uniq:
            metal_atnos_uniq.append(metal_atnos[l])
            coord_num_uniq.append(coord_num[l])
            coord_type_uniq.append(coord_type[l])
            distance_uniq.append(dist_all[l])
            final_metalidx_uniq.append(metal_index[l])
    #oxd_num_uniq, alloxdnum, allatno, vasporcif = struct_oxdnumber(contcarfile, where, metal_atnos_uniq)
    grouplabels_all=[]
    struct_grplabel=[]
    for K,ele in zip(diff_groups,group_labels):
        for i in range(len(K)):
            grouplabels_all.append(ele)
    diff_groupsall=np.hstack(diff_groups)
    for num in metal_atnos_uniq:
        indx=np.where(diff_groupsall==num )
        indx=int(indx[0])
        struct_grplabel.append(grouplabels_all[indx])
    df_structinfo['coordinfo.metal_groupinfo']=[struct_grplabel]
    df_structinfo['coordinfo.metal_atno']=[metal_atnos_uniq]
    df_structinfo['coordinfo.metal_coordnum']=[coord_num_uniq]
    df_structinfo['coordinfo.metal_coordtype']=[coord_type_uniq]
    df_structinfo['coordinfo.metal_coorddistances']=[distance_uniq]
    df_structinfo['coordinfo.fromvasporcif?'] = where
    # if "Present" in water_bonded:
    #     df_structinfo['metal_waterbonded']="Present"
    # else:
    #     df_structinfo['metal_waterbonded']="Not Present"
    # df_structinfo['oxidation_num'] = [oxd_num_uniq]
    # df_structinfo['alloxdnum'] = [alloxdnum]
    # df_structinfo['allatno'] = [allatno]
    # df_structinfo['fromvasporcif?'] = vasporcif
    return df_structinfo,final_metalidx_uniq

def struct_nodeinfo_oxdnumber(contcarfile,where):
    metal_atnos=[]
    grouplabels_all=[]
    struct_grplabel=[]
    struct=Structure.from_file(contcarfile) #read CONTCAR
    runner = OximachineRunner()
    features=runner.run_oximachine(struct)
    metal_index=features['metal_indices']
    at_nos=np.asarray(struct.atomic_numbers)

    if (32 in at_nos) or (51 in at_nos) :
        at_nos=np.asarray(struct.atomic_numbers) #get all atomic numbers
        groups_all=[]
        metal_index=[]
        for i,grp in enumerate(diff_groups):
            check_grp=np.where(np.isin(at_nos,grp))
            while i!=len(diff_groups):
                if check_grp[0].size!=0: #use predefined list for metals, metalloids etc and loop over to find its in a specific group
                    groups_all.append(group_labels[i])
                    metal_index.append(check_grp[0])
                break
        metal_index=np.hstack(metal_index)
    for l in range(len(metal_index)):
        metal_atnos.append(struct.atomic_numbers[metal_index[l]])
    metal_atnos=np.asarray(metal_atnos)
    metal_atnos_oxinum=metal_atnos[(metal_atnos!=32)&(metal_atnos!=51)]
    indxatnos_uniq=np.unique(metal_atnos,return_index=True)[1]
    indxatnos_uniq_oxinum=np.unique(metal_atnos_oxinum,return_index=True)[1]
    metal_atnos_uniq=[metal_atnos[indx] for indx in sorted(indxatnos_uniq)]
    metal_atnos_uniq_oxinum=[metal_atnos_oxinum[indx] for indx in sorted(indxatnos_uniq_oxinum)]
    final_metalidx_uniq1=[metal_index[indx] for indx in sorted(indxatnos_uniq)]
    oxd_num=features['prediction']
    if len(oxd_num)!=0:
        oxd_num_uniq=[oxd_num[indx] for indx in sorted(indxatnos_uniq_oxinum)]
    if len(oxd_num)==0:
        oxd_num_uniq=oxd_num
    for K,ele in zip(diff_groups,group_labels):
        for i in range(len(K)):
            grouplabels_all.append(ele)
    #print ("res",grouplabels_all)
    diff_groupsall=np.hstack(diff_groups)
    for num in metal_atnos_uniq:
        indx=np.where(diff_groupsall==num )
        indx=int(indx[0])
        struct_grplabel.append(grouplabels_all[indx])
    df_structinfo1['metal_groupinfo']=[struct_grplabel]
    df_structinfo1['metal_atno']=[metal_atnos_uniq]
    df_structinfo1['oxidation_num']=[oxd_num_uniq]
    df_structinfo1['alloxdnum']=[features['prediction']]
    df_structinfo1['allatno']=[metal_atnos]
    df_structinfo1['fromvasporcif?']=where
    return df_structinfo1,final_metalidx_uniq1


df_structinfo,final_metalidx_uniq=struct_nodeinfo("./testfiles/EMUDAT_FSR.vasp","vasp")
print (df_structinfo,final_metalidx_uniq)
# df_structinfo1,final_metalidx_uniq1=struct_nodeinfo_oxdnumber("./testfiles/EMUDAT_FSR.vasp","vasp")
# print (df_structinfo1,final_metalidx_uniq1)

