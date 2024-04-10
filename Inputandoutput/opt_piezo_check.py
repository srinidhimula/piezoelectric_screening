from pymatgen.io.vasp.outputs import *
import numpy as np
import pandas as pd
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(suppress=True)

failed = []
failed_duensw=[]
failed_optfinal=[]
failed_piezo=[]

def checkopt_osz():
    os.chdir(grouppath)
    for dirs in os.listdir('.'):
        #os.chdir(dirs)
        if os.path.isdir(dirs):
            if os.path.exists(dirs+'/OSZICAR'):
                linesOZ = open(dirs+'/OSZICAR','r').readlines()
                counts = []
                cutoff = 0
                linesI = open(dirs+'/INCAR','r').readlines()
                for line in linesI:
                    if 'NELM' in line:
                        cutoff = int(line[line.index('=')+1:])
                        break
                if len(linesOZ)>0:
                    for line in linesOZ:
                        try:
                            if line.split()[1].isdigit():
                                counts.append(int(line.split()[1]))
                        except:
                            print(dirs+' HAS BAD OSZ')
                            failed.append(dirs)
                            break
                    if cutoff in counts:
                        if dirs not in failed:
                            failed.append(dirs)
                else:
                    failed.append(dirs)
            else:
                failed.append(dirs)
    print(failed)


def checkopt_outcar(grouppath):
    os.chdir(grouppath)
    for dirs in os.listdir('.'): #initial check in outcar
        #print (dirs)
        #os.chdir(dirs)
        if os.path.isdir(dirs):
            if os.path.exists(dirs+'/OUTCAR'):
                with open(dirs+'/OUTCAR') as file:
                    contents=file.read()
                    search_sentence="reached required accuracy - stopping structural energy minimisation"
                    if search_sentence not in contents:
                        #print ("found")
                    #else:
                        failed.append(dirs)
            else:
                failed.append(dirs)
    for f in failed: # second check from oszicar to see if calc stopped due to nsw reaching limit
        os.chdir(grouppath+'/'+f)
        counts = []
        cutoff = 0
        linesI = open('INCAR','r').readlines()
        for line in linesI:
            if 'NSW' in line:
                cutoff = int(line[line.index('=')+1:])
                break
        if os.path.exists('OSZICAR'):

            linesOZ = open('OSZICAR','r').readlines()
            if len(linesOZ)>0:
                for line in linesOZ:
                    try:
                        if line.split()[0].isdigit():
                            counts.append(int(line.split()[0]))
                    except:
                        print(f+' HAS BAD OSZ')
                        failed_duensw.append(f)
                        break
                if cutoff in counts:
                    if dirs not in failed_duensw:
                        failed_duensw.append(f)
            else:
                failed_duensw.append(f)      

    print("check1,Opt failed while check in outcar",failed)
    print("check2,Opt failed due to nsw while check in oszicar",failed_duensw)
    #failed_optfinal = [item for item in failed if item not in failed_duensw]
    failed_optfinal=failed
    print ("Opt failed final, do not continue further",failed_optfinal)
    return failed_optfinal

def checkpiezo_outcar(grouppath):
    os.chdir(grouppath)
    for dirs in os.listdir('.'): #initial check in outcar
        path_piezo=os.path.join(grouppath,dirs,'piezo_isym')
        #print (dirs)
        #os.chdir(dirs)
        if os.path.isdir(path_piezo):
            if os.path.exists(path_piezo+'/OUTCAR'):
                with open(path_piezo+'/OUTCAR') as file:
                    contents=file.read()
                    search_sentence="PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        (C/m^2)"
                    if search_sentence not in contents:
                        #print ("Piezoelectric calc dint work",dirs)
                        failed_piezo.append(dirs)
            else:
                failed_piezo.append(dirs)
    print ("Piezoelectric calc dint work",failed_piezo)
    return (failed_piezo)

def piezo_e_norm(grouppath,failed_piezo,failed_optfinal,maindir,mofrefcodes): #function to obtain the e matrix and e norm from outcar
    os.chdir(grouppath)
    for dirs in os.listdir('.'):
        if dirs in failed_piezo:
            print ("Piezo didnt finish properly",dirs)
        if dirs not in failed_piezo and dirs not in failed_optfinal:
            path_piezo=os.path.join(grouppath,dirs,'piezo_isym')
            print (path_piezo)
            if os.path.isdir(path_piezo):
                os.chdir(path_piezo)
                os.system("/home/srinidhi/bin/vasptonorm.py OUTCAR ")
                piezooutcar_to_csv(maindir,dirs,mofrefcodes)
                #p=subprocess.run('vasptonorm.py'+' '+'OUTCAR',shell=True)
                #print (p)
                #os.chdir(../../)

def piezooutcar_to_csv(maindir,dirs,mofrefcodes):
    clamp_matrix=np.zeros( (3,6) )
    dynamic_matrix= np.zeros( (3,6) )
    total_e=np.zeros( (3,6) )
    results_matrix=[clamp_matrix,dynamic_matrix,total_e]
    df_clamp= pd.DataFrame()
    df_dynamic=pd.DataFrame()
    df_totale=pd.DataFrame()
    with open("results.txt",'r') as f:
        lines=f.read()
        lines = lines.replace("[", "")
        lines = lines.replace("]", "")
    with open('results.txt','w') as my_file:
        my_file.write(lines)
    lines1=open('results.txt','r').readlines()
    for index,line in enumerate(lines1):
        if "Clamped" in line:
            clamp_num=index
            #print(index, line)
        if "Dynamic" in line:
            dynamic_num=index
        if "Piezo_e" in line:
            total_e_num=index
        if "Norm" in line:
            norm_num=index
    normsplit=lines1[norm_num].split()
    clamped_norm=normsplit[-3]
    dynamic_norm=normsplit[-2]
    total_e_norm=normsplit[-1]
    list_nums=[clamp_num,dynamic_num,total_e_num]
    for l in range(len(list_nums)):
        for k in range(list_nums[l]+1,list_nums[l]+4):
            piezosplit=lines1[k].split()
            for m in range(0,6):
                a=piezosplit[m]
                results_matrix[l][k-(list_nums[l]+1)][m]=a
    df_list=[df_clamp,df_dynamic,df_totale]
    norm_name=["info_ClampedNorm","info_DynamicNorm","info_TotaleNorm"]
    list_norms=[clamped_norm,dynamic_norm,total_e_norm]
    csv_names=["clamped.csv","dynamic.csv","totale.csv"]
    for k,(res_mat,df1) in enumerate(zip(results_matrix,df_list)):
        #print (k,res_mat,df1)
        res_mat=res_mat.reshape(1,-1)
        df1=pd.DataFrame(data=res_mat,columns=["e11","e12","e13","e14","e15","e16","e21","e22","e23","e24","e25","e26","e31","e32","e33","e34","e35","e36"])
        #print (df1)
        df1['info_MOFCode']=dirs
        df1=df1.set_index('info_MOFCode')
        df1[norm_name[k]]=list_norms[k]
        df1["info_PointGroup"]=mofrefcodes[np.where(mofrefcodes[:,0]==dirs),2]
        pos=Poscar.from_file("POSCAR",check_for_POTCAR=False)
        symatoms=[str(m)+str(n) for m,n in zip(pos.site_symbols,pos.natoms)]
        df1["info_Atoms"]='-'.join(symatoms)
        if k==0:
            with open('../INCAR') as f:
                if 'ISPIN = 2' in f.read():
                    print ("ispin exists",dirs)
                    mag=Outcar("../OUTCAR")
                    df1["info_optnetmagmom"]=mag.total_mag
                    mag_piezo=Outcar("./OUTCAR")
                    df1["info_piezonetmagmom"]=mag_piezo.total_mag
                else:
                    df1["info_optnetmagmom"]=0
                    df1["info_piezonetmagmom"]=0
        df1.to_csv(maindir+"/"+csv_names[k],mode='a',header=False)
        #print ("df1",df1)
    #df_totale["info_ClampedNorm"]=clamped_norm
    #df_totale["info_DynamicNorm"]=dynamic_norm
    #print ("clamped",df_clamp)

