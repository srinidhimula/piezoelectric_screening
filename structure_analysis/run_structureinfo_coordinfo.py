from structureinfo import *
from getborncharges import *
qmofpath="/home/srinidhi/qmof_v13/vasp_files/"
maindir="/home/srinidhi/BEC_work/highthroughput/"
mofgroupdir=sys.argv[1].split('/')[-1]
grouppath = os.path.join(maindir,'mof_calcs',mofgroupdir)
os.chdir(grouppath)
for dirs in os.listdir('.'):
    print (dirs)
    path_piezo=os.path.join(grouppath,dirs,'piezo_isym')
    qmofpathdir=os.path.join(qmofpath,dirs)
    if os.path.exists(path_piezo+'/results.txt'):
        bec_allmetalavg=[]
        os.chdir(path_piezo)
        try:
            df_structinfo,final_metalidx_uniq=struct_nodeinfo("POSCAR",where="vasp") #Gives information about coordination around the metal
        except:
            df_structinfo,final_metalidx_uniq=struct_nodeinfo(qmofpathdir+"CONTCAR",where="cif")
        bec_all,bec_errornorm=borncharges_frmoutcar('POSCAR','OUTCAR') #Gets the Born charges from OUTCAR of all atoms and norm of BEC matrix
        if len(bec_all)!=0:
            for m in range(len(final_metalidx_uniq)):
                #bec_metal=bec_all[final_metalidx_uniq[m]])
                bec_metaldiag=bec_all[final_metalidx_uniq[m]].diagonal()
                bec_metalavg=np.average(bec_metaldiag)
                bec_metalavg=round(bec_metalavg,3)
                bec_allmetalavg.append(bec_metalavg)
        #bec_allmetalavg=np.asarray(bec_allmetalavg)
        df_structinfo['BEC_errornorm']=bec_errornorm
        df_structinfo['BEC_metalavg']=[bec_allmetalavg]
        df_structinfo['info_MOFCode']=dirs
        df_structinfo=df_structinfo.set_index('info_MOFCode')
        df_structinfo.to_csv(maindir+"/structure_coordinationinfo_set123.csv",mode='a',header=False)
