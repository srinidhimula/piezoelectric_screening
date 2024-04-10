from opt_piezo_check import *
import csv

def writeheader_first(maindir):

    headerList0=["info_MOFCode","e11","e12","e13","e14","e15","e16","e21","e22","e23","e24","e25","e26","e31","e32","e33","e34","e35","e36","info_ClampedNorm","info_PointGroup","info_Atoms","info_optnetmagmom","info_piezonetmagmom"]
    with open(maindir+"/clamped.csv", 'w') as file0:
        dw = csv.DictWriter(file0, delimiter=',',fieldnames=headerList0)
        dw.writeheader()

    headerList1=["info_MOFCode","e11","e12","e13","e14","e15","e16","e21","e22","e23","e24","e25","e26","e31","e32","e33","e34","e35","e36","info_DynamicNorm","info_PointGroup","info_Atoms","info_optnetmagmom","info_piezonetmagmom"]
    with open(maindir+"/dynamic.csv", 'w') as file1:
        dw = csv.DictWriter(file1, delimiter=',',fieldnames=headerList1)
        dw.writeheader()

    headerList2=["info_MOFCode","e11","e12","e13","e14","e15","e16","e21","e22","e23","e24","e25","e26","e31","e32","e33","e34","e35","e36","info_TotaleNorm","info_PointGroup","info_Atoms","info_optnetmagmom","info_piezonetmagmom"]
    with open(maindir+"/totale.csv", 'w') as file2:
        dw = csv.DictWriter(file2, delimiter=',',fieldnames=headerList2) 
        dw.writeheader()

###intialising essential paths, variables and files
#maindir=os.path.abspath(os.getcwd()) #highthroughput folder
maindir="/home/srinidhi/BEC_work/highthroughput/"
mofgroupdir=sys.argv[1].split('/')[-1]
grouppath = os.path.join(maindir,'mof_calcs',mofgroupdir)
print ("main dir",maindir)
print ("group folder",grouppath)
mofrefcodes=np.loadtxt(maindir+'/mof_lists/'+mofgroupdir,dtype=str)
print (mofrefcodes)

"""calling functions from other scripts"""
#writeheader_first(maindir)
failed_optfinal=checkopt_outcar(grouppath)
failed_piezo=checkpiezo_outcar(grouppath)
print ("Opt not done",failed_optfinal)
print ("Piezo not done",failed_piezo)
piezo_e_norm(grouppath,failed_piezo,failed_optfinal,maindir,mofrefcodes)
