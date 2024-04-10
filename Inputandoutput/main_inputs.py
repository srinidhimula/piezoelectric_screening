from opt_piezo_check import *
from opt_piezo_inputs import *
import csv
###intialising essential paths, variables and files
qmofbase='/home/srinidhi/qmof_v13/vasp_files'
potcar_path='/home/srinidhi/software/VASP/POTCAR/PAW-PBE54/'
maindir=os.path.abspath(os.getcwd()) #highthroughput folder
mofgroupdir=sys.argv[1].split('/')[-1]
grouppath = os.path.join(maindir,'mof_calcs',mofgroupdir)
print ("main dir",maindir)
print ("group folder",grouppath)
#mofrefcodes=np.loadtxt(maindir+'/mof_lists/'+mofgroupdir,dtype=str)
#print (mofrefcodes)

"""calling functions from other scripts"""

def opt_inp():
    mofrefcodes=create_dirs(maindir,grouppath,mofgroupdir)
    INCARopt_edit(grouppath,qmofbase)
    CONTCAR_to_POSCARsym(grouppath,qmofbase,mofrefcodes)
    POTCAR_combined(grouppath,potcar_path)
    kpoints(grouppath,dens=1000)

def piezo_inp():
    failed_optfinal=checkopt_outcar(grouppath)
    piezo_inputs(grouppath,failed_optfinal)
    print ("Opt not done",failed_optfinal)

if sys.argv[2]=="opt":
    print ("Optimisation inputs creating")
    opt_inp()
if sys.argv[2]=="piezo":
    print ("Piezo inputs creating")
    piezo_inp()

