#!/usr/bin/env python3
import sys,os,subprocess
import numpy as np
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifWriter,CifParser
from pymatgen.core.structure import IStructure
from pymatgen.io.vasp.inputs import Kpoints
import shutil
potcars_suffix={'Na':'_pv','K': '_sv','Ca': '_sv','Sc': '_sv','Ti': '_sv','V': '_sv','Cr': '_pv','Mn': '_pv','Ga': '_d','Ge': '_d',
'Rb': '_sv','Sr': '_sv','Y': '_sv','Zr': '_sv','Nb': '_sv','Mo': '_sv','Tc': '_pv','Ru': '_pv','Rh': '_pv','In': '_d','Sn': '_d',
'Cs': '_sv','Ba': '_sv','Pr': '_3','Nd': '_3','Pm': '_3','Sm': '_3','Eu': '_3','Gd': '_3','Tb': '_3','Dy': '_3','Ho': '_3',
'Er': '_3','Tm': '_3','Yb': '_3','Lu': '_3','Hf': '_pv','Ta': '_pv','W': '_sv','Tl': '_d','Pb': '_d','Bi': '_d','Po': '_d',
'At': '_d','Fr': '_sv','Ra': '_sv','Eu':'_3','Yb':'_3','W':'_sv'}

def create_dirs(maindir,grouppath,mofgroupdir):
    #os.chdir('mof_calcs')
    if os.path.isdir(grouppath):
        print ("mofgroupdir exists")
    else:
        os.mkdir(grouppath)
    mofrefcodes=np.loadtxt(maindir+'/mof_lists/'+mofgroupdir,dtype=str)
    print (mofrefcodes)
    for item in mofrefcodes[:,0]:
        mofpath=os.path.join(grouppath,item)
        #print (mofpath)
        if os.path.isdir(mofpath):
            print ("mofdir exists")
        else:
            os.mkdir(mofpath)
    return mofrefcodes

def INCARopt_edit(grouppath,qmofbase):
    os.chdir(grouppath)
    print ("Processing INCAR files")
    for dirs in os.listdir('.'):
        os.chdir(dirs)
        qmof_incarpath=os.path.join(qmofbase,dirs,'INCAR')
        os.system('cp %s INCAR'%qmof_incarpath)
        lines = open('INCAR','r').readlines()
        with open('INCAR','w') as f:
            for line in lines:
                if 'LCHARG' in line:
                    f.write(' LCHARG= .FALSE. \n')
                elif 'LASPH' in line:
                    f.write(' LASPH = .FALSE. \n')
                elif 'LWAVE' in line:
                    f.write(' LWAVE = .FALSE. \n')
                elif 'LAECHG' in line:
                    f.write(' LAECHG = .FALSE. \n')
                elif 'SYMPREC' in line:
                    f.write(' SYMPREC = 1.00E-05 \n')
                elif 'ISYM' in line:
                    f.write(' ISYM = 2 \n')
                elif 'NCORE' in line:
                    f.write(' NCORE = 16 \n')
                elif 'EDIFF' in line:
                    f.write(' EDIFF = 1.00e-07 \n')
                elif 'NSW' in line:
                    f.write(' NSW = 500 \n')
                else:
                    f.write(line)
            f.write('\n')
            f.write('\n KPAR = 4')
            f.write('\n NELMIN = 3')
            f.write('\n IBRION = 2')
            f.write('\n ISIF = 3')
            
        os.chdir('../')
    os.chdir(grouppath)   

def CONTCAR_to_POSCARsym(grouppath,qmofbase,mofrefcodes):
    os.chdir(grouppath)
    print ("Processing POSCAR files")
    for dirs in os.listdir('.'):
        os.chdir(dirs)
        print("entering ",dirs)
        qmof_contcarpath=os.path.join(qmofbase,dirs,'CONTCAR')
        os.system('cp %s POSCAR_P1'%qmof_contcarpath)
        p=Poscar.from_file('POSCAR_P1',check_for_POTCAR=False)
        w =CifWriter(p.structure)
        w.write_file('tmp.cif')
        sgnum=subprocess.check_output([sys.executable,"/home/srinidhi/bin/vasp_inputs/findsym_becwork.py","tmp.cif"])
        sgnum=int(sgnum)
        sgloc=np.where(mofrefcodes[:,0]==dirs)
        sgloc=np.asscalar(sgloc[0])
        parser=CifParser('tmp_'+str(mofrefcodes[sgloc,1])+'.cif')
        print (dirs,mofrefcodes[sgloc,1],sgnum)
        structure=parser.get_structures()[0]
        structure.to(filename="tmp"+str(mofrefcodes[sgloc,1])+".cif", fmt="cif")
        structure.to(filename="POSCAR",fmt="poscar")
        os.chdir('../')
    os.chdir(grouppath)

def POTCAR_combined(grouppath,potcar_path):
    os.chdir(grouppath)
    print ("Processing POTCAR files")
    for dirs in os.listdir('.'):
        os.chdir(dirs)
        print("entering ",dirs)
        content = open('POSCAR','r').readlines()
        atomline=content[5]
        atom_symbols=atomline.split()
        with open("POTCAR","w") as outfile:
            for l in range(len(atom_symbols)):
                if atom_symbols[l] in potcars_suffix.keys():
                    with open(potcar_path+str(atom_symbols[l])+potcars_suffix[atom_symbols[l]]+"/POTCAR") as infile:
                        contents = infile.read()
                        outfile.write(contents)
                else:
                    with open(potcar_path+str(atom_symbols[l])+"/POTCAR") as infile:
                        contents = infile.read()
                        outfile.write(contents)
        #print (atom_symbols)
        os.chdir('../')
    os.chdir(grouppath)

def kpoints(grouppath,dens):
    os.chdir(grouppath)
    print ("Processing KPOINTS files")
    for dirs in os.listdir('.'):
        os.chdir(dirs)
        print("entering ",dirs)
        struct_pymat=IStructure.from_file("POSCAR")
        kpts=Kpoints.automatic_density(struct_pymat,kppa=int(dens),force_gamma=False)
        file=open("KPOINTS","w")
        file.write(str(kpts))
        file.close()
        os.chdir('../')

def INCARpiezo_edit(path_piezo):
    #os.chdir(path_piezo)
    lines = open(path_piezo+'/INCAR','r').readlines()
    with open(path_piezo+'/INCAR','w') as f:
        for line in lines:
            if 'NSW' in line:
                f.write(' #NSW= 400 \n')
            elif 'NCORE' in line:
                f.write(' #NCORE = 32 \n')
            elif 'IBRION' in line:
                f.write(' IBRION = 6 \n')
            elif 'ISIF' in line:
                f.write(' ISIF = 2 \n')
            elif 'NELMIN' in line:
                f.write(' #NELMIN = 3 \n')
            elif 'LWAVE' in line:
                f.write(' LWAVE = .FALSE. \n')
            elif 'MAGMOM' in line:
                f.write(' #MAGMOM \n')
            elif 'LAECHG' in line:
                f.write(' LAECHG = .FALSE. \n')
            #elif 'ALGO' in line:
                #f.write(' ALGO= All \n')
            elif 'NELM' in line:
                f.write(' NELM = 500 \n')
            else:
                f.write(line)
        f.write('\n LEPSILON = .TRUE.')
    #os.chdir(grouppath)


def piezo_inputs(grouppath,failed_optfinal):
    os.chdir(grouppath)
    for dirs in os.listdir('.'):
        if dirs in failed_optfinal:
            print ("Opt didnt finish properly",dirs)
        if dirs not in failed_optfinal:
            path_piezo=os.path.join(grouppath,dirs,'piezo_isym')
            path_opt=os.path.join(grouppath,dirs)
            #os.chdir(dirs)
            #print (dirs)
            if os.path.isdir(path_piezo):
                print ("piezo exists",dirs)
            else:
                print ("Piezo creating in ",dirs)
                os.mkdir(path_piezo)

                src_poscar=path_opt+'/CONTCAR'
                dst_poscar=path_piezo+'/POSCAR'
            
                src_potcar=path_opt+'/POTCAR'
                dst_potcar=path_piezo+'/POTCAR'

                src_kpts=path_opt+'/KPOINTS'
                dst_kpts=path_piezo+'/KPOINTS'

                src_incar=path_opt+'/INCAR'
                dst_incar=path_piezo+'/INCAR'

                shutil.copy(src_poscar,dst_poscar)
                shutil.copy(src_potcar,dst_potcar)
                shutil.copy(src_kpts,dst_kpts)
                shutil.copy(src_incar,dst_incar)
                INCARpiezo_edit(path_piezo)


