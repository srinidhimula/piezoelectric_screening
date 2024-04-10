import sys
from pathlib import Path
import numpy as np
from structureinfo import *
from numpy import linalg as LA

def borncharges_frmoutcar(poscar_path,outcar_path):
    # Read list of elements
    with open(poscar_path) as poscar:
        first_line = next(poscar)
    atom_list = first_line.split()
    atoms = []
    for atom_count in atom_list:
        element = "".join(x for x in atom_count if x.isalpha())
        count = int("".join(x for x in atom_count if x.isdigit()))
        atoms += [element] * count
    # Read BEC data
    with open(outcar_path) as outcar:
        for line in outcar:
            if line ==  " BORN EFFECTIVE CHARGES (including local field effects) (in e, cummulative output)\n":
                #print (line)
                break
        else:
            print("No section with Born effective charges found")
        bec_data = []
        for line in outcar:
            #print (line)
            if line == "\n":
                break
            bec_data.append(line)
    bec_allatoms = []
    if len(bec_data)!=0:
        bec_data = bec_data[1:]
        # Process matrices
        for start in range(1, len(bec_data), 4):
            matrix = np.loadtxt(bec_data[start:start+3])
            bec_allatoms.append(matrix[:, 1:])
        bec_allatoms=np.asarray(bec_allatoms)
        bec_allatoms1=bec_allatoms.reshape(bec_allatoms.shape[0], (bec_allatoms.shape[1]*bec_allatoms.shape[2]))
        sum=np.sum(bec_allatoms1,axis=0)
        sum=sum.reshape(3,3)
        bec_errornorm=LA.norm(sum)
    elif len(bec_data)==0:
        bec_errornorm="BEC not calculated"
    return bec_allatoms,bec_errornorm
