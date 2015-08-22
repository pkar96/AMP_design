import os
import csv
import sys
import copy
import math
import numpy
import string
import random
import shutil
import getopt
import subprocess

###########################
#       FIND MINIMUM      #
###########################

def mini(protein, residue):
     minvalue=100000
     minindex=-1
     residue_vector=[0,0,0]
     for i in range (0,len(protein)):
           for j in range(0,3):
              residue_vector[j]=protein[i][j]-residue[j]
           distance= residue_vector[0]*residue_vector[0]+residue_vector[1]*residue_vector[1]+residue_vector[2]*residue_vector[2]
           if minvalue > distance:
              minvalue=distance
              minindex=i
     return minindex

############################
#  FINDING SECOND MINIMUN  #
############################

def mini2(protein , residue):
     minvalue=100000
     minvalue_2=100000
     minindex=-1
     minindex_2=-1
     residue_vector=[0,0,0]
     for i in range (0,len(protein)):
           for j in range(0,3):
              residue_vector[j]=protein[i][j]-residue[j]
           distance= residue_vector[0]*residue_vector[0]+residue_vector[1]*residue_vector[1]+residue_vector[2]*residue_vector[2]
           if minvalue >= distance:
              minvalue_2=minvalue
              minvalue=distance
              minindex_2=minindex
              minindex=i
     if minindex_2==-1 :
        minindex_2=minindex
     return minindex_2

#############################
#   INDEXING THE COUNTER    #
#############################

def array(res):
    if 'GLY' in res:
        return 0
    if 'ALA' in res:
        return 1
    if 'PRO' in res:
        return 2
    if 'VAL' in res:
        return 3
    if 'LEU' in res:
        return 4
    if 'ILE' in res:
        return 5
    if 'MET' in res:
        return 6
    if 'SER' in res:
        return 7
    if 'THR' in res:
        return 8
    if 'CYS' in res:
        return 9
    if 'ASN' in res:
        return 10
    if 'GLN' in res:
        return 11
    if 'PHE' in res:
        return 12
    if 'TYR' in res:
        return 13
    if 'TRP' in res:
        return 14
    if 'LYS' in res:
        return 15
    if 'ARG' in res:
        return 16
    if 'HIS' in res:
        return 17
    if 'ASP' in res:
        return 18
    if 'GLU' in res:
        return 19
######################################
#  MAIN FREQUENCY COUNTING FUNCTION  #
######################################

def freq_cnt(protein_copy):
 protein=[]
 for i in range(0, len(protein_copy)):
        #extract individual amino acid residues, store as 2D lists:
        #X, Y, Z, residue:
        protein.append([float(protein_copy[i][30:38]), float(protein_copy[i][38:46]), float(protein_copy[i][46:54]), int(protein_copy[i][22:26]) ,protein_copy[i][12:16], protein_copy[i][16:20]])
 '''
 print protein[0][4]+'m'
 print protein[1][4]+'m'
 sys.exit()
 '''
 freq_up=[]
 freq_down=[]
 freq_side=[]

 #initialisng to zero
 for i in range (0, 20):
    freq_up.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    freq_down.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    freq_side.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

 #calculating side frequency

 for i in range(0,len(protein)):
     if ' C  ' == protein[i][4]:
         residue=[protein[i][0], protein[i][1], protein[i][2]]
         protein_base=[]
         base_name=protein[i][5]
         base=array(base_name)
         for j in range(0,len(protein)):
             if ' N  ' == protein[j][4] and protein[j][3]!=i:
                 protein_base.append([float(protein[j][0]),float(protein[j][1]),float(protein[j][2])])
         min_index=mini(protein_base,residue)
         for k in range (0 , len(protein)):
               if protein[k][3]== (min_index+1):
                  break
         residue_name=protein[k][5]
         index=array(residue_name)
         freq_side[base][index]=freq_side[base][index]+1
         freq_side[index][base]=freq_side[index][base]+1

 #calculating up frequency

 for i in range(0,len(protein)):
     if ' H  ' == protein[i][4]:
         residue=[protein[i][0], protein[i][1], protein[i][2]]
         protein_base=[]
         base_name=protein[i][5]
         base=array(base_name)
         for j in range(0,len(protein)):
             if ' O  ' == protein[j][4] and protein[j][3]!=i:
                 protein_base.append([float(protein[j][0]),float(protein[j][1]),float(protein[j][2])])
         min_index=mini(protein_base,residue)
         for k in range (0 , len(protein)):
               if protein[k][3]== (min_index+1):
                  break
         residue_name=protein[k][5]
         index=array(residue_name)
         freq_up[base][index]=freq_up[base][index]+1
         freq_down[index][base]=freq_down[index][base]+1

 return freq_side

##############################
#      MAIN FUNCTION         #
##############################

#usage/error:
if len(sys.argv) != 2:
    print "usage: helix_rotator.py <helix.pdb>"
    sys.exit()

from os import listdir
from os.path import isfile, join
onlyfiles = [ f for f in listdir('/home/prathitha/Desktop/pdb_files') if isfile(join('/home/prathitha/Desktop/pdb_files',f)) ]

protein_copy=[]

#read the input helix into a double-list:
for i in range(0,len(onlyfiles)):
    PDB_file=open('/home/prathitha/Desktop/pdb_files/'+onlyfiles[i])
    file_data=PDB_file.readlines()
    for j in range(0,len(file_data)):
        if file_data[j][0:6] != 'ATOM  ':
            del file_data[j]
    protein_copy.append(file_data)
    
freq=[]
 
#initialising freq to 0
for i in range (0, 20):
    freq.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
#executing all the pdb files one by one
for i in range(0, len(protein_copy)):
   return_value=freq_cnt(protein_copy[i])
   for j in range(0, 20):
    for k in range(0, 20):
      freq[j][k]=freq[j][k]+return_value[j][k]     
for i in range(0,20):
 print freq[i]

'''#extract inputPDB_file, store as list of objects:

protein_copy = []
for i in range(0, len(PDB_data)):
  for j in range(0, len(PDB_data[i])):
     if "ATOM  " in PDB_data[i][j][0:6]:
        protein_copy.append(PDB_data[i][j])
'''

