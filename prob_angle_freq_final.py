#!/usr/bin/python
#script started on: Wed 06 May 2015 04:32:05 PM IST 
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

########################################
#           FUNCTION: Rotator          #
########################################
#perform a rotation of points around an arbitrary axis:
#NOTE: direction of axis vector is from axis[0] to axis[1].
#NOTE: this function works on 3 dimensions or less.
#NOTE: angles must be in radian.
def Rotator(rotate_points, axis, angle):
    
    #turn a right-handed angle into a left-handed angle:
    angle = -angle
    
    #create a rotation matrix based on 'angle' term:
    cos_A = math.cos(angle)
    sin_A = math.sin(angle)
    Ux = axis[0]
    Uy = axis[1]
    Uz = axis[2]
    R = [[cos_A +Ux*Ux*(1-cos_A),
          Ux*Uy*(1-cos_A) -Uz*sin_A,
          Ux*Uz*(1-cos_A) +Uy*sin_A],
          
         [Uy*Ux*(1-cos_A) +Uz*sin_A,
          cos_A +Uy*Uy*(1-cos_A),
          Uy*Uz*(1-cos_A) -Ux*sin_A],
          
         [Uz*Ux*(1-cos_A) -Uy*sin_A,
          Uz*Uy*(1-cos_A) +Ux*sin_A,
          cos_A +Uz*Uz*(1-cos_A)]]
    R = numpy.matrix(R)

    #perform rotations on coordinates:
    rotate_points = rotate_points*R
 
    #output coordinates rotated along a given axis:
    return rotate_points

########################################
#     FUNCTION: Brute-force rotator    #
########################################
def BFR(points, rotation_axis, minimization_axis):
    
    #set minimization axis definition:
    if minimization_axis == [1, 0, 0]:
        minimization_axis = 0
    if minimization_axis == [0, 1, 0]:
        minimization_axis = 1
    if minimization_axis == [0, 0, 1]:
        minimization_axis = 2

    #perform rotation along rotation_axis:
    min_RMSD = 999999
    min_angle = 0
    for angle in range(0,360):
        OP = Rotator(points, rotation_axis, angle*math.pi/180)
        OP = OP.tolist()
        distances_squared = []
        for i in range(0,len(OP)):
            #check Squared-Deviation for points w.r.t minimization_axis:
            distances_squared.append(math.pow(OP[i][minimization_axis],2))
        RMSD = math.sqrt(sum(distances_squared)/len(distances_squared))
        if RMSD < min_RMSD:
            min_RMSD = RMSD
            min_angle = angle

    #print "min angle #1: "+str(min_angle)
    OP = Rotator(points, rotation_axis, min_angle*math.pi/180)
    points = OP.tolist()
    return points

########################################
# FUNCTION: calculate angle bet 2 vec  #
########################################

def Dot_angle(vector_1, vector_2):
    dotproduct = [vector_1[0]*vector_2[0], vector_1[1]*vector_2[1], vector_1[2]*vector_2[2]]
    dotproduct = dotproduct[0]+dotproduct[1]+dotproduct[2]
    length_1 = math.sqrt(math.pow(vector_1[0],2)+math.pow(vector_1[1],2)+math.pow(vector_1[2],2))
    length_2 = math.sqrt(math.pow(vector_2[0],2)+math.pow(vector_2[1],2)+math.pow(vector_2[2],2))
    return float(180)/math.pi*math.acos(float(dotproduct)/float(length_1*length_2))

########################################
# FUNCTION:CONVERT ANGLE IN 360 FORMAT #
########################################

def change_angle(X_angle, Y_angle):
     if Y_angle<90:
        angle=X_angle
     if Y_angle>90:
        angle=360-X_angle
     return angle   
######################################
#  PLOTTIN FREQUENCY MAIN FUNCTION   #
######################################

def plot_freq(protein_copy):
 protein=[]

 for i in range(0, len(protein_copy)):
        #extract individual amino acid residues, store as 2D lists:
        #X, Y, Z, residue:
        protein.append([float(protein_copy[i][30:38]), float(protein_copy[i][38:46]), float(protein_copy[i][46:54]), protein_copy[i][21:26]])

 #move helix to origin:
 mean = [0,0,0]
 for i in range(0, len(protein)):
    for j in range(0,3):
        mean[j] = mean[j] +protein[i][j]

 for i in range(0,3):
    mean[i] = float(mean[i])/len(protein)

 for i in range(0, len(protein)):
    for j in range(0,3):
        protein[i][j] = float(protein[i][j]) -mean[j]

 protein_coordinates = []
 for i in range(0,len(protein)):
    protein_coordinates.append(protein[i][0:3])

 X_axis = [1, 0, 0]
 Y_axis = [0, 1, 0]
 Z_axis = [0, 0, 1]
 #perform rotation along X-axis:
 protein_coordinates = BFR(protein_coordinates, X_axis, Y_axis)

 #perform rotation along Y-axis:
 protein_coordinates = BFR(protein_coordinates, Y_axis, X_axis)
 
 #check if helix is inverted (from Z-coordinates):
 if protein_coordinates[0][2] < protein_coordinates[-1][2]:
    #rotate helix 180-degrees around Y-axis:
    coordinates = Rotator(protein_coordinates, Y_axis, math.pi)
    protein_coordinates=coordinates.tolist()
 #print output as PDB file:
 OP_data = []
 #debugging protein_coordinates
 #print protein_coordinates[0][0]
 #sys.exit() 
 for i in range(0, len(protein_coordinates)):
    coordinate_string = "%8.3f%8.3f%8.3f" % (protein_coordinates[i][0], protein_coordinates[i][1], protein_coordinates[i][2])
 #   print coordinate_string
    OP_string = protein_copy[i][0:30]+coordinate_string+protein_copy[i][55:80]
    OP_data.append(OP_string)

 X_axis = [1, 0, 0]
 Y_axis = [0, 1, 0]
 Z_axis = [0, 0, 1]

 #finding the mean angle of polar residues
 j=0
 angle_sum=0;
 for i in range(0, len(protein_coordinates)):
    if 'ARG' in OP_data[i][16:20] or 'LYS' in OP_data[i][16:20]or 'ASP' in OP_data[i][16:20] or 'GLU' in OP_data[i][16:20] or 'GLN' in OP_data[i][16:20] or 'ASN' in OP_data[i][16:20] or 'HIS' in OP_data[i][16:20] or 'SER' in OP_data[i][16:20] or 'THR' in OP_data[i][16:20] or 'TYR' in OP_data[i][16:20] or 'CYS' in OP_data[i][16:20] or 'MET' in OP_data[i][16:20] or 'TRP' in OP_data[i][16:20] :
       if 'CA' in OP_data[i][13:15]:
          coord=[protein_coordinates[i][0],protein_coordinates[i][1],0]
          X_angle=Dot_angle(coord,X_axis)
          Y_angle=Dot_angle(coord,Y_axis)
          new_angle=change_angle(X_angle,Y_angle)
          angle_sum=new_angle+angle_sum
          j=j+1
 angle_mean = float(angle_sum)/j

 freq_cnt_polar=[]
 freq_cnt_apolar=[]
 
 #initialising freq_cnt to 0
 for i in range(0,360):
    freq_cnt_polar.append(0)
    freq_cnt_apolar.append(0)

 #finding dot_angle and frequency
 for i in range(0, len(protein_coordinates)):
         if 'ARG' in OP_data[i][16:20] or 'LYS' in OP_data[i][16:20]or 'ASP' in OP_data[i][16:20] or 'GLU' in OP_data[i][16:20] or 'GLN' in OP_data[i][16:20] or 'ASN' in OP_data[i][16:20] or 'HIS' in OP_data[i][16:20] or 'SER' in OP_data[i][16:20] or 'THR' in OP_data[i][16:20] or 'TYR' in OP_data[i][16:20] or 'CYS' in OP_data[i][16:20] or 'MET' in OP_data[i][16:20] or 'TRP' in OP_data[i][16:20] :
               if 'CA' in OP_data[i][13:15]:
                    coord=[protein_coordinates[i][0],protein_coordinates[i][1],0]
                    X_angle=Dot_angle(coord,X_axis)
                    Y_angle=Dot_angle(coord,Y_axis)               
                    new_angle=change_angle(X_angle,Y_angle)
                    theta=360-(new_angle-angle_mean)
                    if theta >= 360 :
                         theta=theta-360
                    val=math.trunc(theta)
                    freq_cnt_polar[val]=freq_cnt_polar[val]+1
 for i in range(0, len(protein_coordinates)):
         if 'ALA' in OP_data[i][16:20] or 'ILE' in OP_data[i][16:20]or 'LEU' in OP_data[i][16:20] or 'PHE' in OP_data[i][16:20] or 'VAL' in OP_data[i][16:20] or 'PRO' in OP_data[i][16:20] or 'GLY' in OP_data[i][16:20] :
               if 'CA' in OP_data[i][13:15]:
                    coord=[protein_coordinates[i][0],protein_coordinates[i][1],0]
                    X_angle=Dot_angle(coord,X_axis)
                    Y_angle=Dot_angle(coord,Y_axis)               
                    new_angle=change_angle(X_angle,Y_angle)
                    theta=360-(new_angle-angle_mean)
                    if theta >= 360 :
                         theta=theta-360
                    val=math.trunc(theta)
                    freq_cnt_apolar[val]=freq_cnt_apolar[val]+1
 return freq_cnt_polar, freq_cnt_apolar

########################################
#            MAIN FUNCTION             #
########################################

#usage/error:
if len(sys.argv) != 2:
    print "usage: helix_rotator.py <helix.pdb>"
    sys.exit()

from os import listdir
from os.path import isfile, join
onlyfiles = [ f for f in listdir('/home/prathitha/Desktop/pdb_files') if isfile(join('/home/prathitha/Desktop/pdb_files',f)) ]

#read the input helix into a double-list:
protein_copy=[]
for i in range(0,len(onlyfiles)):
    PDB_file=open('/home/prathitha/Desktop/pdb_files/'+onlyfiles[i])
    file_data=PDB_file.readlines()
    for j in range(0,len(file_data)):
        if file_data[j][0:6] != 'ATOM  ':
            del file_data[j]
    protein_copy.append(file_data)

prob_cnt=[]
freq_cnt_polar=[]
freq_cnt_apolar=[] 
#initialising prob_cnt to 0
for i in range(0,360):
   prob_cnt.append(0.0)
   freq_cnt_polar.append(0)
   freq_cnt_apolar.append(0)
#executing all the pdb files one by one
for i in range(0, len(protein_copy)):
   return_value_polar, return_value_apolar=plot_freq(protein_copy[i])
   for j in range(0, 360):
      freq_cnt_polar[j]=freq_cnt_polar[j]+return_value_polar[j]
      freq_cnt_apolar[j]=freq_cnt_apolar[j]+return_value_apolar[j]
#print freq_cnt_polar, freq_cnt_apolar
for i in range(0,360):
         if (freq_cnt_polar[i]+freq_cnt_apolar[i])==0:
              if i<=90 or i>=270 and i<360 :
                  freq_cnt_polar[i]=freq_cnt_polar[i]+1
              else:
                  freq_cnt_apolar[i]=freq_cnt_apolar[i]+1
         prob_cnt[i]=float(freq_cnt_polar[i])/(float(freq_cnt_polar[i]+freq_cnt_apolar[i]))
for i in range(0,360):
   print prob_cnt[i]
'''
#extract inputPDB_file, store as list of objects:
protein_copy = []
for i in range(0, len(PDB_data)):
  for j in range(0, len(PDB_data[i])):
     if "ATOM  " in PDB_data[i][j][0:6]:
        protein_copy.append(PDB_data[i][j])
'''

