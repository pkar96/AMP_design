from rosetta import *
from toolbox import mutate_residue
rosetta.init()
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

##########################
# TO MAKE A MOVE TO NEW  # 
#    STRUCTURE OR NOT    #
##########################
#completed
def make_move(score,score_best,T):
    delta=score-score_best
    if(delta >0):
        return 1
    else:
        r=random.random()
        if r <= T :
            return 1
    return 0

######################
#  FIND THE MINIMUM  #
######################

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

###################################
# CALCULATE UP AND DOWN FREQUENCY #
###################################

def up_freq(PDB_rows):
  freq_up=[]
  freq_down=[]
  for i in range (0, 20):
    freq_up.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    freq_down.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  for i in range(0,len(PDB_rows)):
     if ' H  ' == PDB_rows[i][0]:
         residue=[PDB_rows[i][3], PDB_rows[i][4], PDB_rows[i][5]]
         protein_base=[]
         base_name=PDB_rows[i][1]
         base=array(base_name)
         for j in range(0,len(PDB_rows)):
             if ' O  ' == PDB_rows[j][0] and PDB_rows[j][2]!=i:
                 protein_base.append([float(PDB_rows[j][3]),float(PDB_rows[j][4]),float(PDB_rows[j][5])])
         min_index=mini(protein_base,residue)
         for k in range (0 , len(PDB_rows)):
               if PDB_rows[k][2]== (min_index+1):
                  break
         residue_name=PDB_rows[k][1]
         index=array(residue_name)
         freq_up[base][index]=freq_up[base][index]+1
         freq_down[index][base]=freq_down[index][base]+1
  return freq_up,freq_down           

####################################
#     CALCULATE SIDE FREQUENCY     #
####################################

def side_freq(PDB_rows):
  freq_side=[]
  for i in range (0, 20):
    freq_side.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  for i in range(0,len(PDB_rows)):
     if ' C  ' == PDB_rows[i][0]:
         residue=[PDB_rows[i][3], PDB_rows[i][4], PDB_rows[i][5]]
         protein_base=[]
         base_name=PDB_rows[i][1]
         base=array(base_name)
         for j in range(0,len(PDB_rows)):
             if ' N  ' == PDB_rows[j][0] and PDB_rows[j][2]!=i:
                 protein_base.append([float(PDB_rows[j][3]),float(PDB_rows[j][4]),float(PDB_rows[j][5])])
         min_index=mini(protein_base,residue)
         for k in range (0 , len(PDB_rows)):
               if PDB_rows[k][2]== (min_index+1):
                  break
         residue_name=PDB_rows[k][1]
         index=array(residue_name)
         freq_side[base][index]=freq_side[base][index]+1
         freq_side[index][base]=freq_side[index][base]+1
  return freq_side

##########################
#   TO CALCULATE SCORE   #
##########################

def score_func(PDB_rows, prob_up, prob_down, prob_side):
    freq_up, freq_down=up_freq(PDB_rows)
    freq_side=side_freq(PDB_rows)     
    #print freq_up[0][0],freq_down[0][0],freq_side[0][0],prob_up[0][0],prob_down[0][0],prob_side[0][0]
    #sys.exit()
    score=0
    for i in range(0, len(freq_up)):
         for j in range(0, len(freq_up[i])):
             score=score+freq_up[i][j]*prob_up[i][j]
                
    for i in range(0, len(freq_down)):
         for j in range(0, len(freq_down[i])):
             score=score+freq_down[i][j]*prob_down[i][j]
    
    for i in range(0, len(freq_side)):
         for j in range(0, len(freq_side[i])):
             score=score+freq_side[i][j]*prob_side[i][j]
    
    return score
 
############################
# PICK UP A RANDOM RESIDUE #
############################
#completed
def random_residue():
     return random.randrange(1,28,1)

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

###########################
#      HELIX ROTATOR      #
###########################

def helix_rotator(PDB_data):
      protein=[]
      for i in range(0, len(PDB_data)):
         if PDB_data[i][0:6] == "ATOM  ":
           #extract individual amino acid residues, store as 2D lists:
           #X, Y, Z, residue:
              protein.append([float(PDB_data[i][30:38]), float(PDB_data[i][38:46]), float(PDB_data[i][46:54]), PDB_data[i][21:26]])

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
      for i in range(0, len(protein_coordinates)):
           coordinate_string = "%8.3f%8.3f%8.3f" % (protein_coordinates[i][0], protein_coordinates[i][1], protein_coordinates[i][2])
           #   print coordinate_string
           OP_string = PDB_data[i][0:30]+coordinate_string+PDB_data[i][55:80]
           OP_data.append(OP_string)
      return protein_coordinates

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

################################
#   CHOOSING POLAR RESIDUES    #
################################

def polar():
    r=random.randrange(0,100,1)/100.00
    if r<= float(5.91/47.89):
       return 'R'
    if r> float(5.91/47.89) and r<= float((5.91+14.41)/47.89):
       return 'K'
    if r> float((5.91+14.41)/47.89) and r<= float((5.91+14.41+2.01)/47.89):
       return 'D'
    if r>float((5.91+14.41+2.01)/47.89) and r<= float((5.91+14.41+2.01+2.65)/47.89):
       return 'E'
    if r> float((5.91+14.41+2.01+2.65)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83)/47.89):
       return 'Q'
    if r> float((5.91+14.41+2.01+2.65+2.83)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71)/47.89):
       return 'N'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63)/47.89):
       return 'H' 
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42)/47.89):
       return 'S'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98)/47.89):
       return 'T'  
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98)/47.89) and r< float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36)/47.89):
       return 'Y'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68)/47.89):
       return 'C'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32)/47.89):
       return 'M'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32)/47.89) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98)/47.89):
       return 'W'

###########################
# CHOOSING APOLAR RESIDUE #
###########################

def apolar():
    r= random.randrange(0,100,1)/100.00
    if r<= float(9.37/51.96):
      return 'A'
    if r>float(9.37/51.96) and r<=float((9.37+7.41)/51.96):  
      return 'I'
    if r>float((9.37+7.41)/51.96) and r<=float((9.37+7.41+10.85)/51.96):
      return 'L'
    if r>float((9.37+7.41+10.85)/51.96) and r<=float((9.37+7.41+10.85+5.23)/51.96):
      return 'F'
    if r>float((9.37+7.41+10.85+5.23)/51.96) and r<=float((9.37+7.41+10.85+5.23+6.53)/51.96):
      return 'V'
    if r>float((9.37+7.41+10.85+5.23+6.53)/51.96) and r<=float((9.37+7.41+10.85+5.23+6.53+2.59)/51.96):
      return 'P'
    if r>float((9.37+7.41+10.85+5.23+6.53+2.59)/51.96) and r<=float((9.37+7.41+10.85+5.23+6.53+2.59+9.98)/51.96):
      return 'G'

######################
# CHANGE THE PROTEIN #
######################

def change_protein(PDB_data, PDB_rows, residue_change, prob_angle):
    protein = []
    for i in range(0, len(PDB_data)):
      if PDB_data[i][0:6] == "ATOM  ":
        #extract individual amino acid residues, store as 2D lists:
        #X, Y, Z, residue:
        protein.append([float(PDB_data[i][30:38]), float(PDB_data[i][38:46]), float(PDB_data[i][46:54]),int(PDB_data[i][22:26])])
    #rotate along the z-axis
    rotated_coordinates=helix_rotator(PDB_data)
    #find the mean of polar residues
    X_axis=[1,0,0]
    Y_axis=[0,1,0]
    j=0
    angle_sum=0
    for i in range(0, len(rotated_coordinates)):
      if 'ARG' in PDB_data[i][16:20] or 'LYS' in PDB_data[i][16:20]or 'ASP' in PDB_data[i][16:20] or 'GLU' in PDB_data[i][16:20] or 'GLN' in PDB_data[i][16:20] or 'ASN' in PDB_data[i][16:20] or 'HIS' in PDB_data[i][16:20] or 'SER' in PDB_data[i][16:20] or 'THR' in PDB_data[i][16:20] or 'TYR' in PDB_data[i][16:20] or 'CYS' in PDB_data[i][16:20] or 'MET' in PDB_data[i][16:20] or 'TRP' in PDB_data[i][16:20] :
        if 'CA' in PDB_data[i][13:15]:
          coord=[rotated_coordinates[i][0],rotated_coordinates[i][1],0]
          X_angle=Dot_angle(coord,X_axis)
          Y_angle=Dot_angle(coord,Y_axis)
          new_angle=change_angle(X_angle,Y_angle)
          angle_sum=new_angle+angle_sum
          j=j+1
    if j==0:
       angle_mean=0
    else:
       angle_mean = float(angle_sum)/j
    #finding the angle at which to enter
    for i in range(0, len(PDB_data)):
           if 'CA' in PDB_data[i][13:15] and residue_change==protein[i][3]:       
                break
    X_axis=[1,0,0]
    Y_axis=[0,1,0]
    coord=[rotated_coordinates[i][0],rotated_coordinates[i][1],0]
    X_angle=Dot_angle(coord,X_axis)
    Y_angle=Dot_angle(coord,Y_axis)               
    new_angle=change_angle(X_angle,Y_angle)
    theta=360-(new_angle-angle_mean)
    if theta >= 360 :
       theta=theta-360
    val=math.trunc(theta)
    r=random.randrange(0,100,1)/100.00 
    if r<=prob_angle[val]:
         new_residue=polar()
    else:
         new_residue=apolar()
    print new_residue, residue_change
    pose=pose_from_pdb("helix.pdb")
    keep_residue=pose.sequence()
    if residue_change==1:
       keep_residue=new_residue+keep_residue[1:]
    else:
       keep_residue=keep_residue[0:residue_change-1]+new_residue+keep_residue[residue_change:]
    pose1=pose_from_sequence(keep_residue,"fa_standard")
    for i in range(1, 27):
       pose1.set_phi(i,-57)
       pose1.set_psi(i,-47)
    os.remove('/media/prathitha/Local Disk1/PyRosetta.Ubuntu-12.04LTS.64Bit.release-r7/new_helix.pdb')
    pose1.dump_pdb("new_helix.pdb")
    print pose_from_pdb("new_helix.pdb")

########################
#  UPDATE TEMPERATURE  #
########################
#completed
def update_temp(T0,Tn,k,N):
     A = (T0-Tn)*float(N+1)/N
     B = T0 -A
     T1 = A/(k+1) +B
     return T1

#########################
#  SIMULATED ANNEALING  #
#########################

def simulated_anneal(PDB_data, PDB_rows, prob_up, prob_down, prob_side, prob_angle):
       init_temp=1
       temp=1
       final_temp=0.000000001
       c=1
       
       #print prob_up[0][0],prob_down[0][0],prob_side[0][0]
       #sys.exit()       
       score_best=score_func(PDB_rows, prob_up, prob_down, prob_side)
       #residue_change=random_residue()
       #change_protein(PDB_data, PDB_rows, residue_change, prob_angle)
       #sys.exit()
       while temp>=final_temp:
          residue_change=random_residue()
          change_protein(PDB_data, PDB_rows, residue_change, prob_angle)
          #read the input helix into a double-list:
          PDB_file = open('/media/prathitha/Local Disk1/PyRosetta.Ubuntu-12.04LTS.64Bit.release-r7/new_helix.pdb') 

          #extract inputPDB_file, store as list of objects:
          new_PDB_data = PDB_file.readlines()
          new_PDB_rows = []
          #print new_PDB_data
          for i in range(0, len(new_PDB_data)):
              if new_PDB_data[i][0:6] == "ATOM  ":
                 #extract individual amino acid residues, store as 2D lists:
                 #X, Y, Z, residue:
                 new_PDB_rows.append([new_PDB_data[i][12:16],new_PDB_data[i][17:20],int(new_PDB_data[i][22:26]),float(new_PDB_data[i][30:38]), float(new_PDB_data[i][38:46]), float(new_PDB_data[i][46:54])])
          new_score=score_func(new_PDB_rows, prob_up, prob_down, prob_side)
          #print new_PDB_data
          #print "score="+str(new_score)
          result=make_move(new_score, score_best, temp)
          print "result=",result
          if result==1:
              PDB_data=new_PDB_data
              PDB_rows=new_PDB_rows
              score_best=new_score
              pose=pose_from_pdb("new_helix.pdb")
              os.remove('/media/prathitha/Local Disk1/PyRosetta.Ubuntu-12.04LTS.64Bit.release-r7/helix.pdb')
              pose.dump_pdb("helix.pdb")
          temp=update_temp(init_temp,final_temp,c,1000)
          print 'score_best=',score_best
          #print "temp=",temp
          c=c+1 
          print "c=",c
       return PDB_data

########################
#    MAIN FUNCTION     #
########################

#getting pdb file
PDB_file = open(sys.argv[1])
PDB_data = PDB_file.readlines()
PDB_rows = []
#print PDB_data[0][22:26]
#sys.exit()
for i in range(0, len(PDB_data)):
   if PDB_data[i][0:6] == "ATOM  ":
      PDB_rows.append([PDB_data[i][12:16],PDB_data[i][17:20],int(PDB_data[i][22:26]),float(PDB_data[i][30:38]), float(PDB_data[i][38:46]), float(PDB_data[i][46:54])])

#gettin freq_cnt file

prob_up_file=open(sys.argv[2])
prob_up_data=prob_up_file.readlines()
prob_up=[]
for i in range(0, len(prob_up_data)):
      prob_up_in=[]
      j=0
      while j<len(prob_up_data[i]):
           prob_up_in.append(int(prob_up_data[i][j:j+3]))
           j=j+3
      prob_up.append(prob_up_in)

prob_down_file=open(sys.argv[3])
prob_down_data=prob_down_file.readlines()
prob_down=[]
for i in range(0, len(prob_down_data)):
      prob_down_in=[]
      j=0
      while j<len(prob_down_data[i]):
           prob_down_in.append(int(prob_down_data[i][j:j+3]))
           j=j+3
      prob_down.append(prob_down_in)

prob_side_file=open(sys.argv[4])
prob_side_data=prob_side_file.readlines()
prob_side=[]
for i in range(0, len(prob_side_data)):
      prob_side_in=[]
      j=0
      while j<len(prob_side_data[i]):
           prob_side_in.append(int(prob_side_data[i][j:j+3]))
           j=j+3
      prob_side.append(prob_side_in)

         
prob_angle_file=open(sys.argv[5])
prob_angle_data=prob_angle_file.readlines()
prob_angle=[]
for i in range(0, len(prob_angle_data)):
        prob_angle.append(float(prob_angle_data[i]))
#print prob_up[0][0],prob_down[0][0],prob_side[0][0]
#sys.exit()
PDB=simulated_anneal(PDB_data, PDB_rows, prob_up, prob_down, prob_side, prob_angle)

for i in range(0, len(PDB)):
      print PDB[i]


