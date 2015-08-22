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

def make_move(score,score_best,T):
    delta=score-score_best
    if(delta >0):
        return 1
    else:
        r=random.random()
        if r <= T :
            return 1
    return 0

##############################
#    SELECTING AN INDEX      #
##############################

def index(residue):
    if 'LYS' in residue:
        return 0
    if 'ARG' in residue or 'ASP' in residue or 'GLU' in residue or 'GLN' in residue or 'ASN' in residue or 'HIS' in residue or 'SER' in residue or 'THR' in residue or 'TYR' in residue or 'CYS' in residue or 'MET' in residue or 'TRP' in residue :
        return 1
    if 'ALA' in residue or 'ILE' in residue or 'LEU' in residue or 'PHE' in residue or 'VAL' in residue or 'PRO' in residue or 'GLY' in residue :
        return 2

#################################
# DISTANCE BETWEEN TWO RESIDUES #
#################################

def distance(protein,i,j):
   dist=math.sqrt(math.pow((protein[i][0]-protein[j][0]),2)+math.pow((protein[i][1]-protein[j][1]),2)+math.pow((protein[i][2]-protein[j][2]),2))
   return dist

######################
#  FREQUENCY COUNT   #
######################

def freq_cnt(PDB_data, r):
  protein=[]
  for i in range(0, len(PDB_data)):
        #extract individual amino acid residues, store as 2D lists:
        #X, Y, Z, residue:
        #print float(PDB_data[0][30:38]),float(PDB_data[0][38:46]),float(PDB_data[0][46:54])
        #sys.exit()
     if PDB_data[i][0:6] == "ATOM  ":
        protein.append([float(PDB_data[i][30:38]), float(PDB_data[i][38:46]), float(PDB_data[i][46:54]), int(PDB_data[i][22:26]) , PDB_data[i][12:16], PDB_data[i][16:20]])
  freq=[]
  #initialisng to zero
  for i in range (0, 3):
    freq.append([0,0,0])
  for i in range(0, len(protein)):
     if 'CA' in protein[i][4]:
       index_x=index(protein[i][5])
       for j in range (0, len(protein)):
        if 'CA' in protein[j][4]:
          if distance(protein,i,j) > (4*r) and distance(protein,i,j) <= (4*r+4) : 
            index_y=index(protein[j][5])
            freq[index_x][index_y]=freq[index_x][index_y]+1
  return freq     

#######################################
# TO CHOOSE THE PROBABILIITY FUNCTION #
#######################################

def choose(i,prob_02,prob_24,prob_46,prob_68):
   if i ==0:
      return prob_02
   if i==1:
      return prob_24
   if i==2:
      return prob_46
   if i==3:
      return prob_68

##########################
#   TO CALCULATE SCORE   #
##########################

def score_func(PDB_data, prob_02,prob_24,prob_46,prob_68):
    score=0.0
    for i in range(0, 4):
        freq=freq_cnt(PDB_data, i)
        prob=choose(i, prob_02,prob_24,prob_46,prob_68)
        for j in range(0, len(prob)):
             for k in range(0, len(prob[j])):
                score=score+freq[j][k]*prob[j][k]
 
    return score
 
############################
# PICK UP A RANDOM RESIDUE #
############################
#completed
def random_residue():
     return random.randrange(1,28,1)

################################
#   CHOOSING       RESIDUES    #
################################

def get_residue():
    r=random.randrange(0,100,1)/100.00
    if r<= float(5.91/100):
       return 'R'
    if r> float(5.91/100) and r<= float((5.91+14.41)/100):
       return 'K'
    if r> float((5.91+14.41)/100) and r<= float((5.91+14.41+2.01)/100):
       return 'D'
    if r>float((5.91+14.41+2.01)/100) and r<= float((5.91+14.41+2.01+2.65)/100):
       return 'E'
    if r> float((5.91+14.41+2.01+2.65)/100) and r<= float((5.91+14.41+2.01+2.65+2.83)/100):
       return 'Q'
    if r> float((5.91+14.41+2.01+2.65+2.83)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71)/100):
       return 'N'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63)/100):
       return 'H' 
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42)/100):
       return 'S'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98)/100):
       return 'T'  
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98)/100) and r< float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36)/100):
       return 'Y'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68)/100):
       return 'C'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32)/100):
       return 'M'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98)/100):
       return 'W'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37)/100):
      return 'A'
    if  r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41)/100): 
      return 'I'
    if  r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85)/100): 
      return 'L'
    if  r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23)/100): 
      return 'F'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23+6.53)/100): 
      return 'V'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23+6.53)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23+6.53+2.59)/100): 
      return 'P'
    if r> float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23+6.53+2.59)/100) and r<= float((5.91+14.41+2.01+2.65+2.83+2.71+2.63+5.42+2.98+1.36+1.68+1.32+1.98+9.37+7.41+10.85+5.23+6.53+2.59+9.98)/100): 
      return 'G'

########################################
#  RETURNING THE PRE DECIDED FREQUENCY #
########################################

def residue_pre_freq(new_residue):
    if new_residue == 'I' or new_residue == 'V' or new_residue == 'L' or new_residue == 'F' or new_residue == 'A' or new_residue == 'G' or new_residue == 'P':
      return 45
    if new_residue=='K':
      return 20
    if new_residue=='C' or  new_residue=='M' or new_residue=='W' or new_residue=='T' or new_residue=='S' or new_residue=='Y' or  new_residue=='Q' or new_residue=='N' or new_residue=='E' or new_residue=='D' or new_residue=='H' or new_residue=='R' :
      return 35
    
########################
# CHANGING THE RESIDUE #
########################

def change_residue(residue_change, new_residue):
    pose=pose_from_pdb("helix.pdb")
    keep_residue=pose.sequence()
    if residue_change==1:
       keep_residue=new_residue+keep_residue[1:]
    else:
       keep_residue=keep_residue[0:residue_change-1]+new_residue+keep_residue[residue_change:]
    return keep_residue

######################################
# RETURNING THE INDEX OF THE RESIDUE #
######################################

def index_single(new_residue):
    if new_residue=='K':
        return 0
    if new_residue=='C' or  new_residue=='M' or new_residue=='W' or new_residue=='T' or new_residue=='S' or new_residue=='Y' or  new_residue=='Q' or new_residue=='N' or new_residue=='E' or new_residue=='D' or new_residue=='H' or new_residue=='R' :
        return 1
    if new_residue == 'I' or new_residue == 'V' or new_residue == 'L' or new_residue == 'F' or new_residue == 'A' or new_residue == 'G' or new_residue == 'P':
        return 2

#####################################
#  FREQUENCY OF THE GIVEN RESIDUE   #
#####################################

def freq_residue(keep_residue,new_residue):
     freq=0
     index_val=index_single(new_residue)
     for i in range(0, len(keep_residue)):
       if index_val==0:
        if keep_residue[i]==new_residue:
          freq=freq+1
       if index_val==1:
        if keep_residue[i]=='C' or keep_residue[i]=='M' or keep_residue[i]=='W' or keep_residue[i]=='T' or keep_residue[i]=='S' or keep_residue[i]=='Y' or keep_residue[i]=='Q' or keep_residue[i]=='N' or keep_residue[i]=='E' or keep_residue[i]=='D' or keep_residue[i]=='H' or  keep_residue[i]=='R' :
          freq=freq+1
        if index_val==2:
         if keep_residue[i]=='I' or keep_residue[i]=='V' or keep_residue[i]=='L' or keep_residue[i]=='F' or keep_residue[i]=='A' or keep_residue[i]=='G' or keep_residue[i]=='P':
          freq=freq+1
     #print freq
     #sys.exit()
     return ((float(freq)/27.0)*100.0)

##########################
# CHECKING THE FREQUENCY #
##########################

def check(keep_residue, new_residue, residue_change):
    if freq_residue(keep_residue,new_residue)< residue_pre_freq(new_residue):
        return keep_residue
    else:
        new_residue=get_residue()    
        keep_residue=change_residue(residue_change, new_residue)
        keep_residue=check(keep_residue, new_residue, residue_change)
        return keep_residue

######################
# CHANGE THE PROTEIN #
######################

def change_protein(PDB_data, PDB_rows, residue_change):
    protein = []
    for i in range(0, len(PDB_data)):
      if PDB_data[i][0:6] == "ATOM  ":
        #extract individual amino acid residues, store as 2D lists:
        #X, Y, Z, residue:
        protein.append([float(PDB_data[i][30:38]), float(PDB_data[i][38:46]), float(PDB_data[i][46:54]),int(PDB_data[i][22:26])])
    new_residue=get_residue()
    #print new_residue, residue_change
    pose=pose_from_pdb("helix.pdb")
    keep_residue=pose.sequence()
    if residue_change==1:
       keep_residue=new_residue+keep_residue[1:]
    else:
       keep_residue=keep_residue[0:residue_change-1]+new_residue+keep_residue[residue_change:]
    keep_residue=check(keep_residue, new_residue, residue_change)
    pose1=pose_from_sequence(keep_residue,"fa_standard")
    for i in range(1, 27):
       pose1.set_phi(i,-57)
       pose1.set_psi(i,-47)
    os.remove('/media/prathitha/Local Disk1/PyRosetta.Ubuntu-12.04LTS.64Bit.release-r7/new_helix.pdb')
    pose1.dump_pdb("new_helix.pdb")
    #print pose_from_pdb("new_helix.pdb")

########################
#  UPDATE TEMPERATURE  #
########################

def update_temp(T0,Tn,i,N, cooling_schedule):
    if cooling_schedule == 0:
        Ti = T0 -i*(T0-Tn)/N
    
    if cooling_schedule == 1:
        Ti = T0*(Tn/T0)**(i/N)
    
    if cooling_schedule == 2:
        A = (T0-Tn)*float(N+1)/N
        B = T0 -A
        Ti = A/(i+1) +B
    
    if cooling_schedule == 3:
        print "warning: cooling_schedule '3' does not work as described"
        print "switching to cooling_schedule '5'"
        cooling_schedule = 5
    
    if cooling_schedule == 4:
        Ti = (T0-Tn)/(1+math.exp(0.01*(i-N/2))) +Tn;
    
    if cooling_schedule == 5:
        Ti = 0.5*(T0 -Tn)*(1+math.cos(i*math.pi/N)) +Tn
    
    if cooling_schedule == 6:
        Ti = 0.5*(T0-Tn)*(1-math.tanh(i*10/N-5)) +Tn;
    
    if cooling_schedule == 7:
        Ti = (T0-Tn)/math.cosh(i*10/N) +Tn;
    
    if cooling_schedule == 8:
        A = (1/N)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i)

    if cooling_schedule == 9:
        A = (1/N**2)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i**2);
    
    return Ti

#########################
#  SIMULATED ANNEALING  #
#########################

def simulated_anneal(PDB_data, PDB_rows, prob_02, prob_24, prob_46, prob_68, cooling_schedule):
       init_temp=1
       temp=1
       final_temp=0.01
       c=1
       
       #print prob_up[0][0],prob_down[0][0],prob_side[0][0]
       #sys.exit()       
       score_best=score_func(PDB_data, prob_02, prob_24, prob_46, prob_68) 
       #residue_change=random_residue()
       #change_protein(PDB_data, PDB_rows, residue_change, prob_angle)
       #sys.exit()
       while temp>final_temp:
          residue_change=random_residue()
          change_protein(PDB_data, PDB_rows, residue_change) 
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
          new_score=score_func(new_PDB_data, prob_02, prob_24, prob_46, prob_68)
          #print new_PDB_data
          #print "score="+str(new_score)
          result=make_move(new_score, score_best, temp) 
          #print "result=",result
          if result==1:
              PDB_data=new_PDB_data
              PDB_rows=new_PDB_rows
              score_best=new_score
              pose=pose_from_pdb("new_helix.pdb")
              os.remove('/media/prathitha/Local Disk1/PyRosetta.Ubuntu-12.04LTS.64Bit.release-r7/helix.pdb')
              pose.dump_pdb("helix.pdb")
          temp=update_temp(init_temp,final_temp,c,1000, cooling_schedule)
          #print 'score_best=',score_best
          print "temp=",temp,"c=",c
          c=c+1  
       return score_best

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

prob_02_file=open(sys.argv[2])
prob_02_data=prob_02_file.readlines()
prob_02=[]
for i in range(0, 3):
      prob_02_in=[]
      for j in range(0, 3):
          prob_02_in.append(float(prob_02_data[3*i+j]))
      prob_02.append(prob_02_in)

prob_24_file=open(sys.argv[3])
prob_24_data=prob_24_file.readlines()
prob_24=[]
for i in range(0, 3):
      prob_24_in=[]
      for j in range(0, 3):
          prob_24_in.append(float(prob_24_data[3*i+j]))
      prob_24.append(prob_24_in)

prob_46_file=open(sys.argv[4])
prob_46_data=prob_46_file.readlines()
prob_46=[]
for i in range(0, 3):
      prob_46_in=[]
      for j in range(0, 3):
          prob_46_in.append(float(prob_46_data[3*i+j]))
      prob_46.append(prob_46_in)
         
prob_68_file=open(sys.argv[5])
prob_68_data=prob_68_file.readlines()
prob_68=[]
for i in range(0, 3):
        prob_68_in=[]
        for j in range(0, 3):
          prob_68_in.append(float(prob_68_data[3*i+j]))
        prob_68.append(prob_68_in)

#for i in range(0,30):
#   cooling_schedule=int(i/3)
score_best=simulated_anneal(PDB_data, PDB_rows, prob_02, prob_24, prob_46, prob_68, 9) 
print "score_best=",score_best
print pose_from_pdb("helix.pdb")
''' pose=pose_from_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAA","fa_standard")
    for i in range(1, 27):
       pose.set_phi(i,-57)
       pose.set_psi(i,-47)
    pose.dump_pdb("helix.pdb")

#for i in range(0, len(PDB)):
#      print PDB[i]
'''

