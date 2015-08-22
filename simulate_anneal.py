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

ch=[['b','b','b','w'],['w','b','w','b'],['b','w','b','w'],['w','b','w','b']
]
def score_func(ch):
     score=0
     for i in range (1 , 3):
        for j in range (1,3):
           if ch[i][j]=='b' and ch[i][j+1]=='b' or ch[i][j]=='w' and ch[i][j+1]=='w':
                  score=score+1 
           if ch[i][j]=='b' and ch[i][j-1]=='b'or ch[i][j]=='w' and ch[i][j+1]=='w':
                  score=score+1 
           if ch[i][j]==ch[i-1][j]:
                  score=score+1
           if ch[i][j]==ch[i+1][j]:
                  score=score+1
     if ch[0][0]==ch[0][1]: 
          score=score+1
     if ch[0][3]==ch[0][2]:
          score=score+1
     if ch[3][0]==ch[3][1]:
          score=score+1
     if ch[3][3]==ch[3][2]:
          score=score+1
     if ch[0][0]==ch[1][0]: 
          score=score+1
     if ch[0][3]==ch[1][3]:
          score=score+1
     if ch[3][0]==ch[2][0]:
          score=score+1
     if ch[3][3]==ch[2][3]:
          score=score+1
     for i in range(1,3):
         if ch[i][0]==ch[i-1][0]:
           score=score+1
         if ch[i][0]==ch[i+1][0]:
           score=score+1
         if ch[i][0]==ch[i][1]:
           score=score+1
         if ch[i][3]==ch[i-1][3]:
           score=score+1
         if ch[i][3]==ch[i+1][3]:
           score=score+1
         if ch[i][3]==ch[i][2]:
           score=score+1
     for i in range(1,3):
         if ch[0][i]==ch[0][i-1]:
           score=score+1
         if ch[0][i]==ch[0][i+1]:
           score=score+1
         if ch[0][i]==ch[1][i]:
           score=score+1
         if ch[3][i]==ch[3][i-1]:
           score=score+1
         if ch[3][i]==ch[3][i+1]:
           score=score+1
         if ch[3][i]==ch[2][i]:
           score=score+1
     return score

def random_coord():
    i=random.randrange(0,4,1)
    j=random.randrange(0,4,1)
    return [i,j]    

def swap(ch, row,column,coord):
    ch=copy.deepcopy(ch)
    r1=coord[0]
    c1=coord[1]
    temp=ch[r1][c1]
    ch[r1][c1]=ch[row][column]
    ch[row][column]=temp
    return ch

def make_move(score,score_best,T):
    delta=score-score_best

    if(delta >0):
        return 1
    else:
        r=random.random()
        if r <= math.exp(delta/T) :
            return 1

def update_temp(T0,Tn,k,N):
     A = (T0-Tn)*float(N+1)/N
     B = T0 -A
     T1 = A/(k+1) +B
     return T1

def simulated_anneal(ch):
       init_temp=1
       temp=1
       final_temp=0.000000001
       i=1
       row=2
       column=2
       score_best=score_func(ch)
       print score_best
       while temp>=final_temp:
           coord=random_coord()
           new_ch=swap(ch,row,column,coord)   
           new_score=score_func(new_ch)
           result=make_move(new_score, score_best, temp)
           if result==1:
               ch=new_ch
               row=coord[0]
               column=coord[1]
               score_best=new_score
           temp=update_temp(init_temp,final_temp,i,1000)
           i=i+1 
       print i, score_best
       return ch

ch=simulated_anneal(ch)
for i in range(0,4):
    print ch[i]
           
    
