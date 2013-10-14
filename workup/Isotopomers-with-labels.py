'''
Created on 29/03/2012
Isotopologue autoscaler: given an input experimental and ab-initio prediction of a normal species, this program will output scaled rotational constants for all common heavy atom isotopologues.
Input coordinate file is of the format:
E# X Y Z
where E# is the atomic number of the atom, and X Y Z are the coordinates in the PA frame.
@author: Cristobal
'''
import subprocess
import os
from easygui import *
import sys
import re
import random
import string
from numpy import *
from numpy import linalg as LA
import pylab
from matplotlib.pyplot import figure, show
import math
from matplotlib import pyplot as plt
import subprocess
import numpy as np
import subprocess
from math import *
from decimal import *

msg = "Rotational constants"
title = "Enter Experimental Rotational Constants"
fieldNames = ["A/MHz","B/MHz","C/MHz"]
exp = []
exp = multenterbox(msg,title,fieldNames)
if exp == None:
    sys.exit(0);

print exp
a = []
b = []
c = []
mass = []
f = open(fileopenbox(msg="Choose Coordinates file"))
for line in f:
    x = line.split()[1]
    a.append(x)
    y = line.split()[2]
    b.append(y)
    z = line.split()[3]
    c.append(z)
    #changes letter atoms to actual masses (you may need to change this to something else depending on your file)
    if line[0]=='1':
        mass.append(1.007825037)
    if line[0]=='8':
        mass.append(15.99491464)
    if line[0]=='7':
        mass.append(14.003074008)
    if line[0]=='6':
        mass.append(12)
    if line[0]=='9':
        mass.append(18.998403)
    if line[0]=='3':
        print('i see yo silicons')
        mass.append(27.9769265325)

print mass
print a
raw3dmatrix=zeros((size(mass),4))#empty matrix for other calcs
for i,row in enumerate(mass):
    #creates matrix for later calcs
    raw3dmatrix[i,0]=row
    raw3dmatrix[i,1]=a[i]
    raw3dmatrix[i,2]=b[i]
    raw3dmatrix[i,3]=c[i]
originalmatrix=raw3dmatrix
def calcabc(dmatrix):
        
    """ Calculates the rotational const from 3d matrix"""     
        
        
    firstmomenttensor = zeros((3,3)) #empty matrix for later tensor
    totalatoms=size(dmatrix[:,0]) #total number of atoms
    temp = zeros((totalatoms,9)) #emtpy matrix to calc all partial sums
       
    for i in range(totalatoms):
        #calcs controbution of each atom to first moment tensor
        temp[i,0]= dmatrix[i,0]*(dmatrix[i,2]**2 + dmatrix[i,3]**2)
        temp[i,1]= -1*(dmatrix[i,0]*dmatrix[i,1]*dmatrix[i,2])
        temp[i,2]= -1*(dmatrix[i,0]*dmatrix[i,1]*dmatrix[i,3])
        temp[i,3]= -1*(dmatrix[i,0]*dmatrix[i,2]*dmatrix[i,1])
        temp[i,4]= dmatrix[i,0]*(dmatrix[i,1]**2 + dmatrix[i,3]**2)
        temp[i,5]= -1*(dmatrix[i,0]*dmatrix[i,2]*dmatrix[i,3])
        temp[i,6]= -1*(dmatrix[i,0]*dmatrix[i,3]*dmatrix[i,1])
        temp[i,7]= -1*(dmatrix[i,0]*dmatrix[i,3]*dmatrix[i,2])
        temp[i,8]= dmatrix[i,0]*(dmatrix[i,1]**2 + dmatrix[i,2]**2)
            
    i=0
    #sums all the controbutions of the atoms into the intertial tensor (3x3 matrix)
    firstmomenttensor[0,0]= sum(temp[:,0])   
    firstmomenttensor[0,1]= sum(temp[:,1])    
    firstmomenttensor[0,2]= sum(temp[:,2])
    firstmomenttensor[1,0]= sum(temp[:,3])
    firstmomenttensor[1,1]= sum(temp[:,4])
    firstmomenttensor[1,2]= sum(temp[:,5])
    firstmomenttensor[2,0]= sum(temp[:,6])
    firstmomenttensor[2,1]= sum(temp[:,7])
    firstmomenttensor[2,2]= sum(temp[:,8])
    
    #calcs eigenvalues and eigenvectors of the intertal tensor...eigenvectors
    #are thrown away and eigenvalues are constants A,B,C
    Iabc,dontcare=LA.eigh(firstmomenttensor)
    
    #converts to MHz from wavenumbers 
    abc= 505379.006/Iabc
    return abc#,raw3dmatrix

def distancefromcenter(dmatrix):
    
    totalmass=sum(dmatrix[:,0]) #total mass of the system
    comx=[]
    comy=[]
    comz=[]
    for row in dmatrix:
        #calcs raw distance from center of mass for each atom
        comx.append(row[0] *row[1])#/totalmass
        comy.append(row[0] *row[2])
        comz.append(row[0] *row[3])
    
    comx=sum(comx)/totalmass 
    
    comy=sum(comy)/totalmass
    
    comz=sum(comz)/totalmass
    
    #shifts the coordiantes to center of mass by value calculated above
    for i,row in enumerate(dmatrix):
        dmatrix[i,1]=dmatrix[i,1]-comx 
        dmatrix[i,2]=dmatrix[i,2]-comy
        dmatrix[i,3]=dmatrix[i,3]-comz
       
    return dmatrix
f.close()

raw3dmatrix=distancefromcenter(raw3dmatrix)

thr_abc = calcabc(raw3dmatrix) # Ab initio Rotational constant 

print thr_abc
factor = [float(exp[0])/thr_abc[0], float(exp[1])/thr_abc[1], float(exp[2])/thr_abc[2]] # Calculates conversion factor exp/thr it is necessary to float exp for some reason
print factor
thr_isotop = []
if ccbox('Would you like rotated coordinates?', 'Isotopomer predictions') ==1:
    rotcoords = open("rotated_output_coords.txt",'w')
    rotcoords.write(array_str(raw3dmatrix))
    rotcoords.close()
    
for i in range(size(originalmatrix[:,0])):
    if raw3dmatrix[i,0]==12:
        tempraw3dmatrix = originalmatrix.copy() #copy the original matrix in every iteration and so we only have monosubstituted species
        tempraw3dmatrix[i,0]=13.003354838      #singly sustitutes a 12C atom
        tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
        c13=calcabc(tempraw3dmatrix)
        temp= []
        temp.append('13C')
        temp.append(c13[0])
        temp.append(c13[1])
        temp.append(c13[2])
        thr_isotop.append(temp)  
    if raw3dmatrix[i,0]==15.99491464:
        tempraw3dmatrix = originalmatrix.copy()
        tempraw3dmatrix[i,0]=17.9991604
        tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
        o18=calcabc(tempraw3dmatrix)
        temp1= []
        temp1.append('18O')
        temp1.append(o18[0])
        temp1.append(o18[1])
        temp1.append(o18[2])
        thr_isotop.append(temp1)      
    if raw3dmatrix[i,0]==14.003074008:
        tempraw3dmatrix = originalmatrix.copy()
        tempraw3dmatrix[i,0]=15.0001088982
        tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
        n15=calcabc(tempraw3dmatrix)
        temp2= []
        temp2.append('15N')
        temp2.append(n15[0])
        temp2.append(n15[1])
        temp2.append(n15[2])
        thr_isotop.append(temp2)
    if raw3dmatrix[i,0]==27.9769265325:
        print('i got here nigga')
        tempraw3dmatrix = originalmatrix.copy()
        tempraw3dmatrix[i,0]=28.976494700
        tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
        si29=calcabc(tempraw3dmatrix)
        temp3= []
        temp3.append('29Si')
        temp3.append(si29[0])
        temp3.append(si29[1])
        temp3.append(si29[2])
        thr_isotop.append(temp3)
    if raw3dmatrix[i,0]==27.9769265325:
        print('i got here nigga2')
        tempraw3dmatrix = originalmatrix.copy()
        tempraw3dmatrix[i,0]=29.97377017
        tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
        si30=calcabc(tempraw3dmatrix)
        temp4= []
        temp4.append('30Si')
        temp4.append(si30[0])
        temp4.append(si30[1])
        temp4.append(si30[2])
        thr_isotop.append(temp4)
  #  if raw3dmatrix[i,0]==1.007825037:
  #      tempraw3dmatrix = originalmatrix.copy()
  #      tempraw3dmatrix[i,0]=2.01410178
  #      tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
  #      h2=calcabc(tempraw3dmatrix)
  #      temp3= []
  #      temp3.append('D')
  #      temp3.append(h2[0])
  #      temp3.append(h2[1])
  #      temp3.append(h2[2])
  #      thr_isotop.append(temp3)
    #if raw3dmatrix[i,0]==34.9689:
       # tempraw3dmatrix = originalmatrix.copy()
       # tempraw3dmatrix[i,0]=36.9659
       # tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
       # cl37=calcabc(tempraw3dmatrix)
       # temp3= []
       # temp3.append('37Cl')
       # temp3.append(cl37[0])
       # temp3.append(cl37[1])
       # temp3.append(cl37[2])
       # thr_isotop.append(temp3)
       # print temp3
#print originalmatrix.copy()
#print tempraw3dmatrix[:,0]
#for row in thr_isotop:
    #print row[0]*factor[0], row[1]*factor[1],row[2]*factor[2]

f = open(filesavebox(msg="Save Isotopomers"),'w') #open("Isotopomers.txt", "w")
f.write("Thr Parent %s %s %s \n" %(thr_abc[0],thr_abc[1],thr_abc[2])) #write thr rotational constant to the file
f.write("Exp Parent %s %s %s \n" %(exp[0],exp[1],exp[2])) #write Exp rotational constant to the file

for row in thr_isotop:
        f.write("     %s   %f %f %f \n" %(row[0],row[1]*factor[0],row[2]*factor[1],row[3]*factor[2])) #for every row
        #f.write("          %f %f %f \n" %(row[1],row[2],row[3]))
f.close()





