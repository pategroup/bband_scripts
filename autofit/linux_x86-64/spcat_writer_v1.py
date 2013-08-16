# -*- coding: utf-8 -*-
"""
Latest version: 11/26/2012
SPCAT INT/VAR FILE GENERATOR & SPECTRUM PREDICTOR
@author: Nathan Seifert
"""
from easygui import *
import string
import math
import subprocess
import sys
from numpy import *
import numpy as np
from matplotlib import pyplot as pp

#BELOW IS THE CONSTRUCTION CODE FOR VAR AND INT FILES
def int_writer(u_A,u_B,u_C, J_min,J_max, inten,Q_rot,freq,temp):#generates SPCAT input file
    input_file = ""
    #print "freq_max=",freq
    input_file += "Molecule \n"
    input_file += "0  91  %s  %s  %s  %s  %s %s  %s\n"%(Q_rot, J_min, J_max,inten,inten,freq, temp)
    input_file += " 001  %s \n" % u_A
    input_file += " 002  %s \n" % u_B
    input_file += " 003  %s \n" % u_C
    fh_int = open("output.int", "w")
    fh_int.write(input_file)
    fh_int.close()

def QrotA(A, B, C, T): #calculates partition function for asymtop
    Qrot = (5.3311*10**(6))*(float(T)**(1.5))*(float(B)*float(C)*float(C))**(-0.5)
    print("Partition function at temperture =" " "+T+"K is "+str(Qrot))
    return Qrot
    
   
def var_writer(A,B,C,DJ,DJK,DK,dJ,dK):# Generates standard FTMW SPCAT input file, no quadrupole
    input_file = ""
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   8  430   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
    input_file +="a   1  1  0  99  0  1  1  1  1  -1   0\n"
    input_file += "           10000  %s""E+000 1.0E-010 \n" %(A)
    input_file += "           20000  %s""E+000 1.0E-010 \n" %(B)
    input_file += "           30000  %s""E+000 1.0E-010 \n" %(C)
    input_file += "             200  %s""E-003 1.0E-010 \n" %(DJ)
    input_file += "            1100  %s""E-003 1.0E-010 \n" %(DJK) 
    input_file += "            2000  %s""E-003 1.0E-010 \n" %(DK)
    input_file += "           40100  %s""E-003 1.0E-010 \n" %(dJ)
    input_file += "           41000  %s""E-003 1.0E-010 \n" %(dK)
    fh_var = open("output.var",'w')
    fh_var.write(input_file)
    fh_var.close()

def var_writerQ(A,B,C,DJ,DJK,DK,dJ,dK,chiaa,chibbcc,spin):# Generates SPCAT input file with quadrupole, no sextics
    input_file = ""
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   10  430   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
    input_file +="a   %s  1  0  99  0  1  1  1  1  -1   0\n" %(spin)
    input_file += "           10000  %s""E+000 1.0E-010 \n" %(A)
    input_file += "           20000  %s""E+000 1.0E-010 \n" %(B)
    input_file += "           30000  %s""E+000 1.0E-010 \n" %(C)
    input_file += "             200  %s""E-003 1.0E-010 \n" %(DJ)
    input_file += "            1100  %s""E-003 1.0E-010 \n" %(DJK) 
    input_file += "            2000  %s""E-003 1.0E-010 \n" %(DK)
    input_file += "           40100  %s""E-003 1.0E-010 \n" %(dJ)
    input_file += "           41000  %s""E-003 1.0E-010 \n" %(dK)
    input_file += "       110010000  %s""E+000 1.0E-010 \n" %(chiaa)
    input_file += "       110040000  %s""E+000 1.0E-010 \n" %(chibbcc)    
    fh_var = open("output.var",'w')
    fh_var.write(input_file)
    fh_var.close()

def var_writermmw(A,B,C,DJ,DJK,DK,dJ,dK,PhiJ,PhiJK,PhiKJ,PhiK,phiJ,phiJK,phiK,LJ,LJJK,LJK,LKKJ,LK,l1,l2,l3,l4,chiaa,chibbcc,spin):# Generates SPCAT input file with sextics, one quadrupole
    input_file = ""
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   26  430   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
    input_file +="a   %s  1  0  99  0  1  1  1  1  -1   0\n" %(spin)
    input_file += "           10000  %s""E+000 1.0E-010 /A    \n" %(A)
    input_file += "           20000  %s""E+000 1.0E-010 /B    \n" %(B)
    input_file += "           30000  %s""E+000 1.0E-010 /C    \n" %(C)
    input_file += "             200  %s""E-003 1.0E-010 /-DJ   \n" %(DJ)
    input_file += "            1100  %s""E-003 1.0E-010 /-DJK  \n" %(DJK) 
    input_file += "            2000  %s""E-003 1.0E-010 /-DK   \n" %(DK)
    input_file += "           40100  %s""E-003 1.0E-010 /-dJ   \n" %(dJ)
    input_file += "           41000  %s""E-003 1.0E-010 /-dK   \n" %(dK)
    input_file += "             300  %s""E-006 1.0E-010 /PhiJ  \n" %(PhiJ)    
    input_file += "            1200  %s""E-006 1.0E-010 /PhiJK \n" %(PhiJK)
    input_file += "            2100  %s""E-006 1.0E-010 /PhiKJ \n" %(PhiKJ)
    input_file += "            3000  %s""E-006 1.0E-010 /PhiK  \n" %(PhiK)
    input_file += "           40200  %s""E-006 1.0E-010 /phiJ  \n" %(phiJ)
    input_file += "           41100  %s""E-006 1.0E-010 /phiJK \n" %(phiJK)
    input_file += "           42000  %s""E-006 1.0E-010 /phiK  \n" %(phiK)
    input_file += "             400  %s""E-009 1.0E-010 /LJ    \n" %(LJ)
    input_file += "            1300  %s""E-009 1.0E-010 /LJJK  \n" %(LJJK)
    input_file += "            2200  %s""E-009 1.0E-010 /LJK   \n" %(LJK)
    input_file += "            3100  %s""E-009 1.0E-010 /LKKJ  \n" %(LKKJ)
    input_file += "            4000  %s""E-009 1.0E-010 /LK    \n" %(LK)
    input_file += "           40300  %s""E-009 1.0E-010 /l1    \n" %(l1)
    input_file += "           50200  %s""E-009 1.0E-010 /l2    \n" %(l2)
    input_file += "           60100  %s""E-009 1.0E-010 /l3    \n" %(l3)
    input_file += "           70000  %s""E-009 1.0E-010 /l4    \n" %(l4)
    input_file += "       110010000  %s""E+000 1.0E-010 /1.5chiaa \n" %(chiaa)
    input_file += "       110040000  %s""E+000 1.0E-010 /0.25chibbcc \n" %(chibbcc)    
    fh_var = open("output.var",'w')
    fh_var.write(input_file)
    fh_var.close()

def prediction(): # Generates predicted spectrum -- code stolen from Cristobal's spectrum generator script
    filtered_cat = []
    msgpr = 'Choose your start frequency (must be greater than 0 MHz) and end frequency (must be less than or equal to the cutoff from the VAR file), and your linewidth)'
    titlepr = 'Prediction spectrum generator'
    choicespr = ['Start freq (MHz)','End freq (MHz)','HWHM Linewidth (MHz)']
    valuespr = []
    valuespr = multenterbox(msgpr,titlepr,choicespr) 
    # Start parsing routine
    start_freq = int(valuespr[0])
    stop_freq = int(valuespr[1])
    x= genfromtxt('output.cat', delimiter=[13,8,8,2,10,3,13,2,2,2,6,2,2,2])
    catx =x[:,0]
    caty =10**x[:,2]
    cat = np.column_stack((catx,caty))
    less_than_test = start_freq < cat[:,0]
    greater_than_test = stop_freq > cat[:,0]
    combined = np.column_stack((less_than_test,greater_than_test))
    tally = np.zeros(np.shape(combined)[0])
    for i,line in enumerate(combined):
        if combined[i,0] == True and combined[i,1] == True:
            tally[i] = 1
        else:
            tally[i] = 0
    counter_start = 0
    for line in tally:
        if line == 0:
            counter_start = counter_start + 1
        else:    
            break
    counter_end = len(tally)-1
    while tally[counter_end] == 0:
        counter_end = counter_end - 1
    filtered_cat = cat[counter_start:counter_end]
    # Generate and plot spectrum
    Test_Spectrum = Generate_Spectrum(start_freq,stop_freq,float(valuespr[2]),filtered_cat)
    np.savetxt('output_spectrum.prn', Test_Spectrum)
    pp.plot(Test_Spectrum[:,0],Test_Spectrum[:,1])
    pp.show()


def Generate_Spectrum(start,stop,linewidth,spectrum): # Routine to generate spectrum
    step = .02
    Nfreq = math.ceil((stop-start)/step)
    spx = np.zeros(np.shape(range(0,int(Nfreq-1)))[0])
    spy = np.zeros(np.shape(range(0,int(Nfreq-1)))[0])
    spectrumx=spectrum[:,0]
    spectrumy=spectrum[:,1]
    sp = np.column_stack((spx,spy))
    spectrum = np.column_stack((spectrumx,spectrumy))
    for n in range(0, int(Nfreq-1)):
        sp[n,0] = start + n*step
        sp[n,1] = 0
    for n in range(0,len(spectrum)-1):
        transint = spectrum[n,1]
        transfreq = spectrum[n,0]
        simstart = transfreq - 4*linewidth
        simstop = transfreq + 4*linewidth                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        nstart = math.floor((simstart-start)/step)
        nstop = math.floor((simstop-start)/step)
        for i in range(int(nstart),int(nstop)):
            sp[i,1] = sp[i,1] + transint*math.exp((((sp[i,0]-transfreq)/linewidth)**2)*(-1*math.log(2)))
            sp[i,1] = np.real(sp[i,1])
    return sp

# INIT: Choose which kind of VAR file you want to create
msg1 = "Choose the desired CAT output." 
title1 = "SPCAT input file generator"
choices = ["Prolate asymmetric, no hyperfine", "Prolate asymmetric, one quadrupolar nucleus", "Prolate asymmetric, with sextics and octics, one quadrupolar nucleus (for mmW)"]    
choice = choicebox(msg1,title1,choices)

title = "SPCAT input file generator"

# Below is the control flow for the choicebox call above
if choice == choices[0]:
    intMsg = "Enter all values below. Set all unused values to 0. Hint: Set Jmax arbitrarily high for a general prediction"
    varMsg = "Enter all values below in MHz unless otherwise specified. Set all unused values to 0. Uses Watson A-reduction."
    intNames = ["uA","uB","uC", "max freq (in GHz)", "temp (K)", "Jmax"]
    intValues = []
    intValues = multenterbox(intMsg, title, intNames)
    varNames = ["A (MHz)","B (MHz)","C (MHz)","-DJ (kHz)","-DJK (kHz)","-DK (kHz)","-dJ (kHz)","-dK (kHz)"]
    varValues = []
    varValues = multenterbox(varMsg, title, varNames)
    int_writer(intValues[0],intValues[1],intValues[2],"00",intValues[5],"-10.0",QrotA(varValues[0],varValues[1],varValues[2],intValues[4]),intValues[3],intValues[4])
    var_writer(varValues[0],varValues[1],varValues[2],varValues[3],varValues[4],varValues[5],varValues[6],varValues[7])

if choice == choices[1]:
    spinMsg = "Enter the spin of the nucleus. Use decimal notation for half-integer spins (e.g. 1.5 for 3/2)."
    spin = enterbox(spinMsg)
    if float(spin) == math.floor(float(spin)):
        spinf = 2*float(spin)+1
        print('Input spin parameter is:'" "+str(spinf))
    else: 
        print('You got here')
        spinf = 2*math.ceil(float(spin))
        print('Input spin parameter is:'" "+str(spinf))
    intMsg = "Enter all values below. Set all unused values to 0. Hint: Set Jmax arbitrarily high for a general prediction"
    varMsg = "Enter all values below in MHz unless otherwise specified. Set all unused values to 0. Uses Watson A-reduction."
    intNames = ["uA","uB","uC", "max freq (in GHz)", "temp (K)", "Jmax"]
    intValues = []
    intValues = multenterbox(intMsg, title, intNames)
    varNames = ["A (MHz)","B (MHz)","C (MHz)","-DJ (kHz)","-DJK (kHz)","-DK (kHz)","-dJ (kHz)","-dK (kHz)","1.5Chi_aa (MHz)", "0.25(Chi_bb-Chi_cc) (MHz)"]
    varValues = []
    varValues = multenterbox(varMsg, title, varNames)
    int_writer(intValues[0],intValues[1],intValues[2],"00",intValues[5],"-10",QrotA(varValues[0],varValues[1],varValues[2],intValues[4]),intValues[3],intValues[4])
    var_writerQ(varValues[0],varValues[1],varValues[2],varValues[3],varValues[4],varValues[5],varValues[6],varValues[7],varValues[8],varValues[9],str(int(spinf)))  
    
if choice == choices[2]:
    spinMsg = "Enter the spin of the nucleus. Use decimal notation for half-integer spins (e.g. 1.5 for 3/2)."
    spin = enterbox(spinMsg)
    if float(spin) == math.floor(float(spin)):
        spinf = 2*float(spin)+1
        print('Input spin parameter is:'" "+str(spinf))
    else: 
        spinf = math.ceil(float(spin))
        spinf = 2*spinf
        print('Input spin parameter is:'" "+str(spinf))   
    intMsg = "Enter all values below. Set all unused values to 0. Hint: Set Jmax arbitrarily high for a general prediction"
    varMsg = "Enter all values below in MHz unless otherwise specified. Set all unused values to 0. Uses Watson A-reduction."
    intNames = ["uA","uB","uC", "max freq (in GHz)", "temp (K)", "Jmax"]
    intValues = []
    intValues = multenterbox(intMsg, title, intNames)
    varNames = ["A (MHz)","B (MHz)","C (MHz)","-DJ (kHz)","-DJK (kHz)","-DK (kHz)","-dJ (kHz)","-dK (kHz)", "PhiJ (Hz)", "PhiJK (Hz)", "PhiKJ (Hz)", "PhiK (Hz)","phiJ (Hz)","phiJK (Hz)","phiK (Hz)","LJ (mHz)","LJJK (mHz)", "LJK (mHz)","LKKJ (mHz)","LK (mHz)","l1 (mHz)","l2 (mHz)", "l3 (mHz)","l4 (mHz)","1.5Chi_aa (MHz)", "0.25(Chi_bb-Chi_cc) (MHz)"]
    varValues = []
    varValues = multenterbox(varMsg, title, varNames)
    int_writer(intValues[0],intValues[1],intValues[2],"00",intValues[5],"-10",QrotA(varValues[0],varValues[1],varValues[2],intValues[4]),intValues[3],intValues[4])
    var_writermmw(varValues[0],varValues[1],varValues[2],varValues[3],varValues[4],varValues[5],varValues[6],varValues[7],varValues[8],varValues[9],varValues[10],varValues[11],varValues[12],varValues[13],varValues[14],varValues[15],varValues[16],varValues[17],varValues[18],varValues[19],varValues[20],varValues[21],varValues[22],varValues[23],varValues[24],varValues[25],str(int(spinf))) 

# Below is the routine to ask the user if they want a predictive spectrum output or not
if ccbox('Would you like a predicted spectrum?', 'SPCAT Prediction Script') ==1:
    a = subprocess.Popen("SPCAT output", stdout=subprocess.PIPE, shell=False)
    a.stdout.read()
    prediction()
    print('Program complete! Spectrum was generated.')
else:
    print('Program complete! No spectrum was generated.')
    sys.exit()


