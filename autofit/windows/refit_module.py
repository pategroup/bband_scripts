import subprocess
import os
from easygui import *
import sys
from multiprocessing import Process
import re
import numpy
import random
import string
import math
from collections import OrderedDict
import shutil
from scipy.interpolate import *

""""
Python Triples Fitter

Written by Ian Finneran, Steve Shipman at NCF

based on previous work by the Pate Lab at UVa

Please comment any changes you make to the code here:

refit module:
-stripping out lots of functionality (using v15 of autofit as a base) to be run from within the GUI instead of individually.
Removed comments describing previous versions of stand-alone autofit because not all of the things described are available here.
This should take a particular fit from a previous autofit run, add new transitions to the fit (passed by the GUI), turn on new
constants to vary (list passed by GUI), and run SPFIT.  Creates files refit.fit, refit.int, refit.par, refit.var, refit.lin.
If the user accepts the fit, these should be copied and moved into a "refits" directory, but in the interest of having all user
interaction go through the GUI, that functionality has been commented out here.


"""
def int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="25.8", temperature="298", flag="default"):#generates SPCAT input file
    input_file = ""
    input_file += "Molecule \n"
    input_file += "0  91  %s  %s  %s  %s  %s %s  %s\n"%(Q_rot, J_min, J_max,inten,inten,freq, temperature)
    input_file += " 001  %s \n" % u_A
    input_file += " 002  %s \n" % u_B
    input_file += " 003  %s \n" % u_C

    if flag == "default":
        fh_int = open("default.int", "w")

    if flag == "refit":
        fh_int = open("refit.int","w")

    fh_int.write(input_file)
    fh_int.close()


def var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag):#generates SPCAT input file

    if main_flow == 'Normal species':
        dA = str(0.02*float(A))  #These are very rough estimates of the uncertainty on the rotational constants.  May need to be considerably refined.
        dB = str(0.01*float(B))
        dC = str(0.01*float(C))
        dDJ = str(0.1*float(DJ))
        dDJK = str(0.1*float(DJK))
        dDK = str(0.1*float(DK))
        ddJ = str(0.3*float(dJ))
        ddK = str(0.3*float(dK))

    elif main_flow == 'Isotopologues':
        dA = str(0.005*float(A))  #Assume 0.5% uncertainties for isotopologues if normal species is known.  Maybe unreasonable.
        dB = str(0.005*float(B))
        dC = str(0.005*float(C))
        dDJ = str(0.005*float(DJ))
        dDJK = str(0.005*float(DJK))
        dDK = str(0.005*float(DK))
        ddJ = str(0.005*float(dJ))
        ddK = str(0.005*float(dK))

    if flag == "uncert":
        fh_var = open("default.var",'w')

    if flag == "refit":
        fh_var = open("refit.var",'w')

    input_file = ""
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   8  430   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
    input_file +="a   1  1  0  99  0  1  1  1  1  -1   0\n"
    input_file += "           10000  %s %s \n" %(A,dA)
    input_file += "           20000  %s %s \n" %(B, dB)
    input_file += "           30000  %s %s \n" %(C, dC)
    input_file += "             200  %s %s \n" %(DJ, dDJ)
    input_file += "            1100  %s %s \n" %(DJK, dDJK)
    input_file += "            2000  %s %s \n" %(DK, dDK)
    input_file += "           40100  %s %s \n" %(dJ, ddJ)
    input_file += "           41000  %s %s \n" %(dK, ddK)
    fh_var.write(input_file)
    fh_var.close()

def par_writer_refit(A,B,C,DJ,DJK,DK,dJ,dK,constants_to_vary): #generates SPFIT par file

    DJ_flag = 0
    DJK_flag = 0
    DK_flag = 0
    dJ_flag = 0
    dK_flag = 0

    for entry in constants_to_vary:
        if entry == 'DJ':
            DJ_flag = 1
        if entry == 'DJK':
            DJK_flag = 1
        if entry == 'DK':
            DK_flag = 1
        if entry == 'dJ':
            dJ_flag = 1
        if entry == 'dK':
            dK_flag = 1

    dA = str(0.5*float(A))  #These allow A, B, and C to vary by 50% and the distortions to vary by a factor of 10.  Too much, too little?
    dB = str(0.5*float(B))
    dC = str(0.5*float(C))
    dDJ = str(abs(10.0*float(DJ)*DJ_flag))
    dDJK = str(abs(10.0*float(DJK)*DJK_flag))
    dDK = str(abs(10.0*float(DK)*DK_flag))
    ddJ = str(abs(10.0*float(dJ)*dJ_flag))
    ddK = str(abs(10.0*float(dK)*dK_flag))

    if float(DJ) == 0 and DJ_flag == 1: # 100 MHz of uncertainty should be good enough for anybody.
        dDJ = '1.00000000E+002'
    if float(DJK) == 0 and DJK_flag == 1:
        dDJK = '1.00000000E+002'
    if float(DK) == 0 and DK_flag == 1:
        dDK = '1.00000000E+002'
    if float(dJ) == 0 and dJ_flag == 1:
        ddJ = '1.00000000E+002'
    if float(dK) == 0 and dK_flag == 1:
        ddK = '1.00000000E+002'
    

    input_file = ""
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   8  500   5    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n" # don't choose more than 497 check transitions or it will crash.
    input_file +="a   1  1  0  50  0  1  1  1  1  -1   0\n"
    input_file += "           10000  %s %s \n" % (A,dA)
    input_file += "           20000  %s %s \n" % (B,dB)
    input_file += "           30000  %s %s \n" % (C,dC)
    input_file += "             200  %s %s \n" % (DJ,dDJ)
    input_file += "            1100  %s %s \n" % (DJK,dDJK)
    input_file += "            2000  %s %s \n" % (DK,dDK)
    input_file += "           40100  %s %s \n" % (dJ,ddJ)
    input_file += "           41000  %s %s \n" % (dK,ddK)
    fh_par = open("refit.par",'w')
    fh_par.write(input_file)
    fh_par.close()

def lin_writer_refit(assignment_list): #writes a lin file for use with SPFIT

    input_file = ""#the next part adds in the three peaks to be fit

    for line in assignment_list:
        input_file += line[1][0:2]+' '+line[1][2:4]+' '+line[1][4:6]+' '+\
                      line[2][0:2]+' '+line[2][2:4]+' '+line[2][4:6]+'                      '+str(line[0])+' '+line[3]+' 1.0000\n'

    fh_lin = open("refit.lin",'w')
    fh_lin.write(input_file)
    fh_lin.close()        

def run_SPCAT_refit(): 
    a = subprocess.Popen("SPCAT refit", stdout=subprocess.PIPE, shell=False)
    a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen

def cat_reader(freq_high,freq_low,flag): #reads output from SPCAT

    if flag == "default":
        fh = open("default.cat")

    if flag == "refit":
        fh = open("refit.cat")

    linelist = []
    for line in fh:
        if line[8:9]==".": 
            freq = line[3:13]
            inten = line[22:29]
            qnum_up = line[55:61]
            qnum_low = line[67:73]
            uncert = line[13:21]
            if float(freq)> freq_low and float(freq)<freq_high:#<<<<<<<<<<<<<<<<<<<<
                linelist.append((inten,freq, qnum_up, qnum_low,uncert))
    linelist.sort()
    fh.close()
    return linelist
    
def trans_freq_refit_reader(peaklist):
    pred_peaks = cat_reader(1000000, 0,flag="refit")
    pred_trans_for_refit = []

    for line in peaklist:
        for peak in pred_peaks:
            if line[2] == peak[2] and line[3] == peak[3]:
                pred_trans_for_refit.append((peak[1],peak[2],peak[3]))

    return pred_trans_for_refit


def match_to_peaklist(pred_trans,peaklist,threshold):

    peak_list_1 = peaklist[0:int(len(peaklist)/4)]#splits peaks into 4 parts to speed up processing
    peak_list_2 = peaklist[int(len(peaklist)/4):int(len(peaklist)/2)]
    peak_list_3 = peaklist[int(len(peaklist)/2):int(len(peaklist)*0.75)]
    peak_list_4 = peaklist[int(len(peaklist)*0.75):len(peaklist)]

    p_1 = float(peak_list_1[0][0])
    p_2 = float(peak_list_1[-1][0])
    p_3 = float(peak_list_2[0][0])
    p_4 = float(peak_list_2[-1][0])
    p_5 = float(peak_list_3[0][0])
    p_6 = float(peak_list_3[-1][0])
    p_7 = float(peak_list_4[0][0])
    p_8 = float(peak_list_4[-1][0])

    best_match_freqs=[]

    for x in range(len(pred_trans)): #matches predicted peaks to peaks in experimental peak list <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        weight = '0.50'
        current_peak = float(pred_trans[x][0])
        
        if current_peak>p_1 and current_peak<p_2:#conditionals to find proper portion of experimental peak list to loop through
            peaks_section = peak_list_1
            peak = p_2
        elif current_peak>p_3 and current_peak<p_4:
            peaks_section = peak_list_2
            peak = p_4
        elif current_peak>p_5 and current_peak<p_6:
            peaks_section = peak_list_3
            peak = p_6
        elif current_peak>p_7 and current_peak<p_8:
            peaks_section = peak_list_4
            peak = p_8
        else:
            peaks_section = peaklist
            peak = p_8
        
        old_omc = 100000.0
        for peak_freq,peak_inten in peaks_section: #find nearest peak in actual spectrum to the given predicted peak
            omc = abs(current_peak-float(peak_freq))
            if omc <= old_omc:
                temp_freq = peak_freq
                old_omc = omc

        if old_omc > threshold: # If the best match is too far off, we don't want to have a bad line in the fit file.  Threshold (in MHz) is a user-determined parameter.
            weight = '0.00'

        best_match_freqs.append((temp_freq,pred_trans[x][1],pred_trans[x][2],weight))

    return best_match_freqs

def refine_fits(job_name,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fit_to_refine,DJ,DJK,DK,dJ,dK,peaklist,main_flow,extra_peaks,constants_to_vary):

    if isotope_ID == "NS only":
        f2=open('input_data_%s.txt'%(job_name),'r')
    else:
        f2=open('input_data_%s_%s.txt'%(job_name,isotope_ID),'r')
    
    refitting_peaks_flag = 0
    refit_peaks_list=[]
    refit_check_peaks_list=[]

    for line in f2:
          
        if line.split() != []:
                                
            if line.split()[0] == "trans_1:" or line.split()[0] == "trans_2:" or line.split()[0] == "trans_3:":
                refitting_peaks_flag = 0
                clean = line[12:53]
                re_split = clean.split("', '")
                tuples = tuple(re_split)
                refit_peaks_list.append(tuples)

            if refitting_peaks_flag == 1:
                clean = line[2:43]
                re_split = clean.split("', '")
                tuples = tuple(re_split)
                refit_check_peaks_list.append(tuples)
            
            if line.split()[0] == "Check":
                refitting_peaks_flag = 1

    f2.close()

    int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="refit")

    result_choice = fit_to_refine.split()[0]
    A_fit = fit_to_refine.split()[6]
    B_fit = fit_to_refine.split()[7]
    C_fit = fit_to_refine.split()[8]

    fit_peaklist = refit_peaks_list + refit_check_peaks_list
    fit_peaklist = list(OrderedDict.fromkeys(fit_peaklist)) # This removes duplicates if some transitions are both check and fit transitions, which might bias SPFIT towards hitting those in particular.

    var_writer(A_fit,B_fit,C_fit,DJ,DJK,DK,dJ,dK,main_flow,flag="refit")
    run_SPCAT_refit()

    if len(extra_peaks) > 0:
        for entry in extra_peaks:
            clean = entry[2:43]
            re_split = clean.split("', '")
            tuples = tuple(re_split)
            extra_peaks_clean.append(tuples)
        extra_peaks = extra_peaks_clean

        unique_extra_peaks = []
        for entry in extra_peaks:
            unique = 1
            for fit_peak in fit_peaklist:
                if (entry[2] == fit_peak[2]) and (entry[3] == fit_peak[3]):
                    unique = 0
            if unique == 1:
                unique_extra_peaks.append(entry)

        fit_peaklist = fit_peaklist + unique_extra_peaks

    else
        pass

    updated_trans = trans_freq_refit_reader(fit_peaklist) # Finds updated predicted frequencies with improved A, B, and C estimates.

    par_writer_refit(A_fit,B_fit,C_fit,DJ,DJK,DK,dJ,dK,constants_to_vary)
    best_matches = match_to_peaklist(updated_trans,peaklist,threshold) # Assigns closest experimental peak frequencies to transitions
    lin_writer_refit(best_matches)

    a = subprocess.Popen("spfit0 refit", stdout=subprocess.PIPE, shell=False)
    a.stdout.read()#used to let SPFIT finish

    #SPFIT_results = open("refit.fit",'r')    # These segments (display and user acceptance of fit and copying of files of accepted fits) should probably be pushed back over to the GUI
    #codebox(msg='SPFIT has finished.',text=SPFIT_results)
    #SPFIT_results.close()

    #fit_decision = buttonbox(msg='Fitting has finished. Would you like to accept the current fit or try again while changing inclusion threshold or allowing different distortions to vary?', choices=('Accept Current Fit','Try Again'))
    #if fit_decision == 'Accept Current Fit':
    #    shutil.copyfile("refit.fit","refits/refit%s.fit"%(str(result_choice)))
    #    shutil.copyfile("refit.int","refits/refit%s.int"%(str(result_choice)))
    #    shutil.copyfile("refit.lin","refits/refit%s.lin"%(str(result_choice)))
    #    shutil.copyfile("refit.par","refits/refit%s.par"%(str(result_choice)))
    #    shutil.copyfile("refit.var","refits/refit%s.var"%(str(result_choice)))


def refit(job_name,main_flow,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fit_to_refine,DJ,DJK,DK,dJ,dK,peaklist,extra_peaks,threshold):
    
    os.chdir(job_name) # Need to figure out what to do here; this may need to be a full path or something.  For now assume that this will work...
    
    if main_flow == 'Normal species':
        #isotope_ID = "NS only" # Allows us to use the same base function (refine_fits) for both the NS only and the isotopologue branches.
        refine_fits(job_name,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fit_to_refine,DJ,DJK,DK,dJ,dK,peaklist,main_flow,extra_peaks,constants_to_vary,threshold)

    if main_flow == 'Isotopologues':
    
        os.chdir(isotope_ID)
        refine_fits(job_name,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fit_to_refine,DJ,DJK,DK,dJ,dK,peaklist,main_flow,extra_peaks,constants_to_vary,threshold)
        os.chdir(os.pardir)

    os.chdir(os.pardir) # Return to starting directory when done





    
