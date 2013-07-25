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

isotopologue module a:
-stripping out lots of functionality (using v15 of autofit as a base) to be run from within the GUI instead of individually.
Removed comments describing previous versions of stand-alone autofit because not all of the things described are available here.
This module takes a processed coordinates file + typical set of rotational constants, check transitions, etc. from the GUI and
performs fitting runs on all singly-substituted isotopologues of selected atoms (also passed by the GUI).  This version includes
multiple potential steps of user interaction because the initial list of fitting transitions may not work for certain isotopologues.

An alternate approach for this module would be to have it process a single isotopologue, with all of the other heavy lifting being
done by the GUI.  That can be done by isotopologue_module_b.


"""
def int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="25.8", temperature="298", flag="default"):#generates SPCAT input file
    input_file = ""
    #print "freq_max=",freq
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
        dA = str(0.005*float(A))  #Assume 1% uncertainties for isotopologues if normal species is known.  Maybe unreasonable.
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
    input_file += "            1100  %s %s \n" %(DJK, dDJK) #need to actually check numbers: SPFIT doesn't read -- as a positive!
    input_file += "            2000  %s %s \n" %(DK, dDK)
    input_file += "           40100  %s %s \n" %(dJ, ddJ)
    input_file += "           41000  %s %s \n" %(dK, ddK)
    fh_var.write(input_file)
    fh_var.close()

def run_SPCAT(): 
    a = subprocess.Popen("SPCAT default", stdout=subprocess.PIPE, shell=False)
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
    
def trans_freq_reader(trans_1, trans_2, trans_3):

    peak_1_freq = 0
    peak_2_freq = 0
    peak_3_freq = 0

    pred_peaks = cat_reader(1000000, 0, flag="default")
    for peak in pred_peaks:
        if trans_1[2] == peak[2] and trans_1[3] == peak[3]:
            peak_1_freq = peak[1]
        if trans_2[2] == peak[2] and trans_2[3] == peak[3]:
            peak_2_freq = peak[1]
        if trans_3[2] == peak[2] and trans_3[3] == peak[3]:
            peak_3_freq = peak[1]
    return peak_1_freq,peak_2_freq,peak_3_freq

def check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high):
    bad_windows = 0
    bad_1 = 0
    bad_2 = 0
    bad_3 = 0

    if (trans_1_center-peak_1_uncertainty) < freq_low:
        bad_windows = 1
        bad_1 = -1
    if (trans_1_center+peak_1_uncertainty) > freq_high:
        bad_windows = 1
        bad_1 = 1
    if (trans_2_center-peak_2_uncertainty) < freq_low:
        bad_windows = 1
        bad_2 = -1
    if (trans_2_center+peak_2_uncertainty) > freq_high:
        bad_windows = 1
        bad_2 = 1
    if (trans_3_center-peak_3_uncertainty) < freq_low:
        bad_windows = 1
        bad_3 = -1
    if (trans_3_center+peak_3_uncertainty) > freq_high:
        bad_windows = 1
        bad_3 = 1
    return bad_windows,bad_1,bad_2,bad_3

def choose_new_transition(full_list):

    msg ="Choose a new transition."
    title = "Microwave Fitting Program"
    choice = choicebox(msg, title, full_list)
        
    choice_clean = []

    clean = choice[2:55]
    re_split = clean.split("', '")
    tuples = tuple(re_split)
    choice_clean.append(tuples)

    trans_out = choice_clean[0]
    trans_out_uncert = float(trans_out[4])

    return trans_out,trans_out_uncert

def final_uncerts(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3):
    
    (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

    while bad_windows ==1:
        while bad_1 == -1:
            bad_wind_decision = buttonbox(msg='The search window for transition 1 extends below the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?', choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_1 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_1,trans_1_uncert) = choose_new_transition(full_list)
                trans_1_center = float(trans_1[1])
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_2 == -1:
            bad_wind_decision = buttonbox(msg='The search window for transition 2 extends below the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?', choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_2 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_2,trans_2_uncert) = choose_new_transition(full_list)
                trans_2_center = float(trans_2[1])
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_3 == -1:
            bad_wind_decision = buttonbox(msg='The search window for transition 3 extends below the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?', choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_3 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_3,trans_3_uncert) = choose_new_transition(full_list)
                trans_3_center = float(trans_3[1])
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_1 == 1:
            bad_wind_decision = buttonbox(msg='The search window for transition 1 extends above the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?', choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_1 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_1,trans_1_uncert) = choose_new_transition(full_list)
                trans_1_center = float(trans_1[1])
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_2 == 1:
            bad_wind_decision = buttonbox(msg='The search window for transition 2 extends above the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?', choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_2 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_2,trans_2_uncert) = choose_new_transition(full_list)
                trans_2_center = float(trans_2[1])
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_3 == 1:
            bad_wind_decision = buttonbox(msg='The search window for transition 3 extends above the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?', choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_3 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_3,trans_3_uncert) = choose_new_transition(full_list)
                trans_3_center = float(trans_3[1])
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        if (trans_1 == trans_2) or (trans_2 == trans_3) or (trans_1 == trans_3):
            same_trans_decision = buttonbox(msg='You do not have three distinct transitions.  Choose new ones manually, or quit?', choices=('Choose manually','Quit'))
            if same_trans_decision == 'Quit':
                quit()
            elif same_trans_decision == 'Choose manually':
                (highest_uncert,trans_1,trans_2,trans_3) = triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag="force manual")
                trans_1_uncert = float(trans_1[4])
                trans_2_uncert = float(trans_2[4])
                trans_3_uncert = float(trans_3[4])
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1_center,trans_1_uncert)))
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2_center,trans_2_uncert)))
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

    return trans_1,trans_2,trans_3,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty


def triples_gen(trans_1_uncert,trans_2_uncert,trans_3_uncert,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,isotopomer_count,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3):

    trans_1_center = float(trans_1[1])
    trans_2_center = float(trans_2[1])
    trans_3_center = float(trans_3[1])

    trans_1_peaks = []
    trans_2_peaks = []
    trans_3_peaks = []

#    if window_decision == 'User-defined, same for each':
#        user_flag = 1
#        same_flag = 1
#    elif window_decision == 'User-defined, different for each':
#        user_flag = 1
#    elif window_decision == 'Three times SPCAT uncertainty':
#        est_unc_flag = 1

    while trans_1_peaks == [] or trans_2_peaks == [] or trans_3_peaks == []: #this loops until there are peaks around each member of the triple
        uncertainty_flag =1

        while uncertainty_flag ==1:            
            (trans_1,trans_2,trans_3,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty) = final_uncerts(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3)
            trans_1_center = float(trans_1[1])
            trans_2_center = float(trans_2[1])
            trans_3_center = float(trans_3[1])

            trans_1_peaks = []
            trans_2_peaks = []
            trans_3_peaks = []

            for freq_p, inten_p in peaklist:
                if abs(float(trans_1_center)-float(freq_p))< peak_1_uncertainty:
                    trans_1_peaks.append((freq_p, inten_p))
                if abs(float(trans_2_center)-float(freq_p))< peak_2_uncertainty: #this bit finds peaks in the real spectrum that are near the predicted peaks
                    trans_2_peaks.append((freq_p, inten_p))
                if abs(float(trans_3_center)-float(freq_p))< peak_3_uncertainty:
                    trans_3_peaks.append((freq_p, inten_p))
            num_of_triples = len(trans_1_peaks)*len(trans_2_peaks)*len(trans_3_peaks) #this tells you how many entries there will be in the all_combo_list
                                    
            # Leave in for now.
            if isotopomer_count == 0:
                decision = buttonbox(msg='There are %s triples in this calculation. Would you like to continue, try new uncertainty, or quit?'%(str(num_of_triples)), choices=('Continue','Quit','new uncertainty'))
            
            if decision == 'Quit':
                quit() # Should change this to return and then do appropriate error handling at calling function.
            elif decision == 'Continue':
                uncertainty_flag = 0
            else:
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1_center,trans_1_uncert)))
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2_center,trans_2_uncert)))
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3_center,trans_3_uncert)))

    
    return trans_1,trans_2,trans_3,trans_1_peaks,trans_2_peaks,trans_3_peaks,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,num_of_triples,decision

def intensity_filter(full_list,peaklist,inten_low,filter_level): # Intensity filter to give more efficient triples searches for isotopologues.

    filtered_peaklist=[]
    filtered_full_list = []
    comparison_level = filter_level*float(inten_low)

    for peak in peaklist:
        if peak[1] >= comparison_level: # Only keep experimental peaks more intense than the lower cutoff.
            filtered_peaklist.append(peak)

    for entry in full_list:
        for peak in filtered_peaklist:
            temp_freq_diff = abs(float(entry[1])-float(peak[0]))
            if temp_freq_diff <= 0.5: # Looking for a strong NS peak within 0.5 MHz of its predicted value.  Too coarse, too tight, OK?
                filtered_full_list.append(entry)
                break # Only need to find one; no need to continue searching through the full experimental peak list after a hit has been found.

    if filter_level == 0:
        filtered_full_list = full_list

    if len(filtered_full_list) < 3:
        print "There aren't enough transitions of appropriate intensity close to predicted positions for an isotopologue search.  Check your NS constants, your scale factor, or your spectral data file."
        return

    return filtered_full_list

def triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag):

    if flag == "undecided":
        triple_style = buttonbox(msg='How do you want to choose fitting transitions?  You can choose from an automatically scored list, or select three transitions manually.', choices=('Automatic scoring','Manual selection'))

    elif flag == "force manual":
        triple_style = 'Manual selection'

    if triple_style == 'Automatic scoring':
        total_check_num = 10 # This is the number of peaks used to generate possible triples, ordered by intensity. 10 = 120 possibilities, 15 = 455 possibilities.
        triples_scores = []
        scaled_triples_scores = []
        max_dependence = 0
        max_RMS = 0
        max_intensity = 0

        for i in range(0,total_check_num-2):
            for j in range(i+1,total_check_num-1):
                for k in range(j+1,total_check_num):
                    boundary_penalty = 0
                    trans_1 = full_list[i]
                    trans_2 = full_list[j]
                    trans_3 = full_list[k]
                    dependence = abs(dependence_test(float(A),float(B),float(C),float(DJ),float(DJK),float(DK),float(dJ),float(dK),trans_1,trans_2,trans_3,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow))
                    worst_RMS = max(float(trans_1[4]),float(trans_2[4]),float(trans_3[4]))
                    RMS_ratio = (worst_RMS/min(float(trans_1[4]),float(trans_2[4]),float(trans_3[4])))
                    RMS_function = RMS_ratio*worst_RMS
                    intensity_avg = abs(float(trans_1[0])+float(trans_2[0])+float(trans_3[0]))/3 # For all sane T and dipoles, intensities are negative, so the greatest sum of abs(intensity) is the smallest set of peaks.
                    
                    if (abs(float(trans_1[1])-freq_low) <= 100) or (abs(float(trans_1[1])-freq_high) <= 100):
                        boundary_penalty += 1
                    if (abs(float(trans_2[1])-freq_low) <= 100) or (abs(float(trans_2[1])-freq_high) <= 100):
                        boundary_penalty += 1
                    if (abs(float(trans_3[1])-freq_low) <= 100) or (abs(float(trans_3[1])-freq_high) <= 100):
                        boundary_penalty += 1
                    triples_scores.append((i,j,k,dependence,RMS_function,intensity_avg,worst_RMS,boundary_penalty))
                    if dependence > max_dependence:
                        max_dependence = dependence
                    if RMS_function > max_RMS:
                        max_RMS = RMS_function
                    if intensity_avg > max_intensity:
                        max_intensity = intensity_avg

        for entry in triples_scores:
            scaled_dep = entry[3]/max_dependence # big dependence is good and should help the score, big values of RMS and intensity are bad and should hurt the score.
            scaled_RMS = 1 - entry[4]/max_RMS
            scaled_inten = 1 - entry[5]/max_intensity
            boundary_penalty = 1 / (1 + entry[7])
            scaled_score = (45*scaled_dep + 45*scaled_RMS + 10*scaled_inten)*boundary_penalty # Reduces score based on how many transitions are close to a spectrum edge.
            stupid_easyGUI_hack = 1 - (scaled_score/100)    # EasyGUI does a sort before displaying choices, presenting them from least to greatest.  This is a dumb way to get it to display triples choices in the order that I want it to, from best to worst.
            scaled_triples_scores.append((entry[0],entry[1],entry[2],entry[3],entry[4],entry[5],scaled_dep,scaled_RMS,scaled_inten,scaled_score,stupid_easyGUI_hack,entry[6]))

        triples_choice_list = []
        for entry in scaled_triples_scores:
            trans_1 = full_list[entry[0]]
            trans_2 = full_list[entry[1]]
            trans_3 = full_list[entry[2]]
            triples_choice_list.append((entry[10],entry[9],entry[11],trans_1[1],trans_1[2],trans_1[3],trans_2[1],trans_2[2],trans_2[3],trans_3[1],trans_3[2],trans_3[3],trans_1[4],trans_2[4],trans_3[4],entry[0],entry[1],entry[2]))

        msg = "Choose a triples combination.  The second number is the triples score, ranging from 0 to 100 with 100 being the best possible.  The third number is the highest uncertainty in MHz from the triple."
        title = "Microwave Fitting Program"
        choice = choicebox(msg,title,triples_choice_list)
        clean_choice = choice[2:-1].split(",")

        highest_uncert = float(clean_choice[2])
        trans_1 = full_list[int(clean_choice[-3])]
        trans_2 = full_list[int(clean_choice[-2])]
        trans_3 = full_list[int(clean_choice[-1])]

    elif triple_style == 'Manual selection':

        user_satisfied = 0

        while user_satisfied == 0:
            msg ="Pick any three transitions.  If more than three are chosen, only the first three will be used."
            title = "Microwave Fitting Program"
            choice = multchoicebox(msg, title, full_list)
        
            choice_clean = []
            for entry in choice:
                clean = entry[2:55]
                re_split = clean.split("', '")
                tuples = tuple(re_split)
                choice_clean.append(tuples)

            if len(choice_clean) >= 3:
                trans_1 = choice_clean[0]
                trans_2 = choice_clean[1]
                trans_3 = choice_clean[2]

                trans_1_uncert = float(trans_1[4])
                trans_2_uncert = float(trans_2[4])
                trans_3_uncert = float(trans_3[4])

                highest_uncert = max(trans_1_uncert,trans_2_uncert,trans_3_uncert)

                dependence = dependence_test(float(A),float(B),float(C),float(DJ),float(DJK),float(DK),float(dJ),float(dK),trans_1,trans_2,trans_3,temperature,freq_high, freq_low,u_A,u_B,u_C,main_flow)
        
                decision1 = buttonbox(msg='These three peaks have a linear dependence of %s and the highest uncertainty in MHz is %s. A linear dependence between -1.0 and 1.0 is probably too low. Would you like to continue, or try three new peaks?'%(str(dependence),highest_uncert), choices=('Continue','Try three new peaks'))

                if decision1 == 'Continue':
                    user_satisfied = 1
                else:
                    pass
            else:
                less_than_3 = buttonbox(msg='You did not choose enough transitions.', choices=('Try again','Quit'))
                if less_than_3 == 'Quit':
                    quit()
                else:
                    pass
    
    return highest_uncert,trans_1,trans_2,trans_3

def dependence_test(A,B,C,DJ,DJK,DK,dJ,dK,trans_1,trans_2,trans_3,T,freq_high, freq_low,u_A,u_B,u_C,main_flow):
    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temperature=T, flag = "default")

    var_writer(A+(2),B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()
    high_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3)

    var_writer(A-(2),B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()
    low_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3)

    dv1A = (float(high_peak_freq[0])-float(low_peak_freq[0]))/4
    dv2A = (float(high_peak_freq[1])-float(low_peak_freq[1]))/4
    dv3A = (float(high_peak_freq[2])-float(low_peak_freq[2]))/4

    var_writer(A,B+(2),C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()
    high_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3)

    var_writer(A,B-(2),C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()
    low_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3)

    dv1B = (float(high_peak_freq[0])-float(low_peak_freq[0]))/4
    dv2B = (float(high_peak_freq[1])-float(low_peak_freq[1]))/4
    dv3B = (float(high_peak_freq[2])-float(low_peak_freq[2]))/4

    var_writer(A,B,C+(2),DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()
    high_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3)

    var_writer(A,B,C-(2),DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()
    low_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3)

    dv1C = (float(high_peak_freq[0])-float(low_peak_freq[0]))/4
    dv2C = (float(high_peak_freq[1])-float(low_peak_freq[1]))/4
    dv3C = (float(high_peak_freq[2])-float(low_peak_freq[2]))/4

    var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")	# This re-runs SPCAT at the initial constants so that other things that read from default.cat are correct after this function is executed.
    run_SPCAT()

    matrix = numpy.array([(dv1A,dv1B,dv1C),(dv2A,dv2B,dv2C),(dv3A,dv3B,dv3C)])
    return numpy.linalg.det(matrix)

def calcabc(dmatrix):
        
    """ Calculates the rotational const from 3d matrix"""     
        
    firstmomenttensor = numpy.zeros((3,3)) #empty matrix for later tensor
    totalatoms=numpy.size(dmatrix[:,0]) #total number of atoms
    temp = numpy.zeros((totalatoms,9)) #emtpy matrix to calc all partial sums
       
    for i in range(totalatoms):
        #calcs contribution of each atom to first moment tensor
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
    #sums all the contributions of the atoms into the intertial tensor (3x3 matrix)
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
    Iabc,dontcare=numpy.linalg.eigh(firstmomenttensor)
    
    #converts to MHz from wavenumbers 
    abc= 505379.006/Iabc
    return abc

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
    
    #shifts the coordinates to center of mass by value calculated above
    for i,row in enumerate(dmatrix):
        dmatrix[i,1]=dmatrix[i,1]-comx 
        dmatrix[i,2]=dmatrix[i,2]-comy
        dmatrix[i,3]=dmatrix[i,3]-comz
       
    return dmatrix

def isotopomers(expt_A,expt_B,expt_C,atoms_to_vary,a,b,c,mass,atom_list): # Uses a file of coordinates (Gaussian output format) + expt'l NS constants to generate A,B,C values of isotopomers

    H_vary_flag = 0
    C_vary_flag = 0
    N_vary_flag = 0
    O_vary_flag = 0
    Si_vary_flag = 0
    S_vary_flag = 0
    Cl_vary_flag = 0
    Br_vary_flag = 0

    raw3dmatrix=numpy.zeros((numpy.size(mass),4))#empty matrix for other calcs
    
    for i,row in enumerate(mass):
        #creates matrix for later calcs
        raw3dmatrix[i,0]=row
        raw3dmatrix[i,1]=a[i]
        raw3dmatrix[i,2]=b[i]
        raw3dmatrix[i,3]=c[i]

    originalmatrix=raw3dmatrix
    raw3dmatrix=distancefromcenter(raw3dmatrix)
    thr_abc = calcabc(raw3dmatrix) # Ab initio Rotational constant 

    factor = [expt_A/thr_abc[0], expt_B/thr_abc[1], expt_C/thr_abc[2]] # Calculates conversion factor exp/thr
    thr_isotop = []

    for entry in atoms_to_vary:
        if entry == 'H':
            H_vary_flag = 1
        if entry == 'C':
            C_vary_flag = 1
        if entry == 'N':
            N_vary_flag = 1
        if entry == 'O':
            O_vary_flag = 1
        if entry == 'Si':
            Si_vary_flag = 1
        if entry == 'S':
            S_vary_flag = 1
        if entry == 'Cl':
            Cl_vary_flag = 1
        if entry == 'Br':
            Br_vary_flag = 1

    H_counter = 1
    C_counter = 1
    N_counter = 1
    O_counter = 1
    Si_counter = 1
    S_counter = 1
    Cl_counter = 1
    Br_counter = 1

    for i in range(numpy.size(originalmatrix[:,0])):

        if raw3dmatrix[i,0]==1.007825037 and H_vary_flag ==1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=2.01410178
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            h2=calcabc(tempraw3dmatrix)
            temp3= []
            temp3.append('D-%s'%(str(H_counter)))
            temp3.append(h2[0])
            temp3.append(h2[1])
            temp3.append(h2[2])
            H_counter += 1
            thr_isotop.append(temp3)

        if raw3dmatrix[i,0]==12 and C_vary_flag ==1:
            tempraw3dmatrix = originalmatrix.copy() #copy the original matrix in every iteration and so we only have monosubstituted species
            tempraw3dmatrix[i,0]=13.003354838      #singly sustitutes a 12C atom
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            c13=calcabc(tempraw3dmatrix)
            temp= []
            temp.append('13C-%s'%(str(C_counter)))
            temp.append(c13[0])
            temp.append(c13[1])
            temp.append(c13[2])
            C_counter += 1
            thr_isotop.append(temp)  

        if raw3dmatrix[i,0]==14.003074008 and N_vary_flag == 1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=15.0001088982
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            n15=calcabc(tempraw3dmatrix)
            temp2= []
            temp2.append('15N-%s'%(str(N_counter)))
            temp2.append(n15[0])
            temp2.append(n15[1])
            temp2.append(n15[2])
            N_counter += 1
            thr_isotop.append(temp2)

        if raw3dmatrix[i,0]==15.99491464 and O_vary_flag == 1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=17.9991604
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            o18=calcabc(tempraw3dmatrix)
            temp1= []
            temp1.append('18O-%s'%(str(O_counter))) # Only does 18O.  17O is not hard to add in if desired...
            temp1.append(o18[0])
            temp1.append(o18[1])
            temp1.append(o18[2])
            O_counter += 1
            thr_isotop.append(temp1)      

        if raw3dmatrix[i,0]==27.9769265325 and Si_vary_flag == 1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=28.976494700
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            si29=calcabc(tempraw3dmatrix)
            temp3= []
            temp3.append('29Si-%s'%(str(Si_counter)))
            temp3.append(si29[0])
            temp3.append(si29[1])
            temp3.append(si29[2])
            thr_isotop.append(temp3)

            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=29.97377017
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            si30=calcabc(tempraw3dmatrix)
            temp4= []
            temp4.append('30Si-%s'%(str(Si_counter)))
            temp4.append(si30[0])
            temp4.append(si30[1])
            temp4.append(si30[2])
            thr_isotop.append(temp4)
            Si_counter += 1

        if raw3dmatrix[i,0]==31.97207100 and S_vary_flag == 1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=33.96786690
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            s34=calcabc(tempraw3dmatrix)
            temp4= []
            temp4.append('34S-%s'%(str(S_counter)))
            temp4.append(s34[0])
            temp4.append(s34[1])
            temp4.append(s34[2])
            thr_isotop.append(temp4)
            S_counter += 1

        if raw3dmatrix[i,0]==34.96885268 and Cl_vary_flag == 1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=36.96590259
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            cl37=calcabc(tempraw3dmatrix)
            temp3= []
            temp3.append('37Cl-%s'%(str(Cl_counter)))
            temp3.append(cl37[0])
            temp3.append(cl37[1])
            temp3.append(cl37[2])
            thr_isotop.append(temp3)
            Cl_counter += 1

        if raw3dmatrix[i,0]==78.9183371 and Br_vary_flag == 1:
            tempraw3dmatrix = originalmatrix.copy()
            tempraw3dmatrix[i,0]=80.9162906
            tempraw3dmatrix=distancefromcenter(tempraw3dmatrix)
            br81=calcabc(tempraw3dmatrix)
            temp3= []
            temp3.append('81Br-%s'%(str(Br_counter)))
            temp3.append(br81[0])
            temp3.append(br81[1])
            temp3.append(br81[2])
            thr_isotop.append(temp3)
            Br_counter += 1

    output_consts=[]
    for row in thr_isotop:
        output_consts.append((row[0],row[1]*factor[0],row[2]*factor[1],row[3]*factor[2]))

    f = open("Isotopomers.txt", "w")
    f.write("Thr Parent %s %s %s \n" %(thr_abc[0],thr_abc[1],thr_abc[2])) #write thr rotational constant to the file
    f.write("Exp Parent %s %s %s \n \n" %(expt_A,expt_B,expt_C)) #write Exp rotational constant to the file
    f.write("Scaled isotopologue constants \n")

    for row in thr_isotop:
        f.write("     %s   %f %f %f \n" %(row[0],row[1]*factor[0],row[2]*factor[1],row[3]*factor[2])) #for every row
    f.close()

    return output_consts


def fit_triples(list_a,list_b,list_c,trans_1,trans_2,trans_3,top_17,peaklist,file_num,A,B,C,DJ,DJK,DK,dJ,dK):
    
    all_combo_file = "all_combo_list%s.txt"%(str(file_num)) 
    
    all_combo_list_file  = open(all_combo_file,"a")
    
    temp_list_a = []
    temp_list_b = []
    temp_list_c = []

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


    pred_ratio = pow(10,max(float(trans_1[0]),float(trans_2[0]),float(trans_3[0]))-min(float(trans_1[0]),float(trans_2[0]),float(trans_3[0])))

    for freq_1,inten_1 in list_a:
      temp_diff = abs(float(trans_1[1])-float(freq_1))
      temp_list_a.append((freq_1,inten_1,temp_diff))

    for freq_2,inten_2 in list_b:
      temp_diff = abs(float(trans_2[1])-float(freq_2))
      temp_list_b.append((freq_2,inten_2,temp_diff))

    for freq_3,inten_3 in list_c:
      temp_diff = abs(float(trans_3[1])-float(freq_3))
      temp_list_c.append((freq_3,inten_3,temp_diff))

    unsorted_full_list = []

    for freq_1,inten_1,diff_1 in temp_list_a:
        for freq_2,inten_2,diff_2 in temp_list_b:
            for freq_3,inten_3,diff_3 in temp_list_c:
                real_ratio = max(float(inten_1),float(inten_2),float(inten_3))/min(float(inten_1),float(inten_2),float(inten_3))
                avg_diff = (diff_1+diff_2+diff_3)/3
                scaled_diff = avg_diff*(abs(real_ratio-pred_ratio)+1) # Freq. difference scaled by deviation of intensity ratio from predicted
                unsorted_full_list.append((freq_1,inten_1,freq_2,inten_2,freq_3,inten_3,scaled_diff))

    sorted_full_list = sorted(unsorted_full_list, key=lambda entry: entry[6])

    for freq_1,inten_1,freq_2,inten_2,freq_3,inten_3,total_diff in sorted_full_list:
        all_combo_list_file.write(str(freq_1)+","+str(inten_1)+","+str(freq_2)+","+str(inten_2)+","+str(freq_3)+","+str(inten_3)+", \n")

    #for freq_1,inten_1,diff_1 in sorted_list_a:#generates all combinations of three peaks from three peaklists
    #    for freq_2,inten_2,diff_2 in sorted_list_b:
    #        for freq_3,inten_3,diff_3 in sorted_list_c:
    #            all_combo_list_file.write(str(freq_1)+","+str(inten_1)+","+str(freq_2)+","+str(inten_2)+","+str(freq_3)+","+str(inten_3)+", \n")
    
    all_combo_list_file.close()
    
    all_combo_list_file = open(all_combo_file)
    
    triples_counter = 0
    output_file = ""
    regular_counter = 0
    #error_counter = 0

    theor_inten_list = []
    for x in range(len(top_17)):
        temp_inten = 10**float(top_17[x][0])
        theor_inten_list.append(temp_inten)

    theor_inten_ratio = (max(theor_inten_list)/min(theor_inten_list))
    theor_inten_avg = (sum(theor_inten_list)/len(theor_inten_list))
    theor_inten_unitless_stdev = numpy.std(theor_inten_list)/theor_inten_avg


    for all_combo_line in all_combo_list_file:
        all_combo_line = all_combo_line.split(",")
        peaks_triple= [(all_combo_line[0],all_combo_line[1]),(all_combo_line[2],all_combo_line[3]),(all_combo_line[4],all_combo_line[5])]
        

        input_file = ""
        input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
        input_file += "   8  500   5    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n" # don't choose more than 497 check transitions or it will crash.
        input_file +="a   1  1  0  50  0  1  1  1  1  -1   0\n"
        input_file += "           10000  %s 1.0E+004 \n" % A
        input_file += "           20000  %s 1.0E+004 \n" % B
        input_file += "           30000  %s 1.0E+004 \n" % C
        input_file += "             200  %s 1.0E-025 \n" % DJ
        input_file += "            1100  %s 1.0E-025 \n" % DJK
        input_file += "            2000  %s 1.0E-025 \n" % DK
        input_file += "           40100  %s 1.0E-025 \n" % dJ
        input_file += "           41000  %s 1.0E-025 \n" % dK
        fh_par = open("default%s.par"%(str(file_num)),'w')
        fh_par.write(input_file)
        fh_par.close()

        input_file = ""#the next part adds in the three peaks to be fit
        input_file += trans_1[2][0:2]+' '+trans_1[2][2:4]+' '+trans_1[2][4:6]+' '+\
                      trans_1[3][0:2]+' '+trans_1[3][2:4]+' '+trans_1[3][4:6]+'                      '+peaks_triple[0][0]+' 0.50 1.0000\n' 
        input_file += trans_2[2][0:2]+' '+trans_2[2][2:4]+' '+trans_2[2][4:6]+' '+\
                      trans_2[3][0:2]+' '+trans_2[3][2:4]+' '+trans_2[3][4:6]+'                      '+peaks_triple[1][0]+' 0.50 1.0000\n' 
        input_file += trans_3[2][0:2]+' '+trans_3[2][2:4]+' '+trans_3[2][4:6]+' '+\
                      trans_3[3][0:2]+' '+trans_3[3][2:4]+' '+trans_3[3][4:6]+'                      '+peaks_triple[2][0]+' 0.50 1.0000\n'
        counter = 0
        for line in top_17:#the hack that adds in the check transitions but doesn't use them in the fit
            input_file += line[2][0:2]+' '+line[2][2:4]+' '+line[2][4:6]+' '+\
                      line[3][0:2]+' '+line[3][2:4]+' '+line[3][4:6]+'                      '+'%s.0'%(str(counter))+' 0.00 1.0000\n'
            counter += 1
        fh_lin = open("default%s.lin"%(str(file_num)), "w")
        fh_lin.write(input_file)
        fh_lin.close()        
        a = subprocess.Popen("SPFIT%s default%s"%(str(file_num),str(file_num)), stdout=subprocess.PIPE, shell=False)
        a.stdout.read()#used to let SPFIT finish

        const_list = []

        fh_var = open("default%s.var"%(str(file_num)))
        for line in fh_var:
            if line[8:13] == "10000":
                temp_A = float(line[15:37])
                const_list.append("%.3f" %temp_A)
            if line[8:13] == "20000":
                temp_B = float(line[15:37])
                const_list.append("%.3f" %temp_B)
            if line[8:13] == "30000":
                temp_C = float(line[15:37])
                const_list.append("%.3f" %temp_C)
        fh_var.close()

        fh_fit = open("default%s.fit"%(str(file_num)))
        file_list = []
        for line in fh_fit:
                file_list.append(line)
        fh_fit.close()

        freq_list = []
        for x in range(len(file_list)):
            if file_list[-x][11:14] == "RMS":
                rms_fit = float(file_list[-x][22:32]) #note - assumes RMS fit error is less than 1 GHz.  Change 22 to 21 if this is a problem.
            if file_list[-x][5:6] == ":" and int(file_list[-x][3:5])>3:
                freq_list.append(file_list[-x][60:71])
            if file_list[-x][40:64]=="EXP.FREQ.  -  CALC.FREQ.":
                break
        read_fit = (const_list[0],const_list[1], const_list[2],freq_list)
        triples_counter +=1
        constants = read_fit[0:3]
        freq_17 = read_fit[3]
        freq_17.reverse()
        A_1 = float(constants[0])
        B_1 = float(constants[1])
        C_1 = float(constants[2])
        omc_list = []
        
        for x in range(len(top_17)): #matches peaks in the top 17 to peaks in experimental peak list <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            real_omc = 1.0
            current_peak = float(freq_17[x])
            if current_peak>p_1 and current_peak<p_2:#conditionals to find proper portion of experimental peak list to loop through
                peaks_section = peak_list_1
                peak = p_2
                regular_counter +=1
            elif current_peak>p_3 and current_peak<p_4:
                peaks_section = peak_list_2
                peak = p_4
                regular_counter +=1
            elif current_peak>p_5 and current_peak<p_6:
                peaks_section = peak_list_3
                peak = p_6
                regular_counter +=1
            elif current_peak>p_7 and current_peak<p_8:
                peaks_section = peak_list_4
                peak = p_8
                regular_counter +=1
            elif current_peak>p_8 or current_peak<p_1:
                peaks_section = peaklist
                peak = p_8
                #error_counter +=1
                real_omc = 0.0# this is the omc if you throw out peaks that go over the edge of the spectrum
            else:
                peaks_section = peaklist
                peak = p_8
                regular_counter +=1
            old_omc = 100000.0
            for peak_freq,peak_inten in peaks_section: #find nearest peak in actual spectrum to the given top 20 peak
                omc = abs(current_peak-float(peak_freq))
                omc_low = abs(current_peak-float(peak))
                if omc>old_omc:
                    omc_low = old_omc
                    omc_inten = temp_inten                    
                    break
                old_omc = omc
                temp_inten = peak_inten
            if real_omc == 1.0:
                real_omc = omc_low
            omc_list.append((omc_low, real_omc, omc_inten))
        omc_avg = [float(omc) for omc, real_omc, omc_inten in omc_list]
        real_omc_avg = [float(real_omc) for omc, real_omc, omc_inten in omc_list]
        omc_inten_scoring = [float(omc_inten) for omc, real_omc, omc_inten in omc_list]
        score = str(len([omc for omc in omc_avg if omc<2.0])) #scores the accuracy of the fit, currently based on a peak being within 2 MHz which may be too coarse
        avg = (sum(omc_avg)/len(omc_avg))+rms_fit
        real_avg = (sum(real_omc_avg)/len(real_omc_avg))+rms_fit

        omc_inten_ratio = (max(omc_inten_scoring)/min(omc_inten_scoring))
        omc_inten_avg = (sum(omc_inten_scoring)/len(omc_inten_scoring))
        omc_inten_unitless_stdev = numpy.std(omc_inten_scoring)/omc_inten_avg

        score_inten_penalty = 1

        if (omc_inten_ratio <= 0.5*theor_inten_ratio) or (omc_inten_ratio >= 1.5*theor_inten_ratio):
            score_inten_penalty += 1

        if (omc_inten_unitless_stdev <= 0.5*theor_inten_unitless_stdev) or (omc_inten_unitless_stdev >= 1.5*theor_inten_unitless_stdev):
            score_inten_penalty += 1
        
        penalized_avg = avg*score_inten_penalty

        if float(A_1)>=float(B_1) and float(B_1)>=float(C_1) and float(C_1)>0:
            if int(score)<10: #makes sorting work properly later
                score = '0'+score  
            output_file += 'score = '+' '+score+' '+"Const = "+str(A_1)+' '+str(B_1)+' '+str(C_1)+' '+"average omc = "+str(avg)+'  '+"avg w/out peaks over edge = "+str(real_avg)+' '+"avg w/ inten penalty = "+str(penalized_avg)+"\n"

            if real_avg <= 0.2: #appends good finds (RMS < 0.2 MHz, ignoring peaks over edge) to interim file for each processor
                interim_output = 'score = '+' '+score+' '+"Const = "+str(A_1)+' '+str(B_1)+' '+str(C_1)+' '+"average omc = "+str(avg)+'  '+"avg w/out peaks over edge = "+str(real_avg)+' '+"avg w/ inten penalty = "+str(penalized_avg)+"\n"
                fh_interim_good = open("interim_good_output%s.txt"%(str(file_num)), "a")
                fh_interim_good.write(interim_output)
                fh_interim_good.close()

            if triples_counter == 100000: #appends to file after every 100000 triples
                fh_final = open("final_output%s.txt"%(str(file_num)), "a")
                fh_final.write(output_file)
                fh_final.close()
                triples_counter = 0
                output_file = ""
    fh_final = open("final_output%s.txt"%(str(file_num)), "a")#writes separate file for each processor
    #print 'out of %s peaks there were %s peaks that werent in the experimental spectrum'%(regular_counter, error_counter) 
    fh_final.write(output_file)
    fh_final.close()
    os.system("sort -r 'final_output%s.txt'>sorted_final_out%s.txt"%(str(file_num),str(file_num)))#sorts output by score
    
def isotopologue_fit(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,processors,inten_high,inten_low,temperature,Jmax,peaklist,freq_low,freq_high,trans_1,trans_2,trans_3,filter_level,atoms_to_vary,a,b,c,mass,atom_list,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,check_peaks_list):

    main_flow = 'Isotopologues'

    a = subprocess.Popen("mkdir %s"%job_name) # Need to do something here with try/except because this may or may not work depending on if we're doing a new fit from scratch or if we're continuing from a previous NS fitting run.
    a.wait()

    freq_uncertainty = 0.0

    for number in range(processors):
        y = subprocess.Popen("cp SPFIT.EXE %s\SPFIT%s.EXE"%(job_name,number), stdout=subprocess.PIPE, shell=True) # Ditto here with try/except.
        y.stdout.read()     

    y = subprocess.Popen("cp SPCAT.EXE %s\SPCAT.EXE"%(job_name), stdout=subprocess.PIPE, shell=True) # Ditto here with try/except.
    y.stdout.read() 
    
    os.chdir(job_name)
    
    int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
    var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()

    full_list = cat_reader(freq_high,freq_low,flag="default")	# New code starts here
    full_list = intensity_filter(full_list,peaklist,inten_low,filter_level)

    trans_1_uncert = float(trans_1[4])
    trans_2_uncert = float(trans_2[4])
    trans_3_uncert = float(trans_3[4])

    isotopomer_count = 0

    calc_isotopomer_ABC = isotopomers(float(A),float(B),float(C),atoms_to_vary,a,b,c,mass,atom_list)
    
    for isotopologue in calc_isotopomer_ABC:
        isotope_ID = isotopologue[0]
        curr_A = isotopologue[1]
        curr_B = isotopologue[2]
        curr_C = isotopologue[3]

        a = subprocess.Popen("mkdir %s"%isotope_ID)
        a.wait()

        for number in range(processors):
            y = subprocess.Popen("cp SPFIT%s.EXE %s\SPFIT%s.EXE"%(number,isotope_ID,number), stdout=subprocess.PIPE, shell=True)
            y.stdout.read()     

        y = subprocess.Popen("cp SPCAT.EXE %s\SPCAT.EXE"%(isotope_ID), stdout=subprocess.PIPE, shell=True)
        y.stdout.read() 

        os.chdir(isotope_ID)

        int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
        var_writer(curr_A,curr_B,curr_C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        run_SPCAT()

        # For now leave in the next bit of code, but we'll ultimately want to move the selection of fitting and check transitions to the GUI and only just do triples fits here.

        # This updates the center frequencies of the fitting transitions to account for the new rotational constants.
        # If this isn't done, then there are issues when trans_1,trans_2,trans_3 are passed over to fit_triples.
        trans_1_center,trans_2_center,trans_3_center = trans_freq_reader(trans_1,trans_2,trans_3)

        trans_1 = (trans_1[0],trans_1_center,trans_1[2],trans_1[3],trans_1[4])
        trans_2 = (trans_2[0],trans_2_center,trans_2[2],trans_2[3],trans_2[4])
        trans_3 = (trans_3[0],trans_3_center,trans_3[2],trans_3[3],trans_3[4])

        trans_1_center = float(trans_1_center)
        trans_2_center = float(trans_2_center)
        trans_3_center = float(trans_3_center)

        if (trans_1_center == 0) or (trans_2_center == 0) or (trans_3_center == 0):
            bad_triples = buttonbox(msg='One of your fitting transitions is no longer within bounds after isotopologue scaling.  Please manually choose new transitions or quit.', choices=('Choose new transitions','Quit'))
            if bad_triples == 'Quit':
                return
            elif bad_triples == 'Choose new transitions':
                (highest_uncert,trans_1,trans_2,trans_3) = triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag="force manual")
                trans_1_uncert = float(trans_1[4])
                trans_2_uncert = float(trans_2[4])
                trans_3_uncert = float(trans_3[4])

        (trans_1,trans_2,trans_3,trans_1_peaks,trans_2_peaks,trans_3_peaks,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,num_of_triples,decision) = triples_gen(trans_1_uncert,trans_2_uncert,trans_3_uncert,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,isotopomer_count,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3)

        top_peaks = check_peaks_list
    
        top_peaks_3cut = []
        for entry in top_peaks:
            if (entry[2] == trans_1[2] and entry[3] == trans_1[3]) or (entry[2] == trans_2[2] and entry[3] == trans_2[3]) or (entry[2] == trans_3[2] and entry[3] == trans_3[3]):
                pass
            else:
                top_peaks_3cut.append(entry)

        job_file = ""
        str(top_peaks_3cut)
        str(trans_1)
        str(trans_2)
        str(trans_3)
        fitting_peaks_str = ""
        for entry in top_peaks_3cut:
            fitting_peaks_str+=str(entry)+"\n"
        
        job_file += "Job Name %s \n Isotope ID %s \n u_A: %s \n u_B: %s \n u_C: %s \n A: %s \n B: %s \n \
C: %s \n DJ: %s \n DJK: %s \n DK: %s \n dJ: %s \n dK: %s \n processors: %s \n freq_high: %s \n freq_low: %s \n \
inten_high: %s \n inten_low: %s \n Temp: %s \n Jmax: %s \n freq_uncertainty: %s \n number of triples: %s \n Check peaks:\n%s \n trans_1: %s \n trans_2: %s \n trans_3: %s "%(job_name,isotope_ID,u_A,u_B,u_C,curr_A,curr_B,curr_C,DJ,DJK,DK,dJ,dK,str(processors),str(freq_high),\
            str(freq_low),str(inten_high),str(inten_low),str(temperature),str(Jmax),str(freq_uncertainty),str(num_of_triples),fitting_peaks_str,str(trans_1),str(trans_2),str(trans_3))
        Job_fh = open("input_data_%s_%s.txt"%(job_name,isotope_ID),"w")
        Job_fh.write(job_file) 
        Job_fh.close()
    
        flag123 = 0
        flag132 = 0
        flag213 = 0
        flag231 = 0
        flag312 = 0
        flag321 = 0

        if (len(trans_1_peaks) >= len(trans_2_peaks)) and (len(trans_1_peaks) >= len(trans_3_peaks)):
            list_a_peaks = trans_1_peaks
            if len(trans_2_peaks) >= len(trans_3_peaks):
                list_b_peaks = trans_2_peaks
                list_c_peaks = trans_3_peaks
                flag123 = 1
            elif len(trans_3_peaks) >= len(trans_2_peaks):
                list_b_peaks = trans_3_peaks
                list_c_peaks = trans_2_peaks
                flag132 = 1
        elif (len(trans_2_peaks) >= len(trans_1_peaks)) and (len(trans_2_peaks) >= len(trans_3_peaks)):
            list_a_peaks = trans_2_peaks
            if len(trans_1_peaks) >= len(trans_3_peaks):
                list_b_peaks = trans_1_peaks
                list_c_peaks = trans_3_peaks
                flag213 = 1
            elif len(trans_3_peaks) >= len(trans_1_peaks):
                list_b_peaks = trans_3_peaks
                list_c_peaks = trans_1_peaks
                flag231 = 1
        elif (len(trans_3_peaks) >= len(trans_1_peaks)) and (len(trans_3_peaks) >= len(trans_2_peaks)):
            list_a_peaks = trans_3_peaks
            if len(trans_1_peaks) >= len(trans_2_peaks):
                list_b_peaks = trans_1_peaks
                list_c_peaks = trans_2_peaks
                flag312 = 1
            elif len(trans_2_peaks) >= len(trans_1_peaks):
                list_b_peaks = trans_2_peaks
                list_c_peaks = trans_1_peaks
                flag321 = 1

        random.shuffle(list_c_peaks)  # Shuffle so that each processor gets a range of values for the third peak, not processor 0 getting only the lowest frequencies.  

        list_c_list = []
        for num in range(processors):
            processors = float(processors)
            num = float(num)
            x = int((num)*(len(list_c_peaks)/processors))
            y = int(len(list_c_peaks)*((num+1)/processors))
            list_c_list.append(list_c_peaks[x:y])
        list_c_list.append("marker")

        if flag123 == 1:
            trans_1_peaks = list_a_peaks
            trans_2_peaks = list_b_peaks
            trans_3_peaks = list_c_list
        if flag213 == 1:
            trans_1_peaks = list_b_peaks
            trans_2_peaks = list_a_peaks
            trans_3_peaks = list_c_list
        if flag312 == 1:
            trans_1_peaks = list_b_peaks
            trans_2_peaks = list_c_list
            trans_3_peaks = list_a_peaks
        if flag132 == 1:
            trans_1_peaks = list_a_peaks
            trans_2_peaks = list_c_list
            trans_3_peaks = list_b_peaks
        if flag231 == 1:
            trans_1_peaks = list_c_list
            trans_2_peaks = list_a_peaks
            trans_3_peaks = list_b_peaks
        if flag321 == 1:
            trans_1_peaks = list_c_list
            trans_2_peaks = list_b_peaks
            trans_3_peaks = list_a_peaks

        processors = int(processors)
        for num in range(processors):
            if trans_1_peaks[-1]=="marker":
                trans_x_peaks = trans_1_peaks[num]
                trans_y_peaks = trans_2_peaks
                trans_z_peaks = trans_3_peaks
        
            if trans_2_peaks[-1]=="marker":
                trans_x_peaks = trans_1_peaks
                trans_y_peaks = trans_2_peaks[num]
                trans_z_peaks = trans_3_peaks

            if trans_3_peaks[-1]=="marker":
                trans_x_peaks = trans_1_peaks
                trans_y_peaks = trans_2_peaks
                trans_z_peaks = trans_3_peaks[num]
        
            vars()["p%s"%str(num)] = Process(target=fit_triples, args=(trans_x_peaks,trans_y_peaks,trans_z_peaks,trans_1,trans_2,trans_3,top_peaks_3cut,peaklist,num,A,B,C,DJ,DJK,DK,dJ,dK))

        for num in range(processors):
            vars()["p%s"%str(num)].start()
        for num in range(processors):
            vars()["p%s"%str(num)].join()
        
        a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 4 -n > sorted_omc_cat_%s.txt'%(isotope_ID), shell=True)
        a.wait()

        a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 6 -n > sorted_inten_omc_cat_%s.txt'%(isotope_ID), shell=True)
        a.wait()

        f = open('sorted_inten_omc_cat_%s.txt'%(isotope_ID),'r')
        fits = []

        for i in range(100):
            r = f.readline()

            if i < 10:
                number = '0'+str(i)
            else:
                number = i

            temp = str(number) + ' ' + r
            fits.append(temp)

        f.close()

        OMC_char_buffer = ""
        for i in (fits):
            OMC_char_buffer += i+''

        f100=open('best100_%s.txt'%(isotope_ID),'w')
        f100.write(OMC_char_buffer)
        f100.close()

        isotopomer_count += 1
        os.chdir(os.pardir)

