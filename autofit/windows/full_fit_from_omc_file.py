from easygui import *
import numpy
from collections import OrderedDict
import subprocess

def run_SPCAT_refit(): 
    a = subprocess.Popen("SPCAT refit", stdout=subprocess.PIPE, shell=False)
    a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen

def int_writer_refit(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="25.8", temp="298"):#generates SPCAT input file
    input_file = ""
    #print "freq_max=",freq
    input_file += "Molecule \n"
    input_file += "0  91  %s  %s  %s  %s  %s %s  %s\n"%(Q_rot, J_min, J_max,inten,inten,freq, temp)
    input_file += " 001  %s \n" % u_A
    input_file += " 002  %s \n" % u_B
    input_file += " 003  %s \n" % u_C
    fh_int = open("refit.int", "w")
    fh_int.write(input_file)
    fh_int.close()

def var_writer_refit(A,B,C,DJ,DJK,DK,dJ,dK):#generates SPCAT input file
    input_file = ""
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   8  430   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
    input_file +="a   1  1  0  99  0  1  1  1  1  -1   0\n"
    input_file += "           10000  %s 1.0E+004 \n" %(A)
    input_file += "           20000  %s 1.0E+004 \n" %(B)
    input_file += "           30000  %s 1.0E+004 \n" %(C)
    input_file += "             200  %s 1.0E-025 \n" %(DJ)
    input_file += "            1100  %s 1.0E-025 \n" %(DJK) #need to actually check numbers: SPFIT doesn't read -- as a positive!
    input_file += "            2000  %s 1.0E-025 \n" %(DK)
    input_file += "           40100  %s 1.0E-025 \n" %(dJ)
    input_file += "           41000  %s 1.0E-025 \n" %(dK)
    fh_var = open("refit.var",'w')
    fh_var.write(input_file)
    fh_var.close()

def par_writer_refit(A,B,C,DJ,DJK,DK,dJ,dK): #generates SPFIT par file
    dA = str(0.5*float(A))  #These allow A, B, and C to vary by 50% and the distortions to vary by a factor of 10.  Too much, too little?
    dB = str(0.5*float(B))
    dC = str(0.5*float(C))
    dDJ = str(abs(10.0*float(DJ)))
    dDJK = str(abs(10.0*float(DJK)))
    dDK = str(abs(10.0*float(DK)))
    ddJ = str(abs(10.0*float(dJ)))
    ddK = str(abs(10.0*float(dK)))
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
    fit_flag = 0.5

    for line in assignment_list:
    	input_file += line[1][0:2]+' '+line[1][2:4]+' '+line[1][4:6]+' '+\
                      line[2][0:2]+' '+line[2][2:4]+' '+line[2][4:6]+'                      '+str(line[0])+' '+line[3]+' 1.0000\n'

    fh_lin = open("refit.lin",'w')
    fh_lin.write(input_file)
    fh_lin.close()        


def cat_refit_reader(freq_high,freq_low): #reads output from SPCAT
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
    pred_peaks = cat_refit_reader(1000000, 0)
    pred_trans_for_refit = []

    for line in peaklist:
        for peak in pred_peaks:
            if line[2] == peak[2] and line[3] == peak[3]:
                pred_trans_for_refit.append((peak[1],peak[2],peak[3]))

    return pred_trans_for_refit

def peakpicker(spectrum,thresh_l,thresh_h):#Code taken from Cristobal's peak-picking script; assumes spectrum is in increasing frequency order
    peaks=[]
    for i in range(1, len(spectrum)-1):
        if spectrum[i,1] > thresh_l and spectrum[i,1] < thresh_h and spectrum[i,1] > spectrum[(i-1),1] and spectrum[i,1] > spectrum[(i+1),1]:
            peaks.append(spectrum[i])

    peakpicks=numpy.zeros((len(peaks),2))
    for i,row in enumerate(peaks):
        peakpicks[i,0]=row[0]
        peakpicks[i,1]=row[1]
    freq_low = spectrum[0,0]
    freq_high = spectrum[-1,0]
    return peakpicks, freq_low, freq_high


def match_to_peaklist(pred_trans,peaklist):

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

        if old_omc > 0.5: # If the best match is too far off, we don't want to have a bad line in the fit file.  0.5 MHz for now, but this should probably be a user-determined parameter.
        	weight = '0.00'

        best_match_freqs.append((temp_freq,pred_trans[x][1],pred_trans[x][2],weight))

    return best_match_freqs



f=open(fileopenbox(msg="Choose a results file (sorted_omc_cat.txt) to refine"))

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

msg = "Choose a previous result for full fitting (with distortions)."
title = "Microwave Fitting Program"
choice = choicebox(msg,title,fits)
A_fit = choice.split()[6]
B_fit = choice.split()[7]
C_fit = choice.split()[8]

f2 = open(fileopenbox(msg="Choose the job input file (input_data_X.txt) from the run corresponding to the results file"))

refitting_peaks_flag = 0
refit_peaks_list=[]
refit_check_peaks_list=[]

for line in f2:
            
    if line.split() != []:
                                
        if line.split()[0] == "u_A:":
            u_A = line.split()[1] 
        if line.split()[0] == "u_B:":
            u_B = line.split()[1] 
        if line.split()[0] == "u_C:":
            u_C = line.split()[1]                                 
        if line.split()[0] == "DJ:":
            DJ_init = line.split()[1]
        if line.split()[0] == "DK:":
            DK_init = line.split()[1]
        if line.split()[0] == "DJK:":
            DJK_init = line.split()[1] 
        if line.split()[0] == "dJ:":
            dJ_init = line.split()[1]
        if line.split()[0] == "dK:":
            dK_init = line.split()[1]
        if line.split()[0] == "freq_high:":
            freq_high = float(line.split()[1]) 
        if line.split()[0] == "inten_high:":
            inten_high = float(line.split()[1]) 
        if line.split()[0] == "inten_low:":
            inten_low = float(line.split()[1]) 
        if line.split()[0] == "Temp:":
            temp = float(line.split()[1])
        if line.split()[0] == "Jmax:":
            Jmax = float(line.split()[1])
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
            #print tuples
            refit_check_peaks_list.append(tuples)
        if line.split()[0] == "Check":
            refitting_peaks_flag = 1

f2.close()

int_writer_refit(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temp=temp)
var_writer_refit(A_fit,B_fit,C_fit,DJ_init,DJK_init,DK_init,dJ_init,dK_init)
run_SPCAT_refit()

fit_peaklist = refit_peaks_list + refit_check_peaks_list
fit_peaklist = list(OrderedDict.fromkeys(fit_peaklist)) # This removes duplicates if some transitions are both check and fit transitions, which might bias SPFIT towards hitting those in particular.

updated_trans = trans_freq_refit_reader(fit_peaklist) # Finds updated predicted frequencies with improved A, B, and C estimates.

fh = numpy.loadtxt(fileopenbox(msg="Choose the experimental spectrum file in two column format: frequency intensity")) #loads full experimental data file, not just list of peaks
(peaklist, freq_low, freq_high) = peakpicker(fh,inten_low,inten_high) # Calls slightly modified version of Cristobal's routine to pick peaks instead of forcing user to do so.

best_matches = match_to_peaklist(updated_trans,peaklist) # Assigns closest experimental peak frequencies to transitions


par_writer_refit(A_fit,B_fit,C_fit,DJ_init,DJK_init,DK_init,dJ_init,dK_init)
lin_writer_refit(best_matches)

a = subprocess.Popen("spfit refit", stdout=subprocess.PIPE, shell=False)
a.stdout.read()#used to let SPFIT finish

