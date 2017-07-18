import sys
from os.path import dirname, join, basename # This is from Pyinstaller issue #1957
try:
    tcl_lib = join(sys._MEIPASS, "lib")
    tcl_new_lib = join(dirname(dirname(tcl_lib)), basename(tcl_lib))
    import shutil
    shutil.copytree(src=tcl_lib, dst=tcl_new_lib)
except:
    pass

import subprocess
import os
from easygui import *
import multiprocessing
import multiprocessing.forking
import timeit
#from multiprocessing import Process
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

version 15c:

-Fixed an issue with reading var files that depended on which version of SPCAT was in use.  Also added a line to fix an 
unusual edge case with omc_inten that would occasionally lead to a crash for particularly poor combinations of predicted
and actual peak positions.

version 15b:

-Added a time estimate to the length of a job based on the number of triples to fit, the estimated time it takes to run
SPCAT 25 times, and the number of cores detected on the user's machine.  Time estimate assumes that the user will run the
job on all available cores.

-Dialogue boxes now give quantum numbers of transitions as well.

-Also fixed a bug in the final_uncerts function and shaved a few seconds off of the cubic_spline routine.

version 15:

-Added an intensity penalty to the scoring of potential fits.  Good fits need to have maximum intensity ratios and unitless
standard deviations of the intensity over all of the check peaks within 50 percent of theoretical values to avoid being
penalized.  Output stored the previous way is still stored in sorted-omc-cat.txt, output sorted the new way is in
sorted-inten-omc-cat.txt.  The "best 100" files are taken from the intensity sorted output.

version-14b:

-Added in manual selection of fitting peaks (rather than only choosing from the automatically-generated list).  Testing of this
revealed underlying bugs in isotopic substitution when transitions could leave the bounds after isotope scaling, requiring new
code to be written to allow the fitting transitions to be updated at various points in the overall algorithm.  The current
situation is still a bit ugly (functions requiring many, many input parameters) - this should be cleaned up in the future.

-Added in a "boundary penalty" to the score of triples in the automatically generated list based on how many are within 100 MHz
of a spectral boundary.

-Provided a user-dependent scale factor to filter acceptable transitions when searching for isotopolgoues.  The factor is used
to require that NS transitions be a certain minimum height (relative to the lower intensity threshold) in the experimental 
spectrum before they can be used as fitting or scoring transitions in the isotopologues.  If the user does not believe that
transitions from the NS isotopologue will be in their spectrum (e.g. pure compound of a synthetically produced isotopically-substituted 
species), they should set this factor to 0.

version-14:

-Merges v13 and v13-isotopes together into a single program.  Some code-tidying has been done to partially de-duplicate things.

-There is now more flexibility on coordinate input for isotopologue searches.  It can handle XYZ file format or either of the
primary formats that come from Gaussian output (5 column or 6 column).  It also shouldn't care if atom_labels are in atomic
number or in chemical symbol ("6" or "C").  Format detection is based on the number of entries in the first row.  If the first 
row contains 1 entry, it assumes XYZ format (which means line 2 is a comment and line 3 is the first line containing information 
about atom type and position).  If the first row contains either 5 or 6 entries, it assumes it is the appopriate format from 
Gaussian output.

version 13-isotopes:

Uses version 13 code as base, with isotopologue constant generator slapped on top.  Later will merge into v14, most likely.
Will run separate triples fitter runs on each isotopologue selected by user.  Right now behavior is that user selects one or
more atom types (list automatically determined from an ab initio structure input), and the triples fitter will examine all
singly-substituted species of those types.  Each fit is saved in a separate sub-directory.

Format of structure file is from Gaussian output (below for EEO):

      1          6           0       -2.670317    0.571148    0.159988
      2          1           0       -2.434542    1.571170   -0.210325
      3          1           0       -2.719870    0.603794    1.250293
      4          1           0       -3.650729    0.279198   -0.226138
      5          6           0       -1.617736   -0.419313   -0.289709
      6          1           0       -1.854133   -1.429010    0.076941
      7          1           0       -1.565139   -0.462200   -1.387392
      8          8           0       -0.364732   -0.006000    0.230899
      9          6           0        0.682634   -0.891655   -0.124503
     10          1           0        0.539846   -1.873119    0.351947
     11          1           0        0.711070   -1.021547   -1.215745
     12          6           0        1.972157   -0.256104    0.344429
     13          1           0        2.816465   -0.900076    0.088216
     14          1           0        1.939463   -0.136824    1.434893
     15          8           0        2.193845    0.987588   -0.298292
     16          1           0        1.384241    1.491448   -0.164781


version 13 features:

-after the fitting has completed, the user can select individual results to further refine by fitting
distortion constants and adding additional transitions.  Any number of the top 100 fits can be refined 
in this fashion.  The fit files thus generated are copied and are available for use with other programs 
(such as AABS).

-User-supplied spectrum is interpolated to 2 kHz resolution with a cubic spline prior to peakpicking.
This currently happens automatically, regardless of actual resolution of user-file.  Should probably check
to make sure that we're not actually downgrading the resolution of the data.

-Some code & comment tidying / compaction.

version 12 features:

-peaks are found from spectral data; user no longer supplies a list of picked peaks first.  High and low 
frequency ranges are automatically determined from the spectrum file.

-only physically reasonable results are saved to file (A >= B >= C, all rotational constants are positive)

-More flexibility with choosing search windows; fixed width for each transition, different for different
transitions, or based on 3x the estimated SPCAT uncertainty.  (The SPCAT uncertainties assume errors of 
2 percent on A, 1 percent on B and C, 10 percent on DJ, DJK, and DK, and 30 percent on dJ, dK.  These may 
need to be revisited.)  Also, this should later be expanded to allow for the selection of isotopic windows
based on ab initio scaling.

-Program alerts user if search windows will exceed bounds of the spectrum and prompts them for new uncertainties
or lets them accept the consequences of their choices.  The user can also elect to quit at this point if they realize
they have a transition far too close to one of the bounds of the spectrum.

version 11 features:

-sorting of triples list to evaluate most promising candidates first

-interim output of good fit results to unsorted file


version 10 features:

-minor bugfixes in displays of numbers in choice boxes and linear-dependence code

-automatic score calculation and ranking of possible triples fitting combinations


version 9 features:

-memory handling


-can use any number of processors, I left the code for my old method of multiprocessing commented out at the bottom
just in case there are bugs with the new version

-input files
if you leave out the "check peaks" or trans_x area it will prompt you for transitions...

example input file:


Job Name conf_I_#2 
 u_A: 0.0 
 u_B: 1.0 
 u_C: 0.0 
 A: 0.98187789E+04 
 B: 0.87351235E+03 
 C: 0.82235840E+03 
 DJ: -0.46322711E-04 
 DJK: 0.80645742E-03 
 DK: -0.23482420E-01 
 dJ: -0.49333549E-05 
 dK: 0.11644082E-03 
 processors: 8 
 freq_high: 18000.0 
 freq_low: 8000.0 
 inten_high: 100000.0 
 inten_low: 0.0 
 Temp: 2.0 
 Jmax: 30.0 
 freq_uncertainty: 900.0 
 number of triples: 7106688 
 Check peaks:
('-6.1163', ' 9360.6137', ' 5 1 4', ' 5 0 5')
('-6.1236', ' 9229.3112', ' 4 1 3', ' 4 0 4')
('-6.1397', ' 9519.9682', ' 6 1 5', ' 6 0 6')
('-6.1689', ' 9125.2491', ' 3 1 2', ' 3 0 3')
('-6.1896', ' 9708.3423', ' 7 1 6', ' 7 0 7')
('-6.2029', '12285.8335', ' 2 1 2', ' 1 0 1')
('-6.2633', ' 9926.8555', ' 8 1 7', ' 8 0 8')
('-6.2672', ' 9047.7751', ' 2 1 1', ' 2 0 2')
 
 trans_1: ('-5.8815', '15499.0588', ' 4 1 4', ' 3 0 3') 
 trans_2: ('-6.0131', '13905.0568', ' 3 1 3', ' 2 0 2') 
 trans_3: ('-6.7334', '15777.1423', ' 6 2 4', ' 7 1 7')





"""
class _Popen(multiprocessing.forking.Popen): #from StackOverflow 24944558 pyinstaller multiprocessing question
    def __init__(self, *args, **kw):
    	if hasattr(sys, 'frozen'):
            # We have to set original _MEIPASS2 value from sys._MEIPASS
            # to get --onefile mode working.
            os.putenv('_MEIPASS2', sys._MEIPASS)
        try:
            super(_Popen, self).__init__(*args, **kw)
        finally:
            if hasattr(sys, 'frozen'):
                # On some platforms (e.g. AIX) 'os.unsetenv()' is not
                # available. In those cases we cannot delete the variable
                # but only set it to the empty string. The bootloader
                # can handle this case.
                if hasattr(os, 'unsetenv'):
                    os.unsetenv('_MEIPASS2')
                else:
                    os.putenv('_MEIPASS2', '')

class Process(multiprocessing.Process):
    _Popen = _Popen

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

def time_estimate():
    s = """\
    int_writer(1,1,1, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temperature=2, flag = "default")
    var_writer(3000,2000,1000,2E-4,2E-4,2E-4,2E-4,2E-4,'Normal species',flag="uncert")
    run_SPCAT()
    """

    outtime = timeit.timeit(stmt=s, number=25, setup="from __main__ import int_writer,var_writer,run_SPCAT")
    scale_factor_one_proc = outtime / 25.0 # Rought estimated time in seconds for a single triple on one core.
    cores_detected = multiprocessing.cpu_count()
    scale_factor = scale_factor_one_proc / cores_detected # Rough estimated time in seconds for a single triple, divided by number of processors.

    return scale_factor


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

def par_writer_refit(A,B,C,DJ,DJK,DK,dJ,dK): #generates SPFIT par file

    constant_list = ['DJ','DJK','DK','dJ','dK']

    DJ_flag = 0
    DJK_flag = 0
    DK_flag = 0
    dJ_flag = 0
    dK_flag = 0

    constants_to_vary = multchoicebox(msg='Choose the distortion constants to vary.  If none are selected, only A, B, and C will be varied in the refined fit.', title='Distortion Constant Choice',choices=constant_list)

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

def cubic_spline(spectrum,new_resolution): # Cubic spline of spectrum to new_resolution; used pre-peak-picking.  Assumes spectrum is already in order of increasing frequency.

    x = spectrum[:,0]
    y = spectrum[:,1]

    old_resolution = (x[-1]-x[0]) / len(spectrum)
    scale_factor = old_resolution / new_resolution

    new_length = int(math.floor(scale_factor*len(spectrum)))

    tck = splrep(x,y,s=0)
    xnew = numpy.arange(x[0],x[-1],new_resolution)
    ynew = splev(xnew,tck,der=0)

    output_spectrum = numpy.column_stack((xnew,ynew))

    return output_spectrum

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

def run_SPCAT(): 
    a = subprocess.Popen("SPCAT default", stdout=subprocess.PIPE, shell=False)
    a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen
 
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
            freq = line[2:13] # Works for frequencies up to 1 THz!
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

def trans_freq_refit_reader(peaklist):
    pred_peaks = cat_reader(1000000, 0,flag="refit")
    pred_trans_for_refit = []

    for line in peaklist:
        for peak in pred_peaks:
            if line[2] == peak[2] and line[3] == peak[3]:
                pred_trans_for_refit.append((peak[1],peak[2],peak[3]))

    return pred_trans_for_refit


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

    threshold = float(enterbox(msg="Enter the OMC threshold in MHz.  Predicted transitions that are not at least this close to an experimental peak will not contribute to the fit."))

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

    trans_1_uncert = float(trans_1[4])
    trans_2_uncert = float(trans_2[4])
    trans_3_uncert = float(trans_3[4])

    while bad_windows ==1:
        while bad_1 == -1:
            bad_wind_decision = buttonbox(msg='The search window for transition 1 (%s - %s) extends below the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?'%(trans_1[2],trans_1[3]), choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_1 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1[2],trans_1[3],trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_1,trans_1_uncert) = choose_new_transition(full_list)
                trans_1_center = float(trans_1[1])
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1[2],trans_1[3],trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_2 == -1:
            bad_wind_decision = buttonbox(msg='The search window for transition 2 (%s - %s) extends below the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?'%(trans_2[2],trans_2[3]), choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_2 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_2[2],trans_2[3],trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_2,trans_2_uncert) = choose_new_transition(full_list)
                trans_2_center = float(trans_2[1])
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2[2],trans_2[3],trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_3 == -1:
            bad_wind_decision = buttonbox(msg='The search window for transition 3 (%s - %s) extends below the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?'%(trans_3[2],trans_3[3]), choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_3 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_3[2],trans_3[3],trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_3,trans_3_uncert) = choose_new_transition(full_list)
                trans_3_center = float(trans_3[1])
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3[2],trans_3[3],trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_1 == 1:
            bad_wind_decision = buttonbox(msg='The search window for transition 1 (%s - %s) extends above the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?'%(trans_1[2],trans_1[3]), choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_1 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_1[2],trans_1[3],trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_1,trans_1_uncert) = choose_new_transition(full_list)
                trans_1_center = float(trans_1[1])
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1[2],trans_1[3],trans_1_center,trans_1_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_2 == 1:
            bad_wind_decision = buttonbox(msg='The search window for transition 2 (%s - %s) extends above the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?'%(trans_2[2],trans_2[3]), choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_2 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_2[2],trans_2[3],trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_2,trans_2_uncert) = choose_new_transition(full_list)
                trans_2_center = float(trans_2[1])
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2[2],trans_2[3],trans_2_center,trans_2_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        while bad_3 == 1:
            bad_wind_decision = buttonbox(msg='The search window for transition 3 (%s - %s) extends above the range of the spectrum.  Would you like to continue, try a new uncertainty for this transition, or choose a new transition?'%(trans_3[2],trans_3[3]), choices=('Continue','New Uncertainty','Choose New Transition'))
            if bad_wind_decision == 'Continue':
                bad_windows = 0
                bad_3 = 0
            elif bad_wind_decision == 'New Uncertainty':
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_3[2],trans_3[3],trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)
            elif bad_wind_decision == 'Choose New Transition':
                (trans_3,trans_3_uncert) = choose_new_transition(full_list)
                trans_3_center = float(trans_3[1])
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3[2],trans_3[3],trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

        if (trans_1 == trans_2) or (trans_2 == trans_3) or (trans_1 == trans_3):
            same_trans_decision = buttonbox(msg='You do not have three distinct transitions.  Choose new ones manually, or quit?', choices=('Choose manually','Quit'))
            if same_trans_decision == 'Quit':
                sys.exit()
            elif same_trans_decision == 'Choose manually':
                (highest_uncert,trans_1,trans_2,trans_3) = triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag="force manual")
                trans_1_uncert = float(trans_1[4])
                trans_2_uncert = float(trans_2[4])
                trans_3_uncert = float(trans_3[4])
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_1[2],trans_1[3],trans_1_center,trans_1_uncert)))
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_2[2],trans_2[3],trans_2_center,trans_2_uncert)))
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 (%s - %s) in MHz.  Its estimated position is %s and its estimated uncertainty is %s MHz."%(trans_3[2],trans_3[3],trans_3_center,trans_3_uncert)))
                (bad_windows,bad_1,bad_2,bad_3)=check_bounds(trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)

    return trans_1,trans_2,trans_3,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty


def triples_gen(window_decision,trans_1_uncert,trans_2_uncert,trans_3_uncert,freq_uncertainty,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,isotopomer_count,decision,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3,scale_factor):

    user_flag = 0
    est_unc_flag = 0
    same_flag = 0

    trans_1_center = float(trans_1[1])
    trans_2_center = float(trans_2[1])
    trans_3_center = float(trans_3[1])

    trans_1_peaks = []
    trans_2_peaks = []
    trans_3_peaks = []

    if window_decision == 'User-defined, same for each':
        user_flag = 1
        same_flag = 1
    elif window_decision == 'User-defined, different for each':
        user_flag = 1
    elif window_decision == 'Three times SPCAT uncertainty':
        est_unc_flag = 1

    while trans_1_peaks == [] or trans_2_peaks == [] or trans_3_peaks == []: #this loops until there are peaks around each member of the triple
        uncertainty_flag =1

        while uncertainty_flag ==1:            
            if freq_uncertainty==0.0 and est_unc_flag ==1 and isotopomer_count == 0:
                peak_1_uncertainty = 3*trans_1_uncert
                peak_2_uncertainty = 3*trans_2_uncert
                peak_3_uncertainty = 3*trans_3_uncert

            if freq_uncertainty==0.0 and user_flag ==1 and same_flag ==1 and isotopomer_count == 0:
                freq_uncertainty = float(enterbox(msg="Enter the frequency uncertainty in MHz.  The largest uncertainty from your fitting peaks is %s MHz."%(highest_uncert)))
                peak_1_uncertainty = freq_uncertainty
                peak_2_uncertainty = freq_uncertainty
                peak_3_uncertainty = freq_uncertainty

            if freq_uncertainty==0.0 and user_flag ==1 and same_flag ==0 and isotopomer_count ==0:
                peak_1_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 1 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_1[2],trans_1[3],trans_1_center,trans_1_uncert)))
                peak_2_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 2 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_2[2],trans_2[3],trans_2_center,trans_2_uncert)))
                peak_3_uncertainty = float(enterbox(msg="Enter the frequency uncertainty for transition 3 (%s - %s) in MHz.  Its estimated position is %s and its uncertainty is %s MHz."%(trans_3[2],trans_3[3],trans_3_center,trans_3_uncert)))
            
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
            
            estimated_time = float(num_of_triples) * scale_factor

            if estimated_time <= 90:
                time_string = str(int(estimated_time)) + " seconds"
            elif estimated_time > 90 and estimated_time <= 5400:
                time_min = math.ceil(estimated_time / 60.0)
                time_string = str(int(time_min)) + " minutes"
            elif estimated_time > 5400 and estimated_time <= 129600:
                time_hours = math.ceil(estimated_time / 3600.0)
                time_string = str(int(time_hours)) + " hours"
            elif estimated_time > 129600:
                time_days = math.ceil(estimated_time / 86400.0)
                time_string = str(int(time_days)) + " days"

            if isotopomer_count == 0:
                decision = buttonbox(msg='There are %s triples in this calculation, which will take roughly %s. Would you like to continue, try new uncertainty, or quit?'%(str(num_of_triples),time_string), choices=('Continue','Quit','new uncertainty'))
            
            if decision == 'Quit':
                sys.exit()
            elif decision == 'Continue':
                uncertainty_flag = 0
            else:
                freq_uncertainty = 0.0
                est_unc_flag = 0
                user_flag = 1
    
    return trans_1,trans_2,trans_3,trans_1_peaks,trans_2_peaks,trans_3_peaks,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,num_of_triples,decision

def refine_fits(job_name,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fits,DJ,DJK,DK,dJ,dK,peaklist,main_flow):

    a = subprocess.Popen("mkdir refits")
    a.wait()

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

    all_fitting_done = 0

    while all_fitting_done == 0:

        add_more_transitions = 0

        msg = "Choose a previous result for full fitting (with distortions)."
        title = "Microwave Fitting Program"
        choice = choicebox(msg,title,fits)

        result_choice = choice.split()[0]
        A_fit = choice.split()[6]
        B_fit = choice.split()[7]
        C_fit = choice.split()[8]

        fit_peaklist = refit_peaks_list + refit_check_peaks_list
        fit_peaklist = list(OrderedDict.fromkeys(fit_peaklist)) # This removes duplicates if some transitions are both check and fit transitions, which might bias SPFIT towards hitting those in particular.

        trans_display = ""
        for entry in fit_peaklist:
            trans_display = trans_display + "%s \n"%(str(entry))

        codebox(msg='These are the transitions that will be used in the refined fit.',text=trans_display)

        more_transitions_decision = buttonbox(msg='Would you like to add additional transitions into the expanded fit?', choices=('Yes','No'))
        if more_transitions_decision == 'Yes':
            add_more_transitions = 1

        var_writer(A_fit,B_fit,C_fit,DJ,DJK,DK,dJ,dK,main_flow,flag="refit")
        run_SPCAT_refit()

        if add_more_transitions == 1:
            extra_peaks = multchoicebox(msg='Choose the transitions to use in the new fit.  Previous fitting and scoring transitions will also be included.  Duplicates will be removed automatically.', title='Result Sorting Setup',choices=cat_reader(freq_high, freq_low,flag="refit"))
            extra_peaks_clean = []
            for entry in extra_peaks:
                #clean = entry[2:43]
                clean = entry[2:44] # Needed for 1 THz freq extension?
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

        updated_trans = trans_freq_refit_reader(fit_peaklist) # Finds updated predicted frequencies with improved A, B, and C estimates.

        fitting_done = 0

        while fitting_done == 0:

            par_writer_refit(A_fit,B_fit,C_fit,DJ,DJK,DK,dJ,dK)
            best_matches = match_to_peaklist(updated_trans,peaklist) # Assigns closest experimental peak frequencies to transitions
            lin_writer_refit(best_matches)

            a = subprocess.Popen("spfit0 refit", stdout=subprocess.PIPE, shell=False)
            a.stdout.read()#used to let SPFIT finish

            SPFIT_results = open("refit.fit",'r')
            codebox(msg='SPFIT has finished.',text=SPFIT_results)
            SPFIT_results.close()

            fit_decision = buttonbox(msg='Fitting has finished. Would you like to accept the current fit or try again while changing inclusion threshold or allowing different distortions to vary?', choices=('Accept Current Fit','Try Again'))
            if fit_decision == 'Accept Current Fit':
                fitting_done = 1
                shutil.copyfile("refit.fit","refits/refit%s.fit"%(str(result_choice)))
                shutil.copyfile("refit.int","refits/refit%s.int"%(str(result_choice)))
                shutil.copyfile("refit.lin","refits/refit%s.lin"%(str(result_choice)))
                shutil.copyfile("refit.par","refits/refit%s.par"%(str(result_choice)))
                shutil.copyfile("refit.var","refits/refit%s.var"%(str(result_choice)))

        all_fit_decision = buttonbox(msg='Would you like to refine another previous result or are you done with this species?', choices=('Refine Another Result','Done with this one!'))
        if all_fit_decision == 'Done with this one!':
            all_fitting_done = 1


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
        sys.exit()

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
                    sys.exit()
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

def read_coords_file():

    a = []
    b = []
    c = []
    an = []
    mass = []
    atom_list = []

    has_H_flag = 0
    has_C_flag = 0
    has_N_flag = 0
    has_O_flag = 0
    has_Si_flag = 0
    has_S_flag = 0
    has_Cl_flag = 0
    has_Br_flag = 0

    f = open(fileopenbox(msg="Choose Coordinates file"))

    type_check = []

    for line in f:
        type_check_temp = line.split()
        type_check.append(type_check_temp)

    f.close()

    if len(type_check[0]) == 1: # XYZ format.  First line is # of atoms, second is a comment.
        for line in type_check:
            if len(line) == 4:
                a.append(line[1])
                b.append(line[2])
                c.append(line[3])
                an.append(line[0])

    if len(type_check[0]) == 5: # One type of Gaussian output.
        for line in type_check:
            a.append(line[2])
            b.append(line[3])
            c.append(line[4])
            an.append(line[1])

    if len(type_check[0]) == 6: # Another type of Gaussian output.
        for line in type_check:
            a.append(line[3])
            b.append(line[4])
            c.append(line[5])
            an.append(line[1])

    for entry in an:  #changes atomic numbers to masses of NS; these are all the atoms in molecules I've taken spectra of.  I need to take data on more types of atoms!
        if entry =='1' or entry == 'H':
            mass.append(1.007825037)
            if has_H_flag == 0:
                has_H_flag = 1
                atom_list.append('H')

        if entry =='6' or entry == 'C':
            mass.append(12)
            if has_C_flag == 0:
                has_C_flag = 1
                atom_list.append('C')

        if entry =='7' or entry == 'N':
            mass.append(14.003074008)
            if has_N_flag == 0:
                has_N_flag = 1
                atom_list.append('N')

        if entry =='8' or entry == 'O':
            mass.append(15.99491464)
            if has_O_flag == 0:
                has_O_flag = 1
                atom_list.append('O')

        if entry =='9' or entry == 'F':
            mass.append(18.998403) # No "has_F_flag" because only one stable isotope.  Change if you are working with molecules from an 18F production facility.

        if entry =='14' or entry == 'Si':
            mass.append(27.9769265325)
            if has_Si_flag == 0:
                has_Si_flag = 1
                atom_list.append('Si')

        if entry =='16' or entry == 'S':
            mass.append(31.97207100)
            if has_S_flag == 0:
                has_S_flag = 1
                atom_list.append('S')

        if entry =='17' or entry == 'Cl':
            mass.append(34.96885268)
            if has_Cl_flag == 0:
                has_Cl_flag = 1
                atom_list.append('Cl')

        if entry =='35' or entry == 'Br':
            mass.append(78.9183371)
            if has_Br_flag == 0:
                has_Br_flag = 1
                atom_list.append('Br')

        if entry =='53' or entry == 'I':
            mass.append(126.904473) # Iodine also only has one stable isotope.

    return a,b,c,mass,atom_list

def isotopomers(expt_A,expt_B,expt_C): # Uses a file of coordinates (Gaussian output format) + expt'l NS constants to generate A,B,C values of isotopomers

    H_vary_flag = 0
    C_vary_flag = 0
    N_vary_flag = 0
    O_vary_flag = 0
    Si_vary_flag = 0
    S_vary_flag = 0
    Cl_vary_flag = 0
    Br_vary_flag = 0

    (a,b,c,mass,atom_list) = read_coords_file()

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

    atoms_to_vary = multchoicebox(msg='Choose atoms to isotopically substitute', title='Isotope Substitution',choices=atom_list)

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
    
    final_omc = []
    triples_counter = 0
    output_file = ""
    regular_counter = 0
    #error_counter = 0



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
            if line.split()[0] == "10000":
                temp_A = float(line.split()[1])
                const_list.append("%.3f" %temp_A)
            if line.split()[0] == "20000":
                temp_B = float(line.split()[1])
                const_list.append("%.3f" %temp_B)
            if line.split()[0] == "30000":
                temp_C = float(line.split()[1])
                const_list.append("%.3f" %temp_C)

        fh_fit = open("default%s.fit"%(str(file_num)))
        file_list = []
        for line in fh_fit:
                file_list.append(line)

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
        
        theor_inten_list = []
        for x in range(len(top_17)):
            temp_inten = 10**float(top_17[x][0])
            theor_inten_list.append(temp_inten)

        theor_inten_ratio = (max(theor_inten_list)/min(theor_inten_list))
        theor_inten_avg = (sum(theor_inten_list)/len(theor_inten_list))
        theor_inten_unitless_stdev = numpy.std(theor_inten_list)/theor_inten_avg

        for x in range(len(top_17)): #matches peaks in the top 17 to peaks in experimental peak list <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            qnum_up = top_17[x][2]
            qnum_low = top_17[x][3]
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
            omc_inten = 100000.0
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
            omc_list.append((omc_low, real_omc, omc_inten))# current_peak,qnum_up,qnum_low)) you can add in this extra output, but its slower
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
    #os.system("sort /R final_output%s.txt > sorted_final_out%s.txt"%(str(file_num),str(file_num)))#sorts output by score
    
if __name__ == '__main__': #multiprocessing imports script as module
    multiprocessing.freeze_support()
    
    u_A = "0.0"#<<<<<<<<<<<<<set default value here
    u_B = "1.0"#<<<<<<<<<<<<<set default value here
    u_C = "0.0"#<<<<<<<<<<<<<set default value here
    A = "9769.6221"#<<<<<<<<<<<<<set default value here 
    B = "868.84665"#<<<<<<<<<<<<<set default value here
    C = "818.51874"#<<<<<<<<<<<<<set default value here, etc...
    DJ = "-0.00004723"
    DJK = "0.0008991"
    DK = "-0.02317"
    dJ = "-0.000005029"
    dK = "-0.000343"
    processors = 6
    #freq_high = float(18000.0)
    #freq_low = float(6000.0)
    inten_high = float('100000.0')
    inten_low = float('4E-004')
    temperature="2"
    Jmax="20"     
        

    main_flow = buttonbox(msg='Are you searching for a normal species spectrum or singly-substituted isotopologues (already have NS experimental constants)?', choices=('Normal species','Isotopologues'))

    scale_factor = time_estimate()

    x = subprocess.Popen("dir", stdout=subprocess.PIPE, shell=True)
    x = x.stdout.read().split()

    file_decision = ynbox(msg="Would you like to use an input file?",title='Input file choice')
        
    if file_decision == True:
        file_decision_flag = 1
       
    elif file_decision == False:
        file_decision_flag = 0
    
    file_flag = 1
    
    while file_flag==1:
        filename = multenterbox("Choose a folder name for the calculation. If the folder already exists, this dialogue will reset.","",["Job Name"] )
        job_name = filename[0]
        marker1 = 0
        for file1 in x:
            if job_name ==file1:
                marker1 = 1
        if marker1 ==0:
            file_flag =0
        
    #a = subprocess.Popen("mkdir %s"%job_name)
    #a.wait()

    if os.path.exists(job_name):
    	job_name = job_name + 'a'

    os.mkdir(job_name)

    fitting_peaks_flag =0
    
    check_peaks_list = []
    fit_peaks_list = []
    freq_uncertainty = 0.0
    if file_decision_flag ==1:

        fh_input = open(fileopenbox(msg="enter the input file"))
        
        for line in fh_input:
            
            if line.split() != []:
                                
                if line.split()[0] == "u_A:":
                    u_A = line.split()[1] 
                if line.split()[0] == "u_B:":
                    u_B = line.split()[1] 
                if line.split()[0] == "u_C:":
                    u_C = line.split()[1] 
                if line.split()[0] == "A:":
                    A = line.split()[1] 
                if line.split()[0] == "B:":
                    B = line.split()[1]
                if line.split()[0] == "C:":
                    C = line.split()[1]
                if line.split()[0] == "DJ:":
                    DJ = line.split()[1]
                if line.split()[0] == "DK:":
                    DK = line.split()[1]
                if line.split()[0] == "DJK:":
                    DJK = line.split()[1] 
                if line.split()[0] == "dJ:":
                    dJ = line.split()[1]
                if line.split()[0] == "dK:":
                    dK = line.split()[1]
                if line.split()[0] == "freq_high:":
                    freq_high = float(line.split()[1])  
                if line.split()[0] == "freq_low:":
                    freq_low = float(line.split()[1])  
                if line.split()[0] == "inten_high:":
                    inten_high = float(line.split()[1]) 
                if line.split()[0] == "inten_low:":
                    inten_low = float(line.split()[1])
                if line.split()[0] == "processors:":
                    processors = int(line.split()[1])
                if line.split()[0] == "Temp:":
                    temperature = float(line.split()[1])
                if line.split()[0] == "Jmax:":
                    Jmax = float(line.split()[1])
                if line.split()[0] == "freq_uncertainty:":
                    freq_uncertainty = float(line.split()[1])                    
                if line.split()[0] == "trans_1:" or line.split()[0] == "trans_2:" or line.split()[0] == "trans_3:":
                    fitting_peaks_flag = 0
                    clean = line[12:53]
                    re_split = clean.split("', '")
                    tuples = tuple(re_split)
                    fit_peaks_list.append(tuples)
                  
                    
                if fitting_peaks_flag == 1:
                    clean = line[2:43]
                    re_split = clean.split("', '")
                    tuples = tuple(re_split)
                    #print tuples
                    check_peaks_list.append(tuples)
                if line.split()[0] == "Check":
                    fitting_peaks_flag = 1
                                    
    if file_decision_flag ==0:
        if main_flow == 'Isotopologues':
            msg = "Enter constants for the normal species isotopologue in MHz. If left blank, they will be set to default values (in parentheses). The inputs will be sent directly to SPCAT, so remember to input -DJ, -DJK, etc. "

        elif main_flow == 'Normal species':
            msg = "Enter your estimated constants in MHz. If left blank, they will be set to default values (in parentheses). The inputs will be sent directly to SPCAT, so remember to input -DJ, -DJK, etc. "
        
        title = "Spectral Fitting Program"
        fieldNames = ["a dipole (%s)"%(u_A),"b dipole (%s)"%(u_B),"c dipole (%s)"%(u_C),"A (%s MHz)"%(A),"B (%s MHz)"%(B),"C (%s MHz)"%(C),"-DJ (%s MHz)"%(DJ),"-DJK (%s MHz)"%(DJK),"-DK (%s MHz)"%(DK),"-dJ (%s MHz)"%(dJ),"-dK (%s MHz)"%(dK),\
                      "number of processors (%s)"%(processors),"upper intensity limit on exp spectrum (%s)"%(inten_high),"lower intensity limit on exp spectrum(%s)"%(inten_low),"temperature (%sK)"%(temperature), "Jmax (%s)"%(Jmax)]

        const = [] 
        const = multenterbox(msg,title, fieldNames) # gui for grabbing data
         
        if const[0]!='':
            u_A = const[0]
        if const[1]!='':
            u_B = const[1] 
        if const[2]!='':
            u_C = const[2]
        if const[3]!='':
            A = const[3]
        if const[4]!='':
            B = const[4]
        if const[5]!='':
            C = const[5]
        if const[6]!='':
            DJ = const[6]
        if const[7]!='':
            DJK = const[7]
        if const[8]!='':
            DK = const[8]
        if const[9]!='':
            dJ = const[9]
        if const[10]!='':
            dK = const[10]
        if const[11] != '':
            processors = int(const[11])
        #if const[12]!='':
        #    freq_high = float(const[12])
        #if const[13]!='':
        #    freq_low = float(const[13])
        if const[12]!='':
            inten_high = float(const[12])
        if const[13]!='':
            inten_low = float(const[13])
        if const[14]!='':
            temperature = const[14]
        if const[15]!='':
            Jmax = const[15]

    fh = numpy.loadtxt(fileopenbox(msg="Enter the spectrum file in two column format: frequency intensity")) #loads full experimental data file, not just list of peaks

    spline_level = integerbox(msg="Enter the resolution to spline the spectrum to, in kHz (minimum of 1 kHz).", title='Spectral Resolution', default=2, lowerbound=1, upperbound=100000)
    spline_level_MHz = float(spline_level)/1000

    spectrum_2kHz = cubic_spline(fh,spline_level_MHz) # Interpolates experimental spectrum to a 2 kHz resolution with a cubic spline.  Gives better peak-pick values.
    (peaklist, freq_low, freq_high) = peakpicker(spectrum_2kHz,inten_low,inten_high) # Calls slightly modified version of Cristobal's routine to pick peaks instead of forcing user to do so.
    
    for number in range(processors):
        y = subprocess.Popen("copy SPFIT.EXE %s\SPFIT%s.EXE"%(job_name,number), stdout=subprocess.PIPE, shell=True)
        y.stdout.read()     

    y = subprocess.Popen("copy SPCAT.EXE %s\SPCAT.EXE"%(job_name), stdout=subprocess.PIPE, shell=True)
    y.stdout.read() 
    
    os.chdir(job_name)
    
    trans_1 = "" 
    trans_2 = ""
    trans_3 = ""
    int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
    var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
    run_SPCAT()

    full_list = cat_reader(freq_high,freq_low,flag="default")	# New code starts here

    if main_flow == 'Isotopologues':

        filter_level = float(enterbox(msg="Isotopologue search transitions are filtered to only include strong NS transitions in the experimental spectrum.  Enter a scale factor (relative to the lower intensity threshold) to use in this filtering.  Enter zero if NS lines are not present."))
        full_list = intensity_filter(full_list,peaklist,inten_low,filter_level)

    (highest_uncert,trans_1,trans_2,trans_3) = triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag="undecided")

    trans_1_uncert = float(trans_1[4])
    trans_2_uncert = float(trans_2[4])
    trans_3_uncert = float(trans_3[4])

    if main_flow == 'Normal species':
        trans_1_center = float(trans_1[1])
        trans_2_center = float(trans_2[1])
        trans_3_center = float(trans_3[1])

        user_flag = 0
        est_unc_flag = 0
        same_flag = 0
        bad_windows = 0

        isotopomer_count = 0 # Allows us to use some of the same functions for both the NS and isotopologue searches

        window_decision = buttonbox(msg='What method would you prefer to use for determining the uncertainty window?', choices=('User-defined, same for each','User-defined, different for each','Three times SPCAT uncertainty'))

        peak_1_uncertainty = 0
        peak_2_uncertainty = 0
        peak_3_uncertainty = 0
        decision = ""

        (trans_1,trans_2,trans_3,trans_1_peaks,trans_2_peaks,trans_3_peaks,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,num_of_triples,decision) = triples_gen(window_decision,trans_1_uncert,trans_2_uncert,trans_3_uncert,freq_uncertainty,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,isotopomer_count,decision,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3,scale_factor)

        if check_peaks_list == []:                                        
            int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
            var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
            run_SPCAT()
            
            choice_method = buttonbox(msg='Choose the method for scoring the results:', title='Scoring Options',choices=\
                            ('Use the 10 most intense peaks (default)','Choose an arbitrary set of peaks from the prediction'))
            if choice_method == 'Use the 10 most intense peaks (default)':
                top_peaks = cat_reader(freq_high, freq_low,flag="default")[0:10] #grab the most intense peaks from the predicted spectrum
            else:
                top_peaks = multchoicebox(msg='Pick the transitions that you want to use for sorting the results', title='Result Sorting Setup',choices=cat_reader(freq_high, freq_low,flag="default"))
                top_peaks_clean = []
                for entry in top_peaks:
                    #clean = entry[2:43]
                    clean = entry[2:44] # Needed for 1 THz freq extension?
                    re_split = clean.split("', '")
                    tuples = tuple(re_split)
                    top_peaks_clean.append(tuples)
                top_peaks = top_peaks_clean
        else:
            top_peaks = check_peaks_list

        top_peaks_3cut = []
        for entry in top_peaks:
            if (entry[2] == trans_1[2] and entry[3] == trans_1[3]) or (entry[2] == trans_2[2] and entry[3] == trans_2[3]) or (entry[2] == trans_3[2] and entry[3] == trans_3[3]):
                pass
            else:
                top_peaks_3cut.append(entry)

        #print num_of_triples

        job_file = ""
        str(top_peaks_3cut)
        str(trans_1)
        str(trans_2)
        str(trans_3)
        fitting_peaks_str = ""
        for entry in top_peaks_3cut:
            fitting_peaks_str+=str(entry)+"\n"
            
        job_file += "Job Name %s \n u_A: %s \n u_B: %s \n u_C: %s \n A: %s \n B: %s \n \
    C: %s \n DJ: %s \n DJK: %s \n DK: %s \n dJ: %s \n dK: %s \n processors: %s \n freq_high: %s \n freq_low: %s \n \
    inten_high: %s \n inten_low: %s \n Temp: %s \n Jmax: %s \n freq_uncertainty: %s \n number of triples: %s \n Check peaks:\n%s \n trans_1: %s \n trans_2: %s \n trans_3: %s "%(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,str(processors),str(freq_high),\
        str(freq_low),str(inten_high),str(inten_low),str(temperature),str(Jmax),str(freq_uncertainty),str(num_of_triples),fitting_peaks_str,str(trans_1),str(trans_2),str(trans_3))
        Job_fh = open("input_data_%s.txt"%(job_name),"w")
        Job_fh.write(job_file) 
        Job_fh.close()

        new_list = []

        new_list = [(len(trans_1_peaks),"trans_1_peaks"),(len(trans_2_peaks),"trans_2_peaks"),(len(trans_3_peaks),"trans_3_peaks")]
        new_list.sort()
            
        list_key = []
        list_a_peaks = [vars()[new_list[0][1]],new_list[0][1]]
        list_b_peaks = [vars()[new_list[1][1]],new_list[1][1]]
        list_c_peaks = vars()[new_list[2][1]]

        random.shuffle(list_c_peaks)  # Shuffle so that each processor gets a range of values for the third peak, not processor 0 getting only the lowest frequencies.  

        list_c_list = []
        for num in range(processors):
            processors = float(processors)
            num = float(num)
            x = int((num)*(len(list_c_peaks)/processors))
            y = int(len(list_c_peaks)*((num+1)/processors))
            list_c_list.append(list_c_peaks[x:y])
        list_c_list.append("marker")
        vars()[new_list[0][1]] = list_a_peaks[0]
        vars()[new_list[1][1]] = list_b_peaks[0]
        vars()[new_list[2][1]] = list_c_list

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
            
        #a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 4 -n > sorted_omc_cat.txt', shell=True)
        #a.wait()
    
        #a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 6 -n > sorted_inten_omc_cat.txt', shell=True)
        #a.wait()

        full_list_for_sort = []

        for i in range(processors):
            fh = open("final_output%s.txt"%str(i))
            for line in fh:
                full_list_for_sort.append(string.split(str(line)))
            fh.close()

        omc_sort = sorted(full_list_for_sort, key=lambda entry: entry[11])
        inten_omc_sort = sorted(full_list_for_sort, key=lambda entry: entry[24])

        fh_omc = open("sorted_omc_cat.txt",'w')
        for item in omc_sort:
            tempstring = ""
            for entry in item:
                tempstring = tempstring + str(entry) + " "

            tempstring = tempstring + "\n"
            fh_omc.write(tempstring)

        fh_omc.close()

        fh_inten_omc = open("sorted_inten_omc_cat.txt",'w')
        for item in inten_omc_sort:
            tempstring = ""
            for entry in item:
                tempstring = tempstring + str(entry) + " "

            tempstring = tempstring + "\n"
            fh_inten_omc.write(tempstring)

        fh_inten_omc.close()

        f = open('sorted_inten_omc_cat.txt','r')
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

        f100=open('best100.txt','w')
        f100.write(OMC_char_buffer)
        f100.close()

        codebox(msg='Fitting routine has finished.  These are the best 100 results, saved in best100.txt.',text=fits)

        refine_fit_decision = buttonbox(msg='Would you like to further refine any of the recent results by allowing distortions to vary?', choices=('Yes!',"No, I'm done!"))
        if refine_fit_decision == "No, I'm done!":
            pass
        else:
            isotope_ID = "NS only" # Allows us to use the same base function (refine_fits) for both the NS only and the isotopologue branches.
            refine_fits(job_name,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fits,DJ,DJK,DK,dJ,dK,peaklist,main_flow)

        continue_decision = buttonbox(msg='Would you like to perform a search for singly-substituted isotopologues?', choices=('Yes','No'))
        if continue_decision == 'No':
            sys.exit()
        elif continue_decision == 'Yes':
            main_flow = 'Isotopologues'
            msg = "Enter constants for the normal species isotopologue in MHz. If left blank, they will be set to default values (in parentheses). The inputs will be sent directly to SPCAT, so remember to input -DJ, -DJK, etc. "

            const = [] 
            const = multenterbox(msg,title, fieldNames) # gui for grabbing data
         
            if const[0]!='':
                u_A = const[0]
            if const[1]!='':
                u_B = const[1] 
            if const[2]!='':
                u_C = const[2]
            if const[3]!='':
                A = const[3]
            if const[4]!='':
                B = const[4]
            if const[5]!='':
                C = const[5]
            if const[6]!='':
                DJ = const[6]
            if const[7]!='':
                DJK = const[7]
            if const[8]!='':
                DK = const[8]
            if const[9]!='':
                dJ = const[9]
            if const[10]!='':
                dK = const[10]
            if const[11] != '':
                processors = int(const[11])
            if const[12]!='':
                inten_high = float(const[12])
            if const[13]!='':
                inten_low = float(const[13])
            if const[14]!='':
                temperature = const[14]
            if const[15]!='':
                Jmax = const[15]

            trans_1 = "" 
            trans_2 = ""
            trans_3 = ""
            int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
            var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
            run_SPCAT()

            full_list = cat_reader(freq_high,freq_low,flag="default")   # New code starts here

            filter_level = float(enterbox(msg="Isotopologue search transitions are filtered to only include strong NS transitions in the experimental spectrum.  Enter a scale factor (relative to the lower intensity threshold) to use in this filtering.  Enter zero if NS lines are not present."))
            full_list = intensity_filter(full_list,peaklist,inten_low,filter_level)

            (highest_uncert,trans_1,trans_2,trans_3) = triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag="undecided")

            trans_1_uncert = float(trans_1[4])
            trans_2_uncert = float(trans_2[4])
            trans_3_uncert = float(trans_3[4])

    if main_flow == 'Isotopologues':
        isotopomer_count = 0

        calc_isotopomer_ABC = isotopomers(float(A),float(B),float(C))
    
        for isotopologue in calc_isotopomer_ABC:
            isotope_ID = isotopologue[0]
            curr_A = isotopologue[1]
            curr_B = isotopologue[2]
            curr_C = isotopologue[3]

            a = subprocess.Popen("mkdir %s"%isotope_ID)
            a.wait()

            for number in range(processors):
                y = subprocess.Popen("copy SPFIT%s.EXE %s\SPFIT%s.EXE"%(number,isotope_ID,number), stdout=subprocess.PIPE, shell=True)
                y.stdout.read()     

            y = subprocess.Popen("copy SPCAT.EXE %s\SPCAT.EXE"%(isotope_ID), stdout=subprocess.PIPE, shell=True)
            y.stdout.read() 

            os.chdir(isotope_ID)

            int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
            var_writer(curr_A,curr_B,curr_C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
            run_SPCAT()

            trans_1_center,trans_2_center,trans_3_center = trans_freq_reader(trans_1,trans_2,trans_3)

            # This updates the center frequencies of the fitting transitions to account for the new rotational constants.
            # If this isn't done, then there are issues when trans_1,trans_2,trans_3 are passed over to fit_triples.
            trans_1 = (trans_1[0],trans_1_center,trans_1[2],trans_1[3],trans_1[4])
            trans_2 = (trans_2[0],trans_2_center,trans_2[2],trans_2[3],trans_2[4])
            trans_3 = (trans_3[0],trans_3_center,trans_3[2],trans_3[3],trans_3[4])

            trans_1_center = float(trans_1_center)
            trans_2_center = float(trans_2_center)
            trans_3_center = float(trans_3_center)

            if (trans_1_center == 0) or (trans_2_center == 0) or (trans_3_center == 0):
                bad_triples = buttonbox(msg='One of your fitting transitions is no longer within bounds after isotopologue scaling.  Please manually choose new transitions or quit.', choices=('Choose new transitions','Quit'))
                if bad_triples == 'Quit':
                    sys.exit()
                elif bad_triples == 'Choose new transitions':
                    (highest_uncert,trans_1,trans_2,trans_3) = triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag="force manual")
                    trans_1_uncert = float(trans_1[4])
                    trans_2_uncert = float(trans_2[4])
                    trans_3_uncert = float(trans_3[4])


        
            user_flag = 0
            est_unc_flag = 0
            same_flag = 0

            if isotopomer_count == 0:
                window_decision = buttonbox(msg='What method would you prefer to use for determining the uncertainty window?', choices=('User-defined, same for each','User-defined, different for each','Three times SPCAT uncertainty'))
                peak_1_uncertainty = 0
                peak_2_uncertainty = 0
                peak_3_uncertainty = 0
                decision = ""

            (trans_1,trans_2,trans_3,trans_1_peaks,trans_2_peaks,trans_3_peaks,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,num_of_triples,decision) = triples_gen(window_decision,trans_1_uncert,trans_2_uncert,trans_3_uncert,freq_uncertainty,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,isotopomer_count,decision,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3,scale_factor)

            if check_peaks_list == []:            
                int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temperature=temperature,flag="default")
                var_writer(curr_A,curr_B,curr_C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
                run_SPCAT()
        
                if isotopomer_count == 0:
                     choice_method = buttonbox(msg='Choose the method for scoring the results:', title='Scoring Options',choices=\
                                      ('Use the 10 most intense peaks (default)','Choose an arbitrary set of peaks from the prediction'))
                    
                if choice_method == 'Use the 10 most intense peaks (default)':
#                    top_peaks = cat_reader(freq_high, freq_low,flag="default")[0:10] #grab the most intense peaks from the predicted spectrum
                    top_peaks_unfiltered = cat_reader(freq_high,freq_low,flag="default")
                    temp_top_peaks = []
                    for entry in top_peaks_unfiltered:
                        for entry2 in full_list:
                            if ((entry[2] == entry2[2]) and (entry[3] == entry2[3])):
                                temp_top_peaks.append(entry)
                                break

                    if len(temp_top_peaks) >= 10:
                        top_peaks = temp_top_peaks[0:10]
                    else:
                        top_peaks = temp_top_peaks

                else:
                    if isotopomer_count == 0:
                        top_peaks_unfiltered = cat_reader(freq_high,freq_low,flag="default")
                        temp_top_peaks = []
                        for entry in top_peaks_unfiltered:
                            for entry2 in full_list:
                                if ((entry[2] == entry2[2]) and (entry[3] == entry2[3])):
                                    temp_top_peaks.append(entry)
                                    break

                        top_peaks = multchoicebox(msg='Pick the transitions that you want to use for sorting the results', title='Result Sorting Setup',choices=temp_top_peaks)
                        top_peaks_clean = []
                        for entry in top_peaks:
                            #clean = entry[2:43]
                            clean = entry[2:44] # Needed for 1 THz freq extension?
                            re_split = clean.split("', '")
                            tuples = tuple(re_split)
                            top_peaks_clean.append(tuples)
                        top_peaks = top_peaks_clean
                        top_original_peaks = top_peaks

                    else:
                        updated_peaks = cat_reader(100000,0,flag="default")
                        updated_top_peaks = []
                        for entry in top_original_peaks:
                            for peak in updated_peaks:
                                if entry[2]==peak[2] and entry[3]==peak[3]:
                                    updated_top_peaks.append(peak)
                        top_peaks = updated_top_peaks
            else:
                top_peaks = check_peaks_list
    
            top_peaks_3cut = []
            for entry in top_peaks:
                if (entry[2] == trans_1[2] and entry[3] == trans_1[3]) or (entry[2] == trans_2[2] and entry[3] == trans_2[3]) or (entry[2] == trans_3[2] and entry[3] == trans_3[3]):
                    pass
                else:
                    top_peaks_3cut.append(entry)

            #print num_of_triples

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

            new_list = []

            new_list = [(len(trans_1_peaks),"trans_1_peaks"),(len(trans_2_peaks),"trans_2_peaks"),(len(trans_3_peaks),"trans_3_peaks")]
            new_list.sort()
    
            list_key = []

            list_a_peaks = [vars()[new_list[0][1]],new_list[0][1]]
            list_b_peaks = [vars()[new_list[1][1]],new_list[1][1]]
            list_c_peaks = vars()[new_list[2][1]]

            random.shuffle(list_c_peaks)  # Shuffle so that each processor gets a range of values for the third peak, not processor 0 getting only the lowest frequencies.  

            list_c_list = []
            for num in range(processors):
                processors = float(processors)
                num = float(num)
                x = int((num)*(len(list_c_peaks)/processors))
                y = int(len(list_c_peaks)*((num+1)/processors))
                list_c_list.append(list_c_peaks[x:y])
            list_c_list.append("marker")
            vars()[new_list[0][1]] = list_a_peaks[0]
            vars()[new_list[1][1]] = list_b_peaks[0]
            vars()[new_list[2][1]] = list_c_list
    
    
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

            isotopomer_count += 1
            os.chdir(os.pardir)

        for isotopologue in calc_isotopomer_ABC:

            isotope_ID = isotopologue[0]

            os.chdir(isotope_ID)

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

            codebox(msg='Fitting routine has finished.  These are the best 100 results for %s, saved in best100_%s.txt.'%(isotope_ID,isotope_ID),text=fits)

            refine_fit_decision = buttonbox(msg='Would you like to further refine any of the recent results by allowing distortions to vary?', choices=('Yes!',"No, I'm done with this isotopologue!"))
            if refine_fit_decision == "No, I'm done with this isotopologue!":
                pass
            else:
                refine_fits(job_name,isotope_ID,u_A,u_B,u_C,Jmax,freq_high,temperature,fits,DJ,DJK,DK,dJ,dK,peaklist,main_flow)

            os.chdir(os.pardir)





    
