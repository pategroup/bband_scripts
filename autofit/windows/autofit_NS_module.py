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

autofit NS module:
-stripping out lots of functionality (using v15 of autofit as a base) to be run from within the GUI instead of individually.
Removed comments describing previous versions of stand-alone autofit because not all of the things described are available here.
This should be able to do a standard fitting run on a normal species case without popping up any additional dialog boxes if
autofit_NS is called with the appropriate arguments.  Currently isotopologue searches and fit refinements aren't included on
completion; those will probably be separated modules later.


"""

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

    all_combo_list_file.close()
    
    all_combo_list_file = open(all_combo_file)
    
    triples_counter = 0
    output_file = ""
    regular_counter = 0

    theor_inten_list = [] # Theoretical intensity ratios and unitless standard deviations for use in scoring fit results below.
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
    fh_final.write(output_file)
    fh_final.close()
    os.system("sort -r 'final_output%s.txt'>sorted_final_out%s.txt"%(str(file_num),str(file_num)))#sorts output by score
    
def autofit_NS(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,freq_high,freq_low,inten_high,inten_low,processors,temperature,Jmax,trans_1,trans_2,trans_3,check_peaks_list,peaklist,trans_1_peaks,trans_2_peaks,trans_3_peaks):

    a = subprocess.Popen("mkdir %s"%job_name) # Need to be able to trust job_name.  Add error handling here later / build it into the GUI.
    a.wait()

    for number in range(processors):
        y = subprocess.Popen("cp SPFIT.EXE %s\SPFIT%s.EXE"%(job_name,number), stdout=subprocess.PIPE, shell=True)
        y.stdout.read()     

    y = subprocess.Popen("cp SPCAT.EXE %s\SPCAT.EXE"%(job_name), stdout=subprocess.PIPE, shell=True)
    y.stdout.read() 
    
    os.chdir(job_name)
    
    num_of_triples = len(trans_1_peaks)*len(trans_2_peaks)*len(trans_3_peaks) #this tells you how many entries there will be in the all_combo_list

    trans_1_uncert = float(trans_1[4])
    trans_2_uncert = float(trans_2[4])
    trans_3_uncert = float(trans_3[4])

    trans_1_center = float(trans_1[1])
    trans_2_center = float(trans_2[1])
    trans_3_center = float(trans_3[1])

    top_peaks_3cut = []
    for entry in check_peaks_list:
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
            
    job_file += "Job Name %s \n u_A: %s \n u_B: %s \n u_C: %s \n A: %s \n B: %s \n C: %s \n DJ: %s \n DJK: %s \n DK: %s \n dJ: %s \n dK: %s \n processors: %s \n freq_high: %s \n freq_low: %s \n inten_high: %s \n inten_low: %s \n Temp: %s \n Jmax: %s \n number of triples: %s \n Check peaks:\n%s \n trans_1: %s \n trans_2: %s \n trans_3: %s "%(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,str(processors),str(freq_high),\
        str(freq_low),str(inten_high),str(inten_low),str(temperature),str(Jmax),str(num_of_triples),fitting_peaks_str,str(trans_1),str(trans_2),str(trans_3))
    
    Job_fh = open("input_data_%s.txt"%(job_name),"w")
    Job_fh.write(job_file) 
    Job_fh.close()

    new_list = []

    new_list = [(len(trans_1_peaks),"trans_1_peaks"),(len(trans_2_peaks),"trans_2_peaks"),(len(trans_3_peaks),"trans_3_peaks")]
    new_list.sort()
            
    #list_key = []
    #list_a_peaks = [globals()[new_list[0][1]],new_list[0][1]]
    #list_b_peaks = [globals()[new_list[1][1]],new_list[1][1]]
    #list_c_peaks = globals()[new_list[2][1]]

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


#    globals()[new_list[0][1]] = list_a_peaks[0]
#    globals()[new_list[1][1]] = list_b_peaks[0]
#    globals()[new_list[2][1]] = list_c_list

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
            
    a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 4 -n > sorted_omc_cat.txt', shell=True)
    a.wait()
    
    a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 6 -n > sorted_inten_omc_cat.txt', shell=True)
    a.wait()

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
    
