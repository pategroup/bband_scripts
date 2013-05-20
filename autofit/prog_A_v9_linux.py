import subprocess
import os
from easygui import *
import sys
from multiprocessing import Process
import re
import numpy
import random
import string
""""
Python Triples Fitter

Written by Ian Finneran, Steve Shipman at NCF

based on previous work by the Pate Lab at UVa

Please comment any changes you make to the code here:










version 9 features:

-memory handling


-can use any number of processors

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
def int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="25.8", temp="298"):#generates SPCAT input file
    input_file = ""
    #print "freq_max=",freq
    input_file += "Molecule \n"
    input_file += "0  91  %s  %s  %s  %s  %s %s  %s\n"%(Q_rot, J_min, J_max,inten,inten,freq, temp)
    input_file += " 001  %s \n" % u_A
    input_file += " 002  %s \n" % u_B
    input_file += " 003  %s \n" % u_C
    fh_int = open("default.int", "w")
    fh_int.write(input_file)
    fh_int.close()
def var_writer(A,B,C,DJ,DJK,DK,dJ,dK):#generates SPCAT input file
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
    fh_var = open("default.var",'w')
    fh_var.write(input_file)
    fh_var.close()

def var_writer_uncert(A,B,C,DJ,DJK,DK,dJ,dK):#generates SPCAT input file
    input_file = ""
    dA = str(0.02*float(A))
    dB = str(0.01*float(B))
    dC = str(0.01*float(C))
    dDJ = str(0.1*float(DJ))
    dDJK = str(0.1*float(DJK))
    dDK = str(0.1*float(DK))
    ddJ = str(0.3*float(dJ))
    ddK = str(0.3*float(dK))
    input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
    input_file += "   8  430   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
    input_file +="a   1  1  0  99  0  1  1  1  1  -1   0\n"
    input_file += "           10000  %s %s \n" %(A,dA)
    input_file += "           20000  %s %s \n" %(B, dB)
    input_file += "           30000  %s %s \n" %(C, dC)
    input_file += "             200  %s %s \n" %(DJ, dDJ)
    input_file += "            1100  %s %s \n" %(DJK, dDJK) #need to actually check numbers: SPFIT doesn't read -- as a positive!
    input_file += "            2000  %s %s \n" %(DK, dDK)
    input_file += "           40100  %s %s \n" %(dJ, dJ)
    input_file += "           41000  %s %s \n" %(dK, dK)
    fh_var = open("default.var",'w')
    fh_var.write(input_file)
    fh_var.close()





def run_SPCAT(): 
    a = subprocess.Popen("./SPCAT default", stdout=subprocess.PIPE, shell=True)
    a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen
 
def cat_reader(freq_high,freq_low): #reads output from SPCAT
    fh = open("default.cat")
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
    
    
    
def trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low):
    pred_peaks = cat_reader(1000000, 0)
    #print "1",trans_1, trans_2, trans_3
    #print "2",pred_peaks 
    for peak in pred_peaks:
        if trans_1[2] == peak[2] and trans_1[3] == peak[3]:
            peak_1_freq = peak[1]

    for peak in pred_peaks:
        if trans_2[2] == peak[2] and trans_2[3] == peak[3]:
            peak_2_freq = peak[1]
        
    for peak in pred_peaks:
        if trans_3[2] == peak[2] and trans_3[3] == peak[3]:
            peak_3_freq = peak[1]
    return peak_1_freq,peak_2_freq,peak_3_freq




def dependence_test(A,B,C,DJ,DJK,DK,dJ,dK,trans_1,trans_2,trans_3,T,freq_high, freq_low,u_A,u_B,u_C):
    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temp=T)
    var_writer(A+(0.01*A),B,C,DJ,DJK,DK,dJ,dK)
    run_SPCAT()
    high_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low)


    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temp=T)
    var_writer(A-(0.01*A),B,C,DJ,DJK,DK,dJ,dK)
    run_SPCAT()
    low_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low)

    dv1A = (float(high_peak_freq[0])-float(low_peak_freq[0]))/(A+2-(A-2))
    dv2A = (float(high_peak_freq[1])-float(low_peak_freq[1]))/(A+2-(A-2))
    dv3A = (float(high_peak_freq[2])-float(low_peak_freq[2]))/(A+(2)-(A-(2)))
    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temp=T)
    var_writer(A,B+(2),C,DJ,DJK,DK,dJ,dK)
    run_SPCAT()
    high_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low)

    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temp=T)
    var_writer(A,B-(2),C,DJ,DJK,DK,dJ,dK)
    run_SPCAT()
    low_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low)

    dv1B = (float(high_peak_freq[0])-float(low_peak_freq[0]))/(B+(2)-(B-(2)))
    dv2B = (float(high_peak_freq[1])-float(low_peak_freq[1]))/(B+(2)-(B-(2)))
    dv3B = (float(high_peak_freq[2])-float(low_peak_freq[2]))/(B+(2)-(B-(2)))

    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temp=T)
    var_writer(A,B,C+(2),DJ,DJK,DK,dJ,dK)
    run_SPCAT()
    high_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low)

    int_writer(u_A,u_B,u_C, J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="100.0", temp=T)
    var_writer(A,B,C-(2),DJ,DJK,DK,dJ,dK)
    run_SPCAT()
    low_peak_freq = trans_freq_reader(trans_1, trans_2, trans_3,freq_high, freq_low)

    dv1C = (float(high_peak_freq[0])-float(low_peak_freq[0]))/(C+(2)-(C-(2)))
    dv2C = (float(high_peak_freq[1])-float(low_peak_freq[1]))/(C+(2)-(C-(2)))
    dv3C = (float(high_peak_freq[2])-float(low_peak_freq[2]))/(C+(2)-(C-(2)))


    matrix = numpy.array([(dv1A,dv1B,dv1C),(dv2A,dv2B,dv2C),(dv3A,dv3B,dv3C)])
    return numpy.linalg.det(matrix)


def fit_triples(list_a,list_b,list_c,trans_1,trans_2,trans_3,top_17,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,file_num,A,B,C,DJ,DJK,DK,dJ,dK):
    
    all_combo_file = "all_combo_list%s.txt"%(str(file_num)) 
    
    all_combo_list_file  = open(all_combo_file,"a")
    


    for freq_1,inten_1 in list_a:#generates all combinations of three peaks from three peaklists
        for freq_2, inten_2 in list_b:
            for freq_3, inten_3 in list_c:
                all_combo_list_file.write(str(freq_1)+","+str(inten_1)+","+str(freq_2)+","+str(inten_2)+","+str(freq_3)+","+str(inten_3)+", \n")
    all_combo_list_file.close()
    
    all_combo_list_file = open(all_combo_file)
    
    
    
    
    
    
        
    
    final_omc = []
    counter = 0
    output_file = ""
    test = 1
    regular_counter = 0
    error_counter = 0
    for all_combo_line in all_combo_list_file:
        all_combo_line = all_combo_line.split(",")
        peaks_triple= [(all_combo_line[0],all_combo_line[1]),(all_combo_line[2],all_combo_line[3]),(all_combo_line[4],all_combo_line[5])]
        
        
        
        
        input_file = ""#the next part adds in the three peaks to be fit
        input_file += trans_1[2][0:2]+' '+trans_1[2][2:4]+' '+trans_1[2][4:6]+' '+\
                      trans_1[3][0:2]+' '+trans_1[3][2:4]+' '+trans_1[3][4:6]+'                      '+peaks_triple[0][0]+' 0.50 1.0000\n' 
        input_file += trans_2[2][0:2]+' '+trans_2[2][2:4]+' '+trans_2[2][4:6]+' '+\
                      trans_2[3][0:2]+' '+trans_2[3][2:4]+' '+trans_2[3][4:6]+'                      '+peaks_triple[1][0]+' 0.50 1.0000\n' 
        input_file += trans_3[2][0:2]+' '+trans_3[2][2:4]+' '+trans_3[2][4:6]+' '+\
                      trans_3[3][0:2]+' '+trans_3[3][2:4]+' '+trans_3[3][4:6]+'                      '+peaks_triple[2][0]+' 0.50 1.0000\n'
        counter = 0
        for line in top_17:#the hack that adds in the 17 peaks but dosen't use them in the fit
            input_file += line[2][0:2]+' '+line[2][2:4]+' '+line[2][4:6]+' '+\
                      line[3][0:2]+' '+line[3][2:4]+' '+line[3][4:6]+'                      '+'%s.0'%(str(counter))+' 0.00 1.0000\n'
            counter += 1
        fh_lin = open("default%s.lin"%(str(file_num)), "w")
        fh_lin.write(input_file)
        fh_lin.close()        
        input_file = ""
        input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
        input_file += "   8  25   5    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n"
        input_file +="a   1  1  0  50  0  1  1  1  1  -1   0\n"
        input_file += "           10000  %s 1.0E+004 \n" % A
        input_file += "           20000  %s 1.0E+004 \n" % B
        input_file += "           30000  %s 1.0E+004 \n" % C
        input_file += "             200  %s 1.0E-025 \n" % DJ
        input_file += "            1100  %s 1.0E-025 \n" % DJK  #need to actually check numbers: SPFIT doesn't read -- as a positive!
        input_file += "            2000  %s 1.0E-025 \n" % DK
        input_file += "           40100  %s 1.0E-025 \n" % dJ
        input_file += "           41000  %s 1.0E-025 \n" % dK
        fh_par = open("default%s.par"%(str(file_num)),'w')
        fh_par.write(input_file)
        fh_par.close()
        a = subprocess.Popen("./SPFIT%s default%s"%(str(file_num),str(file_num)), stdout=subprocess.PIPE, shell=True)
        a.stdout.read()#used to let SPFIT finish
        fh_fit = open("default%s.fit"%(str(file_num)))
        const_list = []
        file_list = []
        for line in fh_fit:
                file_list.append(line)
                if line[13:18] == "10000":
                        const_list.append(line[37:53]) 
                if line[13:18] == "20000":
                        const_list.append(line[37:53])
                if line[13:18] == "30000":
                        const_list.append(line[37:53])
        freq_list = []
        for x in range(len(file_list)):
            if file_list[-x][11:14] == "RMS":
                rms_fit = float(file_list[-x][22:32]) #note - assumes RMS fit error is less than 1 GHz.  Change 22 to 21 if this is a problem.
            if file_list[-x][5:6] == ":" and int(file_list[-x][3:5])>3:
                freq_list.append(file_list[-x][60:71])
            if file_list[-x][40:64]=="EXP.FREQ.  -  CALC.FREQ.":
                break
        read_fit = (const_list[-3],const_list[-2], const_list[-1],freq_list)
        counter +=1
        constants = read_fit[0:3]
        freq_17 = read_fit[3]
        freq_17.reverse()
        A_1 = float(constants[0][:constants[0].find("(")]) #clean up the format of the rotational constants
        B_1 = float(constants[1][:constants[1].find("(")])
        C_1 = float(constants[2][:constants[2].find("(")])
        omc_list = []
        
        for x in range(len(top_17)): #matches peaks in the top 17 to peaks in experimental peak list <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            qnum_up = top_17[x][2]
            qnum_low = top_17[x][3]
            real_omc = 1.0
            freq = float(freq_17[x])
            p_1 = float(peak_list_1[0][0])
            p_2 = float(peak_list_1[-1][0])
            p_3 = float(peak_list_2[0][0])
            p_4 = float(peak_list_2[-1][0])
            p_5 = float(peak_list_3[0][0])
            p_6 = float(peak_list_3[-1][0])
            p_7 = float(peak_list_4[0][0])
            p_8 = float(peak_list_4[-1][0])
            if freq>p_1 and freq<p_2:#conditionals to find proper portion of experimental peak list to loop through
                peaks_section = peak_list_1
                peak = p_2
                regular_counter +=1
            elif freq>p_3 and freq<p_4:
                peaks_section = peak_list_2
                peak = p_4
                regular_counter +=1
            elif freq>p_5 and freq<p_6:
                peaks_section = peak_list_3
                peak = p_6
                regular_counter +=1
            elif freq>p_7 and freq<p_8:
                peaks_section = peak_list_4
                peak = p_8
                regular_counter +=1
            elif freq>p_8 or freq<p_1:
                peaks_section = peaklist
                peak = p_8
                error_counter +=1
                real_omc = 0.0# this is the omc if you throw out peaks that go over the edge of the spectrum
            else:
                peaks_section = peaklist
                peak = p_8
                regular_counter +=1
            current_peak = float(freq_17[x])
            old_omc = 100000.0
            for peak_freq,peak_inten in peaks_section: #find nearest peak in actual spectrum to the given top 20 peak
                omc = abs(current_peak-float(peak_freq))
                omc_low = abs(current_peak-float(peak))
                if omc>old_omc:
                    omc_low = old_omc                    
                    break
                old_omc = omc
            if real_omc == 1.0:
                real_omc = omc_low
            omc_list.append((omc_low, real_omc))# current_peak,qnum_up,qnum_low)) you can add in this extra output, but its slower
        omc_avg = [float(omc) for omc, real_omc in omc_list]
        real_omc_avg = [float(real_omc) for omc, real_omc in omc_list]
        score = str(len([omc for omc in omc_avg if omc<2.0])) #scores the accuracy of the fit
        avg = (sum(omc_avg)/len(omc_avg))+rms_fit
        real_avg = (sum(real_omc_avg)/len(real_omc_avg))+rms_fit
        if int(score)<10: #makes sorting work properly later
            score = '0'+score  
        output_file += 'score = '+' '+score+' '+"Const = "+str(A_1)+' '+str(B_1)+' '+str(C_1)+' '+"average omc = "+str(avg)+'  '+"avg w/out peaks over edge = "+str(real_avg)+"\n"
        if counter == 100000: #appends to file after every 100000 triples
            fh_final = open("final_output%s.txt"%(str(file_num)), "a")
            fh_final.write(output_file)
            fh_final.close()
            counter = 0
            output_file = ""
    fh_final = open("final_output%s.txt"%(str(file_num)), "a")#writes separate file for each processor
    print 'out of %s peaks there were %s peaks that werent in the experimental spectrum'%(regular_counter, error_counter) 
    fh_final.write(output_file)
    fh_final.close()
    os.system("sort -r 'final_output%s.txt'>sorted_final_out%s.txt"%(str(file_num),str(file_num)))#sorts output by score
    
if __name__ == '__main__': #multiprocessing imports script as module
    
  
    
    
    
    
    msg = "Enter all your constants in MHz, if left blank they will be set to default values (in parentheses). The inputs will be sent directly to SPCAT, so remember to input -DJ, -DJK, etc. "
    title = "Spectral Fitting Program"
    fieldNames = ["a dipole (0.0)","b dipole (0.0)","c dipole (0.0)","A (3000.0 MHz)","B (2000.0 MHz)","C (1000.0 MHz)","-DJ (0.0 MHz)","-DJK (0.0 MHz)","-DK (0.0 MHz)","-dJ (0.0 MHz)","-dK (0.0 MHz)",\
                  "number of processors (2)", "upper freqency limit (18000 MHz)", "lower frequency limit (8000 MHz)",\
                  "upper intensity limit on exp spectrum (100000)","lower intensity limit on exp spectrum(0)","temperature (298K)", "Jmax (30)"]
    x = subprocess.Popen("ls", stdout=subprocess.PIPE, shell=True)
    x = x.stdout.read().split()

    file_decision = input("Would you like to use an input file? (1/0)")
    
    
    
    if file_decision == 1:
        file_decision_flag = 1
       
    elif file_decision == 0:
        file_decision_flag = 0
    
    
    
    file_flag = 1
    
    while file_flag==1:
        filename = multenterbox("Choose a folder name for the calculation. If the folder already exists, this dialogue will reset.","",["Job Name"] )
        job_name = filename[0]
        marker1 = 0
        #if re.search("%s"%(str(job_name)), str(x))==None:
        for file1 in x:
            if job_name ==file1:
                marker1 = 1
        if marker1 ==0:
            file_flag =0
    

    
    a = subprocess.Popen("mkdir %s"%job_name,stdout=subprocess.PIPE, shell=True)
    a.wait()

    
    
    
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
                    temp = float(line.split()[1])
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
        const = [] 
        const = multenterbox(msg,title, fieldNames) # gui for grabbing data
         
        if const[0]=='':
            u_A = "0.0"#<<<<<<<<<<<<<set default value here
        else: u_A = const[0]
        if const[1]=='':
            u_B = "1.0"#<<<<<<<<<<<<<set default value here
        else: u_B = const[1] 
        if const[2]=='':
            u_C = "0.0"#<<<<<<<<<<<<<set default value here
        else: u_C = const[2]
        if const[3]=='':
            A = "5399.89204"#<<<<<<<<<<<<<set default value here
        else: A = const[3]
        if const[4]=='':
            B = "1143.24834"#<<<<<<<<<<<<<set default value here
        else: B = const[4]
        if const[5]=='':
            C = "1028.99051"#<<<<<<<<<<<<<set default value here, etc...
        else: C = const[5]
        if const[6]=='':
            DJ = "-0.3083E-003"
        else: DJ = const[6]
        if const[7]=='':
            DJK = "1.649E-003"
        else: DJK = const[7]
        if const[8]=='':
            DK = "-0.01445"
        else: DK = const[8]
        if const[9]=='':
            dJ = "-0.06492E-003"
        else: dJ = const[9]
        if const[10]=='':
            dK = "-1.164E-003"
        else: dK = const[10]
        if const[11] == '':
            processors = 4
        else: processors = int(const[11])
        if const[12]=='':
            freq_high = float(18000.0)
        else: freq_high = float(const[12])
        if const[13]=='':
            freq_low = float(6000.0)
        else: freq_low = float(const[13])
        if const[14]=='':
            inten_high = float('100000.0')
        else: inten_high = float(const[14])
        if const[15]=='':
            inten_low = float('0')
        else: inten_low = float(const[15])
        if const[16]=='':
            temp="2"
        else: temp = const[16]
        if const[17]=='':
            Jmax="15"
        else: Jmax = const[17]
    fh = open(fileopenbox(msg="enter the peaklist file in two column format: frequency intensity")) #gui for choosing peak list from experimental spectrum


    
    
    for number in range(processors):
        y = subprocess.Popen("cp SPFIT %s/SPFIT%s"%(job_name,number), stdout=subprocess.PIPE, shell=True)
        y.stdout.read()     
    y = subprocess.Popen("cp SPCAT %s/SPCAT"%(job_name), stdout=subprocess.PIPE, shell=True)
    y.stdout.read() 
    
    os.chdir(job_name)
    
    trans_1 = "" 
    trans_2 = ""
    trans_3 = ""
    int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temp=temp)
    var_writer_uncert(A,B,C,DJ,DJK,DK,dJ,dK)
    run_SPCAT()

    peaklist = []
    for line in fh:
        a = line.split()
        if float(a[1])>inten_low and float(a[1])<inten_high:
            peaklist.append(tuple(a))
    trans_1_peaks = []
    trans_2_peaks = []
    trans_3_peaks = []
    dependence_flag =1
    new_peaks_flag=0
    while trans_1_peaks == [] or trans_2_peaks == [] or trans_3_peaks == [] and dependence_flag == 1: #this loops until there are peaks around each member of the triple
        if fit_peaks_list == [] or new_peaks_flag==1:
            msg ="Pick any three transitions"
            title = "Microwave Fitting Program"
            choices = cat_reader(freq_high, freq_low)
            choice = multchoicebox(msg, title, choices)
            choice_clean = []
            for entry in choice:
                clean = entry[2:43]
                re_split = clean.split("', '")
                tuples = tuple(re_split)
                choice_clean.append(tuples)
            trans_1 = choice_clean[0]
            trans_2 = choice_clean[1]
            trans_3 = choice_clean[2]
        else:
            trans_1 = fit_peaks_list[0]
            trans_2 = fit_peaks_list[1]
            trans_3 = fit_peaks_list[2]
        dependence = dependence_test(float(A),float(B),float(C),float(DJ),float(DJK),float(DK),float(dJ),float(dK),trans_1,trans_2,trans_3,temp,freq_high, freq_low,u_A,u_B,u_C)
        
        decision1 = buttonbox(msg='These three peaks have a linear dependence of %s. A linear dependence between -1.0 and 1.0 is probably too low. Would you like to continue, or try three new peaks?'%(str(dependence)), choices=('Continue','try three new peaks'))
        
        
        
        
        if decision1 == 'Continue':
            dependence_flag = 0
        else:
            new_peaks_flag = 1
            continue
        
        uncertainty_flag =1
        while uncertainty_flag ==1:
            
            if freq_uncertainty==0.0:
                freq_uncertainty = float(enterbox(msg="Enter the frequency uncertainty in MHz"))
            
            
            
            
            trans_1_peaks = []
            trans_2_peaks = []
            trans_3_peaks = []
            for freq_p, inten_p in peaklist:
                if abs(float(trans_1[1])-float(freq_p))< freq_uncertainty:
                    trans_1_peaks.append((freq_p, inten_p))
            for freq_p, inten_p in peaklist:
                if abs(float(trans_2[1])-float(freq_p))< freq_uncertainty: #this bit finds peaks in the real spectrum that are near the predicted peaks
                    trans_2_peaks.append((freq_p, inten_p))
            for freq_p, inten_p in peaklist:
                if abs(float(trans_3[1])-float(freq_p))< freq_uncertainty:
                    trans_3_peaks.append((freq_p, inten_p))
            num_of_triples = len(trans_1_peaks)*len(trans_2_peaks)*len(trans_3_peaks) #this tells you how many entries there will be in the all_combo_list
            
            
            
            decision = buttonbox(msg='There are %s triples in this calculation. Would you like to continue, try new uncertainty, or quit?'%(str(num_of_triples)), choices=('Continue','Quit','new uncertainty'))
            
            if decision == 'Quit':
                quit()
            elif decision == 'Continue':
                uncertainty_flag = 0
            else:
                freq_uncertainty = 0.0
                pass
        
        if check_peaks_list == []:
            if trans_1_peaks == []:
                msgbox("There are no peaks in the peaklist file near this transition: %s"%(str(trans_1)))#error messages
            elif trans_2_peaks == []:
                msgbox("There are no peaks in the peaklist file near this transition: %s"%(str(trans_2)))
            elif trans_3_peaks == []:
                msgbox("There are no peaks in the peaklist file near this transition: %s"%(str(trans_3)))
        
                        
            int_writer(u_A,u_B,u_C, J_max=Jmax,freq=str((freq_high*.001)), temp=temp)
            var_writer_uncert(A,B,C,DJ,DJK,DK,dJ,dK)
            run_SPCAT()
        
            choice_method = buttonbox(msg='choose the method for scoring the results:', title='Scoring Options',choices=\
                                      ('use the top 10 most intense peaks (default)','choose an arbitrary set of peaks from the prediction'))
            if choice_method == 'use the top 10 most intense peaks (default)':
                top_peaks = cat_reader(freq_high, freq_low)[0:10] #grab the most intense peaks from the predicted spectrum
            else:
                top_peaks = multchoicebox(msg='Pick the transitions that you want to use for sorting the results', title='Result Sorting Setup',choices=cat_reader(freq_high, freq_low))
                top_peaks_clean = []
                for entry in top_peaks:
                    clean = entry[2:43]
                    re_split = clean.split("', '")
                    tuples = tuple(re_split)
                    top_peaks_clean.append(tuples)
                top_peaks = top_peaks_clean
        else:
            top_peaks = check_peaks_list
    
    top_peaks_3cut = []
    for entry in top_peaks:
        if entry != trans_1 and entry != trans_2 and entry != trans_3:
            top_peaks_3cut.append(entry)

    num_of_triples = len(trans_1_peaks)*len(trans_2_peaks)*len(trans_3_peaks) #this tells you how many entries there will be in the all_combo_list
    print num_of_triples

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
    str(freq_low),str(inten_high),str(inten_low),str(temp),str(Jmax),str(freq_uncertainty),str(num_of_triples),fitting_peaks_str,str(trans_1),str(trans_2),str(trans_3))
    Job_fh = open("input_data_%s.txt"%(job_name),"w")
    Job_fh.write(job_file) 
    Job_fh.close()
    all_combo_list = []                                             #which determines how long it will take to execute
    counter = 0
    current_list = "all_combo_list_0.txt"
    
    
    
    #for freq_1,inten_1 in trans_1_peaks:#generates all combinations of three peaks from three peaklists
    #    for freq_2, inten_2 in trans_2_peaks:
    #        for freq_3, inten_3 in trans_3_peaks:
    #            all_combo_list.append(((freq_1,inten_1),(freq_2, inten_2),(freq_3, inten_3)))
    
    

    

    new_list = []

    new_list = [(len(trans_1_peaks),"trans_1_peaks"),(len(trans_2_peaks),"trans_2_peaks"),(len(trans_3_peaks),"trans_3_peaks")]
    new_list.sort()
    
    
    list_key = []

    list_a_peaks = [vars()[new_list[0][1]],new_list[0][1]]
    
    
    
    list_b_peaks = [vars()[new_list[1][1]],new_list[1][1]]
    
    
    
    list_c_peaks = vars()[new_list[2][1]]





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
    
    
        
    
    peak_list_1 = peaklist[0:int(len(peaklist)/4)]#splits peaks into 4 parts to speed up processing
    peak_list_2 = peaklist[int(len(peaklist)/4):int(len(peaklist)/2)]
    peak_list_3 = peaklist[int(len(peaklist)/2):int(len(peaklist)*0.75)]
    peak_list_4 = peaklist[int(len(peaklist)*0.75):len(peaklist)]
    


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
        

        vars()["p%s"%str(num)] = Process(target=fit_triples, args=(trans_x_peaks,trans_y_peaks,trans_z_peaks,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,num,A,B,C,DJ,DJK,DK,dJ,dK))

    for num in range(processors):
        vars()["p%s"%str(num)].start()
    for num in range(processors):
        vars()["p%s"%str(num)].join()
        
    a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 4 -n > sorted_omc_cat.txt', shell=True)
    a.wait()
    
    
    

 
    """
        
    if processors == 2:
        all_combo_list_1 = all_combo_list[0:int(len(all_combo_list)/2)]#splits all combination list into 2 parts for each processor
        all_combo_list_2 = all_combo_list[int(len(all_combo_list)/2):int(len(all_combo_list))]
        p0 = Process(target=fit_triples, args=(all_combo_list_1,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"0",A,B,C,DJ,DJK,DK,dJ,dK))
        p1 = Process(target=fit_triples, args=(all_combo_list_2,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"1",A,B,C,DJ,DJK,DK,dJ,dK))
        p0.start()
        p1.start()
        p0.join()
        p1.join()
    if processors ==1:
        fit_triples(all_combo_list,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"",A,B,C,DJ,DJK,DK,dJ,dK)

    if processors == 4:
        all_combo_list_1 = all_combo_list[0:int(len(all_combo_list)/4)]#splits all combination list into 4 parts for each processor
        all_combo_list_2 = all_combo_list[int(len(all_combo_list)/4):int(len(all_combo_list)/2)]
        all_combo_list_3 = all_combo_list[int(len(all_combo_list)/2):int(0.75*len(all_combo_list))]
        all_combo_list_4 = all_combo_list[int(len(all_combo_list)*0.75):len(all_combo_list)]
        p0 = Process(target=fit_triples, args=(all_combo_list_1,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"0",A,B,C,DJ,DJK,DK,dJ,dK))
        p1 = Process(target=fit_triples, args=(all_combo_list_2,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"1",A,B,C,DJ,DJK,DK,dJ,dK))
        p2 = Process(target=fit_triples, args=(all_combo_list_3,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"2",A,B,C,DJ,DJK,DK,dJ,dK))
        p3 = Process(target=fit_triples, args=(all_combo_list_4,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"3",A,B,C,DJ,DJK,DK,dJ,dK))
        p0.start()
        p1.start()
        p2.start()
        p3.start()        
        p0.join()
        p1.join()
        p2.join()
        p3.join()
        
        

   


    if processors == 8:
        all_combo_list_1 = all_combo_list[0:int(len(all_combo_list)/8)]#splits all combination list into 4 parts for each processor
        all_combo_list_2 = all_combo_list[int(len(all_combo_list)/8):int(len(all_combo_list)/4)]
        all_combo_list_3 = all_combo_list[int(len(all_combo_list)/4):int((3*len(all_combo_list))/8)]
        all_combo_list_4 = all_combo_list[int((3*len(all_combo_list))/8):int(0.5*len(all_combo_list))]
        all_combo_list_5 = all_combo_list[int(len(all_combo_list)*0.5):int((5*len(all_combo_list))/8)]
        all_combo_list_6 = all_combo_list[int((5*len(all_combo_list))/8):int((6*len(all_combo_list))/8)]#splits all combination list into 4 parts for each processor
        all_combo_list_7 = all_combo_list[int((6*len(all_combo_list))/8):int((7*len(all_combo_list))/8)]
        all_combo_list_8 = all_combo_list[int((7*len(all_combo_list))/8):int(len(all_combo_list))]
       
        p0 = Process(target=fit_triples, args=(all_combo_list_1,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"0",A,B,C,DJ,DJK,DK,dJ,dK))
        p1 = Process(target=fit_triples, args=(all_combo_list_2,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"1",A,B,C,DJ,DJK,DK,dJ,dK))
        p2 = Process(target=fit_triples, args=(all_combo_list_3,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"2",A,B,C,DJ,DJK,DK,dJ,dK))
        p3 = Process(target=fit_triples, args=(all_combo_list_4,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"3",A,B,C,DJ,DJK,DK,dJ,dK))
        p4 = Process(target=fit_triples, args=(all_combo_list_5,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"4",A,B,C,DJ,DJK,DK,dJ,dK))
        p5 = Process(target=fit_triples, args=(all_combo_list_6,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"5",A,B,C,DJ,DJK,DK,dJ,dK))
        p6 = Process(target=fit_triples, args=(all_combo_list_7,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"6",A,B,C,DJ,DJK,DK,dJ,dK))
        p7 = Process(target=fit_triples, args=(all_combo_list_8,trans_1,trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,"7",A,B,C,DJ,DJK,DK,dJ,dK))



        p0.start()
        p1.start()
        p2.start()
        p3.start()        
        p4.start()
        p5.start()
        p6.start()
        p7.start()        

        p0.join()
        p1.join()
        p2.join()
        p3.join()
        p4.join()
        p5.join()
        p6.join()
        p7.join()
    """























    """
    entry = 0
    while entry<processors:
        variable_name = "p"+str(entry)
        all_combo_list_part = all_combo_list[int(len(all_combo_list)*(entry/processors)):int(len(all_combo_list)*((entry+1)/processors))]
        vars()[variable_name] = Process(target=fit_triples, args=(all_combo_list_part,trans_1,\
                                                           trans_2,trans_3,top_peaks_3cut,peak_list_1, peak_list_2, peak_list_3, peak_list_4,peaklist,entry,A,B,C,DJ,DJK,DK,dJ,dK))
        vars()[variable_name].start()
        entry+=1
    """


