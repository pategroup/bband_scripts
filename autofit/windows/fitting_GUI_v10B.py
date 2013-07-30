
import sys, os, random

import fitting_GUI_B_v10
import autofit_NS_module
from time import sleep
import subprocess
from multiprocessing import Process
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
import numpy
import math
from scipy.interpolate import *
import os

class AppForm(QMainWindow):
############################<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
############################
############################This next section initializes parameters in the GUI
############################    
############################<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    def __init__(self, parent=None):
        
	QMainWindow.__init__(self, parent)
	self.replot_flag = 0
	self.setWindowTitle('Autofit GUI')
        self.triples_flag = 0
	self.create_menu()
	self.create_main_frame()
	self.create_status_bar()

	self.a_box.setText('3000')
	self.b_box.setText('2000')
	self.c_box.setText('1000')
	self.da_box.setText('1000')
	self.db_box.setText('1000')
	self.dc_box.setText('1000')
	self.dua_box.setText('1')
	self.dub_box.setText('1')
	self.duc_box.setText('1')
	self.J_max_box.setText('20')
	self.DJ_box.setText('0.0')
	self.DJK_box.setText('0.0')
	self.DK_box.setText('0.0')
	self.dJ_box.setText('0.0')
	self.dK_box.setText('0.0')
	self.ua_box.setText('1.0')
	self.ub_box.setText('1.0')
	self.uc_box.setText('1.0')
	self.T_box.setText('5.0')
	self.f_lower_box.setText('0.0')
	self.f_upper_box.setText('30000.0')
	self.I_lower_box.setText('0.01')
	self.I_upper_box.setText('100000.0')
	self.proc_box.setText('8')
	self.skip_box.setText('1')
	self.job_box.setText('job_name')
	self.uncertainty_box_1.setText('100.0')
	self.uncertainty_box_2.setText('100.0')
	self.uncertainty_box_3.setText('100.0')
        self.distortion_cb.toggle()
        self.distortion_cb.toggle()
        self.spectrum_cb.toggle()
        self.spectrum_cb.toggle()
        self.autofit_cb.toggle()
        self.autofit_cb.toggle()
        self.fix_cb.toggle()
        self.fix_cb.toggle()
        self.plot_triples_cb.toggle()
        self.plot_triples_cb.toggle()
	self.a_set  = '3000'
	self.b_set  = '2000'
	self.c_set  = '1000'
	self.autofit_button.setEnabled(False)    
	self.peak_finder_button.setEnabled(False)
	self.initialize()

        
############################<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
############################
############################This section defines all of the relevant functions in the GUI window
############################ Since they are imbedded in a class, they must have the 'self' input
############################    
############################<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



    def var_writer(self,A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag):#generates SPCAT input file
    
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
    def run_SPCAT(self): 
        a = subprocess.Popen("SPCAT default", stdout=subprocess.PIPE, shell=False)
        a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen 

	 
    def cat_reader(self,freq_high,freq_low,flag="default"): #reads output from SPCAT
        
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
    def int_writer(self,u_A,u_B,u_C, J_min="00", J_max='20', inten="-10.0",Q_rot="300000",freq="25.8", temperature="298", flag="default"):#generates SPCAT input file
        input_file = ""
        J_max=int(self.J_max_box.text())
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
  	 
    def check_bounds(self,trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high):
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
 
    def final_uncerts(self,trans_1_center,trans_2_center,trans_3_center,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,u_A,u_B,u_C,main_flow,trans_1,trans_2,trans_3):
        
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
    
    def auto_set_button_push(self):

	self.status_text.setText("OPTIMIZING SELECTED PEAKS FOR AUTOFIT...PLEASE WAIT...")
        self.status_text.setStyleSheet('color: red')
        QApplication.processEvents()
	A = float(unicode(self.a_box.text()))
        B = float(unicode(self.b_box.text()))
        C = float(unicode(self.c_box.text()))
        DJ = float(unicode(self.DJ_box.text()))
        DJK = float(unicode(self.DJK_box.text()))      
        DK = float(unicode(self.DK_box.text()))        
        dJ = float(unicode(self.dJ_box.text()))
        dK = float(unicode(self.dK_box.text()))  
        u_A = float(unicode(self.ua_box.text()))        
        u_B = float(unicode(self.ub_box.text()))
        u_C = float(unicode(self.uc_box.text()))
        freq_high=float(unicode(self.f_upper_box.text()))
        freq_low=float(unicode(self.f_lower_box.text()))
        main_flow = 'Normal species'
	temperature = float(unicode(self.T_box.text()))
	flag = 'Automatic scoring'
	peak_1_uncertainty = float(unicode(self.uncertainty_box_1.text()))
        peak_2_uncertainty = float(unicode(self.uncertainty_box_2.text()))
        peak_3_uncertainty = float(unicode(self.uncertainty_box_3.text()))
	#print fitting_GUI_B_v10.picked_list
	full_list = [(entry[3],entry[0],entry[1],entry[2],entry[4]) for entry in fitting_GUI_B_v10.picked_list] 
	if len(full_list)<3:
	    msg = "You need to choose at least three transitions before doing an automatic setup."
	    
            QMessageBox.about(self, "Helpful hint", msg.strip())
	else:
	     

            
            results = self.triple_selection(full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag)
            
            results_list = []
            uncertainty_list = []
            for line in results:
                #print line[12],line[13],line[14]
                #print line[12],line[13],line[14],uncertainty
                if 0==self.check_bounds(float(line[3]),float(line[6]),float(line[9]),peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,freq_low,freq_high)[0]:
                    results_list.append((line[1],(line[3],line[4],line[5]),(line[6],line[7],line[8]),(line[9],line[10],line[11])))
            
          
            #print len(results_list),uncertainty,
            
            if len(results_list)==0:
              
                msg = "No combinations can work with such a large uncertainty."
	    
                QMessageBox.about(self, "Helpful hint", msg.strip())
            results_list.sort()
            best_option = results_list[-1]

            counter = 0
            self.fit_list = []
            for entry in fitting_GUI_B_v10.picked_list:
                if entry[1]==best_option[1][1] and entry[2]==best_option[1][2]:
                    #print entry[1],best_option[1][1],entry[2],best_option[1][2]
                    fitting_GUI_B_v10.picked_list[counter] = (entry[0],entry[1],entry[2],entry[3],entry[4],'FIT')
                    self.fit_list.append((entry[0],entry[1],entry[2],entry[3],entry[4]))
                elif entry[1]==best_option[2][1] and entry[2]==best_option[2][2]:
                    fitting_GUI_B_v10.picked_list[counter] = (entry[0],entry[1],entry[2],entry[3],entry[4],'FIT')
                    self.fit_list.append((entry[0],entry[1],entry[2],entry[3],entry[4]))
                    #print entry[1],best_option[2][1],entry[2],best_option[2][2]                  
                elif entry[1]==best_option[3][1] and entry[2]==best_option[3][2]:
                    fitting_GUI_B_v10.picked_list[counter] = (entry[0],entry[1],entry[2],entry[3],entry[4],'FIT')
                    #print entry[1],best_option[3][1],entry[2],best_option[3][2]
                    self.fit_list.append((entry[0],entry[1],entry[2],entry[3],entry[4]))
                else: fitting_GUI_B_v10.picked_list[counter] = (entry[0],entry[1],entry[2],entry[3],entry[4],'CHECK')
                counter+=1                  
            self.peaks_button_push()
            window_decision='User-defined, same for each'

            trans_1 = self.fit_list[0]
            trans_2 = self.fit_list[1]
            trans_3 = self.fit_list[2]
            print trans_1,trans_2,trans_3
            #fitting_GUI_B_v10.triples_plt.set_data([],[])
            try:
                peaklist = self.peaklist
                
                (self.trans_1_peaks,self.trans_2_peaks,self.trans_3_peaks,num_of_triples ) = self.triples_gen(peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,trans_1,trans_2,trans_3)
                if num_of_triples!=0:
                    
                    self.num_of_triples_label.setText("%s Triples"%str(num_of_triples))
                    self.num_of_triples_label.setStyleSheet('color: blue')
                else:
                    self.num_of_triples_label.setText("0 Triples")
                            
                    self.num_of_triples_label.setStyleSheet('color: red') 
                    self.num_of_triples = 0
                
            except AttributeError:
                msg = "A list of experimental peaks is needed to generate a list of triples."
	    
                QMessageBox.about(self, "Helpful hint", msg.strip())
            self.plot_triples()
            self.status_text.setText("AUTO SETUP ROUTINE FINISHED")
            self.status_text.setStyleSheet('color: blue')
            self.update_fit_check_status()
    def update_triples_meter(self):
        
        self.fit_list = []
        counter =0
        for entry in fitting_GUI_B_v10.picked_list:
            if entry[5]=='FIT' and counter<3:
                self.fit_list.append((entry[0],entry[1],entry[2],entry[3],entry[4]))
                print "fit here"
                counter+=1 
        freq_high=float(unicode(self.f_upper_box.text()))
        freq_low=float(unicode(self.f_lower_box.text()))

        
        try:
            peaklist = self.peaklist
            
            trans_1 = self.fit_list[0]
            trans_2 = self.fit_list[1]
            trans_3 = self.fit_list[2]
            print 'asd'
            peak_1_uncertainty = float(unicode(self.uncertainty_box_1.text()))
            peak_2_uncertainty = float(unicode(self.uncertainty_box_2.text()))
            peak_3_uncertainty = float(unicode(self.uncertainty_box_3.text()))
            (self.trans_1_peaks,self.trans_2_peaks,self.trans_3_peaks,num_of_triples ) = self.triples_gen(peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,trans_1,trans_2,trans_3)
            self.num_of_triples_label.setText("%s Triples"%str(num_of_triples))
            if num_of_triples!=0:
                
                self.num_of_triples_label.setText("%s Triples"%str(num_of_triples))
                self.num_of_triples_label.setStyleSheet('color: blue')
                self.triples_flag = 1
            else:
                self.num_of_triples_label.setText("0 Triples")
                        
                self.num_of_triples_label.setStyleSheet('color: red') 
                self.num_of_triples = 0
                self.triples_flag = 0
        except AttributeError:
            self.num_of_triples_label.setText("0 Triples")
                    
            self.num_of_triples_label.setStyleSheet('color: red') 
            self.num_of_triples = 0
            self.triples_flag = 0
        except IndexError:
            self.num_of_triples_label.setText("0 Triples")
                    
            self.num_of_triples_label.setStyleSheet('color: red') 
            self.num_of_triples = 0
            self.triples_flag = 0
        self.update_fit_check_status()
    def triples_gen(self,peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,trans_1,trans_2,trans_3):
    

        trans_1_peaks = []
        trans_2_peaks = []
        trans_3_peaks = []
        trans_1_center = trans_1[0]
        print trans_1_center,'trans_1_center' 
        trans_2_center = trans_2[0]
        trans_3_center = trans_3[0]
        for freq_p, inten_p in peaklist:
            print abs(float(trans_1_center)-float(freq_p)),peak_1_uncertainty
            if abs(float(trans_1_center)-float(freq_p))< peak_1_uncertainty:
                trans_1_peaks.append((freq_p, inten_p))
                print 'akdhkashd'
            if abs(float(trans_2_center)-float(freq_p))< peak_2_uncertainty: #this bit finds peaks in the real spectrum that are near the predicted peaks
                print 'akdhkashd1'
                trans_2_peaks.append((freq_p, inten_p))
            if abs(float(trans_3_center)-float(freq_p))< peak_3_uncertainty:
                print 'akdhkashd2'
                trans_3_peaks.append((freq_p, inten_p))
        num_of_triples = len(trans_1_peaks)*len(trans_2_peaks)*len(trans_3_peaks) #this tells you how many entries there will be in the all_combo_list
        return trans_1_peaks,trans_2_peaks,trans_3_peaks,num_of_triples                                             
                                                                            
    def triple_selection(self,full_list,A,B,C,DJ,DJK,DK,dJ,dK,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow,flag):
        #print full_list
        triple_style = flag
        #if flag == "undecided":
        #    triple_style = buttonbox(msg='How do you want to choose fitting transitions?  You can choose from an automatically scored list, or select three transitions manually.', choices=('Automatic scoring','Manual selection'))
    
        #elif flag == "force manual":
        #    triple_style = 'Manual selection'
        
        
        if triple_style == 'Automatic scoring':
            #total_check_num = 10 # This is the number of peaks used to generate possible triples, ordered by intensity. 10 = 120 possibilities, 15 = 455 possibilities.
            if len(full_list)>10:
                total_check_num=10
            else:
                total_check_num = len(full_list)
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
                        dependence = abs(self.dependence_test(float(A),float(B),float(C),float(DJ),float(DJK),float(DK),float(dJ),float(dK),trans_1,trans_2,trans_3,temperature,freq_high,freq_low,u_A,u_B,u_C,main_flow))
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
    
            #msg = "Choose a triples combination.  The second number is the triples score, ranging from 0 to 100 with 100 being the best possible.  The third number is the highest uncertainty in MHz from the triple."
            #title = "Microwave Fitting Program"
            #choice = choicebox(msg,title,triples_choice_list)
            #clean_choice = choice[2:-1].split(",")
    
            #highest_uncert = float(clean_choice[2])
            #trans_1 = full_list[int(clean_choice[-3])]
            #trans_2 = full_list[int(clean_choice[-2])]
            #trans_3 = full_list[int(clean_choice[-1])]
    
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
    
                    dependence = self.dependence_test(float(A),float(B),float(C),float(DJ),float(DJK),float(DK),float(dJ),float(dK),trans_1,trans_2,trans_3,temperature,freq_high, freq_low,u_A,u_B,u_C,main_flow)
            
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
        
        #return highest_uncert,trans_1,trans_2,trans_3
        return triples_choice_list
        
            
    def trans_freq_reader(self,trans_1, trans_2, trans_3):
    
        peak_1_freq = 0
        peak_2_freq = 0
        peak_3_freq = 0
    
        pred_peaks = self.cat_reader(1000000, 0, flag="default")
        for peak in pred_peaks:
            if trans_1[2] == peak[2] and trans_1[3] == peak[3]:
                peak_1_freq = peak[1]
            if trans_2[2] == peak[2] and trans_2[3] == peak[3]:
                peak_2_freq = peak[1]
            if trans_3[2] == peak[2] and trans_3[3] == peak[3]:
                peak_3_freq = peak[1]
        return peak_1_freq,peak_2_freq,peak_3_freq
    def dependence_test(self,A,B,C,DJ,DJK,DK,dJ,dK,trans_1,trans_2,trans_3,T,freq_high, freq_low,u_A,u_B,u_C,main_flow):
        self.int_writer(u_A,u_B,u_C, J_min="00", J_max=int(self.J_max_box.text()), inten="-10.0",Q_rot="300000",freq="100.0", temperature=T, flag = "default")
    
        self.var_writer(A+(2),B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        self.run_SPCAT()
        high_peak_freq = self.trans_freq_reader(trans_1, trans_2, trans_3)
    
        self.var_writer(A-(2),B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        self.run_SPCAT()
        low_peak_freq = self.trans_freq_reader(trans_1, trans_2, trans_3)
    
        dv1A = (float(high_peak_freq[0])-float(low_peak_freq[0]))/4
        dv2A = (float(high_peak_freq[1])-float(low_peak_freq[1]))/4
        dv3A = (float(high_peak_freq[2])-float(low_peak_freq[2]))/4
    
        self.var_writer(A,B+(2),C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        self.run_SPCAT()
        high_peak_freq = self.trans_freq_reader(trans_1, trans_2, trans_3)
    
        self.var_writer(A,B-(2),C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        self.run_SPCAT()
        low_peak_freq = self.trans_freq_reader(trans_1, trans_2, trans_3)
    
        dv1B = (float(high_peak_freq[0])-float(low_peak_freq[0]))/4
        dv2B = (float(high_peak_freq[1])-float(low_peak_freq[1]))/4
        dv3B = (float(high_peak_freq[2])-float(low_peak_freq[2]))/4
    
        self.var_writer(A,B,C+(2),DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        self.run_SPCAT()
        high_peak_freq = self.trans_freq_reader(trans_1, trans_2, trans_3)
    
        self.var_writer(A,B,C-(2),DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")
        self.run_SPCAT()
        low_peak_freq = self.trans_freq_reader(trans_1, trans_2, trans_3)
    
        dv1C = (float(high_peak_freq[0])-float(low_peak_freq[0]))/4
        dv2C = (float(high_peak_freq[1])-float(low_peak_freq[1]))/4
        dv3C = (float(high_peak_freq[2])-float(low_peak_freq[2]))/4
    
        self.var_writer(A,B,C,DJ,DJK,DK,dJ,dK,main_flow,flag="uncert")	# This re-runs SPCAT at the initial constants so that other things that read from default.cat are correct after this function is executed.
        self.run_SPCAT()
    
        matrix = numpy.array([(dv1A,dv1B,dv1C),(dv2A,dv2B,dv2C),(dv3A,dv3B,dv3C)])
        return numpy.linalg.det(matrix)
                                                                             
                                    
    def selectFile_spectrum(self):
        
        self.spectrum_filename = QFileDialog.getOpenFileName()

            
        self.update_exp_spectrum()
        self.peak_finder_button.setEnabled(True)
            
    def update_fit_check_status(self):
        check_counter = 0
        fit_counter = 0
         
        check_flag = 0
        fit_flag = 0
        
        
        
        for line in fitting_GUI_B_v10.picked_list:
            if line[5]=='CHECK':
                check_counter+=1
            if line[5]=='FIT':
                fit_counter+=1            
        
        if check_counter!=0:
            self.check_status_label2.setText(str(check_counter)+" peaks")
            self.check_status_label2.setStyleSheet('color: blue')
            check_flag = 1
        if check_counter==0:
            self.check_status_label2.setText("None Selected") 
            self.check_status_label2.setStyleSheet('color: red')           
        if fit_counter==3:
            self.fitting_status_label2.setText("3 peaks")        
            self.fitting_status_label2.setStyleSheet('color: blue') 
            fit_flag = 1          
        if fit_counter!=3:
            self.fitting_status_label2.setText("%s peaks"%str(fit_counter))        
            self.fitting_status_label2.setStyleSheet('color: red') 
        if self.triples_flag == 1 and check_flag == 1 and fit_flag == 1:
            self.autofit_button.setEnabled(True)  
        else:
            self.autofit_button.setEnabled(False)

    def plot_triples(self):
    
        if self.plot_triples_cb.isChecked()==True:
            fitting_GUI_B_v10.triples_flag = 1
            fitting_GUI_B_v10.ax2.set_ylim([0,self.peaklist.max(0)[1]])    
            freq_high=float(unicode(self.f_upper_box.text()))
            freq_low=float(unicode(self.f_lower_box.text()))
            peak_1_uncertainty = float(unicode(self.uncertainty_box_1.text()))
            peak_2_uncertainty = float(unicode(self.uncertainty_box_2.text()))
            peak_3_uncertainty = float(unicode(self.uncertainty_box_3.text()))
            fitting_GUI_B_v10.peak_1_uncertainty = peak_1_uncertainty
            fitting_GUI_B_v10.peak_2_uncertainty = peak_2_uncertainty
            fitting_GUI_B_v10.peak_3_uncertainty = peak_3_uncertainty
            fitting_GUI_B_v10.peak_list = self.peaklist
            fitting_GUI_B_v10.freq_low = freq_low
            fitting_GUI_B_v10.freq_high = freq_high
            
            
            f_list = []
            i_list = []
            qnum_list = []
            
            for line in fitting_GUI_B_v10.picked_list:
                 
                if line[5]=='FIT':
                    print line
                    f_list.append(float(line[0]))
                    i_list.append(10**float(line[3]))
                    qnum_list.append((line[1],line[2]))
                    
                    
            #fitting_GUI_B_v10.triples_plt.set_data(f_list,i_list)
            #fitting_GUI_B_v10.plt.show()

            #fitting_GUI_B_v10.plt.show()
            try:
                fitting_GUI_B_v10.trans_1 = qnum_list[0]
            except:
                pass
            try:
                fitting_GUI_B_v10.trans_2 = qnum_list[1]
            except: pass
            try:
                fitting_GUI_B_v10.trans_3 = qnum_list[2]
            except:
                pass
            self.update_triples_labels()
            self.update_triples_meter()
            
            fitting_GUI_B_v10.update_plot()
            I_1 = []
            F_1 = []
            I_2 = []
            F_2 = []
            I_3 = []
            F_3 = []
            for line in self.trans_1_peaks:
                
                F_1.append(line[0])
                I_1.append(line[1])
            for line in self.trans_2_peaks:
                F_2.append(line[0])
                I_2.append(line[1])
            for line in self.trans_3_peaks:
                F_3.append(line[0])
                I_3.append(line[1])
            #fitting_GUI_B_v10.trans_1_plt.set_data(F_1,I_1)
            #fitting_GUI_B_v10.trans_2_plt.set_data(F_2,I_2)
            #fitting_GUI_B_v10.trans_3_plt.set_data(F_3,I_3)
            #self.peak_center_cb.toggle()
            #self.peak_center_cb.toggle()

            
            
            
        if self.plot_triples_cb.isChecked()==False:

            fitting_GUI_B_v10.triples_flag = 0   
            fitting_GUI_B_v10.trans_1_plt.set_data([],[])
            fitting_GUI_B_v10.trans_2_plt.set_data([],[])
            fitting_GUI_B_v10.trans_3_plt.set_data([],[])       
            fitting_GUI_B_v10.triples_plt.set_data([],[])
            fitting_GUI_B_v10.plt.show()

            
    def update_exp_spectrum(self):

               
        self.status_text.setText("PLOTTING EXP SPECTRUM...PLEASE WAIT...")
        self.status_text.setStyleSheet('color: red')
        QApplication.processEvents()
        number_of_skips = int(unicode(self.skip_box.text()))
        spectrum_list_F = []
        spectrum_list_I = []
        new_list = []
        try:
            if self.replot_flag==1:
                fh = open(self.spectrum_filename)
                counter =0
                for line in fh:
                    counter+=1
                    if counter%(number_of_skips+1)==0:
                        
                        spectrum_list_F.append(line.split()[0])
                        spectrum_list_I.append(line.split()[1])
                        new_list.append(float(line.split()[1]))
                fitting_GUI_B_v10.exp_plt.set_data(spectrum_list_F,spectrum_list_I)
                
                new_list.sort()
                fitting_GUI_B_v10.ax2.set_ylim([float(new_list[0]),float(new_list[-1])])
                self.status_text.setText("EXP SPECTRUM LOADED")
                self.status_text.setStyleSheet('color: blue')
                self.exp_spec_status_label2.setText(str(self.spectrum_filename))
                self.status_text.setStyleSheet('color: blue')
                QApplication.processEvents()
                fitting_GUI_B_v10.plt.show()
                fh.close()
            self.status_text.setText("EXP SPECTRUM LOADED")
            self.status_text.setStyleSheet('color: blue')
            self.exp_spec_status_label2.setText(str(self.spectrum_filename))
            self.status_text.setStyleSheet('color: blue')
            QApplication.processEvents()
        except AttributeError:
            
            #self.status_text.setText("EXP SPECTRUM ERROR")
            #self.status_text.setStyleSheet('color: yellow')
            pass  
        except IOError:
            pass          
            
            
    def replot_button(self):

        fitting_GUI_B_v10.plt.close()
        self.replot_flag = 1
        A = float(unicode(self.a_box.text()))
        B = float(unicode(self.b_box.text()))
        C = float(unicode(self.c_box.text()))
        dA = float(unicode(self.da_box.text()))
        dB = float(unicode(self.db_box.text()))
        dC = float(unicode(self.dc_box.text()))
        duA = float(unicode(self.dua_box.text()))
        duB = float(unicode(self.dub_box.text()))
        duC = float(unicode(self.duc_box.text()))
        
        DJ = float(unicode(self.DJ_box.text()))
        DJK = float(unicode(self.DJK_box.text()))      
        DK = float(unicode(self.DK_box.text()))        
        dJ = float(unicode(self.dJ_box.text()))
        dK = float(unicode(self.dK_box.text()))  
        ua = float(unicode(self.ua_box.text()))        
        ub = float(unicode(self.ub_box.text()))
        uc = float(unicode(self.uc_box.text()))
        T = float(unicode(self.T_box.text()))
        f_upper=float(unicode(self.f_upper_box.text()))
        f_lower=float(unicode(self.f_lower_box.text()))
        J_max = int(self.J_max_box.text())
        print f_upper,f_lower
        if self.counter != 0:
            picked_list = fitting_GUI_B_v10.picked_list
        
        else:
            picked_list = []
        self.counter+=1
        peaklist = []
        if self.peak_center_cb.isChecked()==True:
            peaklist = self.peaklist
        
        #fitting_GUI_B_v10.initialize(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,J_max,picked_list)  
        flag_1 = 0

         
        fitting_GUI_B_v10.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,duA,duB,duC,f_lower,f_upper,T,J_max,picked_list)
        try:
            a = self.spectrum_filename
            self.update_exp_spectrum()
        except AttributeError:
            fitting_GUI_B_v10.plt.show()    
            
    def initialize(self):
        
        A = float(unicode(self.a_box.text()))
        B = float(unicode(self.b_box.text()))
        C = float(unicode(self.c_box.text()))
        dA = float(unicode(self.da_box.text()))
        dB = float(unicode(self.db_box.text()))
        dC = float(unicode(self.dc_box.text()))
        DJ = float(unicode(self.DJ_box.text()))
        DJK = float(unicode(self.DJK_box.text()))      
        DK = float(unicode(self.DK_box.text()))        
        dJ = float(unicode(self.dJ_box.text()))
        dK = float(unicode(self.dK_box.text()))  
        ua = float(unicode(self.ua_box.text()))        
        ub = float(unicode(self.ub_box.text()))
        uc = float(unicode(self.uc_box.text()))
        T = float(unicode(self.T_box.text()))
        f_upper=float(unicode(self.f_upper_box.text()))
        f_lower=float(unicode(self.f_lower_box.text()))
        J_max = int(self.J_max_box.text())
        print f_upper,f_lower
        if self.counter != 0:
            picked_list = fitting_GUI_B_v10.picked_list
        
        else:
            picked_list = []
        self.counter+=1
        peaklist = []
        if self.peak_center_cb.isChecked()==True:
            peaklist = self.peaklist
        
                
        fitting_GUI_B_v10.initialize(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,J_max,picked_list)
        
        
    def save_input_file(self):
        job_name = unicode(self.job_box.text())
	
	A = float(unicode(self.a_box.text()))
        B = float(unicode(self.b_box.text()))
        C = float(unicode(self.c_box.text()))
        DJ = float(unicode(self.DJ_box.text()))
        DJK = float(unicode(self.DJK_box.text()))      
        DK = float(unicode(self.DK_box.text()))        
        dJ = float(unicode(self.dJ_box.text()))
        dK = float(unicode(self.dK_box.text()))  
        u_A = float(unicode(self.ua_box.text()))        
        u_B = float(unicode(self.ub_box.text()))
        u_C = float(unicode(self.uc_box.text()))
        freq_high=float(unicode(self.f_upper_box.text()))
        freq_low=float(unicode(self.f_lower_box.text()))
        inten_high=float(unicode(self.I_upper_box.text()))
        inten_low=float(unicode(self.I_lower_box.text()))	
	processors = float(unicode(self.proc_box.text()))
	temperature = float(unicode(self.T_box.text()))
	#needed format: (inten,freq, qnum_up, qnum_low,uncert)
	#(entry[0],entry[1],entry[2],entry[3],entry[4])
	#
	check_peaks_str = ""
	for line in fitting_GUI_B_v10.picked_list:
	    if line[5]=='CHECK':
	       temp_line = (line[3],line[0],line[1],line[2],line[4])
	       check_peaks_str+=str(temp_line)+"\n"

	self.fit_list = []
        counter =0
        for entry in fitting_GUI_B_v10.picked_list:
            if entry[5]=='FIT' and counter<3:
                self.fit_list.append((entry[0],entry[1],entry[2],entry[3],entry[4]))
                
                counter+=1 
        freq_high=float(unicode(self.f_upper_box.text()))
        freq_low=float(unicode(self.f_lower_box.text()))

        trans_1 = self.fit_list[0]
        trans_2 = self.fit_list[1]
        trans_3 = self.fit_list[2]
	
	
	trans_1 =  (self.fit_list[0][3],self.fit_list[0][0],self.fit_list[0][1],self.fit_list[0][2],self.fit_list[0][4]) 
	trans_2 = (self.fit_list[1][3],self.fit_list[1][0],self.fit_list[1][1],self.fit_list[1][2],self.fit_list[1][4])
	trans_3 = (self.fit_list[2][3],self.fit_list[2][0],self.fit_list[2][1],self.fit_list[2][2],self.fit_list[2][4])
	Jmax = int(self.J_max_box.text())


        job_file = ""
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save input file', ''))
        job_file += "Job Name %s \n u_A: %s \n u_B: %s \n u_C: %s \n A: %s \n B: %s \n C: %s \n DJ: %s \n DJK: %s \n DK: %s \n dJ: %s \n dK: %s \n processors: %s \n freq_high: %s \n freq_low: %s \n inten_high: %s \n inten_low: %s \n Temp: %s \n Jmax: %s \n number of triples: %s \n Check peaks:\n%s \n trans_1: %s \n trans_2: %s \n trans_3: %s "%(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,str(int(processors)),str(freq_high),\
str(freq_low),str(inten_high),str(inten_low),str(temperature),str(Jmax),str(self.num_of_triples),check_peaks_str,str(trans_1),str(trans_2),str(trans_3))
        
        
   
        Job_fh = open("%s"%(path),"w")
        Job_fh.write(job_file) 
        Job_fh.close()

    def open_input_file(self):
        
        
        path = unicode(QFileDialog.getOpenFileName(self, 
                        'Open input file', ''))
        
        fh_input = open(path)
        
        
        fitting_GUI_B_v10.picked_list = []
        fitting_peaks_flag = 0
        for line in fh_input:
                    
                    if line.split() != []:
                                        
                        if line.split()[0] == "u_A:":
                            u_A = line.split()[1] 
                            self.ua_box.setText(u_A)
                        if line.split()[0] == "u_B:":
                            u_B = line.split()[1]
                            self.ub_box.setText(u_B) 
                        if line.split()[0] == "u_C:":
                            u_C = line.split()[1] 
                            self.uc_box.setText(u_C)
                        if line.split()[0] == "A:":
                            A = line.split()[1] 
                            self.a_box.setText(A)
                        if line.split()[0] == "B:":
                            B = line.split()[1]
                            self.b_box.setText(B)
                        if line.split()[0] == "C:":
                            C = line.split()[1]
                            self.c_box.setText(C)
                        if line.split()[0] == "DJ:":
                            DJ = line.split()[1]
                            self.DJ_box.setText(DJ)
                        if line.split()[0] == "DK:":
                            DK = line.split()[1]
                            self.DK_box.setText(DK)
                        if line.split()[0] == "DJK:":
                            DJK = line.split()[1]
                            self.DJK_box.setText(DJK)
                        if line.split()[0] == "dJ:":
                            dJ = line.split()[1]
                            self.dJ_box.setText(dJ)
                        if line.split()[0] == "dK:":
                            dK = line.split()[1]
                            self.dK_box.setText(dK)
                        if line.split()[0] == "freq_high:":
                            freq_high = line.split()[1] 
                            self.f_upper_box.setText(freq_high) 
                        if line.split()[0] == "freq_low:":
                            freq_low = line.split()[1] 
                            self.f_lower_box.setText(freq_low) 
                        if line.split()[0] == "inten_high:":
                            inten_high = line.split()[1] 
                        if line.split()[0] == "inten_low:":
                            inten_low = line.split()[1]
                        if line.split()[0] == "processors:":
                            processors = float(int(line.split()[1]))
                        if line.split()[0] == "Temp:":
                            temperature = line.split()[1]
                            self.T_box.setText(temperature)
                        if line.split()[0] == "Jmax:":
                            Jmax = line.split()[1]
                            self.J_max_box.setText(Jmax)
                        if line.split()[0] == "freq_uncertainty:":
                            freq_uncertainty = line.split()[1]                    
                        if line.split()[0] == "trans_1:" or line.split()[0] == "trans_2:" or line.split()[0] == "trans_3:":
                            fitting_peaks_flag = 0
                            clean = line[12:65]
                            print clean
                            re_split = clean.split("', '")
                            tuples = tuple(re_split)
                            print tuples
                            fitting_GUI_B_v10.picked_list.append((tuples[1],tuples[2],tuples[3],tuples[0],tuples[4],'FIT'))
                            
                        if fitting_peaks_flag == 1:
                            
                            clean = line[2:55]
                            print clean
                            re_split = clean.split("', '")
                            tuples = tuple(re_split)
                            #print tuples
                            fitting_GUI_B_v10.picked_list.append((tuples[1],tuples[2],tuples[3],tuples[0],tuples[4],'CHECK'))
                        if line.split()[0] == "Check":
                            fitting_peaks_flag = 1
        #self.peaks_button_push()
        fh_input.close()
        self.peaks_button_push()
        for lin in fitting_GUI_B_v10.picked_list:
            print lin
        fitting_GUI_B_v10.update_plot()
        self.update_fit_check_status()
    def status_test(self):
        self.status_text.setText("PEAK FINDING STARTED...PLEASE WAIT...")
        self.status_text.setStyleSheet('color: red')  
        
    def peak_pick_button(self):
        
        self.status_text.setText("PEAK FINDING STARTED...PLEASE WAIT...")
        self.status_text.setStyleSheet('color: red')
        QApplication.processEvents()
        
        fh = numpy.loadtxt(str(self.spectrum_filename))
        
        
        
        

        inten_high=float(unicode(self.I_upper_box.text()))
        inten_low=float(unicode(self.I_lower_box.text()))
        

               
        spectrum_2kHz = self.cubic_spline(fh,0.02) # Interpolates experimental spectrum to a 2 kHz resolution with a cubic spline.  Gives better peak-pick values.
        (self.peaklist, freq_low, freq_high) = self.peakpicker(spectrum_2kHz,inten_low,inten_high) # Calls slightly modified version of Cristobal's routine to pick peaks instead of forcing user to do so.
        self.status_text.setText("PEAK FINDING COMPLETE, %s PEAKS HAVE BEEN FOUND IN THE EXPERIMENTAL SPECTRUM"%str(len(self.peaklist)))
        self.status_text.setStyleSheet('color: blue')
        self.peakfinder_status_label2.setText('%s PEAKS'%str(len(self.peaklist)))
        self.peakfinder_status_label2.setStyleSheet('color: blue')       
        self.update_triples_meter()
        fitting_GUI_B_v10.peak_list = self.peaklist
    def on_about(self):
        msg = """ Stuff about how to use autofit, documentation, etc. here
        """
        QMessageBox.about(self, "About Autofit", msg.strip())
    
    def peakpicker(self,spectrum,thresh_l,thresh_h):#Code taken from Cristobal's peak-picking script; assumes spectrum is in increasing frequency order
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


    def cubic_spline(self,spectrum,new_resolution): # Cubic spline of spectrum to new_resolution; used pre-peak-picking.  Assumes spectrum is already in order of increasing frequency.
    
        x = spectrum[:,0]
        y = spectrum[:,1]
    
        old_resolution = (x[-1]-x[0]) / len(spectrum)
        scale_factor = old_resolution / new_resolution
    
        #new_length = int(math.floor(scale_factor*len(spectrum)))
    
        tck = splrep(x,y,s=0)
        xnew = numpy.arange(x[0],x[-1],new_resolution)
        ynew = splev(xnew,tck,der=0)
    
        output_spectrum = numpy.column_stack((xnew,ynew))
    
        return output_spectrum









    def run_autofit(self):
        self.status_text.setText("AUTOFIT RUN STARTED...")
        self.status_text.setStyleSheet('color: red')
        QApplication.processEvents()
	for line in fitting_GUI_B_v10.picked_list:
	    print line

        
        #self.fit_list 
        
        #self.trans_2_peaks
        #self.trans_3_peaks
        


	job_name = unicode(self.job_box.text())
	
	A = float(unicode(self.a_box.text()))
        B = float(unicode(self.b_box.text()))
        C = float(unicode(self.c_box.text()))
        DJ = float(unicode(self.DJ_box.text()))
        DJK = float(unicode(self.DJK_box.text()))      
        DK = float(unicode(self.DK_box.text()))        
        dJ = float(unicode(self.dJ_box.text()))
        dK = float(unicode(self.dK_box.text()))  
        u_A = float(unicode(self.ua_box.text()))        
        u_B = float(unicode(self.ub_box.text()))
        u_C = float(unicode(self.uc_box.text()))
        freq_high=float(unicode(self.f_upper_box.text()))
        freq_low=float(unicode(self.f_lower_box.text()))
        inten_high=float(unicode(self.I_upper_box.text()))
        inten_low=float(unicode(self.I_lower_box.text()))	
	processors = float(unicode(self.proc_box.text()))
	temperature = float(unicode(self.T_box.text()))
	#needed format: (inten,freq, qnum_up, qnum_low,uncert)
	#(entry[0],entry[1],entry[2],entry[3],entry[4])
	#
	check_peaks_list = []
	for line in fitting_GUI_B_v10.picked_list:
	    if line[5]=='CHECK':
	       temp_line = (line[3],line[0],line[1],line[2],line[4])
	       check_peaks_list.append(temp_line)
	trans_1_peaks = self.trans_1_peaks
	trans_2_peaks = self.trans_2_peaks
	trans_3_peaks = self.trans_3_peaks
	self.fit_list = []
        counter =0
        for entry in fitting_GUI_B_v10.picked_list:
            if entry[5]=='FIT' and counter<3:
                self.fit_list.append((entry[0],entry[1],entry[2],entry[3],entry[4]))
                
                counter+=1 
        freq_high=float(unicode(self.f_upper_box.text()))
        freq_low=float(unicode(self.f_lower_box.text()))

        trans_1 = self.fit_list[0]
        trans_2 = self.fit_list[1]
        trans_3 = self.fit_list[2]
	
	
	trans_1 =  (self.fit_list[0][3],self.fit_list[0][0],self.fit_list[0][1],self.fit_list[0][2],self.fit_list[0][4]) 
	trans_2 = (self.fit_list[1][3],self.fit_list[1][0],self.fit_list[1][1],self.fit_list[1][2],self.fit_list[1][4])
	trans_3 = (self.fit_list[2][3],self.fit_list[2][0],self.fit_list[2][1],self.fit_list[2][2],self.fit_list[2][4])
	Jmax = int(self.J_max_box.text())

	
	
	autofit_NS_module.autofit_NS(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,freq_high,freq_low,inten_high,inten_low,int(processors),temperature,Jmax,trans_1,trans_2,trans_3,check_peaks_list,self.peaklist,trans_1_peaks,trans_2_peaks,trans_3_peaks)
        self.status_text.setText("AUTOFIT RUN COMPLETED")
        self.status_text.setStyleSheet('color: blue')

    def on_item_changed(self):
        print "test"
        print self.model.item(0).checkState()

    def peaks_button_push(self):
        
        self.model = QStandardItemModel()

        x = 0

        for entry in fitting_GUI_B_v10.picked_list:
            new_entry = entry
            item = QStandardItem(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0])+' | '+str(new_entry[3])+' | '+str(new_entry[5]))
 
            # add a checkbox to it
            item.setCheckable(True)

            # Add the item to the model
            self.model.appendRow(item)
            #self.model.appendColumn(item)        
            
            """
            x+=1
            vars()['item%s'%str(x)] = QStandardItem(str(entry)) 

            vars()['item%s'%str(x)].setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) 
            vars()['item%s'%str(x)].setData(QVariant(Qt.Checked), Qt.CheckStateRole) 
            model.appendRow(vars()['item%s'%str(x)]) 
            """
       
        self.lv.setModel(self.model) 
        self.model.itemChanged.connect(self.on_item_changed)
        self.update_triples_labels()
        self.update_fit_check_status()
        #print fitting_GUI_B_v10.picked_list
    
    
    def distortion_toggle(self):
        
        if self.distortion_cb.isChecked()==False:
            self.DJ_box.hide()           
            self.DJK_box.hide()   
            self.DK_box.hide()           
            self.dJ_box.hide()
            self.dK_box.hide()   
            self.DJ_label.hide()           
            self.DJK_label.hide()   
            self.DK_label.hide()           
            self.dJ_label.hide()
            self.dK_label.hide() 
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)
        if self.distortion_cb.isChecked()==True:
            self.DJ_box.show()           
            self.DJK_box.show()   
            self.DK_box.show()           
            self.dJ_box.show()
            self.dK_box.show()   
            self.DJ_label.show()           
            self.DJK_label.show()   
            self.DK_label.show()           
            self.dJ_label.show()
            self.dK_label.show() 
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)     
    def fix_toggle(self):
        if self.fix_cb.isChecked()==False:
            self.a_cb.hide()
            self.b_cb.hide()
            self.c_cb.hide()
            self.DJ_cb.hide()
            self.DJK_cb.hide()
            self.DK_cb.hide()
            self.dJ_cb.hide()
            self.dK_cb.hide()  
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)            
        if self.fix_cb.isChecked()==True:
            self.a_cb.show()
            self.b_cb.show()
            self.c_cb.show()
            self.DJ_cb.show()
            self.DJK_cb.show()
            self.DK_cb.show()
            self.dJ_cb.show()
            self.dK_cb.show() 
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)                                   
    def spectrum_toggle(self):  
        if self.spectrum_cb.isChecked()==False:  
            self.I_upper_label.hide()
            self.I_lower_label.hide()
            self.f_upper_label.hide()
            self.f_lower_label.hide()
            #self.peak_finder_button.hide()
            self.f_lower_box.hide()
            self.f_upper_box.hide()
            self.I_lower_box.hide()
            self.I_upper_box.hide()
            self.skip_box.hide()
            self.skip_label.hide()
            self.peak_center_cb.hide()
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)            
                 
                           
        if self.spectrum_cb.isChecked()==True:     
            self.I_upper_label.show()
            self.I_lower_label.show()
            self.f_upper_label.show()
            self.f_lower_label.show()
            #self.peak_finder_button.show()
            self.f_lower_box.show()
            self.f_upper_box.show()
            self.I_lower_box.show()
            self.I_upper_box.show()
            self.skip_box.show()
            self.skip_label.show()
            self.peak_center_cb.show()
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)     
            
    def delete_button_push(self):
            counter = 0
            new_picked_list = []
            for x in fitting_GUI_B_v10.picked_list:
                if self.model.item(counter).checkState()==0:
                    new_picked_list.append(x)

                else:
                    pass
                counter+=1                                                      
            fitting_GUI_B_v10.picked_list = new_picked_list    
            
               
            self.model = QStandardItemModel()
            
            x = 0
            for entry in fitting_GUI_B_v10.picked_list:
                new_entry = entry
                item = QStandardItem(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0])+' | '+str(new_entry[3])+' | '+str(new_entry[5]))
    
                # add a checkbox to it
                item.setCheckable(True)
            
                # Add the item to the model
                self.model.appendRow(item)
                        
                
                """
                x+=1
                vars()['item%s'%str(x)] = QStandardItem(str(entry)) 
    
                vars()['item%s'%str(x)].setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) 
                vars()['item%s'%str(x)].setData(QVariant(Qt.Checked), Qt.CheckStateRole) 
                model.appendRow(vars()['item%s'%str(x)]) 
                """
        
            self.lv.setModel(self.model) 
            self.model.itemChanged.connect(self.on_item_changed)   
            fitting_GUI_B_v10.update_plot()
            self.update_triples_meter()  
            self.update_fit_check_status()                                                                                                                                                            
    def autofit_toggle(self):                                                  
        if self.autofit_cb.isChecked()==False: 
            self.job_box.hide()
            self.job_label.hide()
            self.proc_label.hide()
            self.proc_box.hide()
	    self.set_fit_button.hide()
	    self.set_check_button.hide()
            self.fix_cb.hide()
            self.uncertainty_box_1.hide()
            self.uncertainty_box_2.hide()
            self.uncertainty_box_3.hide()
            self.uncertainty_label.hide()
            self.uncertainty_label2.hide()
            self.plot_triples_cb.hide()
            self.trans_1_label.hide()
            self.trans_2_label.hide()
            self.trans_3_label.hide()
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)
        if self.autofit_cb.isChecked()==True:        
            self.job_box.show()
            self.job_label.show()
            self.trans_1_label.show()
            self.trans_2_label.show()
            self.trans_3_label.show()             
            self.proc_label.show()
            self.proc_box.show()    
	    self.set_fit_button.show()
	    self.set_check_button.show()
            self.fix_cb.show()
            self.uncertainty_box_1.show()
            self.uncertainty_box_2.show()
            self.uncertainty_box_3.show()
            self.uncertainty_label.show()
            self.uncertainty_label2.show()
            self.plot_triples_cb.show()
            
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)      

    
        
    def open_peaklist_file(self):
        path = unicode(QFileDialog.getOpenFileName(self, 
            'Open peak list file', ''))
        fh = numpy.loadtxt(path)

        self.peaklist = fh
        self.status_text.setText("PEAK LIST FILE LOADED")
        self.status_text.setStyleSheet('color: blue')
        self.peakfinder_status_label2.setText('%s PEAKS'%str(len(self.peaklist)))
        self.peakfinder_status_label2.setStyleSheet('color: blue')   
        self.update_triples_meter()
      
    def save_peaklist_file(self):
        fname = unicode(QFileDialog.getSaveFileName(self, 
            'Save peak list file', ''))
        numpy.savetxt(fname,self.peaklist)
        
    def set_fit_button_push(self):
            counter = 0
            new_picked_list = []
            for x in fitting_GUI_B_v10.picked_list:
                if self.model.item(counter).checkState()==0:
                    new_picked_list.append(x)

                else:
                    new_picked_list.append((x[0],x[1],x[2],x[3],x[4],'FIT'))
                counter+=1                                                      
            fitting_GUI_B_v10.picked_list = new_picked_list    
            
               
            self.model = QStandardItemModel()
            
            x = 0
            for entry in fitting_GUI_B_v10.picked_list:
                new_entry = entry
                item = QStandardItem(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0])+' | '+str(new_entry[3])+' | '+str(new_entry[5]))
    
                # add a checkbox to it
                item.setCheckable(True)
            
                # Add the item to the model
                self.model.appendRow(item)
                        
                
                """
                x+=1
                vars()['item%s'%str(x)] = QStandardItem(str(entry)) 
    
                vars()['item%s'%str(x)].setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) 
                vars()['item%s'%str(x)].setData(QVariant(Qt.Checked), Qt.CheckStateRole) 
                model.appendRow(vars()['item%s'%str(x)]) 
                """
        
            self.lv.setModel(self.model) 
            self.model.itemChanged.connect(self.on_item_changed)   
            fitting_GUI_B_v10.update_plot() 
            self.plot_triples()
            self.update_triples_meter()
            self.update_triples_labels()
            self.update_fit_check_status()
    def update_triples_labels(self):
            counter = 0
            self.trans_1_label.setText('')
            self.trans_2_label.setText('')
            self.trans_3_label.setText('')
            for new_entry in fitting_GUI_B_v10.picked_list:
                if new_entry[5]=='FIT' and counter == 0:
                    self.trans_1_label.setText(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0]))   
                    counter +=1
                    continue
                if new_entry[5]=='FIT' and counter == 1:
                    self.trans_2_label.setText(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0]))                     
                    counter +=1
                    continue
                if new_entry[5]=='FIT' and counter == 2:
                    self.trans_3_label.setText(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0]))     
                    counter +=1
                    continue
            self.update_fit_check_status()
    def plot_current_peaklist(self):
        if self.peak_center_cb.isChecked()==True:
            self.peaklist
            print type(self.peaklist)
            fitting_GUI_B_v10.peak_list_plt.set_data(self.peaklist[:,0],self.peaklist[:,1]) 
            
            fitting_GUI_B_v10.ax2.set_ylim([0,self.peaklist.max(0)[1]])    
            fitting_GUI_B_v10.plt.show()          
        if self.peak_center_cb.isChecked()==False:
            
            fitting_GUI_B_v10.peak_list_plt.set_data([],[]) 
            fitting_GUI_B_v10.plt.show()  
    def set_check_button_push(self):               
            counter = 0
            new_picked_list = []
            for x in fitting_GUI_B_v10.picked_list:
                if self.model.item(counter).checkState()==0:
                    new_picked_list.append(x)

                else:
                    new_picked_list.append((x[0],x[1],x[2],x[3],x[4],'CHECK'))
                counter+=1                                                      
            fitting_GUI_B_v10.picked_list = new_picked_list    
            
               
            self.model = QStandardItemModel()
            
            x = 0
            for entry in fitting_GUI_B_v10.picked_list:
                new_entry = entry
                item = QStandardItem(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0])+' | '+str(new_entry[3])+' | '+str(new_entry[5]))
    
                # add a checkbox to it
                item.setCheckable(True)
            
                # Add the item to the model
                self.model.appendRow(item)
                        
                
                """
                x+=1
                vars()['item%s'%str(x)] = QStandardItem(str(entry)) 
    
                vars()['item%s'%str(x)].setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) 
                vars()['item%s'%str(x)].setData(QVariant(Qt.Checked), Qt.CheckStateRole) 
                model.appendRow(vars()['item%s'%str(x)]) 
                """
        
            self.lv.setModel(self.model) 
            self.model.itemChanged.connect(self.on_item_changed)   
            fitting_GUI_B_v10.update_plot()     
            self.plot_triples()    
            self.update_triples_meter()      
            self.update_triples_labels() 
            self.update_fit_check_status()                          
    def display_peakfind_toggle(self):
        asd                 
    def uncertainty_box_edit(self):
        self.update_triples_meter()   
     
    def read_results_of_autofit(self, text):

        print text       
        text = str(text)
        try:
       	    self.a_box.setText(text.split()[0])
       	    self.b_box.setText(text.split()[1])
            self.c_box.setText(text.split()[2])
        except: pass         

        
	if text=="Display Original Input Constants Constants":
            A = '3000'
            B = '2000'
            C = '1000'
            self.a_box.setText(A)
	    self.b_box.setText(B)
            self.c_box.setText(C)

        """
	self.DJ_box.setText('0.0')
	self.DJK_box.setText('0.0')
	self.DK_box.setText('0.0')
	self.dJ_box.setText('0.0')
	self.dK_box.setText('0.0')
	self.ua_box.setText('1.0')
	self.ub_box.setText('1.0')
	self.uc_box.setText('1.0')
        """                                                                                                                                                                                                          
    def create_main_frame(self):
        self.counter = 0
        self.main_frame = QWidget()
        #list_data = [1,2,3,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    
        

        self.lv = QListView()
        #self.lv = QTreeView()       
        
        self.distortion_cb = QCheckBox('Show Distortions', self)
        #self.uncertainty_cb = QCheckBox('Show Distortion Constants', self)

        self.distortion_cb.stateChanged.connect(self.distortion_toggle)
 
        self.spectrum_cb = QCheckBox('Show Spectral Options', self)
        

        self.spectrum_cb.stateChanged.connect(self.spectrum_toggle)       
                      
        self.autofit_cb = QCheckBox('Show Autofit Options', self)       
        self.autofit_cb.stateChanged.connect(self.autofit_toggle)  
        self.fix_cb = QCheckBox('Manual Fix Const', self)
        self.fix_cb.stateChanged.connect(self.fix_toggle) 
        self.a_cb = QCheckBox('fix', self)
        self.b_cb = QCheckBox('fix', self)
        self.c_cb = QCheckBox('fix', self)
        self.DJ_cb = QCheckBox('fix', self)
        self.DJK_cb = QCheckBox('fix', self)
        self.DK_cb = QCheckBox('fix', self)
        self.dJ_cb = QCheckBox('fix', self)
        self.dK_cb = QCheckBox('fix', self)
        self.peak_center_cb = QCheckBox('Plot Current Peaklist', self)
        self.peak_center_cb.stateChanged.connect(self.plot_current_peaklist)
        
        
                        
        self.a_box = QLineEdit()
	self.b_box = QLineEdit()
	self.c_box = QLineEdit()
        self.da_box = QLineEdit()
	self.db_box = QLineEdit()
	self.dc_box = QLineEdit()
    
        self.DJ_box = QLineEdit()
	self.DJK_box = QLineEdit()
	self.DK_box = QLineEdit()

	self.dJ_box = QLineEdit()
	self.dK_box = QLineEdit()	
	
	self.ua_box = QLineEdit()
	self.ub_box = QLineEdit()
	self.uc_box = QLineEdit()
        
	self.dua_box = QLineEdit()
	self.dub_box = QLineEdit()
	self.duc_box = QLineEdit()
	
        self.T_box = QLineEdit()
        self.J_max_box = QLineEdit()
	self.proc_box = QLineEdit() 
        
        self.f_lower_box = QLineEdit()
        self.f_upper_box = QLineEdit()
        self.I_lower_box = QLineEdit()
        self.I_upper_box = QLineEdit()
        self.job_box = QLineEdit()
	self.uncertainty_box_1 = QLineEdit()
	self.connect(self.uncertainty_box_1,SIGNAL("textEdited(QString)"),self.uncertainty_box_edit)
	self.uncertainty_box_2 = QLineEdit()
	self.connect(self.uncertainty_box_2,SIGNAL("textEdited(QString)"),self.uncertainty_box_edit)
	self.uncertainty_box_3 = QLineEdit()
	self.connect(self.uncertainty_box_3,SIGNAL("textEdited(QString)"),self.uncertainty_box_edit)
        self.plot_triples_cb = QCheckBox('Plot Triples')
        self.connect(self.plot_triples_cb,SIGNAL('stateChanged(int)'), self.plot_triples)
                 
                                    
                                                                        
        self.skip_box = QLineEdit()        
        


        #self.connect(self.textbox2, SIGNAL('editingFinished ()'), self.on_draw)       
        #self.connect(self.textbox3, SIGNAL('editingFinished ()'), self.on_draw) 
        self.draw_button = QPushButton("&Plot")
	self.autofit_button = QPushButton("&Run Autofit")
        self.peaks_button = QPushButton("&Show Selected Peaks")
        #self.spectrum_file_button = QPushButton("&Open Exp Spectrum File")
        self.peak_finder_button = QPushButton("&Run Peakfinder")
        self.connect(self.autofit_button, SIGNAL('clicked()'), self.run_autofit)
        self.connect(self.draw_button, SIGNAL('clicked()'), self.replot_button)
        self.connect(self.peaks_button, SIGNAL('clicked()'), self.peaks_button_push)
        self.connect(self.peak_finder_button, SIGNAL('clicked()'), self.peak_pick_button)
        #self.connect(self.spectrum_file_button, SIGNAL('clicked()'), self.selectFile_spectrum)
        
        self.delete_button = QPushButton("&Delete Transitions")
        self.connect(self.delete_button, SIGNAL('clicked()'), self.delete_button_push)
        
	self.set_fit_button = QPushButton("&Set Fit Transitions")
	self.connect(self.set_fit_button, SIGNAL('clicked()'), self.set_fit_button_push)
	
	self.set_check_button = QPushButton("&Set Check Transitions")
	self.connect(self.set_check_button, SIGNAL('clicked()'),self.set_check_button_push)
	
	self.auto_set_button = QPushButton("&Auto Setup")
        self.connect(self.auto_set_button, SIGNAL('clicked()'), self.auto_set_button_push)
        
        
        

	self.a_label = QLabel('A')
	
	self.b_label = QLabel('B')
	self.c_label = QLabel('C')
	self.DJ_label = QLabel('-DJ')
	self.DJK_label = QLabel('-DJK')
	self.DK_label = QLabel('-DK')
	self.dJ_label = QLabel('-dJ')
	self.dK_label = QLabel('-dK')
	self.ua_label = QLabel('mu a')
	self.ub_label = QLabel('mu b')
	self.uc_label = QLabel('mu c')
	
	self.f_upper_label = QLabel('Upper Frequency Limit')
	self.f_lower_label = QLabel('Lower Frequency Limit')
	self.uncertainty_label = QLabel('Peak Uncertainty')
        self.uncertainty_label2 = QLabel('Peak info')
	self.I_upper_label = QLabel('Upper Intensity Limit')
	self.I_lower_label = QLabel('Lower Intensity Limit')
	self.T_label = QLabel('T')
	self.J_max_label = QLabel('Max J')
	self.proc_label = QLabel('Number of Processors') 
        self.col_label = QLabel('Input Value')
        self.col_label1 = QLabel('Slider Range')
        self.skip_label = QLabel('Skip nth Point:')
        self.job_label = QLabel('Job Name')
        self.num_of_triples_label = QLabel('0 Triples')
        self.num_of_triples_label.setStyleSheet('color: red')
        self.trans_1_label = QLabel('')
        self.trans_2_label = QLabel('')
        self.trans_3_label = QLabel('')
        
        self.exp_spec_status_label1 = QLabel('Exp Spec File:')
        self.exp_spec_status_label2 = QLabel('None Loaded')
        self.exp_spec_status_label2.setStyleSheet('color: blue')
        
        self.peakfinder_status_label1 = QLabel('Peakfinder Status:')
        self.peakfinder_status_label2 = QLabel('None')
        self.peakfinder_status_label2.setStyleSheet('color: red')

        self.fitting_status_label1 = QLabel('Fitting Trans Status:')
        self.fitting_status_label2 = QLabel('None Selected')
        self.fitting_status_label2.setStyleSheet('color: red') 
               
        self.check_status_label1 = QLabel('Check Trans Status:')
        self.check_status_label2 = QLabel('None Selected')
        self.check_status_label2.setStyleSheet('color: red')      
        
        
	hbox = QtGui.QGridLayout()
        #hbox = QHBoxLayout()
        #widget_list = [ a_label, self.textbox,b_label,self.textbox2 ,self.draw_button, self.grid_cb,
        #            slider_label, self.slider]
        #for w in range(len(widget_list)):
        #    hbox.addWidget(widget_list[w])
        #    hbox.setAlignment(w, Qt.AlignVCenter)
        hbox.addWidget(self.distortion_cb,0,2)
        hbox.addWidget(self.spectrum_cb,0,6)
        hbox.addWidget(self.a_label,2,1)
        hbox.addWidget(self.b_label,3,1)
	hbox.addWidget(self.c_label,4,1)
	hbox.addWidget(self.DJ_label,10,1)
	hbox.addWidget(self.DJK_label,11,1)
	hbox.addWidget(self.DK_label,12,1)
	hbox.addWidget(self.dJ_label,13,1)
	hbox.addWidget(self.dK_label,14,1)
	hbox.addWidget(self.uncertainty_label,3,8)
	hbox.addWidget(self.uncertainty_label2,3,7)
	hbox.addWidget(self.uncertainty_box_1,4,8)
	hbox.addWidget(self.uncertainty_box_2,5,8)
	hbox.addWidget(self.uncertainty_box_3,6,8)
	hbox.addWidget(self.trans_1_label,4,7)
	hbox.addWidget(self.trans_2_label,5,7)
	hbox.addWidget(self.trans_3_label,6,7)
	hbox.addWidget(self.ua_label,5,1)
        hbox.addWidget(self.ub_label,6,1)
        hbox.addWidget(self.uc_label,7,1)
        hbox.addWidget(self.T_label,8,1)
        hbox.addWidget(self.f_lower_label,1,5)
        hbox.addWidget(self.f_upper_label,2,5)
        hbox.addWidget(self.I_lower_label,3,5)
        hbox.addWidget(self.I_upper_label,4,5)
        hbox.addWidget(self.job_label,1,7)

        
        hbox.addWidget(self.col_label,1,2)
        hbox.addWidget(self.col_label1,1,3)
        hbox.addWidget(self.skip_label,5,5)


	
        hbox.addWidget(self.a_box,2,2)
        hbox.addWidget(self.b_box,3,2)
	hbox.addWidget(self.c_box,4,2)
        hbox.addWidget(self.da_box,2,3)
        hbox.addWidget(self.db_box,3,3)
	hbox.addWidget(self.dc_box,4,3)
        hbox.addWidget(self.dua_box,5,3)
        hbox.addWidget(self.dub_box,6,3)
	hbox.addWidget(self.duc_box,7,3)
	
	hbox.addWidget(self.DJ_box,10,2)
	hbox.addWidget(self.DJK_box,11,2)
	hbox.addWidget(self.DK_box,12,2)	
	hbox.addWidget(self.dJ_box,13,2)	
	hbox.addWidget(self.dK_box,14,2)	
	
        hbox.addWidget(self.a_cb,2,4)
        hbox.addWidget(self.b_cb,3,4)
        hbox.addWidget(self.c_cb,4,4)
        hbox.addWidget(self.DJ_cb,10,4)
        hbox.addWidget(self.DJK_cb,11,4)
        hbox.addWidget(self.DK_cb,12,4)
        hbox.addWidget(self.dJ_cb,13,4)
        hbox.addWidget(self.dK_cb,14,4)	
        hbox.addWidget(self.fix_cb,7,8)		
        hbox.addWidget(self.peak_center_cb,6,6)	

	
	
	hbox.addWidget(self.ua_box,5,2)
        hbox.addWidget(self.ub_box,6,2)
        hbox.addWidget(self.uc_box,7,2)
        hbox.addWidget(self.T_box,8,2)
        hbox.addWidget(self.J_max_label,9,1)	
        hbox.addWidget(self.J_max_box,9,2)
        hbox.addWidget(self.f_lower_box,1,6)
        hbox.addWidget(self.f_upper_box,2,6)
        hbox.addWidget(self.I_lower_box,3,6)
        hbox.addWidget(self.I_upper_box,4,6)
        
                        
	hbox.addWidget(self.draw_button,0,1)

	#hbox.addWidget(slider_label,2,3)
	hbox.addWidget(self.peak_finder_button,0,5)
	hbox.addWidget(self.autofit_button,0,7)
	
        hbox.addWidget(self.autofit_cb,0,8)
        hbox.addWidget(self.peaks_button,0,9)
        #hbox.addWidget(self.spectrum_file_button,0,5)
	hbox.addWidget(self.proc_label,2,7)
	hbox.addWidget(self.proc_box,2,8)
	hbox.addWidget(self.job_box,1,8)
	hbox.addWidget(self.lv,1,9,8,11)
	hbox.addWidget(self.delete_button,9,9)
	hbox.addWidget(self.set_fit_button,9,10)
	hbox.addWidget(self.set_check_button,9,11)
	
	self.combo = QtGui.QComboBox(self)
	
	self.combo.addItem("Read Autofit Output HERE")
	self.combo.addItem("Display Original Input Constants Constants")
	self.combo.addItem("300 200 100 Score=20, RMS=0.01")
	self.combo.addItem("3200 1300 120 Score=20, RMS=0.01")
	self.combo.addItem("300 200 100 Score=20, RMS=0.01")
	self.combo.addItem("3200 1300 120 Score=20, RMS=0.01")
	self.combo.addItem("300 200 100 Score=20, RMS=0.01")
	self.combo.addItem("3200 1300 120 Score=20, RMS=0.01")
	self.combo.addItem("300 200 100 Score=20, RMS=0.01")
	self.combo.addItem("3200 1300 120 Score=20, RMS=0.01")	
	self.combo.activated[str].connect(self.read_results_of_autofit) 
	hbox.addWidget(self.combo,9,9,9,9)
	
	

    	hbox.addWidget(self.auto_set_button,0,10)
	hbox.addWidget(self.num_of_triples_label,0,11)
	hbox.addWidget(self.skip_box,5,6)
	hbox.addWidget(self.plot_triples_cb,7,7)
	
	hbox.addWidget(self.exp_spec_status_label1,8,5)
        hbox.addWidget(self.exp_spec_status_label2,8,6)

        
        hbox.addWidget(self.peakfinder_status_label1,9,5)
        hbox.addWidget(self.peakfinder_status_label2,9,6)

        hbox.addWidget(self.fitting_status_label1,8,7)
        hbox.addWidget(self.fitting_status_label2,8,8)
               
        hbox.addWidget(self.check_status_label1,9,7)
        hbox.addWidget(self.check_status_label2,9,8)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
        vbox = QVBoxLayout()
        vbox.addLayout(hbox)
        
        
        
        
        
        
        
        self.main_frame.setLayout(vbox)


        
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):
        self.status_text = QLabel("Autofit GUI")
        
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        open_exp_spec = self.create_action("&Open Experimental Spectrum",
            shortcut="Ctrl+O", slot=self.selectFile_spectrum, 
            tip="Open an Experimental Spectrum File")
        
        save_file_action = self.create_action("&Save input file",
            shortcut="Ctrl+S", slot=self.save_input_file, 
            tip="Save an input file")
            
        load_file_action = self.create_action("&Open input file",
            shortcut="Ctrl+I", slot=self.open_input_file, 
            tip="Open an input file")
        load_peaks_action = self.create_action("&Open peak list file",
            shortcut="Ctrl+P", slot=self.open_peaklist_file, 
            tip="Open a peak list file")
        save_peaks_action = self.create_action("&Save peak list file",
            shortcut="Ctrl+X", slot=self.save_peaklist_file, 
            tip="Save a peak list file")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (open_exp_spec,save_file_action,load_file_action,load_peaks_action,save_peaks_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()



if __name__ == "__main__":
    main()
