
import sys, os, random

import fitting_GUI_B_v6
import autofit_NS_module

import subprocess
from multiprocessing import Process
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
import numpy
import math
from scipy.interpolate import *


class AppForm(QMainWindow):
    def __init__(self, parent=None):
	QMainWindow.__init__(self, parent)
	self.setWindowTitle('Autofit GUI')

	self.create_menu()
	self.create_main_frame()
	self.create_status_bar()

	self.a_box.setText('3000')
	self.b_box.setText('2000')
	self.c_box.setText('1000')
	self.da_box.setText('1000')
	self.db_box.setText('1000')
	self.dc_box.setText('1000')
	
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
        self.distortion_cb.toggle()
        self.distortion_cb.toggle()
        self.spectrum_cb.toggle()
        self.spectrum_cb.toggle()
        self.autofit_cb.toggle()
        self.autofit_cb.toggle()
        self.fix_cb.toggle()
        self.fix_cb.toggle()
	self.a_set  = '3000'
	self.b_set  = '2000'
	self.c_set  = '1000'     
	 
        
    def selectFile_spectrum(self):
        self.filename = QFileDialog.getOpenFileName()
        print self.filename






    def replot_button(self):
        
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
        print f_upper,f_lower
        if self.counter != 0:
            picked_list = fitting_GUI_B_v6.picked_list
        
        else:
            picked_list = []
        self.counter+=1
        peaklist = []
        if self.peak_center_cb.isChecked()==True:
            peaklist = self.peaklist
        
        try:
            try:
                
                fitting_GUI_B_v6.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,picked_list,int(unicode(self.skip_box.text())),self.filename,peaklist)
            except AttributeError:
                fitting_GUI_B_v6.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,picked_list,int(unicode(self.skip_box.text())),self.filename)
        except AttributeError:
            fitting_GUI_B_v6.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,picked_list,int(unicode(self.skip_box.text())))


    def save_input_file(self):
               
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', ''))
        

    def open_input_file(self):
        
        path = unicode(QFileDialog.getOpenFileName(self, 
                        'Open file', ''))
        fh_input = open(path)

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
                        if line.split()[0] == "freq_uncertainty:":
                            freq_uncertainty = line.split()[1]                    
                        if line.split()[0] == "trans_1:" or line.split()[0] == "trans_2:" or line.split()[0] == "trans_3:":
                            fitting_peaks_flag = 0
                            clean = line[12:53]
                            re_split = clean.split("', '")
                            tuples = tuple(re_split)
                            #fit_peaks_list.append(tuples)


    def peak_pick_button(self):
        
        fh = fh = numpy.loadtxt(str(self.filename))
        
        
        
        

        inten_high=float(unicode(self.I_upper_box.text()))
        inten_low=float(unicode(self.I_lower_box.text()))
        
        
               
        spectrum_2kHz = self.cubic_spline(fh,0.02) # Interpolates experimental spectrum to a 2 kHz resolution with a cubic spline.  Gives better peak-pick values.
        (self.peaklist, freq_low, freq_high) = self.peakpicker(spectrum_2kHz,inten_low,inten_high) # Calls slightly modified version of Cristobal's routine to pick peaks instead of forcing user to do so.
        print self.peaklist

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
	peaklist = self.peaklist
	autofit_NS_module.autofit_NS(job_name,u_A,u_B,u_C,A,B,C,DJ,DJK,DK,dJ,dK,freq_high,freq_low,inten_high,inten_low,processors,temperature,Jmax,trans_1,trans_2,trans_3,check_peaks_list,peaklist,trans_1_peaks,trans_2_peaks,trans_3_peaks)


    def on_item_changed(self):
        print "test"
        print self.model.item(0).checkState()

    def peaks_button_push(self):
        
        self.model = QStandardItemModel()
        
        x = 0
        for entry in fitting_GUI_B_v6.picked_list:
            new_entry = entry
            item = QStandardItem(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0])+' | '+str(new_entry[3]))
 
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
        
        #print fitting_GUI_B_v6.picked_list

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
            self.peak_finder_button.hide()
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
            self.peak_finder_button.show()
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
            for x in fitting_GUI_B_v6.picked_list:
                if self.model.item(counter).checkState()==0:
                    new_picked_list.append(x)
                else:
                    continue
                counter+=1                                                      
            fitting_GUI_B_v6.picked_list = new_picked_list    
            
               
            self.model = QStandardItemModel()
            
            x = 0
            for entry in fitting_GUI_B_v6.picked_list:
                new_entry = entry
                item = QStandardItem(str(new_entry[1])+' | '+str(new_entry[2])+' | '+str(new_entry[0])+' | '+str(new_entry[3]))
    
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
            fitting_GUI_B_v6.update_plot()
                                                                                                                                                                               
    def autofit_toggle(self):                                                  
        if self.autofit_cb.isChecked()==False: 
            self.job_box.hide()
            self.job_label.hide()
            self.proc_label.hide()
            self.proc_box.hide()
	    self.set_fit_button.hide()
	    self.set_check_button.hide()
            self.fix_cb.hide()
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)
        if self.autofit_cb.isChecked()==True:        
            self.job_box.show()
            self.job_label.show()             
            self.proc_label.show()
            self.proc_box.show()    
	    self.set_fit_button.show()
	    self.set_check_button.show()
            self.fix_cb.show()
            self.layout().setSizeConstraint(QtGui.QLayout.SetFixedSize)      
            
    def display_peakfind_toggle(self):
        asd                                                                                                                                                                                                                                                                      
    def create_main_frame(self):
        self.counter = 0
        self.main_frame = QWidget()
        #list_data = [1,2,3,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

        

        self.lv = QListView()
               
        
        self.distortion_cb = QCheckBox('Show Distortion Constants', self)
        

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
        self.peak_center_cb = QCheckBox('Show Peakfind Results', self)
        
        
                        
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
        self.T_box = QLineEdit()
	self.proc_box = QLineEdit() 
        
        self.f_lower_box = QLineEdit()
        self.f_upper_box = QLineEdit()
        self.I_lower_box = QLineEdit()
        self.I_upper_box = QLineEdit()
        self.job_box = QLineEdit()

        

                  
        self.skip_box = QLineEdit()        
        


        #self.connect(self.textbox2, SIGNAL('editingFinished ()'), self.on_draw)       
        #self.connect(self.textbox3, SIGNAL('editingFinished ()'), self.on_draw) 
        self.draw_button = QPushButton("&Plot")
	self.autofit_button = QPushButton("&Run Autofit")
        self.peaks_button = QPushButton("&Show Selected Peaks")
        self.spectrum_file_button = QPushButton("&Open Exp Spectrum File")
        self.peak_finder_button = QPushButton("&Find Peaks in Exp Spectrum")
        self.connect(self.autofit_button, SIGNAL('clicked()'), self.run_autofit)
        self.connect(self.draw_button, SIGNAL('clicked()'), self.replot_button)
        self.connect(self.peaks_button, SIGNAL('clicked()'), self.peaks_button_push)
        self.connect(self.peak_finder_button, SIGNAL('clicked()'), self.peak_pick_button)
        self.connect(self.spectrum_file_button, SIGNAL('clicked()'), self.selectFile_spectrum)
        
        self.delete_button = QPushButton("&Delete Transitions")
        self.connect(self.delete_button, SIGNAL('clicked()'), self.delete_button_push)
	self.set_fit_button = QPushButton("&Set Fit Transitions")
	self.set_check_button = QPushButton("&Set Check Transitions")
	self.auto_set_button = QPushButton("&Auto Setup")
        
        
        
        
        global a_label
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
	
	
	self.I_upper_label = QLabel('Upper Intensity Limit')
	self.I_lower_label = QLabel('Lower Intensity Limit')
	self.T_label = QLabel('T')
	self.proc_label = QLabel('Number of Processors') 
        self.col_label = QLabel('Input Value')
        self.col_label1 = QLabel('Slider Range')
        self.skip_label = QLabel('Remove nth Point in Exp Spec:')
        self.job_label = QLabel('Job Name')
  
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
	hbox.addWidget(self.DJ_label,9,1)
	hbox.addWidget(self.DJK_label,10,1)
	hbox.addWidget(self.DK_label,11,1)
	hbox.addWidget(self.dJ_label,12,1)
	hbox.addWidget(self.dK_label,13,1)
	
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

	
	hbox.addWidget(self.DJ_box,9,2)
	hbox.addWidget(self.DJK_box,10,2)
	hbox.addWidget(self.DK_box,11,2)	
	hbox.addWidget(self.dJ_box,12,2)	
	hbox.addWidget(self.dK_box,13,2)	
	
        hbox.addWidget(self.a_cb,2,4)
        hbox.addWidget(self.b_cb,3,4)
        hbox.addWidget(self.c_cb,4,4)
        hbox.addWidget(self.DJ_cb,9,4)
        hbox.addWidget(self.DJK_cb,10,4)
        hbox.addWidget(self.DK_cb,11,4)
        hbox.addWidget(self.dJ_cb,12,4)
        hbox.addWidget(self.dK_cb,13,4)	
        hbox.addWidget(self.fix_cb,3,7)		
        hbox.addWidget(self.peak_center_cb,6,6)	

	
	
	hbox.addWidget(self.ua_box,5,2)
        hbox.addWidget(self.ub_box,6,2)
        hbox.addWidget(self.uc_box,7,2)
        hbox.addWidget(self.T_box,8,2)	
        hbox.addWidget(self.f_lower_box,1,6)
        hbox.addWidget(self.f_upper_box,2,6)
        hbox.addWidget(self.I_lower_box,3,6)
        hbox.addWidget(self.I_upper_box,4,6)
        
                        
	hbox.addWidget(self.draw_button,0,1)

	#hbox.addWidget(slider_label,2,3)
	hbox.addWidget(self.peak_finder_button,6,5)
	hbox.addWidget(self.autofit_button,0,7)
        hbox.addWidget(self.autofit_cb,0,8)
        hbox.addWidget(self.peaks_button,0,9)
        hbox.addWidget(self.spectrum_file_button,0,5)
	hbox.addWidget(self.proc_label,2,7)
	hbox.addWidget(self.proc_box,2,8)
	hbox.addWidget(self.job_box,1,8)
	hbox.addWidget(self.lv,1,9,9,11)
	hbox.addWidget(self.delete_button,10,9)
	hbox.addWidget(self.set_fit_button,10,10)
	hbox.addWidget(self.set_check_button,10,11)
	hbox.addWidget(self.auto_set_button,0,10)
	hbox.addWidget(self.skip_box,5,6)
        vbox = QVBoxLayout()
        vbox.addLayout(hbox)
        
        self.main_frame.setLayout(vbox)



        self.setCentralWidget(self.main_frame)
    	
    def create_status_bar(self):
        self.status_text = QLabel("Autofit Setup")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        save_file_action = self.create_action("&Save input file",
            shortcut="Ctrl+S", slot=self.save_input_file, 
            tip="Save an input file")
            
        load_file_action = self.create_action("&Open input file",
            shortcut="Ctrl+O", slot=self.open_input_file, 
            tip="Open an input file")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (save_file_action,load_file_action, None, quit_action))
        
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
