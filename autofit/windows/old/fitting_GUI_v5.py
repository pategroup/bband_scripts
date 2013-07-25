
import sys, os, random

import fitting_GUI_B_v5
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
            picked_list = fitting_GUI_B_v5.picked_list
        
        else:
            picked_list = []
        self.counter+=1
        try:
            try:
                
                fitting_GUI_B_v5.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,picked_list,int(unicode(self.skip_box.text())),self.filename,self.peaklist)
            except AttributeError:
                fitting_GUI_B_v5.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,picked_list,int(unicode(self.skip_box.text())),self.filename)
        except AttributeError:
            fitting_GUI_B_v5.run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,picked_list,int(unicode(self.skip_box.text())))


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
    
        new_length = int(math.floor(scale_factor*len(spectrum)))
    
        tck = splrep(x,y,s=0)
        xnew = numpy.arange(x[0],x[-1],new_resolution)
        ynew = splev(xnew,tck,der=0)
    
        output_spectrum = numpy.zeros((new_length,2))
        for i in range(0, new_length):
            output_spectrum[i,0] = xnew[i]
            output_spectrum[i,1] = ynew[i]
    
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

    def peaks_button_push(self):
        print fitting_GUI_B_v5.picked_list
        




    def create_main_frame(self):
        self.counter = 0
        self.main_frame = QWidget()
        
        self.a_box = QLineEdit()
	self.b_box = QLineEdit()
	self.c_box = QLineEdit()
        self.da_box = QLineEdit()
	self.db_box = QLineEdit()
	self.dc_box = QLineEdit()

        self.DJ_box = QLineEdit()
	self.DJK_box = QLineEdit()
	self.DK_box = QLineEdit()
        self.dDJ_box = QLineEdit()
	self.dDJK_box = QLineEdit()
	self.dDK_box = QLineEdit()
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
        self.peaks_button = QPushButton("&Print Selected Peaks")
        self.spectrum_file_button = QPushButton("&Open Exp Spectrum File")
        self.peak_finder_button = QPushButton("&Find Peaks in Exp Spectrum")
        self.connect(self.autofit_button, SIGNAL('clicked()'), self.run_autofit)
        self.connect(self.draw_button, SIGNAL('clicked()'), self.replot_button)
        self.connect(self.peaks_button, SIGNAL('clicked()'), self.peaks_button_push)
        self.connect(self.peak_finder_button, SIGNAL('clicked()'), self.peak_pick_button)
        self.connect(self.spectrum_file_button, SIGNAL('clicked()'), self.selectFile_spectrum)
        global a_label
	a_label = QLabel('A')
	b_label = QLabel('B')
	c_label = QLabel('C')
	DJ_label = QLabel('-DJ')
	DJK_label = QLabel('-DJK')
	DK_label = QLabel('-DK')
	dJ_label = QLabel('dJ')
	dK_label = QLabel('dK')
	ua_label = QLabel('mu a')
	ub_label = QLabel('mu b')
	uc_label = QLabel('mu c')
	f_upper_label = QLabel('Upper Frequency Limit')
	f_lower_label = QLabel('Lower Frequency Limit')
	
	
	I_upper_label = QLabel('Upper Intensity Limit')
	I_lower_label = QLabel('Lower Intensity Limit')
	T_label = QLabel('T')
	proc_label = QLabel('Number of Processors') 
        col_label = QLabel('Input Value')
        col_label1 = QLabel('Slider Range')
        skip_label = QLabel('Remove nth Point in Exp Spec:')
        job_label = QLabel('Job Name')
  
	hbox = QtGui.QGridLayout()
        #hbox = QHBoxLayout()
        #widget_list = [ a_label, self.textbox,b_label,self.textbox2 ,self.draw_button, self.grid_cb,
        #            slider_label, self.slider]
        #for w in range(len(widget_list)):
        #    hbox.addWidget(widget_list[w])
        #    hbox.setAlignment(w, Qt.AlignVCenter)
        hbox.addWidget(a_label,1,1)
        hbox.addWidget(b_label,2,1)
	hbox.addWidget(c_label,3,1)
	hbox.addWidget(DJ_label,4,1)
	hbox.addWidget(DJK_label,5,1)
	hbox.addWidget(DK_label,6,1)
	hbox.addWidget(dJ_label,7,1)
	hbox.addWidget(dK_label,8,1)
	
	hbox.addWidget(ua_label,9,1)
        hbox.addWidget(ub_label,10,1)
        hbox.addWidget(uc_label,11,1)
        hbox.addWidget(T_label,12,1)
        hbox.addWidget(f_lower_label,13,1)
        hbox.addWidget(f_upper_label,14,1)
        hbox.addWidget(I_lower_label,15,1)
        hbox.addWidget(I_upper_label,16,1)
        hbox.addWidget(job_label,9,4)

        
        hbox.addWidget(col_label,0,2)
        hbox.addWidget(col_label1,0,3)
        hbox.addWidget(skip_label,7,4)


	
        hbox.addWidget(self.a_box,1,2)
        hbox.addWidget(self.b_box,2,2)
	hbox.addWidget(self.c_box,3,2)
        hbox.addWidget(self.da_box,1,3)
        hbox.addWidget(self.db_box,2,3)
	hbox.addWidget(self.dc_box,3,3)

	
	hbox.addWidget(self.DJ_box,4,2)
	hbox.addWidget(self.DJK_box,5,2)
	hbox.addWidget(self.DK_box,6,2)	
	hbox.addWidget(self.dJ_box,7,2)	
	hbox.addWidget(self.dK_box,8,2)	
	
			
	hbox.addWidget(self.dDJ_box,4,3)
	hbox.addWidget(self.dDJK_box,5,3)
	hbox.addWidget(self.dDK_box,6,3)

	
	
	hbox.addWidget(self.ua_box,9,2)
        hbox.addWidget(self.ub_box,10,2)
        hbox.addWidget(self.uc_box,11,2)
        hbox.addWidget(self.T_box,12,2)	
        hbox.addWidget(self.f_lower_box,13,2)
        hbox.addWidget(self.f_upper_box,14,2)
        hbox.addWidget(self.I_lower_box,15,2)
        hbox.addWidget(self.I_upper_box,16,2)
        
                        
	hbox.addWidget(self.draw_button,0,4)

	#hbox.addWidget(slider_label,2,3)
	hbox.addWidget(self.peak_finder_button,1,4)
	hbox.addWidget(self.autofit_button,4,4)
        hbox.addWidget(self.peaks_button,5,4)
        hbox.addWidget(self.spectrum_file_button,6,4)
	hbox.addWidget(proc_label,2,4)
	hbox.addWidget(self.proc_box,3,4)
	hbox.addWidget(self.job_box,10,4)
	
	hbox.addWidget(self.skip_box,8,4)
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
