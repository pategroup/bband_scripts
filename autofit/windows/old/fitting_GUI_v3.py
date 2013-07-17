
import sys, os, random

import fitting_GUI_B_v3
import subprocess
from multiprocessing import Process
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui












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
	self.ua_box.setText('1.0')
	self.ub_box.setText('1.0')
	self.uc_box.setText('1.0')
	
	self.autofit_filename_box.setText('prog_A_GUI_v1.py')
	self.a_set  = '3000'
	self.b_set  = '2000'
	self.c_set  = '1000'      

    def replot_button(self):
        
        A = float(unicode(self.a_box.text()))
        B = float(unicode(self.b_box.text()))
        C = float(unicode(self.c_box.text()))
        dA = float(unicode(self.da_box.text()))
        dB = float(unicode(self.db_box.text()))
        dC = float(unicode(self.dc_box.text()))
        if self.counter != 0:
            picked_list = fitting_GUI_B_v3.picked_list
        
        else:
            picked_list = []
        self.counter+=1
	fitting_GUI_B_v3.run_main(A,B,C,dA,dB,dC,picked_list)



    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_about(self):
        msg = """ A demo of using PyQt with matplotlib:
        
         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        QMessageBox.about(self, "About the demo", msg.strip())
    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        print 'test'
        QMessageBox.information(self, "Click!", msg)


    def run_autofit(self):
	filename = unicode(self.autofit_filename_box.text())
	subprocess.Popen("python ./%s"%filename, shell=True)

    def peaks_button_push(self):
        print fitting_GUI_B_v3.picked_list
		



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

	
	
	self.ua_box = QLineEdit()
	self.ub_box = QLineEdit()
	self.uc_box = QLineEdit()

	self.autofit_filename_box = QLineEdit() 

        self.a_box.setMinimumWidth(100)
	self.b_box.setMinimumWidth(100)

        #self.connect(self.textbox2, SIGNAL('editingFinished ()'), self.on_draw)       
        #self.connect(self.textbox3, SIGNAL('editingFinished ()'), self.on_draw) 
        self.draw_button = QPushButton("&Plot")
	self.autofit_button = QPushButton("&Run Autofit")
        self.peaks_button = QPushButton("&Print Selected Peaks")


        self.connect(self.autofit_button, SIGNAL('clicked()'), self.run_autofit)
        self.connect(self.draw_button, SIGNAL('clicked()'), self.replot_button)
        self.connect(self.peaks_button, SIGNAL('clicked()'), self.peaks_button_push)
        
        
	a_label = QLabel('A')
	b_label = QLabel('B')
	c_label = QLabel('C')
	DJ_label = QLabel('DJ')
	DJK_label = QLabel('DJK')
	DK_label = QLabel('DK')
	ua_label = QLabel('mu a')
	ub_label = QLabel('mu b')
	uc_label = QLabel('mu c')
	autofit_filename_label = QLabel('Autofit Filename') 
        col_label = QLabel('Input Value')
        col_label1 = QLabel('Slider Range')


  
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
	hbox.addWidget(ua_label,7,1)
        hbox.addWidget(ub_label,8,1)
        hbox.addWidget(uc_label,9,1)
        hbox.addWidget(col_label,0,2)
        hbox.addWidget(col_label1,0,3)



	
        hbox.addWidget(self.a_box,1,2)
        hbox.addWidget(self.b_box,2,2)
	hbox.addWidget(self.c_box,3,2)
        hbox.addWidget(self.da_box,1,3)
        hbox.addWidget(self.db_box,2,3)
	hbox.addWidget(self.dc_box,3,3)

	
	hbox.addWidget(self.DJ_box,4,2)
	hbox.addWidget(self.DJK_box,5,2)
	hbox.addWidget(self.DK_box,6,2)	
	hbox.addWidget(self.dDJ_box,4,3)
	hbox.addWidget(self.dDJK_box,5,3)
	hbox.addWidget(self.dDK_box,6,3)
	
	
	
	
	hbox.addWidget(self.ua_box,7,2)
        hbox.addWidget(self.ub_box,8,2)
        hbox.addWidget(self.uc_box,9,2)
	
	hbox.addWidget(self.draw_button,0,4)

	#hbox.addWidget(slider_label,2,3)

	hbox.addWidget(self.autofit_button,4,4)
        hbox.addWidget(self.peaks_button,5,4)
	hbox.addWidget(autofit_filename_label,2,4)
	hbox.addWidget(self.autofit_filename_box,3,4)
        vbox = QVBoxLayout()
        vbox.addLayout(hbox)
        
        self.main_frame.setLayout(vbox)



        self.setCentralWidget(self.main_frame)
    	
    def create_status_bar(self):
        self.status_text = QLabel("Autofit Setup")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_file_action, None, quit_action))
        
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
