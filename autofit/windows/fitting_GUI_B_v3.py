import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import math


def handle_mouse_press( event ):
    if event.dblclick:
        #print "Can haz dblclikz???"
        x = plt.ginput(1)
        #print x
        difference_list = []

        current_peaks = cat_reader()
        for entry in current_peaks:
            inten = entry[1]
            freq = entry[0]
            qnum_up = entry[2]
            qnum_low = entry[3]           
            difference_list.append((math.fabs(float(entry[0])-float(x[0][0])),freq,qnum_up,qnum_low,inten))
        difference_list.sort()
        marker = 0
        new_picked_list = []
        
        for z in picked_list:
            if difference_list[0][2]==z[1] and difference_list[0][3]==z[2]:
                marker = 1
            else:
                new_picked_list.append(z)
        del picked_list  
        global picked_list
        picked_list = new_picked_list
        if marker==0:
                  
            picked_list.append((difference_list[0][1],difference_list[0][2],difference_list[0][3],difference_list[0][4]))        
        data = cat_reader()
        
        s_b = []
        t_b = []
        for x in data:
            #if x[2]==difference_list[0][2] and x[3]==difference_list[0][3] and marker==0:
            #    s_b.append(0.0)
            #    s_b.append(float(10**float(difference_list[0][4]))) 
            #    s_b.append(0.0)
            #    t_b.append(float(difference_list[0][1])-0.0001)
            #    t_b.append(difference_list[0][1])
            #    t_b.append(float(difference_list[0][1])+0.0001)
            for y in picked_list:
                if x[2]==y[1] and x[3]==y[2]:
                    s_b.append(0.0)
                    s_b.append(str(10**float(x[1]))) 
                    s_b.append(0.0)
                    t_b.append(float(x[0])-0.0001)
                    t_b.append(x[0])
                    t_b.append(float(x[0])+0.0001)  
        picked_plt.set_data(t_b,s_b)
        plt.draw()
        #print s_b,t_b
            
            
        
        #print picked_list
def int_writer(u_A='1.0',u_B='0.0',u_C='0.0', J_min="00", J_max="20", inten="-10.0",Q_rot="300000",freq="25.8", temp="298"):#generates SPCAT input file
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


def var_writer_uncert(A,B,C,DJ='0.0',DJK='0.0',DK='0.0',dJ='0.0',dK='0.0'):#generates SPCAT input file
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
    a = subprocess.Popen('SPCAT default', stdout=subprocess.PIPE, shell=False)
    a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen
 
def cat_reader(freq_high=20000,freq_low=0): #reads output from SPCAT
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
                linelist.append((freq,inten, qnum_up, qnum_low,uncert))
    linelist.sort()
    fh.close()
    return linelist


def update(val):
    A_value = str(A_slider.val)
    B_value = str(B_slider.val)
    C_value = str(C_slider.val)  
    int_writer()
    var_writer_uncert(A_value,B_value,C_value)
    run_SPCAT()
    data = cat_reader()
    s = []
    t = []
    s_b = []
    t_b = []
    for x in data:
        s.append(0.0)
        s.append(str(10**float(x[1]))) 
        s.append(0.0)
        t.append(float(x[0])-0.0001)
        t.append(x[0])
        t.append(float(x[0])+0.0001)
        for y in picked_list:
            if x[2]==y[1] and x[3]==y[2]:
                s_b.append(0.0)
                s_b.append(str(10**float(x[1]))) 
                s_b.append(0.0)
                t_b.append(float(x[0])-0.0001)
                t_b.append(x[0])
                t_b.append(float(x[0])+0.0001)            
    #print picked_list


    picked_plt.set_data(t_b,s_b)
    l.set_data(t,s)
    plt.draw()


def reset(event):
    A_slider.reset()
    B_slider.reset()
    C_slider.reset()


def colorfunc(label):
    l.set_color(label)
    plt.draw()


def run_main(A,B,C,dA,dB,dC,last_times_picked_list):
    global picked_list
    picked_list = last_times_picked_list
    plt.close()
    global figure_h
    figure_h = plt.figure()

    int_writer()
    var_writer_uncert(A,B,C)
    run_SPCAT()
    data = cat_reader()
    global t
    global s
    t = []
    s = []
    t_b = []
    s_b = []    
    for x in data:
        s.append(0.0)
        s.append(str(10**float(x[1]))) 
        s.append(0.0)
        t.append(float(x[0])-0.0001)
        t.append(x[0])
        t.append(float(x[0])+0.0001)
        for y in picked_list:
            if x[2]==y[1] and x[3]==y[2]:
                s_b.append(0.0)
                s_b.append(str(10**float(x[1]))) 
                s_b.append(0.0)
                t_b.append(float(x[0])-0.0001)
                t_b.append(x[0])
                t_b.append(float(x[0])+0.0001)  
    ax = figure_h.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    
    a0 = 5
    f0 = 3
    global l
    l, = plt.plot(t,s, lw=2, color='red')
    global picked_plt
    picked_plt, = plt.plot(t_b,s_b,lw=2,color='black')
    #plt.axis([0, 1, -10, 10])
    figure_h.canvas.mpl_connect('button_press_event', handle_mouse_press)
    axcolor = 'lightgoldenrodyellow'
    axA = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    axB  = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axC  = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    global A_slider
    global B_slider
    global C_slider
    
    A_slider = Slider(axA, 'A', A-dA, A+dA, valinit=A)
    B_slider = Slider(axB, 'B', B-dB, B+dB, valinit=B)
    C_slider = Slider(axC, 'C', C-dC, C+dC, valinit=C)
    A_slider.on_changed(update)
    B_slider.on_changed(update)
    C_slider.on_changed(update)
    global button
    global radio
    resetax = plt.axes([0.1, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    button.on_clicked(reset)
    rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
    radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
    radio.on_clicked(colorfunc)
    plt.show()
