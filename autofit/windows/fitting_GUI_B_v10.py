import matplotlib
#matplotlib.use('TkAgg')

import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import math

import pylab as pl

global picked_list
picked_list = []




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import AxesWidget, Slider, Button
from matplotlib.transforms import Affine2D


class VSlider(AxesWidget):
    """
    A slider representing a floating point range

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *hline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None, slidermax=None,
                 dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*

        *valinit*
            The slider initial position

        *label*
            The slider label

        *valfmt*
            Used to format the slider value

        *closedmin* and *closedmax*
            Indicate whether the slider interval is closed

        *slidermin* and *slidermax*
            Used to constrain the value of this slider to the values
            of other sliders.

        additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...)
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin,valinit,0,1, **kwargs)

        self.hline = ax.axhline(valinit,0,1, color='r', lw=1)


        self.valfmt=valfmt
        ax.set_yticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_xticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.02, label, transform=ax.transAxes,
                             verticalalignment='bottom',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.02, valfmt%valinit,
                               transform=ax.transAxes,
                               verticalalignment='top',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active  = False

    def _update(self, event):
        'update the slider position'
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event')
             or (event.name == 'button_press_event' and event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt%val)
        if self.drawon: self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: return
        for cid, func in self.observers.iteritems():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        'remove the observer with connection id *cid*'
        try: del self.observers[cid]
        except KeyError: pass

    def reset(self):
        "reset the slider to the initial value if needed"
        if (self.val != self.valinit):
            self.set_val(self.valinit)











def on_key( event ):
    if event.key=="q" or event.key=="Q":
        check.lines[0][0].set_visible(True)
        check.lines[0][1].set_visible(True)
        check.lines[1][0].set_visible(False)
        check.lines[1][1].set_visible(False)
        check.lines[2][0].set_visible(False)
        check.lines[2][1].set_visible(False)
        check.lines[3][0].set_visible(False)
        check.lines[3][1].set_visible(False)
        check.labels[0].set_text('')
        check.labels[1].set_text('')
        check.labels[2].set_text('')
        check.labels[3].set_text('')
        #print "Can haz dblclikz???"
        x = plt.ginput(1)
        #print x
        difference_list = []
        new_text = ""
        #new_text1 = "Other Nearby Peaks \n"
        current_peaks = cat_reader()
        for entry in current_peaks:
            inten = entry[1]
            freq = entry[0]
            qnum_up = entry[2]
            qnum_low = entry[3]
            uncert = entry[4]
            if math.fabs(float(entry[0])-float(x[0][0]))<100.0 and math.fabs((10**float(entry[1])-float(x[0][1]))/(10**float(entry[1])))<0.2:
                difference_list.append((math.fabs(float(entry[0])-float(x[0][0])),freq,qnum_up,qnum_low,inten,uncert))
        difference_list.sort()
        #print difference_list
        marker = 0
        new_picked_list = []
        try:
            for z in picked_list:
                if difference_list[0][2]==z[1] and difference_list[0][3]==z[2]:
                    new_text+="\n removed: \n "+str(difference_list[0][2])+" | "+str(difference_list[0][3])+' | '+str(difference_list[0][1])+' | '+str(difference_list[0][4])
                    marker = 1
                    check.labels[0].set_text(str(difference_list[0][2])+" | "+str(difference_list[0][3])+' | '+str(difference_list[0][1])+' | '+str(difference_list[0][4]))
                    check.labels[0].set_fontsize(8)
                    if len(difference_list)>1:
                        #new_text1+="\n"+str(difference_list[1][2])+' - '+str(difference_list[1][3])+" "+str(difference_list[1][1])+' '+str(difference_list[1][4])
                        check.labels[1].set_text(str(difference_list[1][2])+' | '+str(difference_list[1][3])+" | "+str(difference_list[1][1])+' | '+str(difference_list[1][4]))
                        check.labels[1].set_fontsize(8)
                        if len(difference_list)>2:
                            #new_text1+="\n"+str(difference_list[2][2])+' - '+str(difference_list[2][3])+" "+str(difference_list[2][1])+' '+str(difference_list[2][4])
                            check.labels[2].set_text(str(difference_list[2][2])+' | '+str(difference_list[2][3])+" | "+str(difference_list[2][1])+' | '+str(difference_list[2][4]))
                            check.labels[2].set_fontsize(8)
                            if len(difference_list)>3:
                                check.labels[3].set_text(str(difference_list[3][2])+' | '+str(difference_list[3][3])+" | "+str(difference_list[3][1])+' | '+str(difference_list[3][4]))
                                check.labels[3].set_fontsize(8)                           
                                #new_text1+="\n"+str(difference_list[3][2])+' - '+str(difference_list[3][3])+" "+str(difference_list[3][1])+' '+str(difference_list[3][4])
                    #if len(difference_list)>1:
                    #   if math.fabs(float(difference_list[0][1])-float(difference_list[1][1]))<0.50:
                    #      picked_list.append((difference_list[1][1],difference_list[1][2],difference_list[1][3],difference_list[1][4]))  
                    #     new_text+="\n removed "+str(difference_list[0][2])+"-"+str(difference_list[0][3])+'   '+str(difference_list[0][4])
                            #print "test1"
                else:
                    new_picked_list.append(z)
            del picked_list  
            global picked_list
            picked_list = new_picked_list
            if marker==0:
                
                check.labels[0].set_text(str(difference_list[0][2])+" | "+str(difference_list[0][3])+' | '+str(difference_list[0][1])+' | '+str(difference_list[0][4]))
                check.labels[0].set_fontsize(8)
                    
                picked_list.append((difference_list[0][1],difference_list[0][2],difference_list[0][3],difference_list[0][4],difference_list[0][5],'CHECK'))  
                new_text+="\n added: \n "+str(difference_list[0][2])+" | "+str(difference_list[0][3])+' | '+str(difference_list[0][1])+' | '+str(difference_list[0][4])
                
                if len(difference_list)>1:
                    #new_text1+="\n"+str(difference_list[1][2])+' - '+str(difference_list[1][3])+" "+str(difference_list[1][1])+' '+str(difference_list[1][4])
                    check.labels[1].set_text(str(difference_list[1][2])+' | '+str(difference_list[1][3])+" | "+str(difference_list[1][1])+' | '+str(difference_list[1][4]))
                    check.labels[1].set_fontsize(8)
                    
                    
                    
                    if len(difference_list)>2:
                        #new_text1+="\n"+str(difference_list[2][2])+' - '+str(difference_list[2][3])+" "+str(difference_list[2][1])+' '+str(difference_list[2][4])
                        check.labels[2].set_text(str(difference_list[2][2])+' | '+str(difference_list[2][3])+" | "+str(difference_list[2][1])+' | '+str(difference_list[2][4]))
                        check.labels[2].set_fontsize(8)
                        if len(difference_list)>3:
                            check.labels[3].set_text(str(difference_list[3][2])+' | '+str(difference_list[3][3])+" | "+str(difference_list[3][1])+' | '+str(difference_list[3][4]))
                            check.labels[3].set_fontsize(8) 
                            #new_text1+="\n"+str(difference_list[3][2])+' - '+str(difference_list[3][3])+" "+str(difference_list[3][1])+' '+str(difference_list[3][4])
                    #if float(difference_list[0][1])-float(difference_list[1][1])<0.010 and difference_list[0][4]==difference_list[1][4]:
                    #   picked_list.append((difference_list[1][1],difference_list[1][2],difference_list[1][3],difference_list[1][4]))    
                    #  new_text+="\n added "+str(difference_list[1][2])+"-"+str(difference_list[1][3])+'   '+str(difference_list[1][4])
            global data
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
            text_box.set_text(new_text)
            text_box.set_fontsize(10)
            #text_box2.set_text(new_text1)
            #picked_plt.set_xlim([f_lower_g,f_upper_g])
            #print new_text
            plt.draw()
            #print s_b,t_b
        except IndexError:
            pass
    
            
        
         #print picked_list
def int_writer(u_A='1.0',u_B='0.0',u_C='0.0', J_min="00", J_max='20', inten="-10.0",Q_rot="300000",freq="25.8", temp="298"):#generates SPCAT input file
    input_file = ""
    J_max = J_max_g
    #print "freq_max=",freq
    input_file += "Molecule \n"
    input_file += "0  91  %s  %s  %s  %s  %s %s  %s\n"%(Q_rot, J_min, J_max,inten,inten,freq, temp)
    input_file += " 001  %s \n" % u_A
    input_file += " 002  %s \n" % u_B
    input_file += " 003  %s \n" % u_C
    fh_int = open("default.int", "w")
    fh_int.write(input_file)
    fh_int.close()

def triples_gen(peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high,trans_1,trans_2,trans_3):


    trans_1_peaks = []
    trans_2_peaks = []
    trans_3_peaks = []
    trans_1_center = trans_1
    #print trans_1_center,'trans_1_center' 
    trans_2_center = trans_2
    trans_3_center = trans_3
    for freq_p, inten_p in peaklist:
        #print abs(float(trans_1_center)-float(freq_p)),peak_1_uncertainty
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
    input_file += "           40100  %s %s \n" %(dJ, ddJ)
    input_file += "           41000  %s %s \n" %(dK, ddK)
    fh_var = open("default.var",'w')
    fh_var.write(input_file)
    fh_var.close()





def run_SPCAT(): 
    a = subprocess.Popen('SPCAT default', stdout=subprocess.PIPE, shell=False)
    a.stdout.read()#seems to be best way to get SPCAT to finish. I tried .wait(), but it outputted everything to screen
 
def cat_reader(): #reads output from SPCAT
    fh = open("default.cat")
    freq_high=float(f_upper_g)
    #print freq_high,freq_low
    
    freq_low=float(f_lower_g)
    linelist = []
    for line in fh:
        if line[8:9]==".": 

            freq = line[3:13]
            inten = line[22:29]
            qnum_up = line[55:61]
            qnum_low = line[67:73]
            uncert = line[13:21]
            #if float(freq)> freq_low and float(freq)<freq_high:#<<<<<<<<<<<<<<<<<<<<
            linelist.append((freq,inten, qnum_up, qnum_low,uncert))
    linelist.sort()
    fh.close()
    return linelist


def update(val):
    A_value = str(A_slider.val)
    B_value = str(B_slider.val)
    C_value = str(C_slider.val)  
    ua_value = str(ua_slider.val) 
    ub_value = str(ub_slider.val) 
    uc_value = str(uc_slider.val) 
    int_writer(ua_value,ub_value,uc_value,temp=T_g,freq=(f_upper_g/1000))
    var_writer_uncert(A_value,B_value,C_value,DJ_g,DJK_g,DK_g,dJ_g,dK_g)
    run_SPCAT()
    



    
    data = cat_reader()
    s = []
    t = []
    s_b = []
    t_b = []
    s_c = []
    t_c = []
    for x in data:
        s.append(0.0)
        s.append(str(10**float(x[1]))) 
        s.append(0.0)
        t.append(float(x[0])-0.0001)
        t.append(x[0])
        t.append(float(x[0])+0.0001)
        if len(trans_1)!=0 and trans_1[0]==x[2] and trans_1[1]==x[3]and triples_flag==1:
            s_c.append(str(10**float(x[1]))) 
            t_c.append(x[0])
            trans_1_center = (x[0])
        if len(trans_2)!=0 and trans_2[0]==x[2] and trans_2[1]==x[3]and triples_flag==1:
            s_c.append(str(10**float(x[1]))) 
            t_c.append(x[0])
            trans_2_center = (x[0])
        if len(trans_3)!=0 and trans_3[0]==x[2] and trans_3[1]==x[3]and triples_flag==1:
            s_c.append(str(10**float(x[1]))) 
            t_c.append(x[0])
            trans_3_center = (x[0])
        for y in picked_list:
            if x[2]==y[1] and x[3]==y[2]:
                s_b.append(0.0)
                s_b.append(str(10**float(x[1]))) 
                s_b.append(0.0)
                t_b.append(float(x[0])-0.0001)
                t_b.append(x[0])
                t_b.append(float(x[0])+0.0001)            
    #print picked_list
    if triples_flag==1:
        
        (trans_1_peaks,trans_2_peaks,trans_3_peaks,num_of_triples ) =triples_gen(peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peak_list,freq_low,freq_high,trans_1_center,trans_2_center,trans_3_center)
        print trans_1_peaks,trans_2_peaks,trans_3_peaks
        I_1 = []
        F_1 = []
        I_2 = []
        F_2 = []
        I_3 = []
        F_3 = []
        
        for line in trans_1_peaks:
            
            F_1.append(line[0])
            I_1.append(line[1])
        for line in trans_2_peaks:
            F_2.append(line[0])
            I_2.append(line[1])
        for line in trans_3_peaks:
            F_3.append(line[0])
            I_3.append(line[1])
        trans_1_plt.set_data(F_1,I_1)
        trans_2_plt.set_data(F_2,I_2)
        trans_3_plt.set_data(F_3,I_3)
        triples_plt.set_data(t_c,s_c)







 
    picked_plt.set_data(t_b,s_b)
    l.set_data(t,s)
    plt.draw()
def update_plot():
    A_value = str(A_slider.val)
    B_value = str(B_slider.val)
    C_value = str(C_slider.val)  
    ua_value = str(ua_slider.val) 
    ub_value = str(ub_slider.val) 
    uc_value = str(uc_slider.val) 
    int_writer(ua_value,ub_value,uc_value,temp=T_g,freq=(f_upper_g/1000))
    var_writer_uncert(A_value,B_value,C_value,DJ_g,DJK_g,DK_g,dJ_g,dK_g)
    run_SPCAT()
    



    
    data = cat_reader()
    s = []
    t = []
    s_b = []
    t_b = []
    s_c = []
    t_c = []
    for x in data:
        s.append(0.0)
        s.append(str(10**float(x[1]))) 
        s.append(0.0)
        t.append(float(x[0])-0.0001)
        t.append(x[0])
        t.append(float(x[0])+0.0001)
        if len(trans_1)!=0 and trans_1[0]==x[2] and trans_1[1]==x[3]and triples_flag==1:
            s_c.append(str(10**float(x[1]))) 
            t_c.append(x[0])
            trans_1_center = (x[0])
        if len(trans_2)!=0 and trans_2[0]==x[2] and trans_2[1]==x[3]and triples_flag==1:
            s_c.append(str(10**float(x[1]))) 
            t_c.append(x[0])
            trans_2_center = (x[0])
        if len(trans_3)!=0 and trans_3[0]==x[2] and trans_3[1]==x[3]and triples_flag==1:
            s_c.append(str(10**float(x[1]))) 
            t_c.append(x[0])
            trans_3_center = (x[0])
        for y in picked_list:
            if x[2]==y[1] and x[3]==y[2]:
                s_b.append(0.0)
                s_b.append(str(10**float(x[1]))) 
                s_b.append(0.0)
                t_b.append(float(x[0])-0.0001)
                t_b.append(x[0])
                t_b.append(float(x[0])+0.0001)            
    #print picked_list
    if triples_flag==1:
        
        (trans_1_peaks,trans_2_peaks,trans_3_peaks,num_of_triples ) =triples_gen(peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peak_list,freq_low,freq_high,trans_1_center,trans_2_center,trans_3_center)
        print trans_1_peaks,trans_2_peaks,trans_3_peaks
        I_1 = []
        F_1 = []
        I_2 = []
        F_2 = []
        I_3 = []
        F_3 = []
        
        for line in trans_1_peaks:
            
            F_1.append(line[0])
            I_1.append(line[1])
        for line in trans_2_peaks:
            F_2.append(line[0])
            I_2.append(line[1])
        for line in trans_3_peaks:
            F_3.append(line[0])
            I_3.append(line[1])
        trans_1_plt.set_data(F_1,I_1)
        trans_2_plt.set_data(F_2,I_2)
        trans_3_plt.set_data(F_3,I_3)
        triples_plt.set_data(t_c,s_c)







 
    picked_plt.set_data(t_b,s_b)
    l.set_data(t,s)
    plt.draw()
def reset(event):
    A_slider.reset()
    B_slider.reset()
    C_slider.reset()
    ua_slider.reset()
    ub_slider.reset()
    uc_slider.reset()

def colorfunc(label):
    l.set_color(label)
    plt.draw()

def func(label):
    try:
        #print 'test'
        check_select = label.split(" | ")
        print label
        new_picked_list = []
        marker = 0
        new_text = ''
        for line in picked_list:
            if line[1]==check_select[0] and line[2]==check_select[1]:
                new_text+="\n removed: \n "+label
                
                marker =1
            else: 
                new_picked_list.append(line)
                
        if marker == 0:
            for line in data:
                if line[2]==check_select[0] and line[3]==check_select[1]:
                    new_picked_list.append((check_select[2],check_select[0],check_select[1],check_select[3],line[4],'CHECK'))
            new_text+="\n added: \n "+label
        global picked_list
        text_box.set_text(new_text)
        text_box.set_fontsize(10)
        picked_list = new_picked_list
        plt.draw()
    except IndexError:
        pass
    update_plot()
def run_main(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,dua,dub,duc,f_lower,f_upper,T,J_max,last_times_picked_list):
    #matplotlib.use('TkAggx')
    
    global DJ_g,DJK_g,DK_g,ua_g,ub_g,uc_g,T_g,dJ_g,dK_g,f_lower_g,f_upper_g,J_max_g,triples_flag
    triples_flag = 0
    J_max_g = J_max
    DJ_g = DJ
    DJK_g = DJK
    DK_g = DK
    ua_g = ua
    ub_g = ub
    uc_g = uc
    T_g=T
    dJ_g = dJ
    dK_g = dK
    f_lower_g = f_lower
    f_upper_g = f_upper
    global trans_1,trans_2,trans_3,peak_list
    trans_1 = []
    trans_2 = []
    trans_3 = []
    
    peak_list = []
    
    
    
    global picked_list
    picked_list = last_times_picked_list
    plt.close()
    global figure_h
    figure_h = plt.figure(figsize=(16.5, 4))
    try:
        figure_h.canvas.manager.window.Move((00,00))
    except AttributeError:
        pass
    #thismanager = pl.get_current_fig_manager()
    #thismanager.window.SetPosition((00, 0))
    #thismanager.window.wm_geometry("+00+0")
    


    int_writer(ua,ub,uc,temp=T_g,freq=(f_upper_g/1000))
    var_writer_uncert(A,B,C,DJ,DJK,DK,dJ_g,dK_g)
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
    ax = figure_h.add_subplot(212)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    pl.subplots_adjust( hspace=0.0,right=0.925 )

    
    
    a0 = 5
    f0 = 3
    global l,triples_plt
    l, = plt.plot(t,s, lw=2, color='red')
    triples_plt, = plt.plot([],[], '*',markersize=15,color='blue')
    ax.set_xlim([f_lower_g,f_upper_g])
    global picked_plt, ax2
    picked_plt, = plt.plot(t_b,s_b,lw=2,color='black')
    #plt.axis([0, 1, -10, 10])
    #figure_h.canvas.mpl_connect('button_press_event', handle_mouse_press)
    ax2 = figure_h.add_subplot(211,sharex=ax)
    ax2.axes.get_xaxis().set_visible(False)

    global peak_list_plt,exp_plt,trans_1_plt,trans_2_plt,trans_3_plt
    peak_list_plt, = ax2.plot([],[],'o',color='red')
    trans_1_plt, = ax2.plot([],[],'s',color='cyan')
    trans_2_plt, = ax2.plot([],[],'s',color='magenta')
    trans_3_plt, = ax2.plot([],[],'s',color='yellow')
    ax2.set_xlim([f_lower_g,f_upper_g])
    global locator   
    locator = ax2.yaxis.get_major_locator() 

    exp_plt, = ax2.plot([],[],lw=2,color='black')




    figure_h.canvas.mpl_connect('key_press_event', on_key)
    axcolor = 'lightgoldenrodyellow'
    axA = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    axB  = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axC  = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    
    axua = plt.axes([.935, 0.32, 0.01, 0.6], axisbg=axcolor)
    axub = plt.axes([.96, 0.32, 0.01, 0.6], axisbg=axcolor)
    axuc = plt.axes([.985, 0.32, 0.01, 0.6], axisbg=axcolor)
    #axub  = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    #axuc  = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    dua_2 = dua
    dub_2 = dub
    duc_2 = duc
    if (ua_g-dua)<0:
       dua_2=ua_g 
    if (ub_g-dub)<0:
       dub_2=ub_g 
    if (uc_g-duc)<0:
       duc_2=uc_g 
    global ua_slider
    ua_slider = VSlider(axua, 'u a', ua_g-dua_2, ua_g+dua, valinit=ua_g)
    ua_slider.on_changed(update)
    global ub_slider
    ub_slider = VSlider(axub, 'u b', ub_g-dub_2, ub_g+dub, valinit=ub_g)
    ub_slider.on_changed(update)
    global uc_slider
    uc_slider = VSlider(axuc, 'u c', uc_g-duc_2, uc_g+duc, valinit=uc_g)
    uc_slider.on_changed(update)
    global A_slider
    global B_slider
    global C_slider
    global rax
    rax = plt.axes([0.0, 0.5, 0.19, 0.4])
    global check
    check = CheckButtons(rax, ('','','',''), (True, False, False,False))
    

    check.on_clicked(func)
    
    
    
    
    A_slider = Slider(axA, 'A', A-dA, A+dA, valinit=A)
    B_slider = Slider(axB, 'B', B-dB, B+dB, valinit=B)
    C_slider = Slider(axC, 'C', C-dC, C+dC, valinit=C)
    
    
    
    A_slider.on_changed(update)
    B_slider.on_changed(update)
    C_slider.on_changed(update)
    global button
    global radio
    resetax = plt.axes([0.1, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset Sliders', color=axcolor, hovercolor='0.975')

    button.on_clicked(reset)
    #rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
    #radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
    #radio.on_clicked(colorfunc)
    global text_box
    #global text_box2
    text_box = plt.text(-1,8, "")
    text_box2 = plt.text(-1,23, "Refine Mouse Selection:                             Select transitions by pushing 'q' and then clicking in the predicted spectrum ")
    #plt.show()

def initialize(A,B,C,dA,dB,dC,DJ,DJK,DK,dJ,dK,ua,ub,uc,f_lower,f_upper,T,J_max,last_times_picked_list):
    #matplotlib.use('TkAggx')
    global DJ_g,DJK_g,DK_g,ua_g,ub_g,uc_g,T_g,dJ_g,dK_g,f_lower_g,f_upper_g,J_max_g
    J_max_g = J_max
    DJ_g = DJ
    DJK_g = DJK
    DK_g = DK
    ua_g = ua
    ub_g = ub
    uc_g = uc
    T_g=T
    dJ_g = dJ
    dK_g = dK
    f_lower_g = f_lower
    f_upper_g = f_upper
    
    
    global picked_list
    picked_list = last_times_picked_list
    plt.close()
    global figure_h
    figure_h = plt.figure(figsize=(16.5, 4))
    try:
        figure_h.canvas.manager.window.Move((00,00))
    except AttributeError:
        pass
    #thismanager = pl.get_current_fig_manager()
    #thismanager.window.SetPosition((00, 0))
    #thismanager.window.wm_geometry("+00+0")
    


    int_writer(ua,ub,uc,temp=T_g,freq=(f_upper_g/1000))
    var_writer_uncert(A,B,C,DJ,DJK,DK,dJ_g,dK_g)
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
    ax = figure_h.add_subplot(212)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    pl.subplots_adjust( hspace=0.0,right=0.97 )

    
    
    a0 = 5
    f0 = 3
    global l,triples_plt
    l, = plt.plot(t,s, lw=2, color='red')
    triples_plt, = plt.plot([],[], '*',markersize=15,color='blue')
    ax.set_xlim([f_lower_g,f_upper_g])
    global picked_plt, ax2
    picked_plt, = plt.plot(t_b,s_b,lw=2,color='black')
    #plt.axis([0, 1, -10, 10])
    #figure_h.canvas.mpl_connect('button_press_event', handle_mouse_press)
    ax2 = figure_h.add_subplot(211,sharex=ax)
    ax2.axes.get_xaxis().set_visible(False)
    
    global peak_list_plt,exp_plt,trans_1_plt,trans_2_plt,trans_3_plt
    peak_list_plt, = ax2.plot([],[],'o',color='red')
    trans_1_plt, = ax2.plot([],[],'s',color='cyan')
    trans_2_plt, = ax2.plot([],[],'s',color='magenta')
    trans_3_plt, = ax2.plot([],[],'s',color='yellow')
    ax2.set_xlim([f_lower_g,f_upper_g])
    global locator   
    locator = ax2.yaxis.get_major_locator() 

    exp_plt, = ax2.plot([],[],lw=2,color='black')

    global peak_1_uncertainty,peak_2_uncertainty,peak_3_uncertainty,peaklist,freq_low,freq_high


    figure_h.canvas.mpl_connect('key_press_event', on_key)
    axcolor = 'lightgoldenrodyellow'
    axA = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    axB  = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axC  = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    
    axua = plt.axes([0.03, 0.22, 0.1, 0.03], axisbg=axcolor)
    axub = plt.axes([0.03, 0.17, 0.1, 0.03], axisbg=axcolor)
    axuc = plt.axes([0.03, 0.12, 0.1, 0.03], axisbg=axcolor)
    #axub  = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    #axuc  = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    global ua_slider
    ua_slider = Slider(axua, 'mu a', ua_g-1, ua_g+1, valinit=ua_g)
    ua_slider.on_changed(update)
    global ub_slider
    ub_slider = Slider(axub, 'mu b', ub_g-1, ub_g+1, valinit=ub_g)
    ub_slider.on_changed(update)
    global uc_slider
    uc_slider = Slider(axuc, 'mu c', uc_g-1, uc_g+1, valinit=uc_g)
    uc_slider.on_changed(update)
    global A_slider
    global B_slider
    global C_slider
    global rax
    rax = plt.axes([0.0, 0.5, 0.19, 0.4])
    global check
    check = CheckButtons(rax, ('','','',''), (True, False, False,False))
    

    check.on_clicked(func)
    
    
    
    
    A_slider = Slider(axA, 'A', A-dA, A+dA, valinit=A)
    B_slider = Slider(axB, 'B', B-dB, B+dB, valinit=B)
    C_slider = Slider(axC, 'C', C-dC, C+dC, valinit=C)
    
    
    
    A_slider.on_changed(update)
    B_slider.on_changed(update)
    C_slider.on_changed(update)
    global button
    global radio
    resetax = plt.axes([0.1, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset Sliders', color=axcolor, hovercolor='0.975')

    button.on_clicked(reset)
    #rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
    #radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
    #radio.on_clicked(colorfunc)
    global text_box
    #global text_box2
    text_box = plt.text(-1,8, "")
    text_box2 = plt.text(-1,23, "Refine Mouse Selection:                             Select transitions by pushing 'q' and then clicking in the predicted spectrum ")
    
