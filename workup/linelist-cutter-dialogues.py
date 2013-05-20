'''
Created on 02/04/2012
Spectrum cutter: Given a linelist and a spectrum, this program will cut all lines from the linelist out of the spectrum with a constant frequency-space cut width (default: 350 kHz)
@author: Cristobal
'''

from numpy import*
import numpy as np
from matplotlib import pyplot as pp
import math
import os
from easygui import *
import sys
import re
import random
import string
from numpy import linalg as LA
import pylab
import math
from matplotlib.pyplot import figure, show
from scipy import *
import pylab
from scipy import pi

spectrum = np.loadtxt(fileopenbox(msg="Choose your spectrum"))
linelist = np.loadtxt(fileopenbox(msg="Choose your linelist"))
linelist = sorted(linelist)
width = 0.6

specstart = spectrum[0,0]
specstep = (spectrum[10000,0]-spectrum[0,0])/10000
def Create_Mask():
    mask = np.zeros(np.shape(spectrum)[0])
    for i,line in enumerate(spectrum):
        mask[i] = 1
    for j,line in enumerate(linelist):
        freq1 = linelist[j]-(width/2)
        freq2 = linelist[j]+(width/2)
        i1 = math.floor((freq1-specstart)/specstep)+1
        i2 = math.floor((freq2-specstart)/specstep)-1
        i1=int(i1)
        i2=int(i2)
        for i in range(i1,i2):
            mask[i] = 0
    return mask

Mask = np.column_stack((spectrum[:,0],Create_Mask()))
#pp.plot(spectrum[:,0],Create_Mask())
applymask = spectrum
for i,line in enumerate(spectrum):
    applymask[i,1] = applymask[i,1]*Mask[i,1]

np.savetxt(filesavebox(msg="Save your Cut Spectrum"),applymask)
print("Idiot check: saved that shit successfully.")
pp.plot(applymask[:,0],applymask[:,1])
pp.show()
