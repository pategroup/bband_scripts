'''
Created on 29/03/2012

@author: Cristobal
'''
import os
from easygui import *
import sys
import re
import random
import string
from numpy import *
from numpy import linalg as LA
import pylab
import math
from matplotlib.pyplot import figure, show
from numpy import *
from scipy import *
import scipy
import scipy.fftpack
import pylab
from scipy import pi
from matplotlib import pyplot as plt

srate=25E9

#timedomain = loadtxt(fileopenbox(msg="Choose your Time Domain File"))
f=open(fileopenbox(msg="Choose your Time Domain File"))
timedomain=[]
for row in f:
    temp=row.split()
    #print row
    timedomain.append(float(temp[size(temp)-1]))

ndata=size(timedomain)

window=kaiser(ndata,9.5)	

for i, row in enumerate(timedomain):
    timedomain[i]=row*window[i]

temp=zeros(ndata*2)
temp[0:ndata]=timedomain
timepad=temp

fft=abs(scipy.fftpack.fft(timepad))/100
freq=(scipy.fftpack.fftfreq(size(timepad),1/srate))/1E6

print size(fft)
print size(freq)
freqs=[]
ft=[]
for i, row in enumerate(freq):
        if row>2000 and row<10000:
                freqs.append(row)
                ft.append(fft[i]) 

out=zeros((size(freqs),2))
#gain=loadtxt(fileopenbox(msg="Choose your Gain Correction Curve"))

for i,row in enumerate(freqs):
    out[i,0]=row
    out[i,1]=ft[i]#*gain[i,1]

f = savetxt(filesavebox(msg="Save FFT"),out)
plt.plot(freqs,ft, 'blue',linewidth=0.8)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (mV)')
plt.title('Spectrum')
plt.grid(True)
plt.show()

