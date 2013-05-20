# -*- coding: utf-8 -*-
'''
Created on 09/04/2012
Outputs a two-column peakpick of an input spectrum given a minimum intensity threshold.
@author: Cristobal
'''
import os
from easygui import *
import sys
import re
import random
import string
from numpy import linalg as LA
import pylab
import math
from numpy import *
from scipy import *
import subprocess
import numpy as np
import subprocess
from math import *

f =loadtxt(fileopenbox(msg="Choose your Spectrum"))
spectrum=f
msg = "Threshold"
title = "Thresold"
fieldNames = [" "]
thr = []
thr = multenterbox(msg,title,fieldNames)
if thr == None:
    sys.exit(0);
print thr
peaks=[]
for i in range(0, len(spectrum)):
    if spectrum[i,1] > float(thr[0]) and spectrum[i,1] > spectrum[(i-1),1] and spectrum[i,1] > spectrum[(i+1),1]:
        peaks.append(spectrum[i])

print len(peaks)
peakpicks=zeros((len(peaks),2))
for i,row in enumerate(peaks):
    peakpicks[i,0]=row[0]
    peakpicks[i,1]=row[1]
savetxt(filesavebox(msg="Save your peakpick file"), peakpicks)
  