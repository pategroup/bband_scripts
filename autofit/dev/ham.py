# Asymmetric top Hamiltonian program
# to replace SPCAT in Autofit
#
#
# Author: Nathan Seifert
# 2013, Licensed under GPLv2

from numpy import *
from numpy import linalg as la
from scipy.optimize import curve_fit
import math
import time
from re import search

# Given an int J_MAX, and three floats corresponding to the distortable rotor rotational constants
# ham will calculate the Hamiltonian matrix for the given J_MAX value.
def ham(J_MAX, rot_A, rot_B, rot_C, DJ=0.0, DJK=0.0, DK=0.0, subdj=0.0,subdk=0.0):
	ham = zeros((2*J_MAX+1,2*J_MAX+1))
	for i in range(-1*(J_MAX),J_MAX+1):
		for j in range(-1*J_MAX,J_MAX+1):
			l =  i + J_MAX
			m = j + J_MAX
			if l == m:
				ham[l,m] = 0.5*(rot_B+rot_C)*(J_MAX*(J_MAX+1))+(rot_A-(rot_B+rot_C)/2.0)*(i**2)-(DJ*J_MAX**2*(J_MAX+1)**2+DJK*i**2*J_MAX*(J_MAX+1)+DK*i**4)
			elif m == l-2:
				ham[l,m] = ((rot_B-rot_C)/4.0-J_MAX*(J_MAX+1)*subdj-(i**2+j**2)*subdk*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i-1)*(i-2))*(J_MAX*(J_MAX+1)-i*(i-1)))
			elif m == l+2:
				ham[l,m] = ((rot_B-rot_C)/4.0-J_MAX*(J_MAX+1)*subdj-(i**2+j**2)*subdk*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i+1)*(i+2))*(J_MAX*(J_MAX+1)-i*(i+1)))
			else:
				ham[l,m] = 0
	return ham
#def ham(J_MAX, rot):
#	ham = zeros((2*J_MAX+1,2*J_MAX+1))
#	for i in range(-1*(J_MAX),J_MAX+1):
#		for j in range(-1*J_MAX,J_MAX+1):
#			l =  i + J_MAX
#			m = j + J_MAX
#			if l == m:
#				ham[l,m] = 0.5*(rot[1]+rot[2])*(J_MAX*(J_MAX+1))+(rot[0]-(rot[1]+rot[2])/2.0)*(i**2)-(rot[3]*J_MAX**2*(J_MAX+1)**2+rot[4]*i**2*J_MAX*(J_MAX+1)+rot[5]*i**4)
#			elif m == l-2:
#				ham[l,m] = ((rot[1]-rot[2])/4.0-J_MAX*(J_MAX+1)*rot[6]-(i**2+j**2)*rot[7]*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i-1)*(i-2))*(J_MAX*(J_MAX+1)-i*(i-1)))
#			elif m == l+2:
#				ham[l,m] = ((rot[1]-rot[2])/4.0-J_MAX*(J_MAX+1)*rot[6]-(i**2+j**2)*rot[7]*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i+1)*(i+2))*(J_MAX*(J_MAX+1)-i*(i+1)))
#			else:
#				ham[l,m] = 0
#	return ham

# STANDARD NOTATION FOR A ROTATIONAL TRANSITION: For a given transition J'Ka'Kc' --> J''Ka''Kc'', 
# these routines will consider the following structure for a set oftransitions:
# NUMPY ARRAY --- [[J'1, Ka'1, Kc'1, J''1, Ka''1, Kc''1],[J'2,Ka'2,Kc'2,J''2,Ka'',Kc''2],....[nth transition]]

# freq_single calculates the frequency for a single transition.
# QN is a numpy array with the row structure shown above for a given rotational transition
def freq_single(QN, A,B,C,DJ,DJK,DK,subdj,subdk):
	set_printoptions(precision=4)
	Jupper = QN[0]
	ind_upper = Jupper - QN[2] + QN[1]
	Jlower = QN[3]
	ind_lower = Jlower - QN[5] + QN[4]

	H_upper = ham(Jupper, A,B,C,DJ,DJK,DK,subdj,subdk)
	E_upper, trash = la.eig(H_upper)
	E_upper = sort(E_upper)

	H_lower = ham(Jlower, A,B,C,DJ,DJK,DK,subdj,subdk)
	E_lower, trash = la.eig(H_lower)
	E_lower = sort(E_lower)

	return float(E_upper[ind_upper]-E_lower[ind_lower])

# Calculates frequencies for all QNs in a single array
def freq_all(QN_array, A,B,C,DJ=0.0,DJK=0.0,DK=0.0,subdj=0.0,subdk=0.0):
	set_printoptions(precision=4)
	output = zeros((shape(QN_array)[0],1))
	for i in range(shape(QN_array)[0]):
		output[i] = freq_single(QN_array[i],A,B,C,DJ,DJK,DK,subdj,subdk)
	return reshape(output,shape(QN_array)[0])

# This calculates the # digit that is nonzero in the standard error of a constant in order to do proper rounding of the constant + two digit std. errors
def digitfind(num):
		return len(search("\.(0*)",str(num)).group(1))

def rms(obs,calc):
    omc_sq = zeros(shape(obs)[0])
    for i,line in enumerate(obs):        
        omc_sq[i] = (obs[i]-calc[i])**2
    n = len(obs)-1
    sigma = sum(sqrt((n**-1)*omc_sq))
    return sigma

def omc(obs, calc):
    omc = zeros(shape(obs)[0])
    for i,line in enumerate(obs):        
        omc[i] = "{0:.2f}".format((obs[i]-calc[i])*1000)
    return omc

def report(popt,covar,trans,freqs):
	
	A_fit = popt[0]
	B_fit = popt[1]
	C_fit = popt[2]
	variance = diagonal(covar)
	std_error = sqrt(variance)
	A_stderr = std_error[0]
	B_stderr = std_error[1]
	C_stderr = std_error[2]
	
	print 'Fit Parameters:'
	
	switch = int(digitfind(A_stderr))
	print('A (MHz): ' + str(round(A_fit,switch+2)) + '(' + "{0:.0f}".format(float(A_stderr)*10**(switch+2)) + ')')
	switch = int(digitfind(B_stderr))
	print('B (MHz): ' + str(round(B_fit,switch+2)) + '(' + "{0:.0f}".format(float(B_stderr)*10**(switch+2)) + ')')
	switch = int(digitfind(C_stderr))
	print('B (MHz): ' + str(round(C_fit,switch+2)) + '(' + "{0:.0f}".format(float(C_stderr)*10**(switch+2)) + ')')
	
	# Calculate RMS errors + OMC 
	calc_freqs = freq_all(trans,A_fit,B_fit,C_fit)
	rms_error = rms(freqs,calc_freqs)
	omcs = omc(freqs,calc_freqs)
	print '\n'
	print '------------ LINELIST ---------'
	print('J Ka Kc  -->  J Ka Kc'+'	 '+'OBS FREQ.'+'    '+'CALC. FREQ'+'   '+'OMC (kHz)')
	
	for i in range(shape(calc_freqs)[0]):
		print(str(trans[i,0])+'  '+str(trans[i,1])+'  '+str(trans[i,2])+'  '+'-->'+'  '+str(trans[i,3])+'  '+str(trans[i,4])+'  '+str(trans[i,5])+'	'+"{:>10.3f}".format(calc_freqs[i])+'    '+"{:>10.3f}".format(freqs[i]))+'  '+"{:>10.2f}".format(omcs[i])
	print 'RMS Error: ' + "{0:.4f}".format(rms_error*10**3) + ' kHz'
	

# Test block
# Hexanal, conformer I

A = 9769.62213
B = 868.846659
C = 818.518746
DJ = 0.000047239
DJK = -0.0006991
DK =  0.023173
subdj = 5.0298E-06
subdk = 0.000343
guess_constants = array([A,B,C])


trans = array([[4,0,4,3,0,3],[5,1,5,4,1,4],[5,0,5,4,0,4],[1,1,0,1,0,1],[2,1,2,1,0,1],[8,2,7,7,2,6],[8,2,6,7,2,5]])
freqs = array([6747.32023,8310.06771,8432.54591,8951.08513,12225.16434,13496.27440,13514.14580])
popt, pcov = curve_fit(freq_all,trans,freqs,guess_constants, sigma=None)
report(popt,pcov, trans,freqs)




# TIME TEST BLOCK
#for i in range(shape(trans)[0]):
#	freq = freq_single(trans[i],A,B,C,DJ,DJK,DK,subdj,subdk)
#	print str(trans[i,0])+" "+str(trans[i,1])+" "+str(trans[i,2])+" ---> "+str(trans[i,3])+" "+str(trans[i,4])+" "+str(trans[i,5])+"	"+str(freq)
#t2 = time.time()
#print str(t2-t1)
