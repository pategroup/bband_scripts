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
import textwrap

# Given an int J_MAX, and three floats corresponding to the distortable rotor rotational constants
# ham will calculate the Hamiltonian matrix for the given J_MAX value.
# dist is an array of Watson Hamiltonian CD parameters, in the following (standard) order:
# [DJ, DJK, DK, dJ, dK]
def ham(J_MAX, rot_A, rot_B, rot_C, DJ,DJK,DK,subdj,subdk):
	ham = zeros((2*J_MAX+1,2*J_MAX+1))
	for i in range(-1*(J_MAX),J_MAX+1):
		for j in range(-1*J_MAX,J_MAX+1):
			l =  i + J_MAX
			m = j + J_MAX
			if l == m:
				ham[l,m] = 0.5*(rot_B+rot_C)*(J_MAX*(J_MAX+1))+(rot_A-(rot_B+rot_C)/2.0)*(i**2)-(DJ*J_MAX**2*(J_MAX+1)**2+DJK*i**2*J_MAX*(J_MAX+1)*DK*i**4)
			elif m == l-2:
				ham[l,m] = ((rot_B-rot_C)/4.0-J_MAX*(J_MAX+1)*subdj-(i**2+j**2)*subdk*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i-1)*(i-2))*(J_MAX*(J_MAX+1)-i*(i-1)))
			elif m == l+2:
				ham[l,m] = ((rot_B-rot_C)/4.0-J_MAX*(J_MAX+1)*subdj-(i**2+j**2)*subdk*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i+1)*(i+2))*(J_MAX*(J_MAX+1)-i*(i+1)))
			else:
				ham[l,m] = 0
	return ham


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
	
# Loop for freq_single for a set of QN numbers QN_array with array structure noted above.
def freq_all(QN_array, A,B,C,DJ,DJK,DK,subdj,subdk):
	output = zeros((shape(QN_array)[0],1))
	for i in range(shape(QN_array)[0]):
		output[i] = freq_single(QN_array[i],A,B,C,DJ,DJK,DK,subdj,subdk)
	return reshape(output,shape(QN_array)[0])

	
# The following string defines the fitting routine that would be appropriate where we want to fit A, B, and C but fix the distortion
# constants to the specified value. 
# Note that the defs for ham(), freq_single(), and freq_all() are identical to those above. 
# However, due to the nature of Scipy's curve_fit Levenburg-Marquardt procedure, you get an error if you try to pass fixed variables
# along with the floated A/B/C set. 
# Therefore, this following routine dynamically creates the function for curve_fit so A&B&C can be passed as floated variables, 
# but the distortion is fixed. 

funcstr = textwrap.dedent('''\
from numpy import *
from numpy import linalg as la
def ham(J_MAX, rot_A, rot_B, rot_C, DJ,DJK,DK,subdj,subdk):
	ham = zeros((2*J_MAX+1,2*J_MAX+1))
	for i in range(-1*(J_MAX),J_MAX+1):
		for j in range(-1*J_MAX,J_MAX+1):
			l =  i + J_MAX
			m = j + J_MAX
			if l == m:
				ham[l,m] = 0.5*(rot_B+rot_C)*(J_MAX*(J_MAX+1))+(rot_A-(rot_B+rot_C)/2.0)*(i**2)-(DJ*J_MAX**2*(J_MAX+1)**2+DJK*i**2*J_MAX*(J_MAX+1)*DK*i**4)
			elif m == l-2:
				ham[l,m] = ((rot_B-rot_C)/4.0-J_MAX*(J_MAX+1)*subdj-(i**2+j**2)*subdk*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i-1)*(i-2))*(J_MAX*(J_MAX+1)-i*(i-1)))
			elif m == l+2:
				ham[l,m] = ((rot_B-rot_C)/4.0-J_MAX*(J_MAX+1)*subdj-(i**2+j**2)*subdk*0.5)*math.sqrt((J_MAX*(J_MAX+1)-(i+1)*(i+2))*(J_MAX*(J_MAX+1)-i*(i+1)))
			else:
				ham[l,m] = 0
	return ham

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
def freq_all(QN_array, A,B,C,{p}):
	output = zeros((shape(QN_array)[0],1))
	for i in range(shape(QN_array)[0]):
		output[i] = freq_single(QN_array[i],A,B,C,DJ,DJK,DK,subdj,subdk)
	return reshape(output,shape(QN_array)[0])
''')

def make_func(**kwargs):
	params = set(('DJ','DJK','DK','subdj','subdk')).difference(kwargs.keys())
	exec funcstr.format(p=','.join(params)) in kwargs
	return kwargs['freq_all']

# ------------ END DYNAMICALLY GENERATED FUNCTIONAL FOR CURVE_FIT ---------------- #
	
# ----------------------- BEGIN PRINTING/STATS ROUTINES ------------------------- #	
	
# This calculates the # digit that is nonzero in the standard error of a constant in order to do proper rounding of the constant + two digit std. errors
def digitfind(num):
		return len(search("\.(0*)",str(num)).group(1))

# Calculates microwave RMS, in MHz
def rms(obs,calc):
    omc_sq = zeros(shape(obs)[0])
    for i,line in enumerate(obs):        
        omc_sq[i] = (obs[i]-calc[i])**2
    n = len(obs)-1
    sigma = sum(sqrt((n**-1)*omc_sq))
    return sigma

# Calculates an array of OMC values for each assigned transition
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
	
	switch = int(digitfind(A_stderr)) # Determines the position after the decimal point of the first nonzero digit
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

# ------------------------ END PRINTING/STATS ROUTINES -----------------=--------- #	


# ----------------------------- TESTING BLOCK ---------------------- #

# Hexanal, conformer I
A = 9769.62213
B = 868.846659
C = 818.518746
d1 = 0.000047239 #DJ
d2 = -0.0006991 #DJK
d3 =  0.023173 #DK
d4 = 5.0298E-06 #dJ
d5 = 0.000343 #dK
#distortion = array([DJ,DJK,DK,subdj,subdk])

guess_constants = array([A,B,C])
trans = array([[4,0,4,3,0,3],[5,1,5,4,1,4],[5,0,5,4,0,4],[1,1,0,1,0,1],[2,1,2,1,0,1],[8,2,7,7,2,6],[8,2,6,7,2,5]])
freqs = array([6747.32023,8310.06771,8432.54591,8951.08513,12225.16434,13496.27440,13514.14580])

# ALLOCATES DYNAMIC FUNCTION FOR FITTING
freq_all = make_func(DJ=d1,DJK=d2,DK=d3,subdj=d4,subdk=d5)

popt, pcov = curve_fit(freq_all,trans,freqs,guess_constants, sigma=None)
report(popt,pcov, trans,freqs)

# ----------------------------- END TESTING BLOCK ---------------------- #



# TIME TEST BLOCK
#for i in range(shape(trans)[0]):
#	freq = freq_single(trans[i],A,B,C,DJ,DJK,DK,subdj,subdk)
#	print str(trans[i,0])+" "+str(trans[i,1])+" "+str(trans[i,2])+" ---> "+str(trans[i,3])+" "+str(trans[i,4])+" "+str(trans[i,5])+"	"+str(freq)
#t2 = time.time()
#print str(t2-t1)
