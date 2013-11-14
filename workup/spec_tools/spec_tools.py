# ======================== SPEC TOOLS =========================== \
# 	        A PYTHON LIBRARY FOR ROTATIONAL SPECTROSCOPY          |
# Developed by, and for, chirped-pulse rotational spectroscopists |
# ----------------------------------------------------------------\
#                 MAINTAINER: Nathan Seifert                      |
#                 e-mail: nas3xf@virginia.edu                     |
#                                                                 |
#                      Brooks Pate Lab                            |
#           Dept. of Chem., University of Virginia                |
# ----------------------------------------------------------------\
#                           LICENSE:                              |
# This module is licensed under the GPLv2 license. For additional |
# information, visit http://www.gnu.org/licenses/gpl-2.0.html     |
# ----------------------------------------------------------------\
#                           PURPOSE:                              |
# This module contains functions that are useful for data workup  |
# and routine analysis of broadband rotational spectra data. This |
# includes FFT routines for free induction decay data, routines   |
# for peakpicking transitions and cutting known transitions out   |
# of a given spectrum, and other useful tools.                    |
# ----------------------------------------------------------------\
#                         AUTHORSHIP:                             |
# The primary maintainer of this module is listed above. Some of  |
# the code used in this library was written by current and former |
# the Pate lab, such as Cristobal Perez and Daniel Zaleski.       |
# ----------------------------------------------------------------\
#                        REQUIREMENTS:                            |
# Currently, this module only requires a few additional libraries | 
# outside of the default Python libraries: numpy, matplotlib.     |
#  MORE TO ADD LATER?                                             |   
#        * This module was developed for Python 2.7.x *           |
# ----------------------------------------------------------------\
#                 TO DO (will be removed gradually)               |
# ARB_PULSE_GEN, SPECTROGRAM                                      |
# ISO_SCALE, SPEC_CUT, SPCAT_WRITER, LAZY_QUAD_CONVERT, GAIN_CORR |
# (ROOM FOR MORE TO ADD LATER)                                    |
#                                                                 |
#                       Specific comments:                        |
# ROT_COORDS:                                                     |
#    Currently rot_coords requires a file input. Adding           | 
# functionality where you can input a 4-column array geometry     |
# directly needs to be added.                                     |
#                                                                 |
# FFT:                                                            |
#	 Currently fft requires a file input. It should be straight-  |
# forward to add functionality where a 1-column FID array can     |
# be inputted and FFTed.                                          |
#                                                                 |
# PEAKPICK:                                                       |
# 	I'm pretty happy about where this is. Again, there needs to   |
# be some functionality where a spectrum as a numpy array in mem  |
# can be pushed to the routine, without having to specify a file- |
# name.                                                           |
# =============================================================== \

# IMPORT BLOCK #
from numpy import *
from scipy import *
from scipy.interpolate import *
import scipy.fftpack as sfft
import matplotlib as mpl
import fileinput
import math
import re, sys
import os
from matplotlib import pyplot as plt
import time
# END IMPORT BLOCK #

# METADATA BLOCK #
__version__ = "1.0_WIP"
__author__ = "Nathan Seifert"
__license__ = "GPLv2"
# END METADATA BLOCK #


# ======================== ROTATE_COORDS ======================= \
# The block below defines functions for the purpose of taking    |
# Cartesian coordinates of an arbitrary orientation and rotates  |
# them into the principal axis. A number of these functions will |
# not be useful for general purpose, but provide support for the | 
# primary functionality.                                         |
# ---------------------------------------------------------------\
#                             USAGE:                             |
#                                                                |
# rot_coords(file COORDS, float mu_A, float mu_B, float mu_C):   |
#                                                                |
# Takes input file COORDS, which should be a four column text    |
# file with the following format for each row entry:             |
#                                                                |   
#     row_i : L_i      X_i   Y_i   Z_i                           |
#                                                                |
# Where L_i is the atomic symbol for the atom in row_i (L_i can  |
# either be numeric or a letter -- e.g. 6 or C for carbon)       |
# and X_i, Y_i, and Z_i are the Cartesian coordinates of that    |
# atom.                                                          |
#                                                                | 
# rot_coords also accepts three floats, representing the         |
# magnitudes of the dipoles in the three Cartesian directions.   |
# If you supply a coordinate file along with the three dipoles,  |
# rot_coords will output an array with the three input dipoles   |
# rotated into the principal axis.                               |
#                                                                |
# 						EXAMPLE USAGE:                           |
#                                                                |
# <import coordinate file FileXYZ>                               | 
# ua = 1.0 #Debye, in FileXYZ frame                              |
# ub = 0.5 #Debye                                                |
# uc = 0.1 #Debye                                                |
# coords, dipoles = spectools.rot_coords(FileXYZ, ua, ub, uc )   |
# ============================================================== |

# =============     BEGIN ROTATE_COORDS BLOCK     ============== #


# _masses takes an input string and matches it with an atomic weight
def _masses(symbol):

	masses_lett = {'H':1.007825037, 'C':12.0, 'O':15.99491464, 'N': 14.003074008, 'Cl': 34.96885268, 'F': 18.99840322, 'Si':27.9769265325, 'S': 31.97207100}
	masses_num = {'1':1.007825037, '6':12.0, '8':15.99491464, '7': 14.003074008, '17': 34.96885268, '9': 18.99840322, '14':27.9769265325, '16': 31.97207100}
	
	if symbol in masses_lett:	
		return masses_lett[symbol]
	if symbol in masses_num:
		return masses_num[symbol]

# comshift(Array matr) shifts coordinates from an input coordinate matrix 
# into a center of mass frame
def _comshift(matr):
	total_mass = sum(matr[:,0])

	# Creates dummy arrays to calc COM shift
	com_x = []
	com_y = []
	com_z = []
	for row in matr:
		com_x.append(row[0] * row[1])
		com_y.append(row[0] * row[2])
		com_z.append(row[0] * row[3])
	com_x = sum(com_x)/total_mass
	com_y = sum(com_y)/total_mass
	com_z = sum(com_z)/total_mass

	# Shifts matr[] into COM frame
	for i, row in enumerate(matr):
		matr[i,1] -= com_x
		matr[i,2] -= com_y 
		matr[i,3] -= com_z

	return matr

# Takes a center of mass-framed coordinate matrix from _comshift,
# calculates the moment of inertia tensor, and diagonalizes the tensor

def _calcabc(matr):
	
	first_mom = zeros((3,3)) # empty first moment tensor
	tot_atom_count = size(matr[:,0]) # counts numbers of row in input coord matrix
	temp = zeros((tot_atom_count,9)) # dummy matrix for later
	
	# Calculates atomic contributions to each tensor element
	for i in range(tot_atom_count):
        # Ixx, Ixy, Ixz  
		temp[i,0] = matr[i,0]*(matr[i,2]**2 + matr[i,3]**2)
		temp[i,1] = -1*(matr[i,0]*matr[i,1]*matr[i,2])
		temp[i,2] = -1*(matr[i,0]*matr[i,1]*matr[i,3])
         
        # Iyx, Iyy, Iyz 
		temp[i,3] = -1*(matr[i,0]*matr[i,2]*matr[i,1]) 
		temp[i,4] = matr[i,0]*(matr[i,1]**2 + matr[i,3]**2)
		temp[i,5] = -1*(matr[i,0]*matr[i,2]*matr[i,3])

		# Izx, Izy, Izz
		temp[i,6] = -1*(matr[i,0]*matr[i,3]*matr[i,1])
		temp[i,7] = -1*(matr[i,0]*matr[i,3]*matr[i,2])
		temp[i,8] = matr[i,0]*(matr[i,1]**2 + matr[i,2]**2)

	# Sums all contributions to form total inertia tensor
	# NOTE: There's gotta be an easy way to merge the above step
	# with the below step...

	first_mom[0,0] = sum(temp[:,0])	
	first_mom[0,1] = sum(temp[:,1])
	first_mom[0,2] = sum(temp[:,2])

	first_mom[1,0] = sum(temp[:,3])
	first_mom[1,1] = sum(temp[:,4])
	first_mom[1,2] = sum(temp[:,5])

	first_mom[2,0] = sum(temp[:,6])
	first_mom[2,1] = sum(temp[:,7])
	first_mom[2,2] = sum(temp[:,8])				

	# Diagonalizes the tensor
	# diag_eigenvecs is the rotation matrix
	# diag_eigenvals is the set of rotational constants (in amu Angstrom**2)
	diag_eigenvals,diag_eigenvecs = linalg.eigh(first_mom)

	# Convert to MHz
	rotconst = 505379.006/diag_eigenvals
    
    # rotconst is an array with three members (A, B, C)
    # diag_eigenvecs is a 3x3 rotation matrix 
	return rotconst, diag_eigenvecs 

# Takes a 4-column coordinate matrix and the rotation matrix output by _calcabc() and returns
# a 3-column coordinate matrix rotated into the principal axis
def _rotate(coords,rotmat):
	coords_mod = delete(coords,s_[:-3],1) # removes atomic symbol column
	coords_mod = dot(coords_mod,rotmat) # matrix transformation into principal axis

	return coords_mod

# Rotates a 3-member array of dipoles into the principal axis
def _rotatedipoles(dipoles,rotmat):
	dipoles = dot(dipoles,rotmat)
	return dipoles



def _write(const,rotmat,syms,rot_coordmat,dipoles,filename):
	set_printoptions(suppress=True)
	filename_out = filename+"_out.geom"
	file = open(filename_out,"wb")
	file.write("Input file name: "+filename+"\n")

	# Writes rotated dipoles
	if dipoles[0] != 0 or dipoles[1] != 0 or dipoles[2] != 0:
		file.write("--------------\n"+"DIPOLES (in Debye)\n"+"--------------\n\n"+"ua: "+str(dipoles[0])+"	"+"ub: "+str(dipoles[1])+"	"+"uc: "+str(dipoles[2])+"\n\n")

	file.write("--------------\n"+"ROTATIONAL CONSTANTS\n"+"--------------\n\n"+"B : "+str(const[0])+" MHz "+"   "+"A : "+str(const[1])+" MHz "+"   "+"C : "+str(const[2])+" MHz "+"   \n\n"+"--------------\n"+"ROTATION MATRIX\n"+"--------------\n\n")

	# Writes rotation matrix
	for i in range(shape(rotmat)[0]):
		for j in range(shape(rotmat)[1]):
			file.write(str(rotmat[i,j])+ "		")
		file.write("\n")
	file.write("\n--------------\n")
	file.write("ROTATED COORDINATES\n"+"--------------\n\n")

	#Writes rotated coordinates with original atomic symbols
	for i in range(size(syms)):
		file.write(str(syms[i]) + "	")
		for j in range(shape(rot_coordmat)[1]):
			file.write(str(rot_coordmat[i,j]) + "		")
		file.write("\n")
	file.close


# Main routine. Coords_matrix is a text filename containing four columns, as specified
# in the usage at the beginning of the ROTATE_COORDS block
def rot_coords(coords_matrix, ua=0.0, ub=0.0, uc=0.0,savetofile=1):

	#set_printoptions(suppress=True)
	symbols = genfromtxt(coords_matrix,dtype=str,usecols=(0))
	coords = genfromtxt(coords_matrix,dtype=float,usecols=(1,2,3))

	# Create merged matrix for the rest of the routine
	coordmat = zeros((size(symbols),4))
	for i in range(size(symbols)):
		coordmat[i,0] = _masses(symbols[i])
		coordmat[i,1] = coords[i,0]
		coordmat[i,2] = coords[i,1]
		coordmat[i,3] = coords[i,2]
	#print coordmat

	constants,rotmatr = _calcabc(_comshift(coordmat))
	rotcoords = _rotate(coordmat,rotmatr)
	dipolestack = zeros(3)

	# dipole rotation -- ua, ub or uc must be nonzero for this to work! Default is all zero.
	if ua != 0.0 or ub != 0.0 or uc != 0.0:
		dipolestack[0] = ua
		dipolestack[1] = ub
		dipolestack[2] = uc
		#print dipolestack
		dipolestack = _rotatedipoles(dipolestack,rotmatr) 
		#print dipolestack

	if savetofile == 1:
		_write(constants,rotmatr,symbols,rotcoords,dipolestack,coords_matrix)
		return constants,rotmatr
	else:
		return constants,rotmatr


# Testing block for rotate_coords #

#filename = input("Enter a coordinate file name: ")
#rot_coords(filename,0.0,0.0,0.0,1)
#set_printoptions(suppress=True)

# =============     END ROTATE_COORDS BLOCK     ================ #


# ========================== FFT =============================== \
# This routine will take a time-domain data file and Fourier     |
# transform it into frequency space. In general, a Kaiser-Bessel |
# window will be applied to the Fourier transform.               |
# FFT() will accept FIDs where the amplitude information is in   |
# the LAST column of the data file. For instance, a 1-column FID |
# will be read correctly, as well as the 5-column format that    |
# is outputted by default on the Tektronix DPO7000-series        |
# oscilloscopes.                                                 |
#----------------------------------------------------------------\
#							USAGE:							     |
# fft(filename, start_freq, end_freq, srate, gaincorr,gain_mat)  |
#                                                                |
# filename = string of the target FID filename                   |
# start_freq, end_freq = frequency bounds in MHz                 |
# srate = sampling rate in samples*Hz, so 50 Gs/s = 50E9         |
# gaincorr = [default is 0] set to 1 if you want to gain correct |
# gain_mat = array with length equivalent to # pts in FT with    |
# gain correction information (GAIN_COR can generate this array) |
# ===============================================================\


# =================     BEGIN FFT BLOCK     ==================== #

def fft(filename, start_freq, end_freq, srate,gaincorr=0,gaincorr_mat=None):
	data = open(filename)
	time_domain = []

	for row in data:
		temp = row.split()
		time_domain.append(float(temp[size(temp)-1]))
	td_size = size(time_domain)
	print "The size of the time domain file is: "+str(td_size)

	# Kaiser-Bessel window with a beta of 9.5 is used.
	# Replace this with another window if you don't want to use it.
	window = kaiser(td_size,9.5)

	for i, row in enumerate(time_domain):
		time_domain[i] = row*window[i]

	temp = zeros(td_size*2)
	temp[0:td_size] = time_domain

	fft = abs(sfft.fft(temp))/100
	freq = (sfft.fftfreq(size(temp),1/srate))/1E6

	out = zeros((size(freq),2))
	for i, row in enumerate(freq):
		if row >= start_freq and row <= end_freq:
			out[i,0] = freq[i]
			out[i,1] = fft[i]
			if gaincorr == 1:
				out[i,1] = fft[i]*gaincorr_mat[i,1]

	return out

# Testing block for fft() #
#filename = input("Input a FID file name")
#spectrum = fft("fid",6000,18000,50E9)
#print shape(spectrum)
#plt.plot(spectrum[:,0],spectrum[:,1],'blue',linewidth=0.8)
#plt.grid(True)
#plt.show()
# =================     END FFT BLOCK     ===================== #


# ========================== PEAKPICK ========================= \
# This routine will take an input spectrum and output an array  |
# of the frequencies and intensities of every transition above  |
# the given intensity threshold.                                |
# There are two subroutines that do essentially the same thing, |
# but one, _ppspline() first splines the spectrum with a given  |
# frequency resolution and then peakpicks the spectrum. This    |
# has the advantage of improving the center-frequency resol-    |
# ution by a significant factor if the splining resolution is   |
# small compared to the point spacing of your spectrum.         |
# ------------------------------------------------------------- \
#							USAGE:                              |
# peakpick(filename,threshold,type,resolution)                  |
#                                                               |
# filename = input filename as string                           |
# threshold = a float indicating the *minimum* intensity a peak |
# can be.                                                       |
# resolution = splining resolution in MHz. Default is 2 kHz.    |
# Routine ignores this if "type" is set to 0.                   |
# type = can be 0 or 1. Default is 0. Choose 1 if you want      |
# a splined peakpick.                                           |
# ============================================================= \

# =================     BEGIN PEAKPICK BLOCK     ==================== #

def peakpick(filename,threshold,type = 0,resolution = 0.002):
	if type == 0:
		return _pp(filename,threshold)
	if type == 1:
		return _ppspline(filename,threshold,resolution)

# Point-spacing resolution-limited peakpick
def _pp(filename, threshold):
	spec = loadtxt(filename)
	peaks = []
	for i in range(0, len(spec)-1):
		if float(spec[i,1]) > threshold and float(spec[i,1]) > float(spec[(i-1),1]) and float(spec[i,1]) > float(spec[(i+1),1]):
			peaks.append(spec[i])
	#print "this shit is this long:"+str(len(peaks))
	picks = zeros((len(peaks),2))
	for i, row in enumerate(peaks):
		picks[i,0] = row[0]
		picks[i,1] = row[1]
	return picks

# Two steps: spline the spectrum completely, then do the same as _pp
# Routine adapted from Steve Shipman's code in autofit
def _ppspline(filename,threshold,resolution):
	spec = loadtxt(filename)

	x = spec[:,0]
	y = spec[:,1]

	old_res = (x[-1]-x[0])/len(spec)
	scale = old_res/resolution

	new_length = int(math.floor(scale*len(spec)))

	tick = splrep(x,y,s=0)
	xnew = arange(x[0],x[-1],resolution)
	ynew = splev(xnew,tick,der=0)

	output = column_stack((xnew,ynew))
	peaks = []
	for i in range(0,len(output)-1):
		if float(output[i,1]) > threshold and float(output[i,1]) > float(output[(i-1),1]) and float(output[i,1]) > float(output[(i+1),1]):
			peaks.append(output[i])
	picks = zeros((len(peaks),2))
	#print "Length of the peakpick is: " + str(len(peaks))
	for i, row in enumerate(peaks):
		picks[i,0] = row[0]
		picks[i,1] = row[1]
	#print "The length of the numpy array peakpick is: "+str(size(picks))
	return picks



#  Test block for peakpick() #
#pick = peakpick("spec",0.01, 1)
#print str(type(pick))

#for i in range(0,10):
#	print str(pick[i,0])+"		"+str(pick[i,1])

# =================     END PEAKPICK BLOCK     ==================== #



# ========================== SPEC_CUT ========================= \
# This routine will take two files -- a 2-column spectrum and a |
# one-column list of frequencies -- and cut out blocks in the   |
# experimental spectrum with center frequencies equivalent to   |
# those in the one-column list with an inputted constant width. |
# For instance, if a target cut transition has frequency 10321  |
# MHz, and the inputted width is 200 kHz, then the output       |
# spectrum will have 10320.8-10321.2 MHz set to 0 intensity.    |
# ------------------------------------------------------------- \

# =================		BEGIN CUT BLOCK    ==================== #

def cut(spectrum_file,linelist_file,width):
	spec = loadtxt(spectrum_file)
	lnlist = loadtxt(linelist_file)

	specstart = spec[0,0]
	step_size = (spec[100,0]-spec[0,0])/100 # This is an arbitrary choice.

	mask = ones(shape(spec)[0])
	for i, line in enumerate(lnlist):
		freq1 = lnlist[i] - (width/2)
		freq2 = lnlist[i] + (width/2)

		i1 = math.floor((freq1-specstart)/step_size)+1
		i2 = math.floor((freq2-specstart)/step_size)-1

		# Sets the bounds of a cut to end/start of spectrum if cut goes over start/stop of spec
		if i1 < spec[0,0]:
			i1 = spec[0,0]
		if i2 > spec[0,0]:
			i2 = spec[-1,0]

		i1 = int(i1)
		i2 = int(i2)

		for i in range(i1,i2):
			mask[i] = 0

	mask_final = column_stack((spec[:,0],mask))
	cut_spec = spec
	for i, line in enumerate(spec):
		cut_spec[i,1] = cut_spec[i,1] * mask_final[i,1]

	return cut_spec

# Test block for cut()
#cutspec = cut("spec","cuts",500)
#plt.plot(cutspec[:,0],cutspec[:,1])
#plt.show()

# =================		END CUT BLOCK    ====================== #


# ========================== ISO_SCALE ======================== \
# This routine will take in a set of experimental rotational    |
# constants and an experimental geometry of the same form as    |
# rot_coords (line 91), and will output scaled rotational       |
# constants for each of the isotopologues.                      |
# -------------------------------------------\\\----------------\
#                             USAGE:                            |
# scale(geomfile,A_exp,B_exp,C_exp,deut_flag)                   |
# where:                                                        |
# - geomfile is a filename corresponding to the geometry used   |
# for scaling                                                   |
# - A/B/C_exp are the experimental rotational constants for     |
# the parent isotopic species (corresp. to geomfile)            |
# - deut_flag takes 0 or 1. If 1, scale() will output all       |
# scaled deuterium isotopologic constants                       |
# ------------------------------------------------------------- \


# ==============	BEGIN ISO_SCALE BLOCK   =================== #


# Returns a list of labels and isotopologue masses for an input mass
def _isomass(mass):

	masses = {1.007825037:['2D',2.0141017778],12.0:['13C',13.0033548378], 14.003074008:['15N',15.0001088982],15.99491464:['18O',17.9991610],34.96885268:['37Cl',36.96590259],27.9769265325:['29Si',28.976494700,'30Si',29.97377017],31.97207100:['34S',33.96786690]}
	
	if mass in masses:
		return masses[mass]
	else:
		return mass

#def _isomass(mass):
#	try:
#		if mass == 1.007825037:
#			return ['2D',2.0141017778]
#		if mass == 12.00000000:
#			return ['13C',13.0033548378]	
#		if mass - 15.9949164 < 0.01:
#			return ['18O',17.9991610]
#		if mass == 14.003074008:
#			return ['15N',15.0001088982]
#		if mass == 34.96885268:
#			return ['37Cl',36.96590259]
#		if mass == 27.9769265325:
#			return ['29Si',28.976494700,'30Si',29.97377017]	
#		if mass == 31.97207100:
#			return ['34S',33.96786690]
#			
#		if isinstance(mass,float):
#			return mass
#		else:
#			print "Messed up mass: " + str(mass)
#			raise ValueError('BLANK_UNKNOWN_SYMBOL')
#	except ValueError:
#		print 'ERROR! An atomic symbol is blank or not recognized.'
#		raise

def scale(geomfile,A_exp=0.0,B_exp=0.0,C_exp=0.0,deut_flag=0):

	symbols = genfromtxt(str(geomfile),dtype=str,usecols=(0))
	coords = genfromtxt(str(geomfile),dtype=float,usecols=(1,2,3))

	# Create merged matrix for the rest of the routine
	coordmat = zeros((size(symbols),4))

	for i in range(size(symbols)):
		coordmat[i,0] = _masses(symbols[i])
		#print coordmat[i,0]
		coordmat[i,1] = coords[i,0]
		coordmat[i,2] = coords[i,1]
		coordmat[i,3] = coords[i,2]

	# Initialize output file
	ns_constants,rotmatr = _calcabc(_comshift(coordmat))
	out = open(str(geomfile)+'_scaled_isotopologues.txt','wb')
	out.write('Scaled isotopologues for geometry file: ' + geomfile + '\r\n')
	out.write('-----------------------------------------------\r\n')
	out.write('Theoretical constants for normal species: '+ str(round(ns_constants[0],3))+' '+str(round(ns_constants[1],3))+' '+str(round(ns_constants[2],3))+"\r\n")
	out.write('Experimental constants for normal species: ' + str(A_exp)+' '+str(B_exp)+' '+str(C_exp)+'\r\n')
	out.write('-----------------------------------------------\r\n')
	
	

	# Main scaling routine
	for i in range(0,shape(coordmat)[0]):
		if coordmat[i,0] == 1.007825037 and deut_flag==0:
			print 'Am I here?'
			continue
		iso_temp = _isomass(coordmat[i,0])
		print iso_temp
		if isinstance(iso_temp,float):
			print 'Am i here?'
			continue
		out.write('\r\n-------------------------------------\r\n')
		out.write('Input atom '+str(i+1)+': '+ symbols[i] + '	'+str(round(coordmat[i,1],3))+' '+str(round(coordmat[i,2],3))+' '+str(round(coordmat[i,3],3))+'\r\n')
		out.write('-------------------------------------\r\n')
		for j in range(0,len(iso_temp)):
			if j % 2 == 0:
				# Check to see if isotope is deuterium and if we want deuterium
				if iso_temp[j] == '2D' and deut_flag == 0:
					continue	
				temp_mat = copy(coordmat)
				label = iso_temp[j]
				temp_mat[i,0] = iso_temp[j+1]
				constants,rotmatr = _calcabc(_comshift(temp_mat))
				print coordmat
				#print constants
				newA = constants[0]*(float(A_exp)/ns_constants[0])
				newB = constants[1]*(float(B_exp)/ns_constants[1])
				newC = constants[2]*(float(C_exp)/ns_constants[2])
				out.write(label+'	'+str(round(newA,3))+' '+str(round(newB,3))+' '+str(round(newC,3))+'\r\n')

		out.write('\r\n')
	out.close()

# ==============	END ISO_SCALE BLOCK   =================== #



# ==============	BEGIN ARB_PUlSE BLOCK   =================== #

# ========================== ARB_PULSE ======================== \
# This routine will generate a chirped pulse for usage on an    |
# arbitrary waveform generator. This was coded primarily for    |
# use on Tektronix AWGs, such as the 7000 series. Functionality |
# on other AWGs is not guaranteed.                              |
# ------------------------------------------------------------- \
#                          USAGE:                               |
# arb_pulse(dict Opt, dict Time, int frames, str out_name)      |
# Input: Opt and Time are dictionaries with the following       |
# format:                                                       |
# Opt = {'Chirp_Start': float start, 'Chirp_Stop': float end,   |
# 'Chirp_Duration': float duration, 'Sample_Rate': float rate}  |
#                                                               |
# Mark = {'DELAY': float delay, 'M1_WIDTH': float width1,       |
# 'M2_WIDTH': 1.0, 'M_PULSE_BUFFER': float pbuff, 'BUFFER':     |
# float buffer, 'PREBUFFER': float prebuff}                     |
#                                                               |
# frames: integer specifying how many chirp frames shall be     |
# made in output pulse                                          |
#                                                               |
# out_name: output file name of choice. Defaults to None.       |
# If none, then arb_pulse merely returns an array with the      |
# pulse. If out_name is a string, then it will save a file      |
# and not return the array.                                     |
#                                                               |
# More information on the variables contained in these dicts    |
# can be found in the code below contained in commented lines.  |
# --------------------------------------------------------------\


# The following diagram shows the layout of the timing variables:
# ==================================================== t increasing --> =================================================================
#
#
# ------------------|                  |---------------|\/\/\/\/\/\/\/\/\/\|---------------|                     |-----------------|
# ------------------|                  |---------------|/\/\/\/\/\/\/\/\/\/|---------------|                     |-----------------|
# <--PREBUFFER---->    <---M1_WID-->     <-M_P_BUF-->   <---CHIRP_DUR-----> <--M_P_BUFF---> <------M2_WIDTH----->                  ^ BUFFER
#  

def _pulse(t, start, stop, duration):
	return math.sin((2*math.pi*(start*10e5)*t)+(2*math.pi*((stop - start)*10e5)*((t**2)/(2*duration*10e-7))))

def _marker(c_on,c_off,s_rate,total_points):
	n_on = int(math.floor(c_on*10e-7*s_rate*10e8))
	n_off = int(math.floor(c_off*10e-7*s_rate*10e8))

	marker = []
	for n in range(0, n_on-1):
		marker.append(0)
	for n in range(n_on-1,n_off-1):
		marker.append(1)
	for n in range(n_off-1,total_points):
		marker.append(0)
	return marker

def _waveform(start,stop,delay,s_rate,duration,total_points):
	N_on = int(math.floor(delay*10e-7*s_rate*10e8)) 
	chirp = []
	for n in range(0, N_on):
		chirp.append(0)
	N_chirp = int(math.ceil(duration*10e-7*s_rate*10e8))
	for n in range(0, N_chirp):
		t = n*(s_rate*10e8)**(-1)
		chirp.append(_pulse(t,start,stop,duration))
	for n in range((N_chirp+N_on),total_points):
		chirp.append(0)
	return chirp


def arb_pulse(Opt, Mark, frames, out_name=None):

	# Initialize options and marker channels

	# These mark the starting and ending frequencies of the chirp, in MHz
	Chirp_Start = float(Opt['Chirp_Start'])
	Chirp_Stop  = float(Opt['Chirp_Stop'])

	# Length of chirp, in microseconds
	Chirp_Duration = float(Opt['Chirp_Duration'])

	# Sampling rate, in GS/s
	Sample_Rate = float(Opt['Sample_Rate'])

	# Chirp delay -- zero padding at beginning of chirp, in microseconds
	Chirp_Delay = float(Mark['DELAY'])

	CH1_ON = Chirp_Delay-float(Mark['M1_WIDTH'])-float(Mark['M_PULSE_BUFFER'])+float(Mark['PREBUFFER'])
	CH1_OFF = Chirp_Delay-float(Mark['M_PULSE_BUFFER'])
	CH2_ON = Chirp_Delay+Chirp_Duration+float(Mark['M_PULSE_BUFFER'])
	CH2_OFF = CH2_ON+float(Mark['M2_WIDTH'])
	# Initialize markers channel timings, in microseconds

	# Buffer time -- extra time zero padding at end of chirp, in microseconds
	BUFFER = float(Mark['BUFFER'])-CH2_OFF

	# Initialize number of points in a single pulse based on the timings initialized previously
	def waveform_time():
		t = Chirp_Delay + Chirp_Duration
		if CH1_OFF > t:
			t = CH1_OFF
		if CH2_OFF > t:
			t = CH2_OFF
		t = t + BUFFER
		return t

	time = waveform_time()
	total_points = int(math.ceil(time*(10**(-6))*Sample_Rate*(10**9)))

	a = _waveform(Chirp_Start,Chirp_Stop,Chirp_Delay,Sample_Rate,Chirp_Duration,total_points)
	b = _marker(CH1_ON,CH1_OFF,Sample_Rate,total_points)
	c = _marker(CH2_ON,CH2_OFF,Sample_Rate,total_points)
	print len(a)

	if out_name == None:
		AWG_data = zeros((frames*len(a),3))
		for n in range(frames):
			for i in range(len(a)):
				AWG_data[n*len(a)+i,0] = a[i]
				AWG_data[n*len(a)+i,1] = b[i]
				AWG_data[n*len(a)+i,2] = c[i]

		return AWG_data

	if isinstance(out_name,str):
		output_filename = str(out_name)
		AWG_data = open(output_filename,'w')
		for n in range(frames):
			for i in range(len(a)):
				AWG_data.write(str(a[i])+'\t'+str(b[i])+'\t'+str(c[i])+'\n')
		AWG_data.close()






