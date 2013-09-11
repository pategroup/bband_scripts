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
# ARB_PULSE_GEN, SPECTROGRAM, PEAKPICK (SPLINE/PT-SP-LIM)         |
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
# =============================================================== \

# IMPORT BLOCK #
from numpy import *
from scipy import *
import scipy.fftpack as sfft
import matplotlib as mpl
import fileinput
import re, sys
import os
from matplotlib import pyplot as plt
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
	try:
		if re.match("[A-Za-z]", symbol): # regex for atomic symbols
			if symbol == 'H':
				return 1.007825037
			if symbol == 'C':
				return 12.00000000
			if symbol == 'O':
				return 15.99491464
			if symbol == 'N':
				return 14.003074008
			if symbol == 'Cl':
				return 34.96885268
			if symbol == 'F':
				return 18.99840322
			if symbol == 'Si':
				return 27.9769265325
			else:
				raise ValueError('BLANK_UNKNOWN_SYMBOL')

		elif re.match("[0-9]", symbol): # regex for atomic numbers
			if symbol == '1':
				return 1.007825037
			if symbol == '6':
				return 12.00000000
			if symbol == '8':
				return 15.99491464
			if symbol == '7':
				return 14.003074008
			if symbol == '17':
				return 34.96885268
			if symbol == '9':
				return 18.99840322
			if symbol == '14':
				return 27.9769265325
			else:
				raise ValueError('BLANK_UNKNOWN_SYMBOL')
	except ValueError:
		print 'ERROR! An atomic symbol is blank or not recognized.'
		raise

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

	file.write("--------------\n"+"ROTATIONAL CONSTANTS\n"+"--------------\n\n"+"A : "+str(const[0])+" MHz "+"   "+"A : "+str(const[1])+" MHz "+"   "+"C : "+str(const[2])+" MHz "+"   \n\n"+"--------------\n"+"ROTATION MATRIX\n"+"--------------\n\n")

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
