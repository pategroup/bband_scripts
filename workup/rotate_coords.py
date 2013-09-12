# Coordinate Rotator & Rotational Constant Output
# --------- For use with Psi4 outputs ----------
# --------- Author: Nathan Seifert -------------
# --------- Licensed under the GPLv2 License ---
# TO DO: grep a Psi4 optimization output and get the coordinates
from numpy import *
from numpy import linalg as la
import fileinput
import re, sys
import os


# --------------- USAGE -------------
# In a terminal, type "python rotate_coords.py <filename>". 
# Currently this does not parse a raw Psi4 output.
# You must supply a coordinate file with the following output:
# row i:  atom L_i   x_i  y_i  z_i
# where L_i is the atomic symbol (C or 6 for carbon, etc)
# NOTE: the rows can be space or tab delimited!



def masses(symbol):
	#print symbol
	# Check to see if symbol string is a number or letter
	if re.match("[A-Za-z]", symbol):
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
	elif re.match("[0-9]+", symbol):
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

# comshift(mat) is the initial step: shifts coordinates from a coordinates matrix of the coordmat form into center of mass frame
def comshift(matr):
	total_mass = sum(matr[:,0])
	
	#Creates dummy matrices to calc COM shift and are thrown away later to re-shift matr[] values
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
	
	# Shifts matr[] to COM frame
	for i, row in enumerate(matr):
		matr[i,1] -= com_x
		matr[i,2] -= com_y
		matr[i,3] -= com_z
	
	return matr
	
# Takes a COM-shifted coordinate matrix, calculates the moment of inertia tensor, and diagonalizes to rotate to principal axis
def calcabc(matr):
	
	first_mom = zeros((3,3)) # Empty first moment tensor
	tot_atom_count = size(matr[:,0]) 
	temp = zeros((tot_atom_count,9))
	
	# Calculates atomic contributions to each tensor element
	for i in range(tot_atom_count):

		temp[i,0]= matr[i,0]*(matr[i,2]**2 + matr[i,3]**2)
		temp[i,1]= -1*(matr[i,0]*matr[i,1]*matr[i,2])
		temp[i,2]= -1*(matr[i,0]*matr[i,1]*matr[i,3])
		temp[i,3]= -1*(matr[i,0]*matr[i,2]*matr[i,1])
		temp[i,4]= matr[i,0]*(matr[i,1]**2 + matr[i,3]**2)
		temp[i,5]= -1*(matr[i,0]*matr[i,2]*matr[i,3])
		temp[i,6]= -1*(matr[i,0]*matr[i,3]*matr[i,1])
		temp[i,7]= -1*(matr[i,0]*matr[i,3]*matr[i,2])
		temp[i,8]= matr[i,0]*(matr[i,1]**2 + matr[i,2]**2)
            
    #sums all the atomic contributions to form total tensor
	#NOTE: There's gotta be an easier way to do this step and the above step in one for loop
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
	diag_eigenvecs,diag_eigenvals = la.eigh(first_mom)
	
	# Goes from amu A**2 to MHz
	rotconst = 505379.006/diag_eigenvecs
	
	return rotconst, diag_eigenvals #Returns a 3-member array with the proper values of A, B and C

def rot_coords(coords,rotmat):
	coords_mod = delete(coords,s_[:-3],1)
	coords_mod = dot(coords_mod,rotmat)
	#masses = coords[:,0].reshape((size(coords[:,0])),1)
	#coords_mod_mass = hstack((masses,coords_mod))
	return coords_mod

def write(const,rotmat,syms,rot_coordmat, filename):
	file = open(filename,"wb")
	file.write("Input file name: "+sys.argv[-1]+"\n")
	file.write("--------------\n"+"ROTATIONAL CONSTANTS\n"+"--------------\n\n"+"A : "+str(const[0])+" MHz "+"   "+"B : "+str(const[1])+" MHz "+"   "+"C : "+str(const[2])+" MHz "+"   \n\n"+"--------------\n"+"ROTATION MATRIX\n"+"--------------\n\n")
	
	#Writes rotation matrix
	for i in range(shape(rotmat)[0]):
		for j in range(shape(rotmat)[1]):
			file.write(str(rotmat[i,j]) + "   ")
		file.write("\n")
	file.write("\n--------------\n")
	file.write("ROTATED COORDINATES\n"+"--------------\n\n")
	
	#Writes rotated coordinates without original atomic symbols
	for i in range(size(syms)):
		file.write(str(syms[i]) + "   ")
		for j in range(shape(rot_coordmat)[1]):
			file.write(str(rot_coordmat[i,j]) + "   ")
		file.write("\n")
	file.close()
				
def main():
	a = []
	b = []
	c = []
	mass = []
	symbols = []
	print "Input file: "+sys.argv[-1]
	print "Output file:"+os.path.splitext(sys.argv[-1])[0]+"_out"+".geom"
	out_file_name = os.path.splitext(sys.argv[-1])[0]+"_out"+".geom"
	print "--------------------"
	for line in fileinput.input(): # Collects X,Y,Z,M into separate arrays from input
		a.append(line.split()[1]) 
		b.append(line.split()[2])
		c.append(line.split()[3])
		symbols.append(line.split()[0])
		mass.append(masses(line.split()[0]))
	coordmat = zeros((size(mass),4)) # Empty matrix for X,Y,Z,M matrix
	for i,row in enumerate(mass):
		coordmat[i,0] = row
		coordmat[i,1] = a[i]
		coordmat[i,2] = b[i]
		coordmat[i,3] = c[i]
	# calcabc(mat) takes a matrix of coordinates and does the principal axis rotation
	constants,rotmatr = calcabc(comshift(coordmat)) 
	
	write(constants,rotmatr,symbols,rot_coords(coordmat,rotmatr),out_file_name)



if __name__ == "__main__":
	set_printoptions(suppress=True) # Turns off scientific notation
	main()
