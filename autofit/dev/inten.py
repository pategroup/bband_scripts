from numpy import *
from numpy import linalg as la
from scipy.optimize	import curve_fit, leastsq
import math
import time
from re import search
import subprocess
import textwrap

def tdipole(J_MAX,M,mu_a,mu_b,mu_c):
	# Make sure M and J_MAX are ints
	J_MAX = int(J_MAX)
	M = int(M)
	
	max = (J_MAX+1)**2
	td_mat = zeros((max,max),dtype='complex')
	
	qn_list = zeros(((J_MAX+1)**2,2))
	for i in range(int(math.fabs(M)),J_MAX+1):
		for k in range(-i,i+1):
			n = i**2+i+k
			qn_list[n,0] = i
			qn_list[n,1] = k
	
	
	for i in range(0, max):
		for j in range(0,max):
			j1 = qn_list[i,0]
			k1 = qn_list[i,1]
			
			j2 = qn_list[j,0]
			k2 = qn_list[j,1]
			
			sm = 0.0
			if j2 == j1 + 1:
				sm = (4*(j1+1)*math.sqrt((2*j1+1)*(2*j1+3)))**(-1.0)*2*math.sqrt((j1+1)**2-M**2)
			elif j2 == j1 and j1 != 0:
				sm = (4*j1*(j1+1))**(-0.50)*2*M
			elif j2 == j1 - 1 and j1 != 0:
				(4*j1*math.sqrt(4*j1**2-1))**(-1.0)*-2*math.sqrt(j1**2-M**2)
			#print str(sm) + ' ' + '('+str(j1)+','+str(j2)+')'
			sma = 0.0
			if k1 == k2 and j2 == j1 + 1 and mu_a != 0:
				sma = mu_a*sm*2*math.sqrt((j1+1)**2-k1**2)
			elif k1 == k2 and j2 == j1 and mu_a != 0:
				sma = mu_a*sm*2*k1
			elif k1 == k2 and j2 == j1 + 1 and mu_a != 0:
				sma == mu_a*sm*-2*math.sqrt(j1**2-k1**2)
			
			smb = 0.0
			if k2 == k1 + 1 and j2 == j1 + 1 and mu_b != 0:
				smb = mu_b*sm*-1*math.sqrt((j1+k1+1)*(j1+k1+2))
			elif k2 == k1 - 1 and j2 == j1 + 1 and mu_b != 0:
				smb = mu_b*sm*1*math.sqrt((j1-k1+1)*(j1+k1+2))
			elif k2 == k1 + 1 and j2 == j1 and mu_b != 0:
				smb = mu_b*sm*math.sqrt(j1*(j1+1)-k1*(k1+1))
			elif k2 == k1 - 1 and j2 == j1 and mu_b != 0:
				smb = mu_b*sm*math.sqrt(j1*(j1+1)-k1*(k1-1))
			if k2 == k1 + 1 and j2 == j1 - 1 and mu_b != 0:
				smb = mu_b*sm*-1*math.sqrt((j1-k1)*(j1-k1-1))
			elif k2 == k1 - 1 and j2 == j1 - 1 and mu_b != 0:
				smb = mu_b*sm*1*math.sqrt((j1+k1)*(j1+k1-1))			
			
			smc = 0.0 
			if k2 == k1 + 1 and j2 == j1 + 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*-1*math.sqrt((j1+k1+1)*(j1+k1+2)))
				#print smc
			elif k2 == k1 - 1 and j2 == j1 + 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*1*math.sqrt((j1-k1+1)*(j1-k1+2))*-1)
				#print smc
			elif k2 == k1 + 1 and j2 == j1 and mu_c != 0:
				smc = complex(0,mu_c*sm*math.sqrt(j1*(j1+1)-k1*(k1+1)))
			elif k2 == k1 - 1 and j2 == j1 and mu_c != 0:
				smc = complex(0,mu_c*sm*math.sqrt(j1*(j1+1)-k1*(k1-1))*-1)
			if k2 == k1 + 1 and j2 == j1 - 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*-1*math.sqrt((j1-k1)*(j1-k1-1)))
			elif k2 == k1 - 1 and j2 == j1 - 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*1*math.sqrt((j1+k1)*(j1+k1-1))*-1)
				#print smc

			sm_sum = sma + smb + smc
			td_mat[i,j] = sm_sum
			
	return td_mat

set_printoptions(precision=4,linewidth=125)
#num = 3+1j
#print  str(num) + ' '+ str(abs(num))
tmom = tdipole(1,0,1,0,1)	
print tmom
