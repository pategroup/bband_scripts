from numpy import *
from numpy import linalg as la
import math

def symtop_qn(J_MAX,M):
	qn_list = zeros(((J_MAX+1)**2,2))
	for i in range(M,J_MAX+1):
		for k in range(-i,i+1):
			n = i**2+i+k
			qn_list[n,0] = i
			qn_list[n,1] = k
	return qn_list

def tdipole(J_MAX,M):
	# Make sure M and J_MAX are ints
	J_MAX = int(J_MAX)
	M = int(M)
	
	max = (J_MAX+1)**2
	td_mat = zeros((max,max),dtype='complex')
	
	qn_list = symtop_qn(J_MAX,M)	
	
	for i in range(0, max):
		for j in range(0,max):
			j1 = qn_list[i,0]
			k1 = qn_list[i,1]
			
			j2 = qn_list[j,0]
			k2 = qn_list[j,1]
			
			sm = 0.0
			if j2 == j1 + 1:
				sm = (4*(j1+1.0)*math.sqrt((2*j1+1.0)*(2*j1+3.0)))**(-1.0)*2*math.sqrt((j1+1)**2-M**2)
			elif j2 == j1 and j1 != 0:
				sm = 1/(4*j1*(j1+1))*2*M
			elif j2 == j1 - 1 and j1 != 0:
				sm = (4.0*j1*math.sqrt(4.0*j1**2-1.0))**(-1.0)*-2*math.sqrt(j1**2-M**2)

			sma = 0.0
			if k1 == k2 and j2 == j1 + 1 and mu_a != 0:
				sma = mu_a*sm*2*math.sqrt((j1+1)**2-k1**2)
				#print 'I see this mu_a j1+1  ' + str(sma)
			elif k1 == k2 and j2 == j1 and mu_a != 0:
				sma = mu_a*sm*2*k1
				#print 'I see this mu_a j1  ' + str(sma)
			elif k1 == k2 and j2 == j1 - 1 and mu_a != 0:
				sma = mu_a*sm*-2.0*math.sqrt(j1**2-k1**2)
							
			smb = 0.0
			if k2 == k1 + 1 and j2 == j1 + 1 and mu_b != 0:
				smb = mu_b*sm*-1*math.sqrt((j1+k1+1)*(j1+k1+2))
			elif k2 == k1 - 1 and j2 == j1 + 1 and mu_b != 0:
				smb = mu_b*sm*1*math.sqrt((j1-k1+1)*(j1-k1+2))
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
				print smc
			elif k2 == k1 - 1 and j2 == j1 + 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*1*math.sqrt((j1-k1+1)*(j1-k1+2))*-1)
				print smc
			elif k2 == k1 + 1 and j2 == j1 and mu_c != 0:
				smc = complex(0,mu_c*sm*math.sqrt(j1*(j1+1)-k1*(k1+1)))
				print smc
			elif k2 == k1 - 1 and j2 == j1 and mu_c != 0:
				smc = complex(0,mu_c*sm*math.sqrt(j1*(j1+1)-k1*(k1-1))*-1)
				print smc
			if k2 == k1 + 1 and j2 == j1 - 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*-1*math.sqrt((j1-k1)*(j1-k1-1)))
				print smc
			elif k2 == k1 - 1 and j2 == j1 - 1 and mu_c != 0:
				smc = complex(0,mu_c*sm*1*math.sqrt((j1+k1)*(j1+k1-1))*-1)
				print smc

			sm_sum = sma + smb + smc
			td_mat[i,j] = sm_sum
	
	output = zeros((max-M**2,max-M**2),dtype='complex')
	for i in range(M**2,max):
		for j in range(M**2,max):
			output[i-M**2,j-M**2] = td_mat[i,j]
	return output, qn_list

def ham(J_MAX, rot_A, rot_B, rot_C, DJ,DJK,DK,subdj,subdk):
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


def asym_top(J_MAX, M):
	max = (J_MAX+1)**2

	temp = zeros((max,max))

	for j in range(abs(M),J_MAX+1):
		ham_j = ham(j, rot_A, rot_B, rot_C, d1,d2,d3,d4,d5)
		for k1 in range(-j,j+1):
			for k2 in range(-j,j+1):
				l = k1+j
				m = k2+j
				lfull = j**2+j+k1
				mfull = j**2+j+k2
				temp[lfull,mfull] = ham_j[l,m]
	output = zeros((max-M**2,max-M**2))
	for i in range(M**2,max):
		for j in range(M**2,max):
			output[i-M**2,j-M**2] = temp[i,j]
	return output

def label_states(energies,vecs,J_MAX,M):
	qn_list = symtop_qn(J_MAX,M)
	labels = zeros((shape(energies)[0],6))
	for i in range(shape(energies)[0]):
		labels[i,0] = energies[i]

		J_avg = 0
		magK_avg = 0
		K_avg = 0
		for m in range(0, shape(energies)[0]):
			J_avg += abs(vecs[m,i])**2*qn_list[m,0]
			magK_avg += abs(vecs[m,i])**2*qn_list[m,1]
			print str(abs(vecs[m,i])**2)+' <--- THIS SHIT FOR N = '+str(i)
			K_avg += vecs[m,i]*qn_list[m,1]
		J = round(J_avg)
		Ka = round(magK_avg)

		labels[i,1] = J
		labels[i,2] = Ka

		if Ka % 2 == 1 and round(K_avg) != 0:
			Kc = J - Ka + 1
		elif Ka % 2 == 1 and round(K_avg) == 0:
			Kc = J - Ka
		elif Ka % 2 == 0 and round(K_avg) != 0:
			Kc = J - Ka + 1
		elif Ka % 2 == 0 and round(K_avg) == 0:
			Kc = J - Ka  
		if Ka == 0:
			Kc = J - Ka

		labels[i,3] = Kc
		labels[i,4] = M
		labels[i,5] = magK_avg

	return labels





def calc_spec(J_MAX,M):
	ham_out = asym_top(J_MAX,M)
	td_out = tdipole(J_MAX,M)

	en, vec = la.eigh(ham_out)
	erow = transpose(en)
	print shape(en)[0]
	labels = label_states(en,vec,J_MAX,M)

	return labels

rot_A = 9769.72213
rot_B = 868.846659
rot_C = 818.518746
d1 = 0.000047239 #DJ
d2 = -0.0008991 #DJK
d3 = 0.023173 #DK
d4 = 5.0298E-06 #dJ
d5 = 0.000343 #dK
mu_a = 1.0
mu_b = 1.0
mu_c = 1.0
set_printoptions(precision=3,linewidth=200,suppress=True)
print calc_spec(2,0)
print symtop_qn(2,0)

#num = 3+1j
#print  str(num) + ' '+ str(abs(num))
#tmom, qnlist = tdipole(3,2,1,1,1)	

#print tmom
#print qnlist
