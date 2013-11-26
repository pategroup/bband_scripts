#import calexcept
import datetime
import math
import subprocess
import os.path
import os
import decimal
import random
import numpy as np 
#import pandas as pn
import time
import struct
import fileinput
from matplotlib import pyplot as pp

class InitializeError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return rept(self.value)

class ExecuteError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return rept(self.value)

class IdiotCheck(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return rept(self.value)

class NotSupportedException(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return rept(self.value)

#   ______     ___       __      .______     _______ .___  ___. 
#  /      |   /   \     |  |     |   _  \   /  _____||   \/   | 
# |  ,----'  /  ^  \    |  |     |  |_)  | |  |  __  |  \  /  | 
# |  |      /  /_\  \   |  |     |   ___/  |  | |_ | |  |\/|  | 
# |  `----./  _____  \  |  `----.|  |      |  |__| | |  |  |  | 
#  \______/__/     \__\ |_______|| _|       \______| |__|  |__|


class calpgm(object):

	# All available parameters for SPCAT

	__PARAMS_RIGID = {'A':10000,'B':20000,'C':30000}
	__PARAMS_DISTORT_A = {'DelJ': 200,'DelJK':1100,'DelK':2000,'delJ':40100,'delK':41000}
	__PARAMS_DISTORT_S = {'DJ': 200, 'DJK':2000,'DK':2000,'d1':40100,'d2':50000}

	__HEX_DISTORT_A = {'PhiJ': 300,'PhiJK':1200,'PhiKJ':2100,'PhiK':3000,'phiJ':40200,'phiJK':41100,'phiK':42000}
	__HEX_DISTORT_S = {'HJ': 300, 'HJK': 1200,'HKJ': 2100,'HK':3000,'h1':40200,'h2':50100,'h3':60000}

	__OCT_DISTORT_A = {'LJ': 400, 'LJJK': 1330, 'LJK': 2200, 'LKKJ': 3100, 'LK': 4000, 'lJ': 40300, 'lJK': 41200, 'lKJ': 42100, 'lK': 43000} 
	__OCT_DISTORT_S = {'LJ': 400, 'LJJK': 1330, 'LJK': 2200, 'LKKJ': 3100, 'LK': 4000, 'l1': 40300, 'l2': 50200, 'l3': 60100, 'l4': 70000}

	__QUAD = {'1_5Xaa': 110010000, '1_5Xbb': 110020000, '1_5Xcc': 110030000, 'Xab': 110610000, 'Xbc': 110210000, 'Xac': 110410000, '0_25Xb-c':110040000}
	
	ALL_PARAMS = dict(__PARAMS_RIGID.items() + __PARAMS_DISTORT_S.items() + __PARAMS_DISTORT_A.items() + __HEX_DISTORT_A.items() + __HEX_DISTORT_S.items() + __OCT_DISTORT_A.items() + __OCT_DISTORT_S.items() + __QUAD.items()) 



	# PARAMETERS

	# Containers for rotational constants, etc
	initial_vals = []
	initial_vals_rigid = [0.0,0.0,0.0]
	current_vals = []
	current_vals_rigid = [0.0,0.0,0.0]

	def spincalc(self,input_spin):
		if input_spin == math.floor(float(input_spin)):
			input_spin = int(2*float(input_spin)+1)
		else:
			input_spin = 2*int(math.ceil(float(input_spin)))

		return input_spin

	def read(self, data):
		if self.current_vals:
			self.current_vals = []
		for i in range(0, len(data)):
			self.current_vals.append([])
			self.current_vals[i].append(data[i][0])
			self.current_vals[i].append(data[i][1])
			try:
				self.current_vals[i].append(data[i][2])
			except IndexError:
				self.current_vals[i].append("1.0E-010")

			if data[i][0] == 'A':
				self.current_vals_rigid[0] = float(data[i][1])
			if data[i][0] == 'B':
				self.current_vals_rigid[1] = float(data[i][1])
			if data[i][0] == 'C':
				self.current_vals_rigid[2] = float(data[i][1])

		if not self.initial_vals:
			self.initial_vals = self.current_vals
			self.initial_vals_rigid = self.current_vals_rigid

	def from_file(self, input_file):
		data = open(input_file, 'r')
		output = []
		for i, line in enumerate(data):
			output.append([])
			output[i].append(line.split()[0])
			output[i].append(line.split()[1])
			try:
				output[i].append(line.split()[2])
			except IndexError:
				output[i].append("1.0E-010")
		#print output
		return output


	def add_spcat_params(self,new_dict):
		self.ALL_PARAMS = dict(ALL_PARAMS.items()+new_dict.items())

	def get_vals(self):
		return self.current_vals

	def get_init_vals(self):
		return self.initial_vals


	def error_message(self,errortype, message, severity=0):
		print '\n\n=========== '+ str(errortype)+ 'ERROR AT: '+datetime.datetime.now().strftime("%a %b %d %I:%M:%S %Y")+'==========='
		print ""
		print 'You screwed up! This is why:'
		print '-----------------'
		print str(message)
		print '-----------------'
		if severity == 0:
			print 'Continue with caution!\n\n'
		if severity == 1:
			print 'Honestly you shouldn\'t have even gotten here! Shame on you!'
		if severity == 2:
			print 'This probably breaks your routine. Check your code!!!'

	# List of kwargs relevant for calpgm (so far)
	# - data: input parameters (list of lists)
	# - name: molecule name (default: "molecule")
	# - max_freq: max frequency for predictions (default: 20 GHz)
	# - dipoles: [ua,ub,uc] in floats
	# - temp : temperature in Kelvin (default 2.0)
	# - spin : nuclear spin (e.g. 1 for nitrogen 14-containing molecules, 1.5 for chlorine, default 0)
	# - reduction: 'a' or 's' (specifies which watson reduction to use, default 'a')
	# - J_min/J_max : min/max J for predictions (0/20 default)
	# - inten: intensity cutoff (log strength, default -10.0)
	def __init__(self,**kwargs):

		self.name = "molecule"
		self.filename = "default"
		self.max_freq = 20.0
		self.dipoles = [1.0,1.0,1.0]
		self.temp = 2
		self.spin = 0
		self.reduction = 'a'
		self.J_min = 0
		self.J_max = 20
		self.inten = -10.0
		self.temp = 2.0

		print 'CALPGM constructor initialized\n'
		#self.spin = self.spincalc(self.spin)
		try:
			for key, value in kwargs.iteritems():

				if key == 'data':
					if isinstance(value,basestring):
						self.read(self.from_file(value))
					else:
						self.read(value)

				elif key == 'filename':
					if isinstance(value,basestring):
						self.filename = value

				elif key == 'max_freq':
					self.max_freq = float(value)

				elif key == 'dipoles':
					self.dipoles = value

				elif key == 'temp':
					self.temp = value

				elif key == 'spin':
					self.spin = self.spincalc(value)

				elif key == 'reduction':
					self.reduction = value

				elif key == 'J_min' or key == 'j_min':
					self.J_min = value

				elif key == 'J_max' or key == 'j_max':
					self.J_max = value

				elif key == 'inten' or key == 'intensity':
					self. inten = value

				elif key == 'new_params':
					try:
						if isinstance(value,dict):
							self.add_params(value)
						else:
							raise InitializeError(str(key)+ ' is not a valid dictionary.')
							pass
					except InitializeError as e:
						self.error_message("InitializeError", e.value,0)


				else:
					raise InitializeError(str(key)+' is not a valid argument of calpgm()')
					pass 


		except InitializeError as e: 
			self.error_message("InitializeError",e.value,0)







#      _______..______     ______      ___   .___________.
#     /       ||   _  \   /      |    /   \  |           |
#    |   (----`|  |_)  | |  ,----'   /  ^  \ `---|  |----`
#     \   \    |   ___/  |  |       /  /_\  \    |  |     
# .----)   |   |  |      |  `----. /  _____  \   |  |     
# |_______/    | _|       \______|/__/     \__\  |__|    



class spcat(calpgm):

	def qrotcalc(self):
	# Calculates pure, low-temp-approx rotational partition function
		A = self.current_vals_rigid[0]
		B = self.current_vals_rigid[1]
		C = self.current_vals_rigid[2]

		return round((5.3311*10**6)*self.temp**(1.5)*(A*B*C)**(-0.5),3)

	def to_var(self):
	# Writes var file
		num_params = len(self.current_vals)

		timestamp = datetime.datetime.now().strftime("%a %b %d %I:%M:%S %Y")

		output  = "%s                                        %s \n" %(self.name,timestamp)
		output += "   %s  999   51    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n" %(str(num_params))
		output += "%s   %s  1  0  99  0  1  1  1  1  -1   0\n" %(self.reduction,str(self.spin))
		for i in range(0, num_params):

			print self.current_vals[i][0]
			key_value = str(self.ALL_PARAMS[self.current_vals[i][0]])
			output += "          %s  %s  %s  /%s \n" %(str(key_value),str(self.current_vals[i][1]),str(self.current_vals[i][2]),str(self.current_vals[i][0]))

		if self.init_var == "":
			self.init_var = output
		self.cur_var = output


	def to_int(self):
	# Writes int file
		 self.pfunc = self.qrotcalc()

		 output  = "%s \n"%(self.name)
		 output += "0  91  %s  %s  %s  %s  %s %s  %s\n"%(str(self.pfunc), str(self.J_min), str(self.J_max),"-10","-10",str(self.max_freq), str(self.temp))
		 output += " 001  %s \n" %(str(self.dipoles[0]))
		 output += " 002  %s \n" %(str(self.dipoles[1]))
		 output += " 003  %s \n" %(str(self.dipoles[2]))

		 if self.init_int == "":
		 	self.init_int = output
		 self.cur_int = output


	def update(self): # Updates var/int in case a parameter is changed (e.g. temp)
		self.to_var()
		self.to_int()

	def get_var(self): # Returns current var file string
		return self.cur_var

	def get_var_init(self): # Returns initial var file string 
		return self.init_var

	def get_int(self): # Returns current int file string
		return self.cur_int

	def get_int_init(self): # Returns initial int file string
		return self.init_int

	def to_file(self,**kwargs): # Writes var and int files to a file, filename specified by self.filename or can be overwritten using filename="new_filename" as an argument.
		if 'filename' in kwargs:
			out_name = kwargs['filename']
		else: 
			out_name = self.filename

		for key, value in kwargs.iteritems():

			if key == 'type' and value == 'var':
					if 'v' in kwargs:
						if kwargs['v'] == 'init' or kwargs['v'] == 'initial' or kwargs['v'] == 'i':
							output = open(out_name+'.var','wb')
							output.write(self.init_var)
							output.close()
						elif kwargs['v'] == 'cur' or kwargs['v'] == 'current' or kwargs['v'] == 'c':
							output = open(out_name+'.var','wb')
							output.write(self.cur_var)
							output.close()
						else:
							output = open(out_name+'.var','wb')
							output.write(self.cur_var)
							output.close()

			if key == 'type' and value == 'int':
					if 'v' in kwargs:
						if kwargs['v'] == 'init' or kwargs['v'] == 'initial' or kwargs['v'] == 'i':
							output = open(out_name+'.int','wb')
							output.write(self.init_int)
							output.close()
						elif kwargs['v'] == 'cur' or kwargs['v'] == 'current' or kwargs['v'] == 'c':
							output = open(out_name+'.int','wb')
							output.write(self.cur_int)
						else:
							output = open(out_name+'.var','wb')
							output.write(self.cur_var)
							output.close()

	def execute(self, **kwargs):
	# Executes SPCAT. Takes in two different optional arguments:
	# filename=" " <--- overrides class filename in case you want to run SPCAT on another file
	# v = 'init' / 'initial' / 'i' <--- runs SPCAT on the initial var/int put into the spcat() object (self.init_int/self.init_var). Good if you want to compare vs new var/int iteration
	# or....
	# v = 'cur' / 'current' / 'c' <--- DEFAULT. runs SPCAT on current var/int in spcat() object, self.cur_var and self.cur_int
	# update = 1 / 0 (default = 0). If set to 1, execute will rerun to_var() and to_int(). This is essential if you update a variable in the class object, such as dipoles or the temperature,
	# and want to get a new prediction. 

		try:
			if 'filename' in kwargs:
				output_name = kwargs['filename']
			else: 
				output_name = self.filename

			if 'v' in kwargs:
				if kwargs['v'] == 'init' or kwargs['v'] == 'initial' or kwargs['v'] == 'i':
					if self.init_int != "" and self.init_var != "":

						if 'update' in kwargs:
							if kwargs['update'] == 1:
								self.to_var()
								self.to_int()
							if kwargs['update'] == 0:
								pass

						self.to_file(type='var',filename=output_name,v='i')
						self.to_file(type='int',filename=output_name,v='i')

						# For *nix systems:
						a = subprocess.Popen("./spcat "+output_name,stdout=subprocess.PIPE,shell=True)
						a.stdout.read()
						

					if self.init_int == "" or self.init_var == "":
						raise ExecuteError('Empty int or var during execute step')

				elif kwargs['v'] == 'cur' or kwargs['v'] == 'current' or kwargs['v'] == 'c':
					if self.cur_var != "" and self.cur_int != "":

						if 'update' in kwargs:
							if kwargs['update'] == 1:
								self.to_var()
								self.to_int()
							if kwargs['update'] == 0:
								pass

						self.to_file(type='var',filename=output_name,v='c')
						self.to_file(type='int',filename=output_name,v='c')
						
						# For *nix systems:
						a = subprocess.Popen("./spcat "+output_name,stdout=subprocess.PIPE,shell=True)
						a.stdout.read()
						

					if self.cur_int == "" or self.cur_var == "":
						raise ExecuteError('Empty int or var during execute step')

		except ExecuteError as e:
			self.error_message("ExecuteError",e.value,1)




	def read_cat(self, **kwargs):
	# Returns a list of lists (cat[i][j]) with the following info extracted from cat file
	# - cat[i][0] : freq  <--- frequency of transition in MHz
	# - cat[i][2] : inten <--- -log(10) of intensity
	#
	# - if system does not contain nuclear quad:
	# - cat[i][3:5] : J_up / Ka_up / Kc_up <--- QNs of upper state
	# - cat[i][6:8] : J_down / Ka_down / Kc_down <--- QNs of lower state
	# - cat[i][9] : if component=1 (see below), this will either be 'a' or 'b' or 'c'

	# - if system does contain nuclear quad:
	# - cat[i][3:6] : J_up / Ka_up / Kc_up / F_up
	# - cat[i][7:10]: J_down / Ka_down / Kc_down / F_down
	# - cat[i][1] : uncert <--- line uncertainties (0 unless you set uncertainties in constants in var file)
	# - cat[i][11] : if components=1 (see below), this will either be 'a' or 'b' or 'c'

	# Possible input arguments:
	# min_freq = float /GHz <--- default 0.0, sets minimium frequency for appending from cat file to output
	# max_freq = float /GHz <--- default 20.0, sets maximum frequency (overrides self.max_freq)
	# min_inten = float <--- default -10.0, overrides self.inten
	# max_inten = float <--- default 0, can't be any larger than this. Probably no reason to change it.
	# component = 0/1 <--- default 0; if 1, read_cat() will append additional row to output with either 'a', 'b', or 'c', specifying the dipole component associated with this transition

	# pretty = 0/1 <--- default 0. If 1, read_cat() returns a Pandas DataFrame instead with the following columns:
	# 'freq' / 'uncert' / 'inten' / 'J_up' / 'Ka_up' / 'Kc_up' / 'J_down' / 'Ka_down' / 'Kc_down'

		# Checks to see if user wants to filter cat by frequency
		if "min_freq" in kwargs:
			min_freq = kwargs['min_freq']
		else:
			min_freq = 0.0

		if "max_freq" in kwargs:
			max_freq = kwargs['max_freq']
		else:
			max_freq = self.max_freq


		# Checks to see if user wants to filter cat by intensity	
		try: 
		
			if "min_inten" in kwargs:
				if min_inten <= 0.0: # log scale idiot check
					min_inten = kwargs['min_inten']
			else: 
				min_inten = self.inten

			if "max_inten" in kwargs:
				if max_inten > min_inten and max_inten <= 0.0: # Idiot check
					max_inten = kwargs['max_inten']
				else:
					max_inten = 0.0
					raise IdiotCheck('Check your intensities!')
			else:
				max_inten = 0.0

		except IdiotCheck as e:
			self.error_message("IdiotCheck",e.value,1)

		try:
			if not os.path.isfile(self.filename+".cat"):
				raise IOError
			else:
				pass 

		except IOError:
			msg = 'Can\'t find a CAT file with the filename:' + self.filename+".cat\n"
			msg += 'SPCAT has not been executed yet for this object. Please run execute() before you read_cat().'
			self.error_message("IOError",msg,2)

		if self.spin == 1:
			names = ['freq','uncert','inten','J_up',"Ka_up","Kc_up","J_down","Ka_down","Kc_down"]
		if self.spin != 1:
			names = ['freq','uncert','inten','J_up',"Ka_up","Kc_up","F_up","J_down","Ka_down","Kc_down","F_down"]

		cat_file = []
		f = open(self.filename+".cat")
		for line in f:
			#print line
			if float(line[3:13]) > min_freq*1000 and float(line[3:13]) < max_freq*1000:
				#print line[13:21]
				if float(line[22:29]) < max_inten and float(line[22:29]) > min_inten:
					#print 'Got here too!'
					if self.spin == 1:
						if 'component' in kwargs and kwargs['component'] == 1:
							if (int(line[57:59])-int(line[69:71])) % 2 == 0 and (int(line[59:61])-int(line[71:73])) % 2 == 1:
								cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[67:69]),int(line[69:71]),int(line[71:73]),'a'])

							if (int(line[57:59])-int(line[69:71])) % 2 == 1 and (int(line[59:61])-int(line[71:73])) % 2 == 1:
								cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[67:69]),int(line[69:71]),int(line[71:73]),'b'])

							if (int(line[57:59])-int(line[69:71])) % 2 == 1 and (int(line[59:61])-int(line[71:73])) % 2 == 0:
								cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[67:69]),int(line[69:71]),int(line[71:73]),'c'])

						else:
							cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[67:69]),int(line[69:71]),int(line[71:73])])
					if self.spin != 1:
						# Routine checks to see if a/b/c-type  
						if 'component' in kwargs and kwargs['component'] == 1:
							if (int(line[57:59])-int(line[69:71])) % 2 == 0 and (int(line[59:61])-int(line[71:73])) % 2 == 1:
								cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[61:64]),int(line[67:69]),int(line[69:71]),int(line[71:73]),int(line[73:75]),'a'])

							if (int(line[57:59])-int(line[69:71])) % 2 == 1 and (int(line[59:61])-int(line[71:73])) % 2 == 1:
								cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[61:64]),int(line[67:69]),int(line[69:71]),int(line[71:73]),int(line[73:75]),'b'])

							if (int(line[57:59])-int(line[69:71])) % 2 == 1 and (int(line[59:61])-int(line[71:73])) % 2 == 0:
								cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[61:64]),int(line[67:69]),int(line[69:71]),int(line[71:73]),int(line[73:75]),'c'])
						else:
							cat_file.append([float(line[3:13]),float(line[13:21]),float(line[22:29]),int(line[55:57]),int(line[57:59]),int(line[59:61]),int(line[61:64]),int(line[67:69]),int(line[69:71]),int(line[71:73]),int(line[73:75])])

						# Writes 
		#if "pretty" in kwargs:
			#if kwargs['pretty'] == 1:
				#names = ['freq','uncert','inten','J_up',"Ka_up","Kc_up","J_down","Ka_down","Kc_down"]
				#df = pn.DataFrame(np.array(cat_file),columns=names)
				#return df

		return cat_file



	def cat_filter(self,catfile, **kwargs): 
	# Filters cat file, generally returned by read_cat(), based on user-supplied specifications. 
	# These specifications include J ranges, Ka ranges, frequency, intensity, and uncert.
	# Since a "cat file" is not a variable of the spcat() object, you will need to supply it a cat_file object, e.g. from read_cat(). It will
	# accept the list of lists format (pretty=0, default), but at the moment NOT the Pandas output (the author hasn't quite figured out how to do this yet!!!)


	# Flags are stored in the dictionary below, and have a three member list, flag[i,j,k,l]; i = 0/1 (on/off), j = corresponding column in catfile, k = inputted filter value, l = 1/-1 (-1 if lower bound, 1 if upper bound -- for boolean filters)
		flags = {'freq_up':[0,0,0.0,1],'freq_min':[0,0,0.0,-1], 'J_max':[0,3,0,1],'J_min':[0,3,0,-1],'min_inten':[0,2,0.0,-1],'max_inten':[0,2,0.0,1],'uncert':[0,1,0.0,1],'Ka_up_max':[0,4,0,1],'Ka_up_min':[0,4,0,-1],'Ka_down_max':[0,7,0,1],'Ka_down_min':[0,7,0,-1],'Kc_up_max':[0,5,0,1],'Kc_down_max':[0,8,0,1],'Kc_up_min':[0,5,0,-1],'Kc_down_min':[0,8,0,-1], 'component':[0,9,0,0]}
		
		# Shifts Ka/Kc labels due to addition of F quantum numbers
		if self.spin != 1:
			flags['Ka_down_max'][1] = 8
			flags['Ka_down_min'][1] = 8
			flags['Kc_down_max'][1] = 9
			flags['Kc_down_min'][1] = 9
			flags['component'][1] = 11

	#	try:
	#		if isinstance(catfile,pn.DataFrame):
	#			raise NotSupportedException('Author has been too lazy to implement Pandas support for cat_filter(). Make sure read_cat is not taking in pretty=1 as an argument.')

	#	except NotSupportedException as e: # Because I haven't figured out how to do filters on DataFrames. Will come soon!!!! 
	#		self.error_message("NotSupportedException",e.value,2)

		temp_dict = {}
		for key in flags: # Sets filter flags based on kwargs input
			if key in kwargs:
				flags[key][0] = 1
				flags[key][2] = kwargs[key]
			if flags[key][0] == 1:
				temp_dict[key] = flags[key]



		#Filter routine
		output = []

		for i in range(0,len(catfile)):
			filter_val = False
			for key in temp_dict:
				if temp_dict[key][3] == 1:
					index = int(temp_dict[key][1])
					if float(catfile[i][index]) > temp_dict[key][2]:
						filter_val = True

				if temp_dict[key][3] == -1:
					index = int(temp_dict[key][1])
					if float(catfile[i][index]) < temp_dict[key][2]:
						filter_val = True

				if temp_dict[key][3] == 0: # Filter by dipole type
					index = int(temp_dict[key][1])
					# Check to see if dipole filter is string -- e.g. only one dipole component selected
					if isinstance(temp_dict[key][2],basestring):
						#print str(catfile[i][index]) + '<--- CATFILE //// DICTIONARY CALL --->' + str(temp_dict[key][2]) + ' INDEX ON CATFILE CALL ---> ' + str(index)
						if catfile[i][index] != temp_dict[key][2]:
							filter_val = True
					else:
						if not catfile[i][index] in temp_dict[key][2]:
							filter_val = True
						
			if not filter_val:
				output.append(catfile[i])

		return output

	# Searches through an input cat file for a given transition
	# catfile is a list of lists, e.g. an output from cat_filter or read_cat()
	# transition is a list of quantum numbers for the target transition of the form: [J', Ka', Kc', J'', Ka'', Kc'']
	# for molecules with a single nuclear quad, transition will be of the form: [J', Ka', Kc', F', J'', Ka'', Kc'', F'']
	def line_search(self,catfile,transition,**kwargs):

		# Make sure quantum numbers are int
		for i in range(0,len(transition)):
			transition[i] = int(transition[i])

		# Filter input cat file to limit search
		to_search = self.cat_filter(catfile,J_max = transition[0], J_min = transition[0], Ka_up_max = transition[1], Ka_up_min = transition[1])
		print "Filter results: "
		for i in range(0, len(to_search)):
			print to_search[i]

		
		for i in range(0, len(to_search)):

				j = 0
				found = True
				
				while j < len(transition) and found:
					if transition[j] != to_search[i][j+3]:
						found = False
					j += 1

				if not found:
					pass
				if found:
					return to_search[i]


	def predict_spec(self,catfile,**kwargs):
		if 'min_freq' in kwargs:
			min_freq = kwargs['min_freq']
		else:
			min_freq = 0.0

		if 'max_freq' in kwargs:
			max_freq = kwargs['max_freq']
		else:
			max_freq = self.max_freq*1000

		if 'linewidth' in kwargs:
			linewidth = kwargs['linewidth']
		else:
			linewidth = 0.03

		if 'step_size' in kwargs:
			step_size = kwargs['step_size']
		else:
			step_size = 0.02


		num_points = math.ceil((max_freq-min_freq)/step_size)

		spectrum = np.column_stack(([catfile[i][0] for i in range(0,len(catfile))],[10**(catfile[i][2]) for i in range(0,len(catfile))]))


		
		t1 = time.time()
		#sp_x = 
		sp = np.column_stack(([min_freq+n*step_size for n in range(0,int(num_points-1))],np.zeros(np.shape(range(0,int(num_points-1)))[0])))
		#for n in range(0, int(num_points-1)):
		#	sp[n,0] = min_freq + n*step_size
		#	sp[n,1] = 0
		t2 = time.time()
		print 'Takes this long to build the spectrum inputs: ' + str(t2-t1)

		for n in range(0, len(spectrum)-1):
			if spectrum[n,0] < min_freq:
				continue 
			if spectrum[n,0] > max_freq:
				break

			else:
				temp_int = spectrum[n,1]
				temp_freq = spectrum[n,0]

				sim_start = temp_freq - 4*linewidth
				sim_stop = temp_freq + 4*linewidth

				nstart = math.floor((sim_start-min_freq)/step_size)
				nstop = math.floor((sim_stop-min_freq)/step_size)

				for i in range(int(nstart),int(nstop)):
					sp[i,1] = sp[i,1] + temp_int*math.exp((((sp[i,0]-temp_freq)/linewidth)**2)*(-1*math.log(2)))
					sp[i,1] = np.real(sp[i,1])

		return sp

		


	def __init__(self, **kwargs):

		self.init_var = ""
		self.init_int = ""

		self.cur_var = ""
		self.cur_int = ""


		self.pfunc = 1.0 # Partition function for the object stored as a float


		super(spcat,self).__init__(**kwargs)
		if 'data' in kwargs:
			if kwargs['data'] != "" and kwargs['data'] != None:
				self.to_var()
				self.to_int()
		pass


#      _______..______    _______  __  .___________.
#     /       ||   _  \  |   ____||  | |           |
#    |   (----`|  |_)  | |  |__   |  | `---|  |----`
#     \   \    |   ___/  |   __|  |  |     |  |     
# .----)   |   |  |      |  |     |  |     |  |     
# |_______/    | _|      |__|     |__|     |__|   


class spfit(spcat):


	def read_fit(self, filename, **kwargs):

		self.errors = [0.0,0.0]

		# Booleans that will cut off reading of file once we get through final step of fit file

		do_RMS = True
		do_vars = False
		do_lines = False

		# Reverses parameter dictionary. This comprehension works in 2.7 and up, but will not work in any version of Python lower.
		inv_params = {v:k for k,v in self.ALL_PARAMS.iteritems()}


		# Iterate through fit file backwards!
		for line in reversed(open(filename).readlines()):
			#print line

			if do_RMS == True:
			# Check for microwave RMS (unitless)
				if line.split()[0] == "END":
					#print "Microwave RMS error: " + line.split()[8]
					self.errors[0] = float(line.split()[8])

				if line.split()[0] == 'MICROWAVE' and line.split()[1] == 'RMS':
					#print 'Got here'
					self.errors[1] = float(line.split()[3])
					do_RMS = False
					do_vars = True

			is_param_line = False
			# Pull parameters
			if do_vars == True:
				try:
					float(line.split()[1])
					is_param_line = True
				except ValueError:
					pass

				if is_param_line:

					if int(line.split()[1]) in inv_params:

						# This will pull out the fit value of the parameter, sans error
						temp_val = ""
						temp_exponent = ""
						temp_error = ""

						temp_line = ""
						for i in range(2,len(line.split())-1):
							temp_line += line.split()[i]

						got_neg = False
						got_val = False
						got_error = False
						got_exponent = False

						for i,c in enumerate(temp_line):
							if c == "-" and got_neg == False and i == 0:
								temp_val += c
								got_neg = True

							# get error
							if c == "(" and not got_error:
								got_val = True
								k = i+1
								while k < len(temp_line):
									if temp_line[k] != ")":
										temp_error += temp_line[k]
										k += 1
									if temp_line[k] == ")":
										got_error = True
										k = len(temp_line)+1

							# get exponent, if there is one
							if c == "E" and not got_exponent:
								got_exponent = True
								got_neg = True # In case val is positive, makes sure not to pull  negative from "change this iteration" column
								temp_exponent += temp_line[i+1:i+4]


							elif not got_val:
								if c == ".":
									temp_val += c

								try:
									if c != "-":
										float(c)
										temp_val += c
								except ValueError:
									pass

							#elif got_neg and got_error and got_val:
							#	break



						num_digits = decimal.Decimal(temp_val).as_tuple().exponent

						if not temp_exponent:
							temp_exponent = "0"
						temp_error = float(temp_error)*(10**(num_digits+int(temp_exponent)))
						temp_error = round(temp_error, -1*num_digits+10)
						temp_val = float(temp_val)*10**(int(temp_exponent))
						temp_val = round(temp_val,-1*num_digits+10)


						self.fit_vars_cur.append([line.split()[1],float(temp_val),temp_error])
						

				# Let's do the linelist now!

				# Make sure we're done reading the parameters		
				if not is_param_line and line.split()[0] == 'NEW' and line.split()[1] == 'PARAMETER':
					do_vars = False

			# Now we need to start reading lines
			if  line.split()[0] == 'NORMALIZED':
				#print 'GOt here'
				do_lines = True

				# Linelist line is organized as so:
				# #LINE: J K K J K K EXP_FREQ CALC_FREQ OMC <junk>
				# #LINE is an int, so we gotta check and then cut off once we hit #LINE = 1
			is_line = False
			if do_lines:
				try:
					int(line.split()[0][:-1]) # line.split()[0] for a line is N: so we gotta kick off the : from the string
					is_line = True
				except ValueError:
					pass
				#self.spin = 1
				#print 'Spin is: '+ str(self.spin)
				if is_line and self.spin == 1:
					self.linelist.append([int(line.split()[1]),int(line.split()[2]), int(line.split()[3]),int(line.split()[4]),int(line.split()[5]),int(line.split()[6]),float(line.split()[7]),float(line.split()[8]),float(line.split()[9])])
				if is_line and self.spin > 1:
					self.linelist.append([int(line.split()[1]),int(line.split()[2]), int(line.split()[3]),int(line.split()[4]),int(line.split()[5]),int(line.split()[6]),int(line.split()[7]),int(line.split()[8]),float(line.split()[9]),float(line.split()[10]),float(line.split()[11])])
					#print "J K K F UPPER: "+ line.split()[1] + " " + line.split()[2] + " " + line.split()[3] + line.split()[4] + " / J K K LOWER: "+ line.split()[5] + " " + line.split()[6] + " " + line.split()[7] + " " + line.split()[8] + " / EXP: " + line.split()[9] + " / PRED: " + line.split()[10] + " / OMC: " + line.split()[11]

			# Stops read_fit since we're done at this point.
			if line.split()[0] == "EXP.FREQ.":
				break 
				#pass

		return self.fit_vars_cur

	def __init__(self,**kwargs):
		self.linelist = []
		self.fit_vars_cur = []
		self.errors = [0.0,0.0]
		self.spin = 0

		super(spfit, self).__init__(**kwargs)

		#self.read_fit()

		


# EXAMPLE BLOCK


# example = spcat(data='data') # Creates instance of SPCAT object. Everything in this API will be controlled through the SPCAT object
							   #spcat.__init__ creates a new var and int object based on the input parameters; in this case 
							   #everything is default (defaults set in calpgm() class), and it points to the file 'data' to get rotational constants

# time1 = time.time()
# example.execute(v='c') # Runs SPCAT on current data & parameters, as above call sets
# time2 = time.time()
# print "It took: "+ str(round((time2-time1)*1000,3)) + " miliseconds to initialize and run\n"

# time1 = time.time()
# outputcat = example.read_cat(pretty=0) # Reads and returns cat file as a list of lists
# time2 = time.time()
# print "It took: "+ str(round((time2-time1)*1000,3)) + " miliseconds to run the reader\n"

# time1 = time.time()
# filtered = example.cat_filter(outputcat, freq_up=10000,freq_min=2000,Ka_up_max=0) # Filters cat file with frequency range 2-10 GHz, and Ka can only be 0
# time2 = time.time()

# print "It took: "+ str(round((time2-time1)*1000,3)) + " miliseconds to run the filter\n"
# for i in filtered:
	# print i
# print "--------"

# filtered = example.cat_filter(outputcat,freq_up=10000,freq_min=2000,Ka_up_max=0,J_max=3) # Does the above, but with J no greater than 3.
# for i in filtered:
	# print i

# example.dipoles = [1.0,0,0] # Update the dipole moments to be just a-type

# example.execute(v='c',update=1) # Reruns SPCAT with update=1, which uses the above dipole changes in calculating the CAT file

# filtered = example.cat_filter(example.read_cat(),freq_up=10000,freq_min=2000,Ka_up_max=0,J_max=3) # Filters the updated CAT file
# print "---------"
# for i in filtered:
	# print i

# example.execute(v='i') # Reruns SPCAT with initial parameters (all dipoles default). Each spcat object will always store the initial var/int file. They can
					     # be changed, of course, by pointing to example.init_var and example.init_int to a new var/int string (such as that from example.get_var())
# filtered = example.cat_filter(example.read_cat(),freq_up=10000,freq_min=2000,Ka_up_max=0,J_max=3)  # Filters the "initial" CAT file
# print "---------"
# for i in filtered:
	# print i

#For instance:
# example.init_var = example.get_var()
# example.init_int = example.get_int()
#Replaces initial var and int in example object with the updated dipoles

# example.execute(v='c')
# filtered = example.cat_filter(example.read_cat(),freq_up=10000,freq_min=2000,Ka_up_max=0,J_max=3) # Results in the same as the third cat_filter call in this example block
# print "---------"
# for i in filtered:
	# print i

#Easy!
