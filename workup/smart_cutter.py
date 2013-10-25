from numpy import *
from scipy.interpolate import *
import os
import sys
import time
# <----- Program Flow ----> 
# 1. Input spectrum
#    a. 2 Input parameters: Approximate noise threshold; Intensity cutoff for peaks
# 			(peak cutoff should be as near as possible to the baseline, or at significantly low enough to not screw up the 
#	noise analysis step
#
# 2. Spline to good resolution -- 2 kHz?
# 3. Peakpick on spline up to intensity cutoff, add measure of freq difference between peakpick and previous peakpicked line
# 4. Input cut target list and match to peakpick
# 5. Measure linewidth 
# 6. Cut line 
# 7. Build up the mask, apply
# 8. Respline the resulting spectrum back down to a lower resolution and save.


def pp(spec, threshold):
	peaks = []
	counters = []
	for i in range(0, len(spec)-1):
		if float(spec[i,1]) > threshold and float(spec[i,1]) > float(spec[(i-1),1]) and float(spec[i,1]) > float(spec[(i+1),1]):
			peaks.append(spec[i])
			counters.append(float(i))


	#print "this shit is this long:"+str(len(peaks))
	picks = zeros((len(peaks),4))
	for i, row in enumerate(peaks):
		picks[i,0] = row[0]
		picks[i,1] = row[1]
		if i == 0:
			picks[i,2] = counters[i]
			picks[i,3] = 0
		else:
			# i,3 is placeholder for linewidth
			picks[i,2] = counters[i]
			picks[i,3] = 0
	return picks

def spline(spec,new_resolution):
    x = spec[:,0]
    y = spec[:,1]

    old_resolution = (x[-1]-x[0]) / len(spec)
    scale_factor = old_resolution / new_resolution

    new_length = int(math.floor(scale_factor*len(spec)))

    tck = splrep(x,y,s=0)
    xnew = arange(x[0],x[-1],new_resolution)
    ynew = splev(xnew,tck,der=0)

    output_spectrum = column_stack((xnew,ynew))

    return output_spectrum

def noise(spec,peakpick, peakpick_threshold):
	thesum = 0.0
	n = 0
	listit = []
	for i in range(0,len(peakpick)-1):
		counter_1 = int(width(spec,peakpick[i],peakpick_threshold)[1])
		counter_2 = int(width(spec,peakpick[i+1],peakpick_threshold)[0])
		if spec[counter_2,0] > spec[counter_1,0]:
			for j in range(counter_1,counter_2):
					thesum += spec[j,1]
					listit.append(spec[j,1])
					n += 1
	avg = thesum/n
	sigma = 0.0
	diff = 0.0
	for val in listit:
		diff += (val-avg)**2
	sigma = math.sqrt(diff/n)
	return math.sqrt(2.0)*math.sqrt(sigma**2+avg**2)




def width(spec_splined,transition,noise_level):
	counter = int(transition[2])
	c_lower = 0
	c_upper = 0
	#print "initial intensity is: " + str(spec_splined[counter,1])
	while spec_splined[counter-c_lower,1] >= noise_level:
		if (spec_splined[counter-c_lower-1,1]>spec_splined[counter-c_lower,1] and spec_splined[counter-c_lower+1,1]>spec_splined[counter-c_lower,1]):
			break
		c_lower += 1

	while spec_splined[counter+c_upper,1] >= noise_level:
		if (spec_splined[counter+c_upper-1,1]>spec_splined[counter+c_upper,1] and spec_splined[counter+c_upper+1,1]>spec_splined[counter+c_upper,1]):
			break
		c_upper += 1
	
	return array([int(counter-c_lower),int(counter+c_upper)])

def search(peakpick,key,lo=0,hi=None):
	if hi is None:
		hi = len(peakpick)
	while lo < hi:
		mid = (lo+hi)/2
		midval = peakpick[mid,0]
		if midval - key <= -0.05: #30 kHz tolerance default
			lo = mid+1
		elif midval - key >= 0.05: #30 kHz tolerance default
			hi = mid
		else:
			return [peakpick[mid,0],peakpick[mid,1],peakpick[mid,2]]
	return [key,-1.0]

def cut_it(spec,peaks,noise_level):
	
	n = 0
	mask = ones(shape(spec)[0])
	for i in range(0,shape(peaks)[0]):
		foo = peaks[i]
		bounds = None
		if foo[1] != -1.0:
			bounds = width(spec,foo,noise_level)
		if foo[1] == -1.0:
			n += 1
		#print "for peak with freq " + str(peaks[i,0]) + " there is a width of "+ str((bounds[1]-bounds[0])*0.002) + "MHz"
		if bounds != None:
			for i in range(int(bounds[0]),int(bounds[1])+1):
				mask[i] = 0
	cut_spec = spec
	for i in range(0,len(spec)):
		cut_spec[i,1] = cut_spec[i,1]*mask[i]
	if n != 0:
		choice = str(raw_input("\n"+"There are "+ str(n)+" peaks that were not found in the search. Would you like to see a list? (yes/no or 1/0): "))
		yes = set(['yes','y','ye','','1'])
		no = set(['no','n','0'])
		choice = choice.lower()
		if choice in yes:
			for i in range(0,shape(peaks)[0]):
				if peaks[i][1] == -1.0:
					print "Peak #"+str(i)+ "  FREQ: "+str(peaks[i][0])
		choice = str(raw_input("\n"+"Would you like to save the cut spectrum anyway? If not, program will restart so you can try again: "))
		if choice in yes:
			return cut_spec
		if choice in no:
			return None
	return cut_spec



def cut(spec,peakpick,cut_list,noise_level):

	# Correlate cut_list with peakpick
	peaks = zeros((len(cut_list),3))

	for i in range(0,len(cut_list)):
		search_peak = cut_list[i]
		peak = search(peakpick,search_peak,0,shape(peakpick)[0])
		if peak[1] == -1.0:
			peaks[i,0] = peak[0]
			peaks[i,1] = -1.0
		else:
			peaks[i,0] = peak[0]
			peaks[i,1] = peak[1]
			peaks[i,2] = peak[2]

	#for i in range(0,len(cut_list)):
	#	print "Input freq: " + str(cut_list[i]) + "    Output freq: " + str(peaks[i])

	return cut_it(spec,peaks,noise_level)



def main(threshold=None):
	spec_file_name = str(sys.argv[-2])
	cut_list_file_name = str(sys.argv[-1])
	
	print "============= SMART CUTTER ==========="
	print "	Author: Nathan Seifert"
	print "	E-mail: nas3xf@virginia.edu"
	print "======================================"+"\n\n"
	print "Loading spectrum file: " + spec_file_name
	t1 = time.time()
	spec = genfromtxt(spec_file_name,dtype=float,delimiter="	")
	t2 = time.time()
	print "Loaded spectrum file. Time it took to load: "+str(round(float(t2-t1),3))+" seconds"+"\n"
	#time.sleep(1.0) 
	original_res = 1000*spec[1,0]-1000*spec[0,0]
	print "Spectrum has "+str(shape(spec)[0])+" data points and a frequency resolution of: "+str(1000*spec[1,0]-1000*spec[0,0])+" kHz"
	print "Intensities are used in the unit of the spectrum file."+"\n"
	time.sleep(1)

	print "\nEnter the following parameter."
	print "NOTE: Make sure it's as close to the baseline as possible to account for all observable transitions!"
	if threshold == None:
		peakpick_threshold = float(input("ENTER GUESS NOISE LEVEL: "))
	else: 
		peakpick_threshold = threshold
	

	spline_res = 0.0020
	print "\n"+"Now splining the data to a resolution of "+str(spline_res*1000)+" kHz"
	t1 = time.time()
	data = spline(spec,spline_res)
	t2 = time.time()
	print "Splining complete. It took: "+str(round(float(t2-t1),3))+" seconds and the resulting spectrum has a resolution of: "+str(1000*data[1,0]-1000*data[0,0])+" kHz"
	
	time.sleep(1.0)
	#print "\n"+ "Writing splined spectrum..."
	#spline_output = open(os.path.splitext(sys.argv[-2])[0]+"_splined.dat","wb")
	#t1 = time.time()
	#for line in data:
	#	spline_output.write(str(line[0])+"	"+str(line[1])+"\n")
	#t2 = time.time()
	#spline_output.close()
	#print "Wrote splined spectrum in "+str(round(float(t2-t1),3))+" seconds"


	print "\n"+"Running peakpick..."
	peakpick = pp(data,peakpick_threshold)
	print "Number of picked transitions: " + str(len(peakpick))

	time.sleep(1.0)
	print "\n"+"Running noise analysis..."
	avg_noise = noise(data,peakpick, peakpick_threshold)
	print "Noise analysis completed. P2P mean noise: "+ str(round(avg_noise*1000,3))+" uV"+"\n"
	time.sleep(0.5)

	cut_list = genfromtxt(cut_list_file_name,dtype=float)
	print "Cut list has " + str(len(cut_list)) +" entries"
	time.sleep(0.5)
	print "Cutting spectrum and resplining to original resolution...."
	t1 = time.time()
	cut_spec = cut(data,peakpick,cut_list,avg_noise)
	if cut_spec != None:
		cut_spec = spline(cut_spec,original_res/1000)
		new_res_test = 1000*cut_spec[1,0]-1000*spec[0,0]
		t2 = time.time()
		print "Success! Cutting and resplining spectrum took " + str(round(t2-t1,3))+" seconds and has a frequency resolution of "+ str(round(new_res_test,3))+ "kHz" 
	elif cut_spec == None:
		test = str(raw_input("\nDo you want to use the same intensity threshold you used previously? "))
		yes = set(['yes','y','ye','','1'])
		no = set(['no','n','0'])
		print "\n============================ "
		print "    RESTARTING PROGRAM....    "
		print " ============================ \n\n\n"
		time.sleep(3.0)
		if test in yes:
			main(threshold=peakpick_threshold)
		if test in no:
			main()


	time.sleep(1.0)

	print "Saving cut spectrum as... "  + os.path.splitext(sys.argv[-2])[0]+"_cut.dat"
	cut_output = open(os.path.splitext(sys.argv[-2])[0]+"_cut.dat","wb")
	for i in range(0,len(cut_spec)):
		cut_output.write(str(cut_spec[i,0])+"	"+str(cut_spec[i,1])+"\r\n")
	cut_output.close()
	print "\nFinished writing spectrum! Thanks for using the program!!"


	
	
		
	




	#
	# Create spectrum array
	#with genfromtxt(spec_file_name,dtype=float) as spec:
	#	for line in spec:
	#		print line
			

if __name__ == "__main__": 
	set_printoptions(suppress=True,precision=5)
	main()