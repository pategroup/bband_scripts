import spec_tools
import numpy as np
import matplotlib.pyplot as plt



#options =  {'Chirp_Start': 2000.0, 'Chirp_Stop': 2100.0, 'Chirp_Duration': 1.0, 'Sample_Rate': 24}
#markers = {'DELAY': 1.5, 'CH1_ON': 0.1, 'CH1_OFF': 1.1, 'CH2_ON': 2.8, 'CH2_OFF': 3.8, 'BUFFER': 52.0-3.8}

options= {'Chirp_Start': 2000.0, 'Chirp_Stop': 2100.0, 'Chirp_Duration': 0.25, 'Sample_Rate': 24}
timings = {'DELAY': 1.5, 'M1_WIDTH': 1.0, 'M2_WIDTH':  1.0, 'M_PULSE_BUFFER': 0.5,'BUFFER': 52.0}

frames = 8

durations = [0.2,0.3]

# for i in range(0,len(durations)):
	# options.update({'Chirp_Duration':durations[i]}) 
	# chirp = spec_tools.arb_pulse(options,timings,frames)

	# duration_ns = durations[i]*1000

	# filename = str(round(float(options['Chirp_Start']/1000),3))+'_'+str(round(float(options['Chirp_Stop']/1000),3))+'_GHz_chirp_'+str(int(duration_ns))+'ns_length_string_of_'+str(frames)+'.txt'
	# output = open(filename,'w')
	# for i in range(len(chirp)):
		# output.write(str(chirp[i,0])+'\t'+str(chirp[i,1])+'\t'+str(chirp[i,2])+'\n')
	# output.close()

pulse = spec_tools.arb_pulse(options,timings, 2)

x = range(0,len(pulse))
#print x[0]
y1 = pulse[:,0]
y2 = pulse[:,1]
y3 = pulse[:,2]

plt.plot(x,y1,x,y2,'b',x,y3,'r')
plt.show()

#print np.shape(x)
#print np.shape(y)



