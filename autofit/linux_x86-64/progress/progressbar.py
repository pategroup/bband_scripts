import time
import sys

t1 = time.time()
t2 = time.time()
veloc = 0
count = 0

# This works nicely!!!
def update_progress(progress):
	global count
	if count == 0:
		global t1
		t1 = time.time() # Gets time at beginning of speed calc
	count += 1
	if count == 10:
		global veloc
		global t2
		t2 = time.time()
		veloc = 10.0/(t2-t1)
		count = 0
	print type(progress)	
	print type(procs)
	print type(((progress*procs)/(int(input_num)/10)))
	print type((progress*procs)/float(input_num)*100)
	sys.stdout.write('\r'+str(progress*procs)+'/'+str(input_num)+' :: '+'[{0}] {1}%'.format('#'*((progress*procs)/(int(input_num)/10)),str((progress*procs)/float(input_num)*100))+' :: '+str(int(veloc*procs))+' Hz') # using print() prints new lines
	sys.stdout.flush()

def run(numtrip):
	for i in range(int(numtrip)):
		time.sleep(0.01)
		update_progress(i)


input_num = raw_input("Input a number of triples:")
procs = int(raw_input("Input number of procs"))
run(input_num)
