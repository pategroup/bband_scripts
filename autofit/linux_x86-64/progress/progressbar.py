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

	sys.stdout.write('\r'+str(progress)+'/'+str(input_num)+' :: '+'[{0}] {1}%'.format('#'*(progress/(int(input_num)/10)),progress/float(input_num)*100)+' :: '+str(int(veloc))+' Hz') # using print() prints new lines
	sys.stdout.flush()

def run(numtrip):
	for i in range(int(numtrip)):
		time.sleep(0.01)
		update_progress(i)


input_num = raw_input("Input a number of triples:")
run(input_num)

