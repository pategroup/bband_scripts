import time
import sys

t1 = time.time()

def update_progress(progress):
	sys.stdout.write('\r[{0}] {1}%'.format('#'*(progress/10),progress))  # using print() prints new lines
	sys.stdout.flush()

def run(numtrip):
	for i in range(int(numtrip)):
		time.sleep(0.1)
		update_progress(i)

input_num = raw_input("Input a number of triples:")
run(input_num)

