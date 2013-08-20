import time
import sys
import curses 

def progress():
	curses.initscr()
	win = curses.newwin(3,32,14,10)
	win.border(0)
	win.addstr(1,1,"Progress ")

	pos = 10
	for i in range(15):
		win.addstr(1,pos,".")
		win.refresh()
		time.sleep(0.1)
		pos += 1
	win.addstr(1,26,"Done!")
	win.refresh()
	time.sleep(0.5)

progress()