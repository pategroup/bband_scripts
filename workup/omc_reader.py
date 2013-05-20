from easygui import *
from numpy import *
f=open('sorted_omc_cat.txt','r')
fits = []
for i in range (1000):
    r=f.readline()
    fits.append(r)
    #print fits
f.close()
data = ""
for i in (fits):
    data+= i +''

f2=open("1000fits.txt","w")
f2.write(data)
f2.close()



