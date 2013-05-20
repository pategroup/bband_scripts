====================README FOR WORKUP SCRIPTS============

This folder contains scripts useful for working up broadband 
data as well as basic analysis of spectra. 

For the time being, there are six scripts available in this folder,
with their descriptions listed below.

All scripts are licensed under the GPLv2, as stated in the 
primary repository README file.

==========================================================

1) fft-zero-padding-with heading.py
This script takes a one-column FID file (a list of amplitudes, 
time values are implicit from the sample rate chosen), FTs the FID
and plots the resultant FT.
Currently, sample rate and frequency boundries must be edited manually.
You can edit the sample rate by changing variable srate (line 25) to
the requested sample rate. You can change the frequency boundries
for the FT in line 54.

2) omc_reader.py
Put this script in an autofit result folder (with the sorted_omc_cat.txt) and run it. It will take the top 1000 fits
from the sorted file and create a new file called 1000fits.txt

3) isotopomers-with-labels.py
Given an input geometry and experimental rotational constants for the normal species, this program will output scaled rotational constants for a number of common single isotopologues. In order to do
double isotopologues or calculate the constants for an uncommon isotopologue, you will have to manually edit the program yourself (the flow is fairly straightforward to figure out). 
KNOWN ISSUE: For a time being, you cannot define an atom with an element number larger than 9. Inputting a set of coordinates with atomic numbers > 9 will cause the program to crash. Currently the program has silicons, for instance, as element 3.

4) peakpicker.py
This script takes in an input FT and a intensity threshold and output
a two-column list of frequencies and intensities for every line above the input threshold.

5) linelist-cutter-dialogues.py
This script takes in an input FT and a one-column frequency linelist, and cuts all lines in the linelist out of the spectrum. In order to change the frequency window for each cut, change variable width on line 28. Default is 600 kHz centered at the frequency in the linelist.

6) SPCAT_writer-v1.py
This script will take in an input of constants for various types of species (symmetric top, asymmetric top, asymmetric top with quadrupole, etc) and generate .int/.var files suitable for use with SPCAT. If SPCAT is contained in the folder with this script, you can also have it predict, write and plot predicted spectra for the generated .int/.var fines. Currently this functionality only works in Windows (it calls SPCAT.exe), though it is easy to fix this for an SPCAT binary in linux (replace the SPCAT.exe calls with ./spcat).


