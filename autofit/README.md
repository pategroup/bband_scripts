==== Autofit README ====

1) Authorship
Primary authorship of this program is associated with Steve Shipman, of New College of Florida
and Ian Finneran, now at Caltech. The original idea for this program originated with 
the Pate group at UVa, but the implentation stored here in Python was primarily developed by
Steve and Ian. 

2) Purpose
This program is a brute-force triples fitter, primary for use with assigning unknown
speices in broadband rotational spectra. The program flow is summarized as follows:

	a) Autofit takes in a 2-column linelist, the first column being frequencies
and the second column being intensities. 
        b) An input of "guess" rotational constants and parameters (such as dipole moments)
are made, and this input array is passed to SPCAT to predict transition intensities and frequencies.
        c) From the prediction output, three transitions (the fixed transition triple) are chosen to fit,
ideally three transitions with sufficient linear independence (this is checked in-program after the triple is selected)
        d) An additional set of N user-chosen lines are also input, which in turn leads to each "triples" fit being an N+3 line fit;
the use of these additional N lines are for the purpose of ranking the triples fits
        e) Autofit then checks every set of possible transitions within the input frequency error window with the N+3 rotational 
transition parameters, and outputs them as a sorted text file.

3) Requirements
This program is coded exclusively in Python 2.7. Required third-party libraries include Easygui and Numpy, 
which are easily found online or in your favorite Linux repository.

4) Notes about this git repository
Currently, the only version available on this repository is the Linux version, which has been tested
and is working in Ubuntu 12.10 LTS, 64-bit. The SPCAT/SPFIT binaries included in this repository are
32-bit binaries, so this should work fine in 32-bit Linux, as well.

Porting to Windows is easy; any call for the SPCAT/SPFIT binaries can be replaced with Windows 
executable calls (replace "./SPCAT default" with "SPCAT.exe default"). Use of this script in
Windows will require the use of a .bash shell, such as Cygwin, as the default directory handling is done using bash commands.
It will NOT work using cmd.exe.

5) Benchmarking
Currently, the program is highly CPU-dependent, and results vary dramatically on the desktop you chose to use.
With a Core i7-3770S CPU @ 3.10 GHz, using 8 cores (4x2 w/ multithreading), an effective triples fitting rate
of ~250 Hz is achieved, so a 10,000,000 triples fit run (a standard scan size) will take around 11 hours.

6) Support
Please direct all questions regarding this program to the current administrator of the git repository,
Nathan Seifert. His e-mail is nas3xf[at]virginia.edu.


 

 
