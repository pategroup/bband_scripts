Autofit
=======

Authorship
----------
Primary authorship of this program is associated with Steve Shipman, of New College of Florida
and Ian Finneran, now at Caltech. The original idea for this program originated with 
the Pate group at UVa, but the implentation stored here in Python was primarily developed by
Steve and Ian. 

Purpose
-------
This program is a brute-force triples fitter, primary for use with assigning unknown
speices in broadband rotational spectra. The program flow is summarized as follows:<br /><br />
		1) Autofit takes in a 2-column linelist, the first column being frequencies
and the second column being intensities. <br />
        2) An input of "guess" rotational constants and parameters (such as dipole moments)
are made, and this input array is passed to SPCAT to predict transition intensities and frequencies.<br />
        3) From the prediction output, three transitions (the fixed transition triple) are chosen to fit,
ideally three transitions with sufficient linear independence (this is checked in-program after the triple is selected)<br />
        4) An additional set of N user-chosen lines are also input, which in turn leads to each "triples" fit being an N+3 line fit;
the use of these additional N lines are for the purpose of ranking the triples fits<br />
        5) Autofit then checks every set of possible transitions within the input frequency error window with the N+3 rotational 
transition parameters, and outputs them as a sorted text file.<br />

Requirements
------------
This program is coded exclusively in Python 2.7. Required third-party libraries include Easygui and Numpy, 
which are easily found online or in your favorite Linux repository. Current work with a new autofit gui
requires PyQT 4.

Contents of this Repository
---------------------------
In the /autofit root folder, you can find a Linux-compatible autofit script and SPCAT/SPFIT binaries.
The SPX binaries were compiled using 32-bit gcc, so make sure you have valid Linux 32-bit libraries
(for instance, ia32libs in Ubuntu) if you're running a 64-bit distribution. Otherwise the SP(X) binaries
will not run. The version in the root, _v9, is old, and does not contain many of the new features
available in the most recent versions of autofit.

It is fairly trivial to port a new version of Autofit to Linux (the subprocess calls need to be changed
in order to support Linux executables), but this problem has been on the backburner so at the moment
we will leave it up to the interested user to do it themselves (feel free to contact those listed under Support
if you have any questions)

The most recent versions, tested & working in Windows, can be found in the /windows subdirectory. There are 
a number of files in this folder, but they can be categorized as follows:

- prog_A_vX.py : This is the main, standalone Autofit program. <br />
- fitting_GUI_Vx.py/fitting_GUI_B_Vx.py : This is an experimental graphical spectrum plotting/prediction program that is currently
in development. It will be integrated with an autofitting routine. At the moment it is unstable but works as intended. It requires
SPCAT and SPFIT to be in the directory with it, as well as the scripts listed below. At the moment, it only works in Windows,
but it should be possible to port it over to Linux easily since the GUI is developed exclusively in PyQT. It does not require prog_A_vX.py to run.<br />

- refit_module.py / autofit_NS_module.py / isotopologue_module_b.py : These are stripped down autofitting modules for use in integrating into another
program, such as the fitting GUI mentioned above. They do not work alone as intended in prog_A_vX.py; rather they offer another program a library of
functions in order to implement an autofitting routine. The function of these libraries are as follows:<br />
	1) autofit_NS_module.py : provides basic autofitting routines for the purpose of finding singular species, e.g. a parent species in a set of isotopologues.<br />
	2) isotopologue_module_b.py : provides routines for searching for isotopologues of a parent species, given the experimental constants of the parent species and its ab initio geometry 
	(for the purposes of scaling predicted isotopologue rotational constants)<br />
	3) refit_module.py: provides routines for taking an autofit A/B/C triple and fitting it with SPFIT (with distortion and additional lines, for instance).<br />
- isotopologue_module_a.py: This is a stand-alone routine for fitting isotopologues given a parent species experimental constants/geometry, as with module_b. This is instead standalone;
however, this has been integrated into the latest version of prog_A_vX.py, so there is really no reason to use this file anymore for practical use (just use prog_A_vX.py)

Benchmarking
------------
Currently, the program is highly CPU-dependent, and results vary dramatically on the desktop you chose to use.
With a Core i7-3770S CPU @ 3.10 GHz, using 8 cores (4x2 w/ multithreading), an effective triples fitting rate
of ~250 Hz is achieved, so a 10,000,000 triples fit run (a standard scan size) will take around 11 hours.

Support
-------
Please direct all questions regarding this program to the current administrator of the git repository,
Nathan Seifert. His e-mail is nas3xf[at]virginia.edu.


 

 
