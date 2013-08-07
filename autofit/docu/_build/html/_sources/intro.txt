.. _intro:

Introduction to Autofit
***********************

Welcome to the Autofit documentation! 

Authorship
==========
The original implementation and idea for Autofit was developed in the lab of Brooks Pate at the University of Virginia. However,
the modern Python-based implementation has primarily been developed by Steve Shipman, of New College of Florida, and
Ian Finneran, a former undergraduate student of Steve and now a graduate student in Geoff Blake's group at Caltech.
General maintenance, support, marketing, etc -- as well as tender loving care outside of programming -- has been primarily done by the Pate group at UVa.


Motivation
==========

Scientific
----------
In short, Autofit is an automated triples fitting tool for the purpose of identifying unknown target spectra in a broadband rotational spectrum.
Since the implementation of high speed digital electronics has become routine in broadband microwave spectroscopy, new tools were required in
order to decouple spectra of low abundance molecules out of broadband scans that are increasingly becoming more sensitive and spectrally dense.
Some of the most dense spectra taken in the Pate lab have dynamic ranges of over 20000:1, and with line densities of over 1 MHz\ :sup:`-1`\ . These kinds of 
line densities make it unfeasable to manually fit spectra using traditional visual cues, especially for species such as isotopologues which are often
components of the dense "weeds" seen in deep-averaged rotational spectra. 

Computational
-------------
Originally, Autofit was developed as a script for PTC's Mathcad. Although Mathcad is easy to use, it is sadly proprietary, not used by the majority 
of physical chemists, and is rather unstable (and extremely slow!) for computation with large data sets. Python is a natural choice for applications such as autofit,
since it is well supported by the scientific community, it is open source (and the majority of useful scientific and mathematical libraries are as well), it is lightweight
and completely portable, and the language itself has a remarkably small learning curve for scientists with respect to other languages such as FORTRAN and C. Additionally,
it is relatively easy to scale computations across multiple CPUs or CPU cores using Python, which is a necessity for making Autofit an efficient tool. 

What Autofit Can Do for You!
============================
As stated in the motivation, Autofit is an automated triples fitter for broadband rotational spectra. As a program, it merely serves as a generator and processor of input and output files
for Herb Pickett's legendary SPCAT/SPFIT program package, which can be found  `on the JPL website. <http://spec.jpl.nasa.gov/>`_
Since Autofit is built exclusively in Python, the only functional requirements are compatible binaries for SPCAT and SPFIT and a Python interpreter for whatever platform you want to run Autofit on. 
Here is, in a jiffy, what Autofit does:

* Given a set of ab-initio or guess rotational constants and an experimental spectrum, Autofit will fit a A/B/C rotational constant triplet to every possible set of 3 lines possible within an error window.

	* For example, if you want to fit the 2\ :sub:`02`\ - 1\ :sub:`01`\ transition, which is predicted to be at 10 GHz, with an error window of 300 MHz, Autofit will check every experimental line between 9.7 and 10.3 GHz. 
	* In this example (and in general), Autofit will do this for the three transitions you choose, and then output a sorted list of A/B/C triples ranked by goodness of fit (which is calculated by fitting an additional set of transitions to build up an effective fit)

* For a given A/B/C triple output by Autofit, Autofit can take that triple and refine the fit with distortion and/or additional transitions in the fit.
* For a given set of experimental rotational constants for a parent species, and an ab-initio structure for the parent species, Autofit can search for isotopologues separately and automatically by scaling predicted constants via the ab initio geometry.
* In principle, Autofit can scale to an arbitrary number of CPUs or cores. It does this by dividing the total list of triples to check into N segments, where N is the number of CPUs/cores being used in the calculation. 

	* As it stands, Autofit can fit around 35-50 triples a second per core on a typical last-generation Intel CPU (Ivy Bridge). So for a 4-core (8 core via Hyperthreading) system, approximately 250 triples a second can be fit. A typical Autofit run of 10 million triples therefore takes around 10-12 hours.

* Although this is currently in development, Autofit will soon be integrated into a GUI where choosing transitions to check for Autofit and analyzing Autofit results will be possible, therefore removing the requirement of using an additional program (such as JB95 or PGOPHER) for checking Autofit results by hand.

What Autofit Can't Do for You!
==============================
Rotational spectroscopists, wipe the sweat off your brow -- Autofit is not here to make your job obsolete. Rather, Autofit has been designed to speed up the purely mechanical (and often mind-numbingly painful) process of fitting spectra. 
For instance:

* Autofit does not fit nuclear quadrupole hyperfine. If you run Autofit on a molecule, such as one containing nitrogen, you will often find multiple fits that are roughly close to each other in RMS error. This is due to the fact that Autofit finds the same rigid rotor spectrum, but fits to different resolved hyperfine components.
* Internal rotation is not considered at all. For fitting a rigid rotor spectrum with Autofit, you will usually get at least two fits: one with the A state fit, which will fit to the rigid rotor model closely, and a higher energy E state fit. More than likely you will also get fits where A and E transitions are mixed, albeit at a higher RMS error.
* Autofit does NOT fit distortion! You can enter in distortion at the beginning for SPCAT to consider, but it will be treated as a static variable in the predictions. There is functionality to take an Autofit result and then go and fit distortion by adding additional lines to the fit, but this is effectively an automated tool for doing "final" fits in SPFIT. None of it is done automatically during the autofitting process.

And most importantly, though this point is becoming less important as we continue to refine the program to make it easier to use:

* Autofit does not making fitting rigid rotor spectra a mindless activity. Some knowledge, or at least common sense, about what rotational transitions are good choices as triples fitting parameters is required. Mindlessly choosing transitions from the predicted list will almost never give you good results. In the Pate group, we have a saying that has become true over and over again: **Autofit never fails you; only you fail Autofit!**

However, fear not -- there are a number of helpful tips and tricks available in this documentation that will ease your Autofit journey. For a trained rotational spectroscopist, use your best judgement -- what combinations of lines do you typically use to get a good rough fit on A, B and C?

Licensing
=========
As it stands currently (as of August 2013), Autofit is licensed with the GPLv2 license. This is temporary, as it is likely we will move to a more liberal license 
in the future. You can find the actual terms of the GPLv2 license `here. <http://www.gnu.org/licenses/gpl-2.0.html>`_ 

Support
=======
Currently, maintenance and support is managed by Nathan Seifert of the Pate group. He is happy to take any questions you might have on the program, both in terms
of scientific and computational support. You can contact him at **nas3xf[at]virginia[dot]edu**. 