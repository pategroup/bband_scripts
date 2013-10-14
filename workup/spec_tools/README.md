spec_tools README
=================

Introduction
------------
This is a work-in-progess module for Python 2.7.x 
that provides standard functionality for data work-up
and analysis pertinent to broadband Fourier transform
microwave spectroscopy.

Currently, the module is standalone and still a work in
progress. The features listed below are NOT necessarily
the features available in the current version in this 
Github repository.

Authorship
----------
This module is being developed by Nathan Seifert,
of the Pate lab at University of Virginia. Nathan
provides the primary support and maintenance for this
module. However, some of the original code for some
of the functions available in spec_tools were originally
developed in Python by other members of the Pate lab,
such as Cristobal Perez (UVa) and Daniel Zaleski (currently
at Newcastle University).

License
-------
This module is licensed under the GPLv2 license. For additional
inforamtion, visit http://www.gnu.org/licenses/gpl-2.0.html

Features
--------
Planned features for this module include:
- Generation of chirped pulses
- A Fourier transform routine for time-domain data
- Rotation of molecular geometries and dipoles into the principal axis
- Spectrum peakpicking, as well as spectrum cutting from a given linelist
- Scaling and prediction of isotopologue rotational constants
- Quadrupole constant converter between (1.5aa/.25bb-cc <--> aa/bb/cc)
- Creation of gain correction curves for normalizing noise floors in FT spectra
- Generator of input files for SPCAT, as well as predicted spectrum generation
- Functionality to generate spectrogram for the purpose of analyzing pulse purity data
- And more!



