Getting Started with Autofit
****************************

Okay, so hopefully the introduction hasn't turned you off from the idea of fitting spectra automatically (though if you're reading this I'm sure you're still quite interested!). One of the motivations of writing this documentation
is not only to teach a new user how to use Autofit, but also to provide a detailed repository of information and verification for those who are still skeptical of the method. Rotational spectroscopy is uniquely beautiful amongst
the spectroscopic techniques, and the amount of rich chemical information we can gather from a simple microwave spectrum is often incredible. However, we as a scientific community have gotten to the point where the absolute amount of information we can gather from a single spectrum is beyond the capabilities of efficient analysis using manual fitting techniques. 

This is where Autofit comes in. We have been using this program in the Pate lab for a couple of years now, and we have attributed much of our speed, efficiency and success to intelligent use of Autofit. This would be an apt time to say that "the proof is in the pudding," 
but it is more instructive for us to educate and have Autofit prove to you its incredible potential as a essential tool for broadband rotational spectroscopy, hopefully alongside Pickett's CALPGM suite and the vast majority of tools on Zbiginew Kisiel's `PROSPE. <http://info.ifpan.edu.pl/~kisiel/prospe.htm>`_ 

Prerequisites for Autofit
=========================

* Python 2.7.5 -- Currently, as standard with most scientific applications in Python, Python 3.3 is NOT supported, though in reality this is due to the required additional libraries listed below not necessarily supporting 3.3.

	* Python 2.7.5 for Windows can be downloaded `here. <http://www.python.org/getit/>`_
	
* Numpy / Scipy 
	
	* Numpy is available `here. <https://pypi.python.org/pypi/numpy>`_
	* Scipy is also available `here. <http://sourceforge.net/projects/scipy/files/scipy/0.12.0/scipy-0.12.0-win32-superpack-python2.7.exe>`_

* Matplotlib, available `here. <http://matplotlib.org/>`_
* Easygui, available `here. <http://easygui.sourceforge.net/>`_

**For the purposes of creating/moving files, Autofit uses *bash shell* commands. Therefore, Autofit will *not* work in the standard Windows command prompt. You need to use a bash shell emulator.** We recommend either of the two listed below:

* Cygwin, available `here. <http://www.cygwin.com/>`_ NOTE: There are a lot of options in the Cygwin install, but Autofit will run out of the box in a default Cygwin installation.
* The mingw32 distribution included in Git for Windows, available `here. <http://git-scm.com/downloads>`_ This is a nice package in case you want to keep up to date with our Github repository for the latest updates to Autofit. 

If you want to use the experimental GUI programs currently being developed, you will also need PyQT4, binaries of which can be obtained `here. <http://www.riverbankcomputing.com/software/pyqt/download>`_

* SPCAT / SPFIT. There are 64-bit binaries available on the `JPL repository <http://spec.jpl.nasa.gov/ftp/pub/calpgm/>`_ but they pop up these annoying windows that will spam your computer when you run Autofit (in fact they'll pop up each time you fit a triple, so 100-250 times a second!). **Instead, we recommend using the 32-bit binaries available on our Github repository, listed in the next section.**

Obtaining Autofit
=================
You can find the latest version of Autofit and the experimental GUI at our `Github repository. <https://github.com/pategroup/bband_scripts/tree/master/autofit>`_ 

If you're comfortable with using Git, you can clone our repository by entering in the command::
	
	git clone https://github.com/pategroup/bband_scripts.git
	
Otherwise, here are direct links to the most useful programs in the repository for the average user:

* Autofit, stand-alone: `Current version: v15 <https://raw.github.com/pategroup/bband_scripts/master/autofit/windows/prog_A_v15.py>`_
* SPCAT.exe available `here. <https://github.com/pategroup/bband_scripts/raw/master/autofit/windows/SPCAT.EXE>`_
* SPFIT.exe available `here. <https://github.com/pategroup/bband_scripts/raw/master/autofit/windows/SPFIT.EXE>`_ 

There are many other useful scripts available in the repository, including some basic workup scripts for broadband spectra and the new experimental fitting GUI/Autofit analysis program. Feel free to browse the repository and see if anything catches your eye (the README files in the root directory and in the autofit directory give details about what each script in the repository is).

In addition, there are SPCAT/SPFIT binaries that are tested (and working) in x86 (32-bit only!) Linux, as well as a tested but older version of autofit (v9) for Linux available in the `/autofit root of the Github repository. <https://github.com/pategroup/bband_scripts/tree/master/autofit>`_ If you are running 64-bit Linux, these binaries will NOT work unless you have installed the proper libraries (e.g. ia32libs in Ubuntu) for running 32-bit binaries in 64-bit Linux!!!

In terms of best practices:

* Include SPCAT.exe / SPFIT.exe in a new folder with the Autofit program. Autofit will create new directories for each job you start with it, so make sure you have permissions to create directories.

* When you start a new job in Autofit, Autofit requires that the job name (and hence the directory name) is unique. So if you already have a directory in your autofit folder called "/foobar" from some old run and you want to start a new job called foobar, make sure to delete the old /foobar folder.

* Autofit is completely CPU dependent and not all that reliant on RAM availability, so if you have a computer that likes to overheat when you work the CPU hard, either avoid running Autofit on all available cores or find a better cooling solution.

Notes on "Out of the Box" Performance with Autofit
==================================================
* Running Autofit in an environment such as IPython or Spyder is probably doomed for failure. Always run in a simple bash shell by using the command::
	
	python prog_A_vX.py

* Since ALL of the testing done on Autofit has been on either Windows 7 or Ubuntu/Debian, functionality for Autofit on OS X is unknown. There have been some preliminary suggestions that it does NOT work out of the box on OS X. This could be due to the fact that the OS X terminal is not bash, but rather tsch (this might have changed, but it would be prudent to make sure). Changing the terminal to bash could alleviate these issues.
Also there are no official binaries of Easygui for OS X, though it might be possible to use the Linux tar in lieu of this. Note: this is pure speculation!

* We find that the best out of the box performance is gained by installing a very simple Python environment and avoiding use of distributions such as Enthought. It could very well work just fine on Enthought or Anaconda, but installing the packages listed above to create a very minimal Python environment tends to work 100% of the time. Again, feel free to e-mail us if you have issues with this.
