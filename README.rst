TelFit Code
===========

This code will fit the telluric absorption spectrum in data, using
LBLRTM. More details can be found in the `documentation`_ and examples, but the
installation and run procedure are outlined below. If you use this code,
please cite `my paper`_.

Installation
------------
**New for 2021:** Python 3 support!

As of version `v1.4.0` TelFit now only runs on Python 3.  Python 2 support has been dropped, though it is technically still possible to install old Python 2 versions pre-`v1.4.0`.  The code has been tested on Python 3.7, and may work with other minor versions. Code installation requires the following packages:

  - matplotlib
  - **numpy v1.6 or greater**
  - scipy v0.13 or greater
  - astropy v0.2 or greater
  - lockfile
  - pysynphot v0.7 or greater
  - fortranformat
  - **cython**
  - **requests**

The bolded entries are required *before* installation, so make sure you get them from pip, apt-get/yum, or conda (depending on your OS and python distribution). The setup script will attempt to install the rest if you don't have them, but I suggest doing it yourself just to make sure nothing goes wrong. Once you have the dependencies, simply type

.. code:: bash

    pip install TelFit

to install TelFit. It may take a while, as it needs to build the LBLRTM code and some of its standard input files.

Running TelFit
--------------

To run TelFit, you should create a script like in the examples. The key
parts of the script are the inputs to the TelluricFitter class. You
should:

-  Initialize fitter: fitter = TelluricFitter()
-  Define variables to fit: must provide a dictionary where the key is
   the name of the variable, and the value is the initial guess value
   for that variable. Example: fitter.FitVariable({“ch4”: 1.6, “h2o”:
   45.0})
-  Edit values of constant parameters: similar to FitVariable, but the
   variables given here will not be fit. Useful for settings things like
   the telescope pointing angle, temperature, and pressure, which will
   be very well-known. Example: fitter.AdjustValue({“angle”: 50.6})
-  Set bounds on fitted variables (fitter.SetBounds): Give a dictionary
   where the key is the name of the variable, and the value is a list of
   size 2 of the form [lower\_bound, upper\_bound]
-  Import data (fitter.ImportData): Copy data as a class variable. Must
   be given as a DataStructures.xypoint instance
-  Perform the fit: (fitter.Fit): Returns a DataStructures.xypoint
   instance of the model. The x-values in the returned array are the
   same as the data.

.. _my paper: http://adsabs.harvard.edu/abs/2014AJ....148...53G
.. _documentation:  http://telfit.readthedocs.org/en/latest/
