Installation
=============

This code requires the following packages:

  - matplotlib
  - **numpy v1.6 or greater**
  - scipy v0.13 or greater
  - astropy v0.2 or greater
  - lockfile
  - pysynphot v0.7 or greater
  - fortranformat
  - **cython**
  - **requests**

The bolded entries are required *before* installation, so make sure you get them from pip, apt-get/yum, or conda (depending on your OS and linux distribution). The setup script will attempt to install the rest if you don't have them, but I suggest doing it yourself just to make sure nothing goes wrong. Once you have the dependencies, simply type

.. code:: bash

    pip install TelFit

to install TelFit. It may take a while, as it needs to build the LBLRTM code and some of its standard input files.