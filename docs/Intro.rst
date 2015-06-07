Introduction to Telluric Modeling with TelFit
=============================================

Do you have spectra observed from Earth? Do you hate needing to spend precious telescope time just to observe a telluric standard star? If so, TelFit might work for you! TelFit is a Python code written specifically to model and fit the telluric absorption spectrum in astronomical data. It is essentially a wrapper around the `LBLRTM`_ FORTRAN code, and provides an easy to use, object-oriented interface to the code. TelFit vastly simplifies the process of gnerating a telluric model, and can efficiently remove the telluric absorption lines in high-resolution optical or near-infrared spectra. 





.. _LBLRTM: http://rtweb.aer.com/lblrtm.html