Simple Model Fitting
=====================

Alright making a model is pretty simple, but you probably want to go further and find the *best* model for your data. This is a more complicated thing, so requires a bit more code. The rest of this tutorial will assume that you have imported a few modules:

.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from telfit import TelluricFitter, DataStructures

The first thing we need to do is make an instance of the TelluricFitter class, and set some basic information about your observatory.

.. code:: python

    fitter = TelluricFitter()
    observatory = {"latitude": 30.5,
                   "altitude": 2.5}
    fitter.SetObservatory(observatory)

You just need to give the TelluricFitter class a latitude and altitude. For your convenience, there are a few common observatories already loaded into TelFit, so you can also set the observatory with

.. code:: python

    fitter.SetObservatory(obsname)

where "obsname" is one of the following:

  - McDonald
  - CTIO
  - Mauna Kea
  - La Silla
  - Paranal

You can easily add your own observatory if you wish by editing the source code for the SetObservatory method. Alright, now that we have initialized the fitter class we need to tell it how what variables to fit and give their initial values. In this example, we will just fit the relative humidity and oxygen mixing ratio. You should be able to get an initial value for the humidity from your fits header or archived weather data; if not, guess with something like 50%. We will also set some constant values that are important.

.. code:: python

    fitter.FitVariable({'h2o': humidity_guess,
                        'o2': 2.12e5})   # This is a good initial guess for the O2 abundance

    fitter.AdjustValue({"angle": angle,
                        "pressure": pressure,
                        "resolution": resolution,
                        "wavestart": data.x[0] - 20.0,
                        "waveend": data.x[-1] + 20.0})

You will notice that the input for both of the methods is a dictionary with parameter name and value. For the FitVariable method, the value is an initial guess. For the AdjustValue method, the value you give it is fixed (with the exception of detector resolution). Now let's set some bounds on the fitted parameters as well. You can set bounds on any parameter, but it will be ignored if it is a constant.

.. code:: python
    
    fitter.SetBounds({'h2o': [1.0, 99.0],
                      'o2': [5e4, 1e6],
                      'resolution': 50000., 70000.})

You will notice that we set bounds on the 'resolution' parameter as well. The detector resolution is always fit, and you should always give it bounds that are pretty close to the true value.

Finally, we can perform the fit. This will take quite some time to run, so go get some coffee if you are following along...

.. code:: python

    model = fitter.Fit(data=data)

There are lots of options to Fit to control how it does everything. Check the API if the defaults aren't working for you. The data object you pass to Fit must be an object of type DataStructures.xypoint, which you can create with the following code (assuming you have your data in numpy arrays called 'wave', 'flux', and 'continuum'). In old versions of telfit, wave had to be in nanometers; now, it will accept any astropy units compatible with nanometers.

.. code:: python
    
    data = DataStructures.xypoint(x=wave, y=flux, cont=continuum)

The continuum does not need to be great, since TelFit estimates that as well after dividing out the telluric model.