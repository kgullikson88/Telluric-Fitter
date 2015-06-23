Generating a Telluric Model with TelFit
=======================================

Making a telluric model used to be hard. With TelFit, it is just a few lines of code. 

.. code:: python

    from telfit import Modeler

    # Set the start and end wavelength, in nm
    wavestart = 500.0
    waveend = 900.0

    # Make the model
    modeler = Modeler()
    model = modeler.MakeModel(humidity=50.0, 
                              lowfreq=1e7/waveend, 
                              highfreq=1e7/wavestart)
