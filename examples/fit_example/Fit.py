"""

    This file is part of the TelFit program.

    TelFit is free software: you can redistribute it and/or modify
    it under the terms of the MIT license.

    TelFit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the MIT license
    along with TelFit.  If not, see <http://opensource.org/licenses/MIT>.


"""


import numpy
import matplotlib.pyplot as plt
import TelluricFitter
import DataStructures
from astropy import units, constants


if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter(debug=False)

  #Set the observatory location with a keyword
  fitter.SetObservatory("McDonald")

  #Can also give a python dictionary like below
  #observatory = {"latitude": 30.5,
  #               "altitude": 2.5}
  #fitter.SetObservatory(observatory)

  #Read in the spectrum and store as an xypoint
  x,y,c,e = numpy.loadtxt("Spectrum.dat", unpack=True)
  data = DataStructures.xypoint(x=x, y=y, cont=c, err=e)
  

  #The following are parameters that should be in the fits header for an observation
  angle = 37.4    #Zenith distance
  pressure = 796.8    #Pressure, in hPa
  humidity = 54.0    #Percent humidity, at the observatory altitude
  temperature = 292.5   #Temperature in Kelvin
  resolution = 60000.0    #Resolution lambda/delta-lambda


  #Define variables to be fit, and give initial guesses
  #Let temperature float, but we will tightly constrain it with bounds (it should be known well)
  fitter.FitVariable({"h2o": humidity, 
                      "o2": 2.12e5,
                      "temperature": temperature})

  #Adjust parameters that will not be fit, but are important
  fitter.AdjustValue({"angle": angle,
                      "pressure": pressure,
                      "resolution": resolution,
                      "wavestart": data.x[0] - 20.0,
                      "waveend": data.x[-1] + 20.0})

  #Set bounds on the variables being fit (resolution is always fit inside each loop)
  fitter.SetBounds({"h2o": [1.0, 99.0],
                    "o2": [5e4, 1e6],
                    "temperature": [temperature-5.0, temperature+5.0],
                    "resolution": [resolution/2.0, resolution*2.0]})

  #Import the data to be fit to the fitter
  fitter.ImportData(data)

  #Perform the fit. See the documentation for all the options for fit
  model = fitter.Fit(resolution_fit_mode="gauss", adjust_wave="model")

  #Get the improved continuum from the fitter
  data.cont = fitter.data.cont
  
  #Plot the data along with the best fit
  plt.figure(1)
  plt.plot(data.x, data.y/data.cont, 'k-', label="Data")
  plt.plot(model.x, model.y, 'r-', label="Telluric Model")
  plt.legend(loc='best')
  plt.title("Comparison of the Data to the Telluric Model")

  #Also make a plot of the corrected data
  plt.figure(2)
  plt.plot(data.x, data.y / model.y, 'k-')
  plt.title("Telluric-Corrected Data")
  plt.show()
  
  #Save the best-fit model spectrum
  numpy.savetxt("Model.dat", numpy.transpose((model.x, model.y)), fmt="%.8g")


