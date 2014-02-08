"""
This module provides the 'TelluricFitter' class, used
to fit the telluric lines in data.

Usage:
  - Initialize fitter: fitter = TelluricFitter()
  - Define variables to fit: must provide a dictionary where
      the key is the name of the variable, and the value is
      the initial guess value for that variable.
      Example: fitter.FitVariable({"ch4": 1.6, "h2o": 45.0})
  - Edit values of constant parameters: similar to FitVariable,
      but the variables given here will not be fit. Useful for 
      settings things like the telescope pointing angle, temperature,
      and pressure, which will be very well-known.
      Example: fitter.AdjustValue({"angle": 50.6})
  - Set bounds on fitted variables (fitter.SetBounds): Give a dictionary
      where the key is the name of the variable, and the value is
      a list of size 2 of the form [lower_bound, upper_bound]
  - Import data (fitter.ImportData): Copy data as a class variable.
      Must be given as a DataStructures.xypoint instance
  - Perform the fit: (fitter.Fit):
      Returns a DataStructures.xypoint instance of the model. The 
      x-values in the returned array are the same as the data.
   - Optional: retrieve a new version of the data, which is 
      wavelength-calibrated using the telluric lines and with
      a potentially better continuum fit using
      data2 = fitter.data  


      

    This file is part of the TelFit program.

    TelFit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TelFit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TelFit.  If not, see <http://www.gnu.org/licenses/>.


"""


import matplotlib.pyplot as plt
import numpy
import sys
import os
import subprocess
import scipy
from scipy.interpolate import  UnivariateSpline
from scipy.optimize import leastsq, minimize, fminbound
from scipy.linalg import svd, diagsvd
from scipy import mat
import MakeModel
import DataStructures
import FittingUtilities



class TelluricFitter:
  def __init__(self, debug=False, debug_level=2):
    #Set up parameters
    self.parnames = ["pressure", "temperature", "angle", "resolution", "wavestart", "waveend",
                     "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "no",
                     "so2", "no2", "nh3", "hno3"]
    self.const_pars = [795.0, 273.0, 45.0, 50000.0, 2200.0, 2400.0,
                       50.0, 368.5, 3.9e-2, 0.32, 0.14, 1.8, 2.1e5, 1.1e-19,
                       1e-4, 1e-4, 1e-4, 5.6e-4]
    self.bounds = [[0.0, 1e30] for par in self.parnames]  #Basically just making sure everything is > 0
    self.fitting = [False]*len(self.parnames)
    
    #Latitude and altitude (to nearest km) of the observatory
    #  Defaults are for McDonald Observatory
    self.observatory = {"latitude": 30.6,
                        "altitude": 2.0}
    self.data = None
    self.resolution_bounds = [10000.0, 100000.0]

    homedir = os.environ['HOME']
    self.resolution_fit_mode="SVD"
    self.fit_primary = False
    self.adjust_wave = "model"
    self.first_iteration=True
    self.continuum_fit_order = 7
    self.wavelength_fit_order = 3
    self.debug = debug
    self.debug_level = debug_level   #Number from 1-5, with 5 being the most verbose
    self.Modeler = MakeModel.Modeler(debug=self.debug)
    self.parvals = [[] for i in range(len(self.parnames))]
    self.chisq_vals = []
    self.ignore = []
    self.shift = 0   #The wavelength shift to make the model and data align

    #Just open and close chisq_summary, to clear anything already there
    outfile = open("chisq_summary.dat", "w")
    outfile.close()


    
  """
    Display the value of each of the parameters, and show whether it is being fit or not
  """
  def DisplayVariables(self, fitonly=False):
    print "%.15s\tValue\t\tFitting?\tBounds" %("Parameter".ljust(15))
    print "-------------\t-----\t\t-----\t\t-----"
    for i in range(len(self.parnames)):
      if (fitonly and self.fitting[i]) or not fitonly:
        if len(self.bounds[i]) == 2:
          print "%.15s\t%.5E\t%s\t\t%g - %g" %(self.parnames[i].ljust(15), self.const_pars[i], self.fitting[i], self.bounds[i][0], self.bounds[i][1])
        else:
          print "%.15s\t%.5g\t\t%s" %(self.parnames[i].ljust(15), self.const_pars[i], self.fitting[i])


  """
    Add one or more variables to the list being fit. vardict must be a dictionary
      where the key is the parameter name and the value is the value of that parameter.
  """
  def FitVariable(self, vardict):
    for par in vardict.keys():
      try:
        idx = self.parnames.index(par)
        self.const_pars[idx] = vardict[par]
        self.fitting[idx] = True
      except ValueError:
        print "Error! Bad parameter name given. Currently available are: "
        self.DisplayVariables()
        raise ValueError

        
  """
    Similar to FitVariable, but this just adjusts the value of a constant parameter
  """
  def AdjustValue(self, vardict):
    for par in vardict.keys():
      try:
        idx = self.parnames.index(par)
        self.const_pars[idx] = vardict[par]
      except ValueError:
        print "Error! Bad parameter name given. Currently available are: "
        self.DisplayVariables()
        raise ValueError

  """
    Similar to FitVariable, but it sets bounds on the variable. This can technically
      be done for any variable, but is only useful to set bounds for those variables
      being fit (and detector resolution)
  """
  def SetBounds(self, bounddict):
    for par in bounddict.keys():
      try:
        idx = self.parnames.index(par)
        self.bounds[idx] = bounddict[par]
        if par == "resolution":
          self.resolution_bounds = bounddict[par]
      except ValueError:
        print "Error! Bad parameter name given. Currently available are: "
        self.DisplayVariables()
        raise ValueError


  """
    Set the observatory. Can either give a dictionary with the latitude and altitude,
      or give the name of the observatory. Some names are hard-coded in here.
  """
  def SetObservatory(self, observatory):
    if type(observatory) == str:
      if observatory.lower() == "ctio":
        self.observatory["latitude"] = -30.6
        self.observatory["altitude"] = 2.2
      if observatory.lower() == "la silla":
        self.observatory["latitude"] = -29.3
        self.observatory["altitude"] = 2.4
      if observatory.lower() == "paranal":
        self.observatory["latitude"] = -24.6
        self.observatory["altitude"] = 2.6
      if observatory.lower() == "mauna kea":
        self.observatory["latitude"] = 19.8
        self.observatory["altitude"] = 4.2
      if observatory.lower() == "mcdonald":
        self.observatory["latitude"] = 30.7
        self.observatory["altitude"] = 2.1
        
    elif type(observatory) == dict:
      if "latitude" in observatory.keys() and "altitude" in observatory.keys():
        self.observatory = observatory
      else:
        print "Error! Wrong keys in observatory dictionary! Keys must be"
        print "'latitude' and 'altitude'. Yours are: ", observatory.keys()
        raise KeyError
    else:
      raise ValueError("Error! Unrecognized input to TelluricFitter.SetObservatory()")
    

  """
    Function for the user to give the data. The data should be in the form of
      a DataStructures.xypoint structure.
  """
  def ImportData(self, data):
    if not isinstance(data, DataStructures.xypoint):
      raise TypeError( "ImportData Error! Given data is not a DataStructures.xypoint structure!" )
    self.data = data.copy()
    return



  """
    Edits the atmosphere profile for a given parameter. This is just a wrapper
      for the MakeModel.Modeler method.
  """
  def EditAtmosphereProfile(self, profilename, profile_height, profile_value):
    self.Modeler.EditProfile(profilename, profile_height, profile_value)
    
  
  """
    Tells the fitter to ignore certain regions of the spectrum
      in the chi-squared calculation. Useful for stellar or interstellar
      lines.
  """
  def IgnoreRegions(self, region):
    if not isinstance(region, list) or len(region) == 0:
      raise TypeError("Must give a non-empty list to TelluricFitter.IgnoreRegions")
    
    if isinstance(region[0], list):
      #The user gave a list of lists. Append each one to self.ignore
      for r in region:
        self.ignore.append(r)
    elif isinstance(region[0], int):
      #The user gave a single region. Append to self.ignore
      self.ignore.append(region)
    else:
      raise TypeError("Unrecognized variable type for region given in TelluricFitter.IgnoreRegions")
    
    return



  """
    Finally, the main fitting function. Before calling this, the user MUST
      1: call FitVariable at least once, specifying which variables will be fit
      2: import data into the class using the ImportData method.
    resolution_fit_mode controls which function is used to estimate the resolution.
      3: Set resolution bounds (any other bounds are optional)
      "SVD" is for singlular value decomposition, while "gauss" is for convolving with a gaussian
      (and fitting the width of the guassian to give the best fit)

    continuum_fit_mode controls how the continuum is fit in the data. Choices are 'polynomial' and 'smooth'

    fit_primary determines whether an iterative smoothing is applied to the data to approximate the primary star (only works for primary stars with broad lines)
    
    return_resolution controls whether the best-fit resolution is returned to the user or not.

    adjust_wave can be set to either 'data' or 'model'. To wavelength calibrate the data to the telluric lines, set to 'data'. If you think the wavelength calibration is good on the data (such as Th-Ar lines in the optical), then set to 'model' Note that currently, the vacuum --> air conversion for the telluric model is done in a very approximate sense, so adjusting the data wavelengths may introduce a small (~1 km/s) offset from what it should be.
  """
  def Fit(self, resolution_fit_mode="SVD", fit_primary=False, return_resolution=False, adjust_wave="model", continuum_fit_order=7, wavelength_fit_order=3):

    self.resolution_fit_mode=resolution_fit_mode
    self.fit_primary = fit_primary
    self.adjust_wave = adjust_wave
    self.continuum_fit_order = continuum_fit_order
    self.wavelength_fit_order = wavelength_fit_order
    self.return_resolution=return_resolution

    #Check to make sure the user gave data to fit
    if self.data == None:
      print "\n\nError! Must supply data to fit\n\n!"
      return

    #Make sure resolution bounds are given (resolution is always fit)
    idx = self.parnames.index("resolution")
    if len(self.bounds[idx]) < 2 and self.resolution_fit_mode != "SVD":
      print "Must give resolution bounds!"
      inp = raw_input("Enter the lowest and highest possible resolution, separated by a space: ")
      self.resolution_bounds = [float(inp.split()[0]), float(inp.split()[1])]


    #Make fitpars array
    fitpars = [self.const_pars[i] for i in range(len(self.parnames)) if self.fitting[i] ]
    if len(fitpars) < 1:
      print "\n\nError! Must fit at least one variable!\n\n"
      return
    

    #Set up the fitting logfile and logging arrays
    self.parvals = [[] for i in range(len(self.parnames))]
    self.chisq_vals = []
    outfile = open("chisq_summary.dat", "a")
    outfile.write("\n\n\n\n")
    for i in range(len(self.parnames)):
      if self.fitting[i]:
        outfile.write("%s\t" %self.parnames[i])
    outfile.write("X^2\n")
    outfile.close()


    #Perform the fit
    self.first_iteration = True
    errfcn = lambda pars: numpy.sum(self.FitErrorFunction(pars))
    bounds = [self.bounds[i] for i in range(len(self.parnames)) if self.fitting[i]]
    optdict = {"rhobeg": [1,5,1000.0]}
    optdict = {"eps": 5}
    #output = minimize(errfcn, fitpars, method='SLSQP', bounds=bounds, options=optdict, tol=0.001)
    
    #print "Fitting done with output message: \n%s" %output.message
    #fitpars = output.x
    fitpars, success = leastsq(self.FitErrorFunction, fitpars, diag=1.0/numpy.array(fitpars))

    #Finally, return the best-fit model
    if self.fit_primary:
      return self.GenerateModel(fitpars, separate_primary=True, return_resolution=return_resolution)
    else:
      return self.GenerateModel(fitpars, return_resolution=return_resolution)
    


  """
    The error function for the fitter. This should never be called directly!
  """
  def FitErrorFunction(self, fitpars):
    if self.return_resolution:
      model, resolution = self.GenerateModel(fitpars, return_resolution=True)
    else:
      model = self.GenerateModel(fitpars)
    outfile = open("chisq_summary.dat", 'a')
    weights = 1.0 / self.data.err**2

    #Find the regions to use (ignoring the parts that were defined as bad)
    good = numpy.arange(self.data.x.size, dtype=numpy.int32)
    for region in self.ignore:
      x0 = min(region)
      x1 = max(region)
      tmp1 = [self.data.x[i] in self.data.x[good] for i in range(self.data.x.size)]
      tmp2 = numpy.logical_or(self.data.x<x0, self.data.x>x1)
      good = numpy.where(numpy.logical_and(tmp1, tmp2))[0]

    
    return_array = (self.data.y - self.data.cont*model.y)[good]**2 * weights[good]
    #Evaluate bound conditions and output the parameter value to the logfile.
    fit_idx = 0
    for i in range(len(self.bounds)):
      if self.fitting[i]:
        if len(self.bounds[i]) == 2:
          return_array += FittingUtilities.bound(self.bounds[i], fitpars[fit_idx])
        outfile.write("%.12g\t" %fitpars[fit_idx])
        self.parvals[i].append(fitpars[fit_idx])
        fit_idx += 1
      elif len(self.bounds[i]) == 2 and self.parnames[i] != "resolution":
        return_array += FittingUtilities.bound(self.bounds[i], self.const_pars[i])
    outfile.write("%g\n" %(numpy.sum(return_array)/float(weights.size)))
    
    self.chisq_vals.append(numpy.sum(return_array)/float(weights.size))
    print "X^2 = ", numpy.sum(return_array)/float(weights.size)
    outfile.close()
    
    return return_array




  """
  This function does the actual work of generating a model with the given parameters,
    fitting the continuum, making sure the model and data are well aligned in
    wavelength, and fitting the detector resolution

  """
  def GenerateModel(self, pars, nofit=False, separate_primary=False, return_resolution=False):
    data = self.data
    #Update self.const_pars to include the new values in fitpars
    #  I know, it's confusing that const_pars holds some non-constant parameters...
    fit_idx = 0
    for i in range(len(self.parnames)):
      if self.fitting[i]:
        self.const_pars[i] = pars[fit_idx]
        fit_idx += 1
    self.DisplayVariables(fitonly=True)
    
    #Extract parameters from pars and const_pars. They will have variable
    #  names set from self.parnames
    fit_idx = 0
    for i in range(len(self.parnames)):
      #Assign to local variables by the parameter name
      if self.fitting[i]:
        exec("%s = %g" %(self.parnames[i], pars[fit_idx]))
        fit_idx += 1
      else:
        exec("%s = %g" %(self.parnames[i], self.const_pars[i]))
      
      #Make sure everything is within its bounds
      if len(self.bounds[i]) > 0:
        lower = self.bounds[i][0]
        upper = self.bounds[i][1]
        exec("%s = %g if %s < %g else %s" %(self.parnames[i], lower, self.parnames[i], lower, self.parnames[i]))
        exec("%s = %g if %s > %g else %s" %(self.parnames[i], upper, self.parnames[i], upper, self.parnames[i]))

    
    wavenum_start = 1e7/waveend
    wavenum_end = 1e7/wavestart
    lat = self.observatory["latitude"]
    alt = self.observatory["altitude"]
    

    #Generate the model:
    model = self.Modeler.MakeModel(pressure, temperature, wavenum_start, wavenum_end, angle, h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, lat=lat, alt=alt, wavegrid=None, resolution=None)
    
    #Shift the x-axis, using the shift from previous iterations
    if self.adjust_wave == "data":
      data.x += self.shift
    elif self.adjust_wave == "model":
      model.x -= self.shift

    #Save each model if debugging
    if self.debug and self.debug_level >= 5:
      model_name = "Models/transmission"+"-%.2f" %pressure + "-%.2f" %temperature + "-%.1f" %h2o + "-%.1f" %angle + "-%.2f" %(co2) + "-%.2f" %(o3*100) + "-%.2f" %ch4 + "-%.2f" %(co*10)
      numpy.savetxt(model_name, numpy.transpose((model.x, model.y)), fmt="%.8f")
      
    #Interpolate to constant wavelength spacing
    xgrid = numpy.linspace(model.x[0], model.x[-1], model.x.size)
    model = FittingUtilities.RebinData(model, xgrid)

    #Use nofit if you want a model with reduced resolution. Probably easier
    #  to go through MakeModel directly though...
    if data == None or nofit:
      return FittingUtilities.ReduceResolution(model, resolution)

    

    model_original = model.copy()
  
    #Reduce to initial guess resolution
    if "SVD" in self.resolution_fit_mode and not self.first_iteration:
      broadening_fcn = self.broadstuff[0]
      xarr = self.broadstuff[1]
      Model = UnivariateSpline(model_original.x, model_original.y, s=0)
      model_new = Model(xarr)
      model = DataStructures.xypoint(x=xarr)
      Broadened = UnivariateSpline(xarr, numpy.convolve(model_new, broadening_fcn, mode="same"),s=0)
      model.y = Broadened(model.x)
      model = FittingUtilities.RebinData(model, data.x)
      
    elif "gauss" in self.resolution_fit_mode or self.first_iteration:
      if (resolution - 10 < self.resolution_bounds[0] or resolution+10 > self.resolution_bounds[1]):
        resolution = numpy.mean(self.resolution_bounds)
      model = FittingUtilities.ReduceResolution(model.copy(), resolution)
      model = FittingUtilities.RebinData(model.copy(), data.x.copy())

    else:
      sys.exit("Error! Unrecognized resolution fit mode: %s" %self.resolution_fit_mode)

    
    #Shift the data (or model) by a constant offset. This gets the wavelength calibration close
    shift = FittingUtilities.CCImprove(data, model, tol=0.1)
    if self.adjust_wave == "data":
      data.x += shift
    elif self.adjust_wave == "model":
      model_original.x -= shift
    else:
      sys.exit("Error! adjust_wave parameter set to invalid value: %s" %self.adjust_wave)
    self.shift += shift

      
    #Need to reduce resolution to initial guess again if the model was shifted
    if self.adjust_wave == "model":
      if "SVD" in self.resolution_fit_mode and not self.first_iteration:
        Model = UnivariateSpline(model_original.x, model_original.y, s=0)
        model_new = Model(xarr)
        model = DataStructures.xypoint(x=xarr)
        Broadened = UnivariateSpline(xarr, numpy.convolve(model_new, broadening_fcn, mode="same"),s=0)
        model.y = Broadened(model.x)
        model = FittingUtilities.RebinData(model, data.x)
      
      elif "gauss" in self.resolution_fit_mode or self.first_iteration:
        model = FittingUtilities.ReduceResolution(model_original.copy(), resolution)
        model = FittingUtilities.RebinData(model.copy(), data.x.copy())

    
    resid = data.y/model.y
    nans = numpy.isnan(resid)
    resid[nans] = data.cont[nans]
    """
    if self.fit_primary:
      data2 = data.copy()
      data2.y /= model.y
      primary_star = data.copy()
      primary_star.y = FittingUtilities.savitzky_golay(data2.y, 91, 4)
      primary_star.y /= primary_star.y.mean()
      PRIMARY_STAR = UnivariateSpline(primary_star.x, primary_star.y, s=0)

      model2 = model.copy()
      model2.y *= primary_star.y
      resid /= primary_star.y
    """
      
    #As the model gets better, the continuum will be less affected by
    #  telluric lines, and so will get better
    data.cont = FittingUtilities.Continuum(data.x, resid, fitorder=self.continuum_fit_order, lowreject=2, highreject=10)
    
    if self.fit_primary:
      primary_star = data.copy()
      primary_star.y = FittingUtilities.savitzky_golay(resid/data.cont, 61, 4)
      data.cont /= primary_star.y



    
    if self.debug and self.debug_level >= 4:
      print "Saving data and model arrays right before fitting the wavelength"
      print "  and resolution to Debug_Output1.log"
      numpy.savetxt("Debug_Output1.log", numpy.transpose((data.x, data.y, data.cont, model.x, model.y)))

      
    #Fine-tune the wavelength calibration by fitting the location of several telluric lines
    #if self.fit_primary:
    if 1>2:
      modelfcn, mean = self.FitWavelength(data, model2.copy(), fitorder=self.wavelength_fit_order)
    else:
      modelfcn, mean = self.FitWavelength(data, model.copy(), fitorder=self.wavelength_fit_order)
      
    if self.adjust_wave == "data":
      test = modelfcn(data.x - mean)
      xdiff = [test[j] - test[j-1] for j in range(1, len(test)-1)]
      if min(xdiff) > 0 and numpy.max(test - data.x) < 0.1:
        print "Adjusting data wavelengths by at most %.8f" %numpy.max(test - model.x)
        data.x = test.copy()
      else:
        print "Warning! Wavelength calibration did not succeed!"
    elif self.adjust_wave == "model":
      test = modelfcn(model_original.x - mean)
      test2 = modelfcn(model.x - mean)
      xdiff = [test[j] - test[j-1] for j in range(1, len(test)-1)]
      if min(xdiff) > 0 and numpy.max(test2 - model.x) < 0.1:
        model.x = test2.copy()
        model_original.x = test.copy()
        print "Adjusting model wavelengths by at most %.8f" %numpy.max(test2 - model.x)
      else:
        print "Warning! Wavelength calibration did not succeed!"
    else:
      sys.exit("Error! adjust_wave set to an invalid value: %s" %self.adjust_wave)

    if self.debug and self.debug_level >= 4:
      print "Saving data and model arrays after fitting the wavelength"
      print "  and before the resolution fit to Debug_Output2.log"
      numpy.savetxt("Debug_Output2.log", numpy.transpose((data.x, data.y, data.cont, model.x, model.y)))

      
    #Fit instrumental resolution
    done = False
    while not done:
      done = True
      if "SVD" in self.resolution_fit_mode:
        #if self.fit_primary:
	if 1>2:
          model2 = model_original.copy()
          prim = PRIMARY_STAR(model2.x)
          prim[prim < 0.0] = 0.0
          prim[prim > 10.0] = 10.0
          model2.y *= prim
          model, self.broadstuff = self.Broaden(data.copy(), model2, full_output=True)
        else:
          model, self.broadstuff = self.Broaden(data.copy(), model_original.copy(), full_output=True)
      elif "gauss" in self.resolution_fit_mode:
        #if self.fit_primary:
	if 1>2:
	  model2 = model_original.copy()
          prim = PRIMARY_STAR(model2.x)
          prim[prim < 0.0] = 0.0
          prim[prim > 10.0] = 10.0
          model2.y *= prim
          model, resolution = self.FitResolution(data.copy(), model2.y, resolution)
        else:
          model, resolution = self.FitResolution(data.copy(), model_original.copy(), resolution)
      else:
        done = False
        print "Resolution fit mode set to an invalid value: %s" %self.resolution_fit_mode
        self.resolution_fit_mode = raw_input("Enter a valid mode (SVD or guass): ")
    
    
    self.data = data
    self.first_iteration = False
    if separate_primary:
      if return_resolution:
        return primary_star, model, resolution
      else:
        return primary_star, model
    elif return_resolution:
      return model, resolution
    else:
      return model


  
  #Wavelength-fitting function that just shifts lines, instead of fitting them to gaussians
  def WavelengthErrorFunction(self, shift, data, model):
    modelfcn = UnivariateSpline(model.x, model.y, s=0)
    weight = 1e9 * numpy.ones(data.x.size)
    weight[data.y > 0] = 1.0/numpy.sqrt(data.y[data.y > 0])
    weight[weight < 0.01] = 0.0
    newmodel = modelfcn(model.x + float(shift))
    if shift < 0:
      newmodel[model.x - float(shift) < model.x[0]] = 0
    else:
      newmodel[model.x - float(shift) > model.x[-1]] = 0
    returnvec = (data.y - newmodel)**2*weight
    return returnvec


  #Gaussian absorption line
  def GaussianFitFunction(self, x,params):
    cont = params[0]
    depth = params[1]
    mu = params[2]
    sig = params[3]
    return cont - depth*numpy.exp(-(x-mu)**2/(2*sig**2))


  #Returns the residuals between the fit from above and the actual values
  def GaussianErrorFunction(self, params, x, y):
    return self.GaussianFitFunction(x,params) - y
  
  """
    Function to fine-tune the wavelength solution of a generated model
      It does so looking for telluric lines in both the
      data and the telluric model. For each line, it finds the shift needed
      to make them line up, and then fits a function to that fit over the
      full wavelength range of the data.

    Wavelength calibration MUST already be very close for this algorithm
      to succeed!
  """
  def FitWavelength(self, data_original, telluric, tol=0.05, oversampling=4, fitorder=3):
    print "Fitting Wavelength"
    old = []
    new = []
    #Find lines in the telluric model
    linelist = FittingUtilities.FindLines(telluric, debug=self.debug, tol=0.995)
    if len(linelist) < fitorder:
      fit = lambda x: x
      mean = 0.0
      return fit, mean
    linelist = telluric.x[linelist]
    
    if self.debug and self.debug_level >= 5:
      logfilename = "FitWavelength.log"
      print "Outputting data and telluric model to %s" %logfilename
      numpy.savetxt(logfilename, numpy.transpose((data_original.x, data_original.y, data_original.cont, data_original.err)), fmt="%.8f")
      infile = open(logfilename, "a")
      infile.write("\n\n\n\n\n")
      numpy.savetxt(infile, numpy.transpose((telluric.x, telluric.y)), fmt="%.8f")
      infile.close()

    #Interpolate to finer spacing
    xgrid = numpy.linspace(data_original.x[0], data_original.x[-1], data_original.x.size*oversampling)
    data = FittingUtilities.RebinData(data_original, xgrid)
    model = FittingUtilities.RebinData(telluric, xgrid)
  
    #Begin loop over the lines
    numlines = 0
    for line in linelist:
      if line-tol > data.x[0] and line+tol < data.x[-1]:
        numlines += 1
        #Find line in the model
        left = numpy.searchsorted(model.x, line - tol)
        right = numpy.searchsorted(model.x, line + tol)
        minindex = model.y[left:right].argmin() + left

        mean = model.x[minindex]
        left2 = numpy.searchsorted(model.x, mean - tol*2)
        right2 = numpy.searchsorted(model.x, mean + tol*2)

        argmodel = DataStructures.xypoint(right2 - left2)
        argmodel.x = numpy.copy(model.x[left2:right2])
        argmodel.y = numpy.copy(model.y[left2:right2])

        #Do the same for the data
        left = numpy.searchsorted(data.x, line - tol)
        right = numpy.searchsorted(data.x, line + tol)
        minindex = data.y[left:right].argmin() + left

        mean = data.x[minindex]

        argdata = DataStructures.xypoint(right2 - left2)
        argdata.x = numpy.copy(data.x[left2:right2])
        argdata.y = numpy.copy(data.y[left2:right2]/data.cont[left2:right2])

        #Fit argdata to gaussian to find the actual line location:
        cont = 1.0
        depth = cont - argdata.y[argdata.y.size/2]
        mu = argdata.x[argdata.x.size/2]
        sig = 0.025
        params = [cont, depth, mu, sig]
        params,success = leastsq(self.GaussianErrorFunction, params, args=(argdata.x, argdata.y))
        mean = params[2]
        
        #Do a cross-correlation, to get the wavelength solution close
        ycorr = scipy.correlate(argdata.y-1.0, argmodel.y-1.0, mode="full")
        xcorr = numpy.arange(ycorr.size)
        maxindex = ycorr.argmax()
        lags = xcorr - (argdata.x.size-1)
        distancePerLag = (argdata.x[-1] - argdata.x[0])/float(argdata.x.size)
        offsets = -lags*distancePerLag
        shift = offsets[maxindex]

        #Now, fit the shift more precisely
        shift, success = leastsq(self.WavelengthErrorFunction, shift, args=(argdata, argmodel))
        
        if self.debug and self.debug_level >= 3:
          print argdata.x[0], argdata.x[-1], argdata.x.size
          print "wave: ", mean, "\tshift: ", shift, "\tsuccess = ", success
          plt.figure(1)
          plt.plot(model.x[left:right]-shift, model.y[left:right])
          plt.plot(argmodel.x, argmodel.y)
          plt.plot(argdata.x, argdata.y)
          
        if (success < 5):
          old.append(mean)
          if self.adjust_wave == "data":
            new.append(mean + float(shift))
          elif self.adjust_wave == "model":
            new.append(mean - float(shift))
          else:
            sys.exit("Error! adjust_wave set to an invalid value: %s" %self.adjust_wave)

            
    if self.debug and self.debug_level >= 3:
      plt.figure(2)
      plt.plot(old, new, 'ro')
      plt.title("Fitted Line shifts")
      plt.xlabel("Old Wavelength")
      plt.ylabel("New Wavelength")

    numlines = len(old)
    print "Found %i lines in this order" %numlines
    fit = lambda x: x
    mean = 0.0
    if numlines < fitorder:
      return fit, mean
    
    #Check if there is a large gap between the telluric lines and the end of the order (can cause the fit to go crazy)
    keepfirst = False
    keeplast = False
    if min(old) - data.x[0] > 0.5:
      old.insert(0, data.x[0])
      new.insert(0, data.x[0])
      keepfirst = True
    if data.x[-1] - max(old) > 0.5:
      old.append(data.x[-1])
      new.append(data.x[-1])
      keeplast = True
      
    #Iteratively fit to a cubic with sigma-clipping
    done = False
    while not done and len(old) >= fitorder:
      done = True
      mean = numpy.mean(old)
      fit = numpy.poly1d(numpy.polyfit(old - mean, new, fitorder))
      residuals = fit(old - mean) - new
      std = numpy.std(residuals)
      badindices = numpy.where(numpy.logical_or(residuals > 2*std, residuals < -2*std))[0]
      for badindex in badindices[::-1]:
        if (badindex == 0 and keepfirst) or (badindex == len(residuals)-1 and keeplast):
          continue
        del old[badindex]
        del new[badindex]
        done = False
    if self.debug and self.debug_level >= 3:
      plt.figure(3)
      plt.plot(old, fit(old - mean) - new, 'ro')
      plt.title("Residuals")
      plt.xlabel("Wavelength")
      plt.ylabel("Delta-lambda")
      plt.figure(4)
      plt.plot(data_original.x, data_original.y/data_original.cont)
      plt.plot(telluric.x, telluric.y)
      plt.plot(fit(telluric.x-mean), telluric.y)
      plt.show()
      
    return fit, mean



  """
    Fits the instrumental resolution with a Gaussian
  """
  def FitResolution(self, data, model, resolution=75000.0):
    ####resolution is the initial guess####

    print "Fitting Resolution"

    #Subsample the model to speed this part up (it doesn't affect the accuracy much)
    xgrid = numpy.linspace(model.x[0], model.x[-1], model.size()/5)
    newmodel = FittingUtilities.RebinData(model, xgrid)

    ResolutionFitErrorBrute = lambda resolution, data, model: numpy.sum(self.ResolutionFitError(resolution, data, model))
    
    resolution = fminbound(ResolutionFitErrorBrute, self.resolution_bounds[0], self.resolution_bounds[1], xtol=1, args=(data,newmodel))
    
    print "Optimal resolution found at R = ", float(resolution)
    newmodel = FittingUtilities.ReduceResolution(newmodel, float(resolution))
    return FittingUtilities.RebinData(newmodel, data.x), float(resolution)

  
  """
    This function gets called by scipy.optimize.fminbound in FitResolution (above)
  """
  def ResolutionFitError(self, resolution, data, model):
    resolution = max(1000.0, float(int(float(resolution) + 0.5)))
    if self.debug and self.debug_level >= 5:
      print "Saving inputs for R = ", resolution
      print " to Debug_ResFit.log and Debug_ResFit2.log"
      numpy.savetxt("Debug_ResFit.log", numpy.transpose((data.x, data.y, data.cont)))
      numpy.savetxt("Debug_Resfit2.log", numpy.transpose((model.x, model.y)))
    newmodel = FittingUtilities.ReduceResolution(model, resolution, extend=False)
    newmodel = FittingUtilities.RebinData(newmodel, data.x)
    weights = 1.0/data.err**2
    returnvec = (data.y - data.cont*newmodel.y)**2*weights + FittingUtilities.bound(self.resolution_bounds, resolution)
    if self.debug:
      print "Resolution-fitting X^2 = ", numpy.sum(returnvec)/float(weights.size), "at R = ", resolution
    if numpy.isnan(numpy.sum(returnvec**2)):
      print "Error! NaN found in ResolutionFitError!"
      outfile=open("ResolutionFitError.log", "a")
      outfile.write("#Error attempting R = %g\n" %(resolution))
      numpy.savetxt(outfile, numpy.transpose((data.x, data.y, data.cont, newmodel.x, newmodel.y)), fmt="%.10g")
      outfile.write("\n\n\n\n")
      numpy.savetxt(outfile, numpy.transpose((model.x, model.y)), fmt="%.10g")
      outfile.write("\n\n\n\n")
      outfile.close()
      raise ValueError
    return returnvec

  


  """
    -Fits the broadening profile using singular value decomposition
    -oversampling is the oversampling factor to use before doing the SVD
    -m is the size of the broadening function, in oversampled units
    -dimension is the number of eigenvalues to keep in the broadening function. (Keeping too many starts fitting noise)

    -NOTE: This function works well when there are strong telluric lines and a flat continuum.
           If there are weak telluric lines, it's hard to not fit noise.
           If the continuum is not very flat (i.e. from the spectrum of the actual
             object you are trying to telluric correct), the broadening function
             can become multiply-peaked and oscillatory. Use with care!
  """
  def Broaden(self, data, model, oversampling = 5, m = 101, dimension = 20, full_output=False):
    n = data.x.size*oversampling
    
    #n must be even, and m must be odd!
    if n%2 != 0:
      n += 1
    if m%2 == 0:
      m += 1
  
    #resample data
    Spectrum = UnivariateSpline(data.x, data.y/data.cont, s=0)
    Model = UnivariateSpline(model.x, model.y, s=0)
    xnew = numpy.linspace(data.x[0], data.x[-1], n)
    ynew = Spectrum(xnew)
    model_new = FittingUtilities.RebinData(model, xnew).y

    #Make 'design matrix'
    design = numpy.zeros((n-m,m))
    for j in range(m):
      for i in range(m/2,n-m/2-1):
        design[i-m/2,j] = model_new[i-j+m/2]
    design = mat(design)
    
    #Do Singular Value Decomposition
    try:
      U,W,V_t = svd(design, full_matrices=False)
    except numpy.linalg.linalg.LinAlgError:
      outfilename = "SVD_Error.log"
      outfile = open(outfilename, "a")
      numpy.savetxt(outfile, numpy.transpose((data.x, data.y, data.cont)))
      outfile.write("\n\n\n\n\n")
      numpy.savetxt(outfile, numpy.transpose((model.x, model.y, model.cont)))
      outfile.write("\n\n\n\n\n")
      outfile.close()
      sys.exit("SVD did not converge! Outputting data to %s" %outfilename)
      
    #Invert matrices:
    #   U, V are orthonormal, so inversion is just their transposes
    #   W is a diagonal matrix, so its inverse is 1/W
    W1 = 1.0/W
    U_t = numpy.transpose(U)
    V = numpy.transpose(V_t)
  
    #Remove the smaller values of W
    W1[dimension:] = 0
    W2 = diagsvd(W1,m,m)
    
    #Solve for the broadening function
    spec = numpy.transpose(mat(ynew[m/2:n-m/2-1]))
    temp = numpy.dot(U_t, spec)
    temp = numpy.dot(W2,temp)
    Broadening = numpy.dot(V,temp)
    
    #Make Broadening function a 1d array
    spacing = xnew[2] - xnew[1]
    xnew = numpy.arange(model.x[0], model.x[-1], spacing)
    model_new = Model(xnew)
    Broadening = numpy.array(Broadening)[...,0]
    
    #Ensure that the broadening function is appropriate:
    maxindex = Broadening.argmax()
    if maxindex > m/2.0 + m/10.0 or maxindex < m/2.0 - m/10.0:
      #The maximum should be in the middle because we already did wavelength calibration!
      outfilename = "SVD_Error2.log"
      numpy.savetxt(outfilename, numpy.transpose((Broadening, )) )
      print "Warning! SVD Broadening function peaked at the wrong location! See SVD_Error2.log for the broadening function"
      
      idx = self.parnames.index("resolution")
      resolution = self.const_pars[idx]
      model = FittingUtilities.ReduceResolution(model, resolution)
      
      #Make broadening function from the gaussian
      centralwavelength = (data.x[0] + data.x[-1])/2.0
      FWHM = centralwavelength/resolution;
      sigma = FWHM/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))
      left = 0
      right = numpy.searchsorted(xnew, 10*sigma)
      x = numpy.arange(0,10*sigma, xnew[1] - xnew[0])
      gaussian = numpy.exp(-(x-5*sigma)**2/(2*sigma**2))
      return FittingUtilities.RebinData(model, data.x), [gaussian/gaussian.sum(), xnew]
      
    elif numpy.mean(Broadening[maxindex-int(m/10.0):maxindex+int(m/10.0)]) < 3* numpy.mean(Broadening[int(m/5.0):]):
      outfilename = "SVD_Error2.log"
      numpy.savetxt(outfilename, numpy.transpose((Broadening, )) )
      print "Warning! SVD Broadening function is not strongly peaked! See SVD_Error2.log for the broadening function"
      
      idx = self.parnames.index("resolution")
      resolution = self.const_pars[idx]
      model = FittingUtilities.ReduceResolution(model, resolution)
      
      #Make broadening function from the gaussian
      centralwavelength = (data.x[0] + data.x[-1])/2.0
      FWHM = centralwavelength/resolution;
      sigma = FWHM/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))
      left = 0
      right = numpy.searchsorted(xnew, 10*sigma)
      x = numpy.arange(0,10*sigma, xnew[1] - xnew[0])
      gaussian = numpy.exp(-(x-5*sigma)**2/(2*sigma**2))
      return FittingUtilities.RebinData(model, data.x), [gaussian/gaussian.sum(), xnew]
    
    #If we get here, the broadening function looks okay.
    #Convolve the model with the broadening function
    model = DataStructures.xypoint(x=xnew)
    Broadened = UnivariateSpline(xnew, numpy.convolve(model_new,Broadening, mode="same"),s=0)
    model.y = Broadened(model.x)
    
    #Fit the broadening function to a gaussian
    params = [0.0, -Broadening[maxindex], maxindex, 10.0]
    params,success = leastsq(self.GaussianErrorFunction, params, args=(numpy.arange(Broadening.size), Broadening))
    sigma = params[3] * (xnew[1] - xnew[0]) 
    FWHM = sigma * 2.0*numpy.sqrt(2.0*numpy.log(2.0))
    resolution = numpy.median(data.x) / FWHM
    idx = self.parnames.index("resolution")
    self.const_pars[idx] = resolution
    
    print "Approximate resolution = %g" %resolution
    
    x2 = numpy.arange(Broadening.size)

    if full_output:
      return FittingUtilities.RebinData(model, data.x), [Broadening, xnew]
    else:
      return FittingUtilities.RebinData(model, data.x)

