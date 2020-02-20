"""
    A set of useful functions that I use often while fitting. Most, but not all,
      of these functions are used in the TelFit program. Note that this is
      Cython code, which needs to be compiled!



    This file is part of the TelFit program.

    TelFit is free software: you can redistribute it and/or modify
    it under the terms of the MIT license.

    TelFit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the MIT license
    along with TelFit.  If not, see <http://opensource.org/licenses/MIT>.
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.interpolate import UnivariateSpline as smoother
from scipy.signal import argrelmin, fftconvolve
import matplotlib.pyplot as plt
from pysynphot.observation import Observation
from pysynphot.spectrum import ArraySourceSpectrum, ArraySpectralElement

import DataStructures


cimport numpy as np
cimport cython
from libc.math cimport exp, log, sqrt
import os
from astropy import units as u

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


#Define bounding functions:
# lower bound:            lbound(boundary_value, parameter)
# upper bound:            ubound(boundary_value, parameter)
# lower and upper bounds: bound([low, high], parameter)
# fixed parameter:        fixed(fixed_value, parameter)
lbound = lambda p, x: 1e2*np.sqrt(p-x) + 1e-3*(p-x) if (x<p) else 0
ubound = lambda p, x: 1e2*np.sqrt(x-p) + 1e-3*(x-p) if (x>p) else 0
bound  = lambda p, x: lbound(p[0],x) + ubound(p[1],x)
fixed  = lambda p, x: bound((p,p), x)



def ensure_dir(f):
  """
    Ensure that a directory exists. Create if it doesn't
  """
  d = os.path.dirname(f)
  if not os.path.exists(d):
    os.makedirs(d)



def CCImprove(data, model, be_safe=True, tol=0.2, debug=False):
  """
  Improve the wavelength solution by a constant shift
  if be_safe: it will not allow the solution to change by more than tol
  tol: the largest allowable shift (in nm), if be_safe == True
  data, model: xypoint instances of the data and model, respectively
    The model size MUST be >= the data size!
  """
  correction = data.size() + (model.size() - np.searchsorted(model.x, data.x[-1]))
  #correction = data.y.size + float(np.searchsorted(model.x, data.x[0]))/2.0 - 1
  ycorr = np.correlate(data.y/data.cont-1.0, model.y/model.cont-1.0, mode="full")
  xcorr = np.arange(ycorr.size)
  lags = xcorr - correction
  lags = xcorr - xcorr.size / 2
  distancePerLag = (data.x[-1] - data.x[0])/(float(data.x.size) - 1.0)
  offsets = -lags*distancePerLag
  offsets = offsets[::-1]
  ycorr = ycorr[::-1]

  if be_safe:
    left = np.searchsorted(offsets, -tol)
    right = np.searchsorted(offsets, tol)
  else:
    left, right = 0, ycorr.size

  maxindex = ycorr[left:right].argmax() + left

  if debug:
    return offsets[maxindex], DataStructures.xypoint(x=offsets, y=ycorr)
  else:
    return offsets[maxindex]


def Continuum(x, y, fitorder=3, lowreject=2, highreject=4, numiter=10000, function="poly"):
  """
  This function fits the continuum spectrum by iteratively removing
  points too far from the mean. The defaults work well for removing telluric lines in
  spectra of reasonable S/N ratio.

  x/y are assumed to by np arrays holding the spectrum to be fit.

  Note: Only the 'poly' function has been tested! Change to spline with extreme caution!
  """
  done = False
  x2 = np.copy(x)
  y2 = np.copy(y)
  iteration = 0
  while not done and iteration < numiter:
    iteration += 1
    done = True
    if function == "poly":
      fit = np.poly1d(np.polyfit(x2 - x2.mean(), y2, fitorder))
    elif function == "spline":
      fit = smoother(x2, y2, s=fitorder)
    residuals = y2 - fit(x2 - x2.mean())
    mean = np.mean(residuals)
    std = np.std(residuals)
    badpoints = np.where(np.logical_or((residuals - mean) < -lowreject*std, residuals - mean > highreject*std))[0]
    if badpoints.size > 0 and x2.size - badpoints.size > 5*fitorder:
      done = False
      x2 = np.delete(x2, badpoints)
      y2 = np.delete(y2, badpoints)
  return fit(x - x2.mean())







def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
    Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1]/m.sum(), y, mode='valid')


def Iterative_SV(y, window_size, order, lowreject=3, highreject=3, numiters=100, expand=0, deriv=0, rate=1):
  """
  Iterative version of the savitzky-golay smoothing function.
  See the documentation for savitzky_golay for details
  """
  done = False
  iteration = 0
  while not done and iteration < numiters:
    iteration += 1
    done = True
    smoothed = savitzky_golay(y, window_size, order, deriv, rate)

    reduced = y/smoothed
    sigma = np.std(reduced)
    mean = np.mean(reduced)
    badindices = np.where(np.logical_or((reduced - mean)/sigma < -lowreject, (reduced - mean)/sigma > highreject))[0]

    # Now, expand the outliers by 'expand' pixels on either
    exclude  = []
    for outlier in badindices:
      for i in range(max(0, outlier - expand), min(outlier+expand+1, reduced.size)):
        exclude.append(i)

    #Remove duplicates from 'exclude'
    badindices = []
    [badindices.append(i) for i in exclude if not i in badindices]
    badindices = np.array(badindices)
    if badindices.size > 0:
      done = False
      y[badindices] = smoothed[badindices]
  return smoothed





def FindLines(spectrum, tol=0.99, linespacing = 0.01, debug=False):
  """
  Function to find the spectral lines, given a model spectrum
  spectrum:        An xypoint instance with the model (Must be linearly spaced!)
  tol:             The line strength needed to count the line
                      (0 is a strong line, 1 is weak)
  linespacing:     The minimum spacing (nm) between two consecutive lines.
                      FindLines will choose the strongest line if there are
                      several too close.
  """
  #First, convert the inputs into inputs for scipy's argrelmin
  xspacing = float(max(spectrum.x) - min(spectrum.x))/float(spectrum.size())
  N = int( linespacing / xspacing + 0.5)
  lines = list(argrelmin(spectrum.y / spectrum.cont, order=N)[0])

  #Check for lines that are too weak.
  for i in range(len(lines)-1, -1, -1):
    idx = lines[i]
    xval = spectrum.x[idx]
    yval = spectrum.y[idx] / spectrum.cont[idx]
    if yval > tol:
      lines.pop(i)
    elif debug:
      plt.plot([xval, xval], [yval-0.01, yval-0.03], 'r-')

  if debug:
      plt.plot(spectrum.x, spectrum.y / spectrum.cont, 'k-')
      plt.title("Lines found in FittingUtilities.FindLines")
      plt.xlabel("Wavelength (nm)")
      plt.ylabel("Flux")
      plt.show()
  return np.array(lines)






def RebinData(data, xgrid, synphot=True):
  """
  This function rebins (x,y) data onto the grid given by the array xgrid
    It is designed to rebin to a courser wavelength grid, but can also
    interpolate to a finer grid.
  if synphot=True, it uses pySynphot for the rebinning, which conserves flux
    Otherwise, it just interpolates the data and continuum which is faster
    but could cause problems.
  """
  if synphot:
    #synphot chokes with astropy units, so remove them before proceeding
    data, xunits, yunits = data.strip_units()
    if isinstance(xgrid, u.Quantity):
        xgrid = xgrid.value

    #Do the actual rebinning
    newdata = DataStructures.xypoint(x=xgrid)
    newdata.y = rebin_spec(data.x, data.y, xgrid)
    newdata.cont = rebin_spec(data.x, data.cont, xgrid)
    newdata.err = rebin_spec(data.x, data.err, xgrid) #Unlikely to be correct!

    #pysynphot has edge effect issues on the first and last index.
    firstindex = np.argmin(np.abs(data.x - xgrid[0]))
    lastindex = np.argmin(np.abs(data.x - xgrid[-1]))
    newdata.y[0] = data.y[firstindex]
    newdata.y[-1] = data.y[lastindex]
    newdata.cont[0] = data.cont[firstindex]
    newdata.cont[-1] = data.cont[lastindex]
    newdata.err[0] = data.err[firstindex]
    newdata.err[-1] = data.err[lastindex]

    #Re-apply units
    newdata.x *= xunits
    newdata.y *= yunits
    newdata.err *= yunits
    newdata.cont *= yunits
    return newdata

  else:
    data_spacing = data.x[1] - data.x[0]
    grid_spacing = xgrid[1] - xgrid[0]
    newdata = DataStructures.xypoint(x=xgrid)
    if grid_spacing < 2.0*data_spacing:
      # Interpolate
      Model = spline(data.x, data.y)
      Continuum = spline(data.x, data.cont)
      newdata.y = Model(newdata.x)
      newdata.cont = Continuum(newdata.x)

    else:
      # Add up the pixels to rebin (actually re-binning).
      left = np.searchsorted(data.x, (3*xgrid[0]-xgrid[1])/2.0)
      for i in range(xgrid.size-1):
        right = np.searchsorted(data.x, (xgrid[i]+xgrid[i+1])/2.0)
        newdata.y[i] = np.mean(data.y[left:right])
        newdata.cont[i] = np.mean(data.cont[left:right])
        left = right
      right = np.searchsorted(data.x, (3*xgrid[-1]-xgrid[-2])/2.0)
      newdata.y[xgrid.size-1] = np.mean(data.y[left:right])

    return newdata



def rebin_spec(wave, specin, wavnew):
  """
  This takes in an x and y array, as well as the desired new x-array,
    and outputs the re-binned y array.
  """
  spec = ArraySourceSpectrum(wave=wave, flux=specin)
  f = np.ones(len(wave))
  filt = ArraySpectralElement(wave, f)
  obs = Observation(spec, filt, binset=wavnew, force='taper')

  return obs.binflux



def ReduceResolution(data,resolution, extend=True):
  """
  This function reduces the resolution by convolving with a gaussian kernel
  It assumes constant x-spacing!

  data is an xypoint instance containing the spectrum to smooth
  resolution: a float with the detector resolution (lam/dlam)

  WARNING! If the wavelength range is very large, it will give incorrect
    results because it uses a single Gaussian kernel for the whole range.
    This works just fine for small wavelength ranges such as in CRIRES,
    or with echelle spectra. If you have a large wavelength range, use
    ReduceResolution2!!
  """
  #Make the gaussian kernel
  centralwavelength = (data.x[0] + data.x[-1])/2.0
  xspacing = data.x[1] - data.x[0]
  FWHM = centralwavelength/resolution;
  sigma = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))
  left = 0
  right = np.searchsorted(data.x, 10*sigma)
  x = np.arange(0,10*sigma, xspacing)
  gaussian = np.exp(-(x-5*sigma)**2/(2*sigma**2))

  #Extend the xy axes to avoid edge-effects, if desired
  if extend:

    before = data.y[-gaussian.size//2+1:]
    after = data.y[:gaussian.size//2]
    extended = np.r_[before, data.y, after]

    first = data.x[0] - float(int(gaussian.size//2.0+0.5))*xspacing
    last = data.x[-1] + float(int(gaussian.size//2.0+0.5))*xspacing
    x2 = np.linspace(first, last, extended.size)

    conv_mode = "valid"

  else:
    extended = data.y.copy()
    x2 = data.x.copy()
    conv_mode = "same"

  newdata = data.copy()
  newdata.y = fftconvolve(extended, gaussian/gaussian.sum(), mode=conv_mode)

  return newdata


"""
  The following is np code to quickly convolve a spectrum
    for ReduceResolution2. DO NOT CALL THIS DIRECTLY!
"""
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[DTYPE_t, ndim=1] convolve(np.ndarray[DTYPE_t, ndim=1] x,
                                             np.ndarray[DTYPE_t, ndim=1] y,
                                             np.ndarray[DTYPE_t, ndim=1] output,

                                             double R,
                                             double nsig):
  cdef int i, n, start, end, length
  cdef double dx, sigma, total, conv, g, x0
  cdef np.ndarray[DTYPE_t, ndim=1] sig

  dx = x[1] - x[0]    #Assumes constant x-spacing!

  #Determine the edges
  sig = x/(2.0*R*sqrt(2.0*log(2.0)))
  n1 = np.searchsorted((x-x[0])/sig, nsig)
  n2 = np.searchsorted((x-x[x.size-1])/sig, -nsig)

  #Convolution outer loop
  for n in range(n1, n2):
    sigma = sig[n]
    length = int(sigma/dx * nsig + 0.5)
    x0 = x[n]
    total = 0.0
    conv = 0.0

    #Inner loop
    for i in range(-length, length+1):
      g = exp(-(x[n+i]-x0)**2 / (2.0*sigma**2))
      total += g
      conv += g*y[n+i]
    output[n] = conv/total
  return output




def ReduceResolution2(data,resolution, extend=True, nsig=5):
  """
  This is a more accurate version of ReduceResolution.
    It is also a bit slower, so don't use if you don't
    have to!
  """
  sig1 = data.x[0]/(2.0*resolution*np.sqrt(2.0*np.log(2.0)))
  sig2 = data.x[-1]/(2.0*resolution*np.sqrt(2.0*np.log(2.0)))
  dx = data.x[1] - data.x[0]
  n1 = int(sig1*(nsig+1)/dx + 0.5)
  n2 = int(sig2*(nsig+1)/dx + 0.5)

  if extend:
    #Extend array to try to remove edge effects (do so circularly)
    before = data.y[-n1:]
    after = data.y[:n2]
    #extended = np.append(np.append(before, data.y), after)
    extended = np.r_[before, data.y, after]

    first = data.x[0] - n1*dx
    last = data.x[-1] + n2*dx
    x2 = np.linspace(first, last, extended.size)
    convolved = np.ones(extended.size)
    convolved = convolve(x2, extended, convolved, resolution, nsig)
    convolved = convolved[n1:convolved.size-n2]

  else:
    extended = data.y.copy()
    x2 = data.x.copy()
    convolved = np.ones(extended.size)
    convolved = convolve(x2, extended, convolved, resolution, nsig)


  newdata = data.copy()
  newdata.y = convolved
  return newdata



def ReduceResolutionFFT(data, resolution, extend=True, loglinear=True, nsig=5):
    """
    Uses an fft to quickly reduce the resolution of a spectrum to the desired value
    :param data: An xypoint structure with the data to be convolved
    :param resolution: The value of the resolution dlam/lam
    :keyword extend: Flag to decide whether to extend the array before convolving to reduce edge effects
    :keyword loglinear: If true, the spectrum is already in log-linear spacing. If not, we have to do it here...
    :keyword nsig: The number of sigma to go on either side of the peak.
    """
    if not loglinear:
        x0 = np.log10(data.x[0])
        x1 = np.log10(data.x[-1])
        xgrid = np.logspace(x0, x1, data.size())
        data = RebinData(data, xgrid)

    # Make the broadening kernel.
    sigma = 1.0 / (2.0*resolution*np.sqrt(2*np.log(2.0)))
    dx = np.log(data.x[100]/data.x[99])
    d_logx = np.arange(0.0, nsig*sigma, dx)
    d_logx = np.r_[-d_logx[::-1][:-1], d_logx]
    B = np.exp(-0.5*(d_logx/sigma)**2)
    B /= B.sum()  #Normalize

    # Extend the data, if requested
    if extend:
        n = int(nsig*sigma/dx + 1)
        before = data.y[-n:]
        after = data.y[:n]
        y = np.r_[before, data.y, after]
    else:
        y = data.y

    # Do the convolution
    conv = fftconvolve(y, B, mode='same')

    # Shorten the final array if we extended it
    if extend:
        conv = conv[before.size:-after.size]

    newdata = data.copy()
    newdata.y = conv
    return newdata








def ReduceResolutionAndRebinData(data,resolution,xgrid):
  """
  Just a convenince fcn which combines two of the above functions
  """
  data = ReduceResolution(data,resolution)
  return RebinData(data,xgrid)
