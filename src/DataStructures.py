"""
    This file provides the xypoint data structure, which I use to store spectra for TelFit
    

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


import numpy
import sys


class xypoint:
  def __init__(self, size=100, x=None, y=None, cont=None, err=None):
    if x != None:
      size = x.size
    if y != None:
      size = y.size
    if cont != None:
      size = cont.size
    if err != None:
      size = err.size
      
    if x == None:
      self.x = numpy.arange(size)
    else:
      self.x = x.copy()
    if y == None:
      self.y = numpy.zeros(size)
    else:
      self.y = y.copy()
    if cont == None:
      self.cont = numpy.ones(size)
    else:
      self.cont = cont.copy()
    if err == None:
      self.err = numpy.ones(self.y.size)*1e9
      self.err[self.y > 0] = numpy.sqrt(self.y[self.y > 0])
    else:
      self.err = err.copy()
      self.err[self.err <=0] = 1e9    #Making sure we don't divide by zero
      
  def copy(self):
      copy = xypoint(self.x.size)
      copy.x = self.x.copy()
      copy.y = self.y.copy()
      copy.cont = self.cont.copy()
      copy.err = self.err.copy()
      return copy
  def size(self):
      return self.x.size
  def output(self, outfilename):
    numpy.savetxt(outfilename, numpy.transpose((self.x, self.y, self.cont, self.err)) )
  def __getitem__(self, index):
    if isinstance(index, slice):
      start = max(0, index.start)
      stop = min(index.stop, self.size())
      if index.step == None:
        step = 1
      else:
        step = index.step
      x = numpy.array([self.x[i] for i in range(start, stop, step)])
      y = numpy.array([self.y[i] for i in range(start, stop, step)])
      cont = numpy.array([self.cont[i] for i in range(start, stop, step)])
      err = numpy.array([self.err[i] for i in range(start, stop, step)])
      return xypoint(x=x, y=y, cont=cont, err=err)
    elif isinstance(index, list):
      x = self.x[index]
      y = self.y[index]
      cont = self.cont[index]
      err = self.err[index]
      return xypoint(x=x, y=y, cont=cont, err=err)
    else:
      return [self.x[index], self.y[index], self.cont[index], self.err[index]]
  def __len__(self):
    return self.size()


  
  
"""  
Function to combine a list of xypoints into a single
  xypoint. Useful for combining several orders/chips
  or for coadding spectra

Warning! This function is basically un-tested! It is NOT used
in TelFit!
  
  ***Optional keywords***
  snr: the spectra will be weighted by the signal-to-noise ratio
       before adding
  xspacing: the x-spacing in the final array
  numpoints: the number of points in the final array. If neither
             numpoints nor xspacing is given, the x-spacing in the
             final array will be determined by averaging the spacing
             in each of the xypoints.
  interp_order: the interpolation order. Default is cubic
"""
def CombineXYpoints(xypts, snr=None, xspacing=None, numpoints=None, interp_order=3):
  from scipy.interpolate import InterpolatedUnivariateSpline
  if snr == None or type(snr) != list:
    snr = [1.0]*len(xypts)
    
  #Find the maximum range of the x data:
  first = 1e30
  last = -1
  xspacing2 = 0.0
  for xypt in xypts:
    if xypt.x[0] < first:
      first = xypt.x[0]
    if xypt.x[-1] > last:
      last = xypt.x[-1]
    xspacing2 += (xypt.x[-1] - xypt.x[0]) / float(xypt.size() - 1)
  
  if xspacing == None and numpoints == None:
    xspacing = xspacing2 / float(len(xypts))
  if numpoints == None:
    if xspacing == None:
      xspacing = xspacing2 / float(len(xypts))
    numpoints = (last - first)/xspacing
  x = numpy.linspace(first, last, numpoints)
  
  full_array = xypoint(x=x, y=numpy.ones(x.size))
  numvals = numpy.ones(x.size)  #The number of arrays each x point is in
  normalization = 0.0
  for xypt in xypts:
    interpolator = InterpolatedUnivariateSpline(xypt.x, xypt.y/xypt.cont, k=interp_order)
    left = numpy.searchsorted(full_array.x, xypt.x[0])
    right = numpy.searchsorted(full_array.x, xypt.x[-1], side='right')
    if right < xypt.size():
      right += 1
    numvals[left:right] += 1.0
    full_array.y[left:right] += interpolator(full_array.x[left:right])
    
  full_array.y[numvals > 0] /= numvals[numvals > 0]  
  return full_array
  
  
  













