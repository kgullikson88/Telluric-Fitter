"""
    This file provides the xypoint data structure, which I use to store spectra for TelFit
    

    This file is part of the TelFit program.

    TelFit is free software: you can redistribute it and/or modify
    it under the terms of the MIT license.

    TelFit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the MIT license
    along with TelFit.  If not, see <http://opensource.org/licenses/MIT>.

"""
from __future__ import print_function, division

import numpy as np
from astropy import units as u


class xypoint:
    def __init__(self, size=100, x=None, y=None, cont=None, err=None, array=None):
        """
          Generate an xypoint instance. Can do one of the following three things:
          1: Give a size only. This will generate an xypoint of a given size.
             You should fill it with your data in your own script
          2: Give x,y,cont,and err arrays, as np arrays. This will just
             copy the arrays into the xypoint structure
          3: Give a multidimensional array. For now at least, it must have
             shape = (size,4). That is the shape returned by np.loadtxt
             with unpack=False.
        """
        if array != None:
            self.x = array[:, 0]
            self.y = array[:, 1]
            self.cont = array[:, 2]
            self.err = array[:, 3]
            return

        if x is not None:
            size = x.size
        if y is not None:
            size = y.size
        if cont is not None:
            size = cont.size
        if err is not None:
            size = err.size

        if x is None:
            self.x = np.arange(size, dtype='float64')
        else:
            self.x = x.copy()
        if y is None:
            self.y = np.zeros(size)
        else:
            self.y = y.copy()
        if cont is None:
            self.cont = np.ones(size)
        else:
            self.cont = cont.copy()
        if err is None:
            self.err = np.ones(self.y.size) * 1e9
            self.err[self.y > 0] = np.sqrt(self.y[self.y > 0])
        else:
            self.err = err.copy()
            self.err[self.err <= 0] = 1e9  # Making sure we don't divide by zero

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
        np.savetxt(outfilename, np.transpose((self.x, self.y, self.cont, self.err)))

    def __getitem__(self, index):
        if isinstance(index, slice):
            start = max(0, index.start)
            stop = min(index.stop, self.size())
            if index.step == None:
                step = 1
            else:
                step = index.step
            x = np.array([self.x[i] for i in range(start, stop, step)])
            y = np.array([self.y[i] for i in range(start, stop, step)])
            cont = np.array([self.cont[i] for i in range(start, stop, step)])
            err = np.array([self.err[i] for i in range(start, stop, step)])
            return xypoint(x=x, y=y, cont=cont, err=err)
        elif isinstance(index, list) or isinstance(index, tuple) or isinstance(index, np.ndarray):
            x = self.x[index]
            y = self.y[index]
            cont = self.cont[index]
            err = self.err[index]
            return xypoint(x=x, y=y, cont=cont, err=err)
        else:
            return [self.x[index], self.y[index], self.cont[index], self.err[index]]

    def __len__(self):
        return self.size()

    def toarray(self, norm=False):
        """
        Turns the data structure into a multidimensional np array
        If norm == True, it will have shape self.size(), 2 and the y
          axis will be divided by the continuum axis
        Otherwise, it will have shape self.size(), 4
        """
        if norm:
            return np.array((self.x, self.y / self.cont)).transpose()
        else:
            return np.array((self.x, self.y, self.cont, self.err)).transpose()

    def strip_units(self):
        """
        Strips units from an xypoint structure.

        Returns:
          A copy of the xypoint with no units
          The x-units
          the y-units
        """
        xunits = self.x.unit if isinstance(self.x, u.Quantity) else 1.0
        yunits = self.y.unit if isinstance(self.y, u.Quantity) else 1.0
        x = self.x.value if isinstance(self.x, u.Quantity) else self.x
        y = self.y.value if isinstance(self.y, u.Quantity) else self.y
        err = self.err.value if isinstance(self.err, u.Quantity) else self.err
        cont = self.cont.value if isinstance(self.cont, u.Quantity) else self.cont
        return xypoint(x=x, y=y, cont=cont, err=err), xunits, yunits





def CombineXYpoints(xypts, snr=None, xspacing=None, numpoints=None, interp_order=3):
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
    from scipy.interpolate import InterpolatedUnivariateSpline

    if snr is None or type(snr) != list:
        snr = [1.0] * len(xypts)

    # Find the maximum range of the x data:
    first = np.min([o.x[0] for o in xypts])
    last = np.max([o.x[-1] for o in xypts])
    avg_spacing = np.mean([(o.x[-1] - o.x[0]) / float(o.size() - 1) for o in xypts])

    if xspacing is None and numpoints is None:
        xspacing = avg_spacing
    if numpoints is None:
        if xspacing is None:
            xspacing = avg_spacing
        numpoints = (last - first) / xspacing
    x = np.linspace(first, last, numpoints)

    full_array = xypoint(x=x, y=np.zeros(x.size))
    numvals = np.zeros(x.size, dtype=np.float)  #The number of arrays each x point is in
    normalization = 0.0
    for xypt in xypts:
        interpolator = InterpolatedUnivariateSpline(xypt.x, xypt.y / xypt.cont, k=interp_order)
        left = np.searchsorted(full_array.x, xypt.x[0])
        right = np.searchsorted(full_array.x, xypt.x[-1], side='right')
        if right < xypt.size():
            right += 1
        numvals[left:right] += 1.0
        full_array.y[left:right] += interpolator(full_array.x[left:right])
    full_array.y[numvals > 0] /= numvals[numvals > 0]
    return full_array
  
  
  













