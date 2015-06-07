import sys
import os
from scipy.interpolate import InterpolatedUnivariateSpline as interp

import numpy
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy import units, constants
from astropy import units, constants
from telfit import TelluricFitter, Modeler, DataStructures, FittingUtilities
import MakeModel


# Define regions to ignore in the fit
badregions = [[588.98, 589.037],  # Na D line 1
              [589.567, 589.632],  #Na D line 2
              [627.4, 629.0],  # Very deep O2 band
              [686.4, 690.7]]  # Very deep O2 band


def FindOrderNums(orders, wavelengths):
    """
      Given a list of xypoint orders and
      another list of wavelengths, this
      finds the order numbers with the
      requested wavelengths
    """
    nums = []
    for wave in wavelengths:
        for i, order in enumerate(orders):
            if order.x[0] < wave and order.x[-1] > wave:
                nums.append(i)
                break
    return nums


def GetProfile(filename):
    """
      This function reads in the GDAS profile, sorts the data
        by the atmosphere height, and converts the values
        to the appropriate units
      Returns: pressure, height, temperature, and relative humidity
               at several layers, interpolated to the observation time
    """
    # Get the interpolated values
    P, Z, T, D = numpy.loadtxt(filename, usecols=(0, 1, 2, 3), unpack=True)

    #Sort by height
    sorter = numpy.argsort(Z)
    Z = Z[sorter]
    P = P[sorter]
    T = T[sorter]
    D = D[sorter]

    #Convert dew point temperature to ppmv
    Pw = numpy.zeros(D.size)
    for i, dewpoint in enumerate(D):
        Pw[i] = MakeModel.VaporPressure(dewpoint + 273.15)
    h2o = Pw / (P - Pw) * 1e6

    #Convert height to km, and temperature to kelvin
    Z /= 1000.0
    T += 273.15
    return Z, P, T, h2o


def OutputFitsFileExtensions(column_dicts, template, outfilename, mode="append", headers_info=[]):
    """
    Function to output a fits file
    column_dict is a dictionary where the key is the name of the column
       and the value is a numpy array with the data. Example of a column
       would be the wavelength or flux at each pixel
    template is the filename of the template fits file. The header will
       be taken from this file and used as the main header
    mode determines how the outputted file is made. Append will just add
       a fits extension to the existing file (and then save it as outfilename)
       "new" mode will create a new fits file.
       header_info takes a list of lists. Each sub-list should have size 2 where the first element is the name of the new keyword, and the second element is the corresponding value. A 3rd element may be added as a comment
    """

    # Get header from template. Use this in the new file
    if mode == "new":
        header = pyfits.getheader(template)

    if not isinstance(column_dicts, list):
        column_dicts = [column_dicts, ]
    if len(headers_info) < len(column_dicts):
        for i in range(len(column_dicts) - len(headers_info)):
            headers_info.append([])

    # Generate the hdu list
    if mode == "append":
        hdulist = pyfits.open(template)
    elif mode == "new":
        header = pyfits.getheader(template)
        pri_hdu = pyfits.PrimaryHDU(header=header)
        hdulist = pyfits.HDUList([pri_hdu, ])

    # Make a fits binary table with the column data
    for i in range(len(column_dicts)):
        column_dict = column_dicts[i]
        header_info = headers_info[i]
        columns = []
        for key in column_dict.keys():
            columns.append(pyfits.Column(name=key, format="D", array=column_dict[key]))
        cols = pyfits.ColDefs(columns)
        tablehdu = pyfits.new_table(cols)

        #Add keywords to extension header
        num_keywords = len(header_info)
        header = tablehdu.header
        for i in range(num_keywords):
            info = header_info[i]
            if len(info) > 2:
                header.set(info[0], info[1], info[2])
            elif len(info) == 2:
                header.set(info[0], info[1])

        hdulist.append(tablehdu)

    #Output to file
    hdulist.writeto(outfilename, clobber=True, output_verify='ignore')
    hdulist.close()


if __name__ == "__main__":
    # Initialize fitter
    fitter = TelluricFitter()
    fitter.SetObservatory("CTIO")

    # Read in the fits file using astropy
    # This is not guaranteed to work!
    fname = "fitsfile.fits"
    outfilename = "Corrected.fits"
    orders = []
    hdulist = pyfits.open(fname)
    header = hdulist[0].header
    for hdu in hdulist[1:]:
        data = hdu.data
        order = DataStructures.xypoint(x=data['wavelength'],
                                       y=data['flux'],
                                       cont=data['continuum'],
                                       err=data['error'])
        orders.append(order)
    hdulist.close()

    # Read in the appropriate parameters:
    angle = float(header['ZD'])
    resolution = 80000.0
    humidity = header['OUTHUM']
    pressure = header['OUTPRESS']
    temperature = header['OUTTEMP'] + 273.15  # Convert to Kelvin

    # Adjust the atmosphere profile using GDAS data
    height, Pres, Temp, h2o = GetProfile("GDAS_Atmosphere.dat")
    fitter.EditAtmosphereProfile("Temperature", height, Temp)
    fitter.EditAtmosphereProfile("Pressure", height, Pres)
    fitter.EditAtmosphereProfile("H2O", height, h2o)


    # Adjust fitter values
    fitter.FitVariable({"h2o": humidity})  # We will fit humidity
    fitter.AdjustValue({"angle": angle,  # Leave these constant
                        "pressure": pressure,
                        "temperature": temperature,
                        "resolution": resolution,
                        "o2": 2.12e5})
    fitter.SetBounds({"h2o": [1.0, 99.0],
                      "temperature": [temperature - 10, temperature + 10],
                      "o2": [5e4, 1e6],
                      "resolution": [70000, 90000]})

    #Ignore the interstellar sodium D lines and parts of the O2 bands
    fitter.IgnoreRegions(badregions)


    ### ---------------------------------------------
    ###              Main Algorithm - Humidity fit
    ### ---------------------------------------------

    # Determine the H2O abundance first by fitting several orders
    resolution = []
    h2o = []
    o2 = []
    waveshifts = []
    wave0 = []
    weights = []
    for i in FindOrderNums(orders, [595, 650, 717, 726]):
        print "\n***************************\nFitting order %i: " % (i)
        order = orders[i]

        # Set wavelength region for the fit
        fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                            "waveend": order.x[-1] + 20.0})

        # Determine the continuum (refined by the fitting algorithm)
        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)

        # Perform the fit
        primary, model, R = fitter.Fit(data=order.copy(),
                                       #Can give data here as well as separately like in the simpler fit example
                                       resolution_fit_mode="gauss",
                                       fit_source=True,  #Fit the broad stellar lines too
                                       return_resolution=True,  #Return the best fit resolution
                                       adjust_wave="model",  #Assume the data is well-calibrated
                                       wavelength_fit_order=3)
        # Save the best-fit values
        resolution.append(R)
        waveshifts.append(fitter.shift)
        wave0.append(fitter.data.x.mean())
        h2o.append(fitter.GetValue("h2o"))

        #Weight by the chi-squared value, but also weight strong lines more
        weights.append((1.0 - min(model.y)) / fitter.chisq_vals[-1])

    # Determine the average humidity (weighted)
    humidity = numpy.sum(numpy.array(h2o) * numpy.array(weights)) / numpy.sum(weights)
    fitter.AdjustValue({"h2o": humidity})  #This turns h2o off as a fitted variable!


    ### ---------------------------------------------
    ###              Main Algorithm - O2 fit
    ### ---------------------------------------------

    # Now, determine the O2 abundance
    fitter.FitVariable({"o2": 2.12e5})
    for i in FindOrderNums(orders, [630, 690]):
        order = orders[i]
        fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                            "waveend": order.x[-1] + 20.0})
        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
        primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
        primary, model, R = fitter.Fit(data=order.copy(),
                                       resolution_fit_mode="gauss",
                                       fit_source=True,
                                       return_resolution=True,
                                       adjust_wave="model",
                                       wavelength_fit_order=3)
        resolution.append(R)
        waveshifts.append(fitter.shift)
        wave0.append(fitter.data.x.mean())
        o2.append(fitter.GetValue("o2"))
        weights.append((1.0 - min(model.y)) / fitter.chisq_vals[-1])

    # Determine the average of the other parameter values
    chi2 = numpy.array(weights)
    o2 = numpy.array(o2)
    resolution = numpy.array(resolution)
    waveshifts = numpy.array(waveshifts)
    wave0 = numpy.array(wave0)
    velshifts = waveshifts / wave0 * constants.c.cgs.value * units.cm.to(units.km)
    vel = numpy.sum(velshifts * chi2) / numpy.sum(chi2)
    o2 = numpy.sum(o2 * chi2[-2:]) / numpy.sum(chi2[-2:])
    resolution = numpy.sum(resolution[:-2] * chi2[:-2]) / numpy.sum(chi2[:-2])



    ### ---------------------------------------------
    ###              Main Algorithm - Application
    ### ---------------------------------------------

    #Prepare for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Finally, apply these parameters to all orders in the data
    for i, order in enumerate(orders):
        print "\n\nGenerating model for order %i of %i\n" % (i, len(orders))
        fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                            "waveend": order.x[-1] + 20.0,
                            "o2": o2,
                            "h2o": humidity,
                            "resolution": resolution})
        fitter.resolution_fit_mode = "gauss"

        #Get the full fitpars array (kind of a hack...)
        fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]

        #Fit the continuum
        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)

        #Import the data (We are not calling 'Fit' anymore)
        fitter.ImportData(order)

        #Generate the model using the best-fit parameters
        primary, model = fitter.GenerateModel(fitpars,
                                              separate_source=True,
                                              return_resolution=False)

        # Get the data back from the fitter.
        # (It has a better continuum)
        data = fitter.data

        #If there are very few deep telluric lines, the wavelength calibration can be off
        #Set the wavelength calibration automatically
        if min(model.y) > 0.98:
            wave0 = order.x.mean()
            fitter.shift = vel / (constants.c.cgs.value * units.cm.to(units.km)) * wave0
            model = fitter.GenerateModel(fitpars, separate_source=False, nofit=True)
            model.x /= (1.0 + vel / (constants.c.cgs.value * units.cm.to(units.km)))
            model = FittingUtilities.RebinData(model, order.x)
            data = order.copy()
            data.cont = FittingUtilities.Continuum(data.x, data.y, fitorder=3, lowreject=2, highreject=5)

        # Set up data structures for OutputFitsFile
        columns = {"wavelength": data.x,
                   "flux": data.y,
                   "continuum": data.cont,
                   "error": data.err,
                   "model": model.y,
                   "primary": primary.y}

        header_info = []

        if i == 0:
            OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info, ], mode="new")
        else:
            OutputFitsFileExtensions(columns, outfilename, outfilename, headers_info=[header_info, ], mode="append")

        #Plot
        ax.plot(data.x, data.y / data.cont, 'k-')
        ax.plot(model.x, model.y, 'r-')
    plt.show()
      