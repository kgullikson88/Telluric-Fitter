"""

    This file provides the MakeModel class. It is what directly interfaces
      with LBLRTM to make the telluric model. You can call this function 
      from the bash shell using the command 'python MakeModel.py' to generate
      a model transmission spectrum. The input settings for the model can be
      adjusted at the bottom of this file (after the line that reads
      'if __name__ == "__main__":')


    This file is part of the TelFit program.

    TelFit is free software: you can redistribute it and/or modify
    it under the terms of the MIT license.

    TelFit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the MIT license
    along with TelFit.  If not, see <http://opensource.org/licenses/MIT>.


"""

import sys
import subprocess
import os
from collections import defaultdict
import struct
import warnings
import time
import FittingUtilities
import copy

import scipy.interpolate
import lockfile
import numpy as np
from astropy import units

import DataStructures
import MakeTape5




# Dectionary giving the number to molecule name for LBLRTM
MoleculeNumbers = {1: "H2O",
                   2: "CO2",
                   3: "O3",
                   4: "N2O",
                   5: "CO",
                   6: "CH4",
                   7: "O2",
                   8: "NO",
                   9: "SO2",
                   10: "NO2",
                   11: "NH3",
                   12: "HNO3",
                   13: "OH",
                   14: "HF",
                   15: "HCl",
                   16: "HBr",
                   17: "HI",
                   18: "ClO",
                   19: "OCS",
                   20: "H2CO",
                   21: "HOCl",
                   22: "N2",
                   23: "HCN",
                   24: "CH3Cl",
                   25: "H2O2",
                   26: "C2H2",
                   27: "C2H6",
                   28: "PH3",
                   29: "COF2",
                   30: "SF6",
                   31: "H2S",
                   32: "HCOOH",
                   33: "HO2",
                   34: "O",
                   35: "ClONO2",
                   36: "NO+",
                   37: "HOBr",
                   38: "C2H4",
                   39: "CH3OH"}

"""
This is the main code to generate a telluric absorption spectrum.
The pressure, temperature, etc... can be adjusted all the way 
on the bottom of this file.
"""


class Modeler:
    def __init__(self, debug=False,
                 NumRunDirs=4,
                 TelluricModelingDirRoot=os.environ["TELLURICMODELING"],
                 nmolecules=12):

        Atmosphere = defaultdict(list)
        indices = {}
        self.debug = debug
        if not TelluricModelingDirRoot.endswith("/"):
            TelluricModelingDirRoot = TelluricModelingDirRoot + "/"
        self.TelluricModelingDirRoot = TelluricModelingDirRoot

        #Determine working directories
        self.FindWorkingDirectory()
        TelluricModelingDir = self.TelluricModelingDir
        ModelDir = self.ModelDir

        #Read in MIPAS atmosphere profile for upper atmosphere layers
        if debug:
            print "Generating new atmosphere profile"
        filename = TelluricModelingDir + "MIPAS_atmosphere_profile"
        infile = open(filename)
        lines = infile.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith("*") and "END" not in line:
                if line.find("HGT") > 0 and line.find("[") > 0:
                    numlevels = int(lines[i - 1].split()[0])
                    indices["Z"] = i
                elif line.find("PRE") > 0 and line.find("[") > 0:
                    indices["P"] = i
                elif line.find("TEM") > 0 and line.find("[") > 0:
                    indices["T"] = i
                else:
                    molecule = line.split("*")[-1].split()[0]
                    indices[molecule] = i
        infile.close()

        #Determine the number of lines that follow each header
        levelsperline = 5.0
        linespersection = int(numlevels / levelsperline + 0.9)

        #Fill atmosphere structure
        layers = []
        for j in range(indices["Z"] + 1, indices["Z"] + 1 + linespersection):
            line = lines[j]
            levels = line.split()
            [layers.append(float(level)) for level in levels]
        #Pressure
        for j in range(linespersection):
            line = lines[j + indices["P"] + 1]
            levels = line.split()
            for i, level in enumerate(levels):
                Atmosphere[layers[int(j * levelsperline + i)]].append(float(level))
        #Temperature
        for j in range(linespersection):
            line = lines[j + indices["T"] + 1]
            levels = line.split()
            for i, level in enumerate(levels):
                Atmosphere[layers[int(j * levelsperline + i)]].append(float(level))
                Atmosphere[layers[int(j * levelsperline + i)]].append([])
        #Abundances
        for k in range(1, nmolecules + 1):
            for j in range(linespersection):
                line = lines[j + indices[MoleculeNumbers[k]] + 1]
                levels = line.split()
                for i, level in enumerate(levels):
                    Atmosphere[layers[int(j * levelsperline + i)]][2].append(float(level))


        #Save some local variables as class variables, so that other functions can use them
        self.Atmosphere = Atmosphere
        self.layers = layers
        self.nmolecules = nmolecules

        self.Cleanup()
        return


    def EditProfile(self, profilename, profile_height, profile_value):
        """
        This function will take a np array as a profile, and stitch it into the
          MIPAS atmosphere profile read in in __init__

        -profilename:  A string with the name of the profile to edit.
                       Should be either 'pressure', 'temperature', or
                       one of the molecules given in the MakeModel.MoleculeNumbers
                       dictionary
        -profile_height:  A np array with the height in the atmosphere (in km)
        -profile_value:   A np array with the value of the profile parameter at
                          each height given in profile_height.
        """
        #Translate the (string) name of the profile to a number
        profilenum = -1
        if profilename.lower() == "pressure":
            profilenum = 0
        elif profilename.lower() == "temperature":
            profilenum = 1
        else:
            profilenum = 2
            molnum = -1
            for n in MoleculeNumbers:
                if MoleculeNumbers[n] == profilename:
                    molnum = n - 1
            if molnum < 0:
                print "Error! Profilename given in Modeler.EditProfile is invalid!"
                print "You gave: ", profilename
                print "Valid options are: "
                for n in MoleculeNumbers:
                    print MoleculeNumbers[n]
                raise ValueError
        if profilenum < 0:
            print "Error! Profilename given in Modeler.EditProfile is invalid!"
            print "You gave: ", profilename
            print "Valid options are either a molecule name, 'pressure', or 'temperature'"
            raise ValueError

        #Okay, they gave a valid profile name. Lets grab some class variables
        Atmosphere = self.Atmosphere
        layers = np.array(self.layers)

        #Create a list of the relevant molecule profile
        mipas = []
        for layer in layers:
            if profilenum < 2:
                mipas.append(Atmosphere[layer][profilenum])
            else:
                mipas.append(Atmosphere[layer][profilenum][molnum])

        #Check if the profile value is near 0 for the last few entries.
        #This has happened with GDAS profiles that report the dew point
        #near absolute zero when it doesn't get a measurement
        while abs(profile_value[-1]) < 1:
            profile_height = profile_height[:-1]
            profile_value = profile_value[:-1]

        #We now need to adjust the MIPAS profile to avoid discontinuities
        #Do so differently for layers below the profile given as compared to
        #  those above it
        profile_fcn = scipy.interpolate.interp1d(profile_height, profile_value, kind='linear')
        left = np.searchsorted(layers, profile_height[0])
        right = np.searchsorted(layers, profile_height[-1]) - 1
        newprofile = list(mipas)
        newprofile[:left] -= (mipas[left] - profile_fcn(layers[left])) * np.exp(-(layers[left] - layers[:left]))
        newprofile[right:] -= (mipas[right] - profile_fcn(layers[right])) * np.exp(-(layers[right:] - layers[right]))
        newprofile[left:right] = profile_fcn(layers[left:right])


        #Now, put the newprofile array into Atmosphere
        for i, layer in enumerate(layers):
            if profilenum < 2:
                Atmosphere[layer][profilenum] = newprofile[i]
            else:
                Atmosphere[layer][profilenum][molnum] = newprofile[i]
        self.Atmosphere = Atmosphere


    def Cleanup(self):
        """
          Release the lock on the directory. This can be called on its own, but
          should never need to be.
        """
        lock = self.lock
        #Unlock directory
        try:
            lock.release()
        except lockfile.NotLocked:
            warnings.warn("The model directory was somehow unlocked prematurely!")
        return


    def FindWorkingDirectory(self):
        """
        Find a run directory to work in. This is necessary so that you
          can run several instances of MakeModel (or TelFit) at once.
          Should not need to be called on directly by the user.
        """
        #Determine output filename
        TelluricModelingDirRoot = self.TelluricModelingDirRoot
        found = False
        possible_rundirs = [d for d in os.listdir(self.TelluricModelingDirRoot) if
                            d.startswith('rundir') and "." not in d]
        while not found:
            for test in possible_rundirs:
                test = "%s%s" % (TelluricModelingDirRoot, test)
                if not test.endswith("/"):
                    test = test + "/"

                lock = lockfile.FileLock(test)
                if not lock.is_locked():
                    TelluricModelingDir = test
                    ModelDir = "%sOutputModels/" % TelluricModelingDir
                    lock.acquire()
                    found = True
                    break
            if not found:
                print "Un-locked directory not found! Waiting 10 seconds..."
                time.sleep(10)
        if self.debug:
            print "Telluric Modeling Directory: %s" % TelluricModelingDir
            print "Model Directory: %s" % ModelDir

        self.TelluricModelingDir = TelluricModelingDir
        self.ModelDir = ModelDir
        self.lock = lock


    def MakeModel(self, pressure=795.0, temperature=283.0, lowfreq=4000, highfreq=4600, angle=45.0, humidity=50.0,
                  co2=368.5, o3=3.9e-2, n2o=0.32, co=0.14, ch4=1.8, o2=2.1e5, no=1.1e-19, so2=1e-4, no2=1e-4, nh3=1e-4,
                  hno3=5.6e-4, lat=30.6, alt=2.1, wavegrid=None, resolution=None, save=False, libfile=None):
        """
        Here is the important function! All of the variables have default values,
          which you will want to override for any realistic use.
        Arguments/Units:
          Pressure:               Pressure at telescope altitude (hPa)
          Temperature:            Temperature at telescope altitude (Kelvin)
          lowfreq, highfreq:      low and high wavenumber (cm^-1)
          angle:                  The zenith angle of the telescope (degrees)
          humidity:               percent humidity
          co2 - hno3 abundances:  ppmv concentration
          lat:                    The latitude of the observatory (degrees)
          alt:                    The altitude of the observatory above sea level (km)
          wavegrid:               If given, the model will be resampled to this grid.
                                  Should be a np array
          resolution:             If given, it will reduce the resolution by convolving
                                  with a gaussian of appropriate width. Should be a float
                                  with R=lam/dlam
          save:                   If true, the generated model is saved. The filename will be
                                  printed to the screen.
          libfile:                Useful if generating a telluric library. The filename of the
                                  saved file will be written to this filename. Should be a string
                                  variable. Ignored if save==False

        """

        self.FindWorkingDirectory()

        #Make the class variables local
        TelluricModelingDir = self.TelluricModelingDir
        debug = self.debug
        lock = self.lock
        layers = np.array(self.layers)
        ModelDir = self.ModelDir

        #Make a deep copy of atmosphere so that I don't repeatedly modify it
        Atmosphere = copy.deepcopy(self.Atmosphere)


        #Convert from relative humidity to concentration (ppm)
        Psat = VaporPressure(temperature)
        Pw = Psat * humidity / 100.0
        h2o = Pw / (pressure - Pw) * 1e6


        #Start by scaling the abundances from those at 'alt' km
        #  (linearly interpolate)
        keys = sorted(Atmosphere.keys())
        lower = max(0, np.searchsorted(keys, alt) - 1)
        upper = min(lower + 1, len(keys) - 1)
        if lower == upper:
            raise ZeroDivisionError(
                "Observatory altitude of %g results in the surrounding layers being the same!" % alt)
        scale_values = list(Atmosphere[lower])
        scale_values[2] = list(Atmosphere[lower][2])
        scale_values[0] = (Atmosphere[upper][0] - Atmosphere[lower][0]) / (keys[upper] - keys[lower]) * (
        alt - keys[lower]) + Atmosphere[lower][0]
        scale_values[1] = (Atmosphere[upper][1] - Atmosphere[lower][1]) / (keys[upper] - keys[lower]) * (
        alt - keys[lower]) + Atmosphere[lower][1]
        for mol in range(len(scale_values[2])):
            scale_values[2][mol] = (Atmosphere[upper][2][mol] - Atmosphere[lower][2][mol]) / (
            keys[upper] - keys[lower]) * (alt - keys[lower]) + Atmosphere[lower][2][mol]

        #Do the actual scaling
        pressure_scalefactor = (scale_values[0] - pressure) * np.exp(-(layers - alt) ** 2 / (2.0 * 10.0 ** 2))
        temperature_scalefactor = (scale_values[1] - temperature) * np.exp(-(layers - alt) ** 2 / (2.0 * 10.0 ** 2))
        for i, layer in enumerate(layers):
            Atmosphere[layer][0] -= pressure_scalefactor[i]
            Atmosphere[layer][1] -= temperature_scalefactor[i]
            Atmosphere[layer][2][0] *= h2o / scale_values[2][0]
            Atmosphere[layer][2][1] *= co2 / scale_values[2][1]
            Atmosphere[layer][2][2] *= o3 / scale_values[2][2]
            Atmosphere[layer][2][3] *= n2o / scale_values[2][3]
            Atmosphere[layer][2][4] *= co / scale_values[2][4]
            Atmosphere[layer][2][5] *= ch4 / scale_values[2][5]
            Atmosphere[layer][2][6] *= o2 / scale_values[2][6]
            Atmosphere[layer][2][7] *= no / scale_values[2][7]
            Atmosphere[layer][2][8] *= so2 / scale_values[2][8]
            Atmosphere[layer][2][9] *= no2 / scale_values[2][9]
            Atmosphere[layer][2][10] *= nh3 / scale_values[2][10]
            Atmosphere[layer][2][11] *= hno3 / scale_values[2][11]


        #Now, Read in the ParameterFile and edit the necessary parameters
        parameters = MakeTape5.ReadParFile(parameterfile=TelluricModelingDir + "ParameterFile")
        parameters[48] = "%.1f" % lat
        parameters[49] = "%.1f" % alt
        parameters[51] = "%.5f" % angle
        parameters[17] = lowfreq
        freq, transmission = np.array([]), np.array([])

        #Need to run lblrtm several times if the wavelength range is too large.
        maxdiff = 1999.9
        if (highfreq - lowfreq > maxdiff):
            while lowfreq + maxdiff <= highfreq:
                parameters[18] = lowfreq + maxdiff

                MakeTape5.WriteTape5(parameters, output=TelluricModelingDir + "TAPE5", atmosphere=Atmosphere)

                #Run lblrtm
                cmd = "cd " + TelluricModelingDir + ";sh runlblrtm_v3.sh"
                try:
                    command = subprocess.check_call(cmd, shell=True)
                except subprocess.CalledProcessError:
                    print "Error: Command '%s' failed in directory %s" % (cmd, TelluricModelingDir)
                    sys.exit()

                #Read in TAPE12, which is the output of LBLRTM
                freq, transmission = self.ReadTAPE12(TelluricModelingDir, appendto=(freq, transmission))
                lowfreq = lowfreq + 2000.00001
                parameters[17] = lowfreq

        parameters[18] = highfreq
        MakeTape5.WriteTape5(parameters, output=TelluricModelingDir + "TAPE5", atmosphere=Atmosphere)

        #Run lblrtm for the last time
        cmd = "cd " + TelluricModelingDir + ";sh runlblrtm_v3.sh"
        try:
            command = subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "Error: Command '%s' failed in directory %s" % (cmd, TelluricModelingDir)
            sys.exit()

        #Read in TAPE12, which is the output of LBLRTM
        freq, transmission = self.ReadTAPE12(TelluricModelingDir, appendto=(freq, transmission))

        #Convert from frequency to wavelength units
        wavelength = units.cm.to(units.nm) / freq

        #Correct for index of refraction of air (use IAU standard conversion from
        #  Morton, D. C. 1991, ApJS, 77, 119
        wave_A = wavelength * units.nm.to(units.angstrom)  #Wavelength in angstroms
        n = 1.0 + 2.735182e-4 + 131.4182 / wave_A ** 2 + 2.76249e8 / wave_A ** 4
        wavelength /= n

        if save:
            #Output filename
            model_name = ModelDir + "transmission" + "-%.2f" % pressure + "-%.2f" % temperature + "-%.1f" % humidity + "-%.1f" % angle + "-%.2f" % (
            co2) + "-%.2f" % (o3 * 100) + "-%.2f" % ch4 + "-%.2f" % (co * 10)
            print "All done! Output Transmission spectrum is located in the file below:"
            print model_name
            np.savetxt(model_name, np.transpose((wavelength[::-1], transmission[::-1])), fmt="%.8g")
            if libfile != None:
                infile = open(libfile, "a")
                infile.write(model_name + "\n")
                infile.close()

        self.Cleanup()  #Un-lock the working directory

        if wavegrid != None:
            model = DataStructures.xypoint(x=wavelength[::-1], y=transmission[::-1])
            return FittingUtilities.RebinData(model, wavegrid)

        return DataStructures.xypoint(x=wavelength[::-1], y=transmission[::-1])


    def ReadTAPE12(self, directory, filename="TAPE12_ex", appendto=None):
        """
        Here is a function to read in the binary output of lblrtm, and convert
          it into arrays of frequency and transmission.
        Warning! Some values are hard-coded in for single precision calculations.
          You MUST compile lblrtm as single precision or this won't work!
        Not meant to be called directly by the user.
        """
        debug = self.debug
        if not directory.endswith("/"):
            directory = directory + "/"
        infile = open("%s%s" % (directory, filename), 'rb')
        content = infile.read()
        infile.close()

        offset = 1068
        size = struct.calcsize('=ddfl')
        pv1, pv2, pdv, numpoints = struct.unpack('=ddfl', content[offset:offset + size])
        v1 = pv1
        v2 = pv2
        dv = pdv
        if debug:
            print 'info: ', pv1, pv2, pdv, numpoints
        npts = numpoints
        spectrum = []
        while numpoints > 0:
            offset += size + struct.calcsize("=4f")
            size = struct.calcsize("=%if" % numpoints)
            temp1 = struct.unpack("=%if" % numpoints, content[offset:offset + size])
            offset += size
            temp2 = struct.unpack("=%if" % numpoints, content[offset:offset + size])
            npts += numpoints
            junk = [spectrum.append(temp2[i]) for i in range(numpoints)]

            offset += size + 8
            size = struct.calcsize('=ddfl')
            if len(content) > offset + size:
                pv1, pv2, pdv, numpoints = struct.unpack('=ddfl', content[offset:offset + size])
                v2 = pv2
            else:
                break

        v = np.arange(v1, v2, dv)
        spectrum = np.array(spectrum)
        if v.size < spectrum.size:
            v = np.r_[v, v2 + dv]
        if debug:
            print "v, spec size: ", v.size, spectrum.size

        if appendto != None and appendto[0].size > 0:
            old_v, old_spectrum = appendto[0], appendto[1]
            #Check for overlap (there shouldn't be any)
            last_v = old_v[-1]
            firstindex = np.searchsorted(v, last_v)
            v = np.r_[old_v, v[firstindex:]]
            spectrum = np.r_[old_spectrum, spectrum[firstindex:]]

        return v, spectrum


def VaporPressure(T):
    """
      This function uses equations and constants from
      http://www.vaisala.com/Vaisala%20Documents/Application%20notes/Humidity_Conversion_Formulas_B210973EN-F.pdf
      to determine the vapor pressure at the given temperature

      T must be a float with the temperature (or dew point) in Kelvin
    """
    #Constants
    c1, c2, c3, c4, c5, c6 = -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.8022502
    a0, a1 = -13.928169, 34.707823
    if T > 273.15:
        theta = 1.0 - T / 647.096
        Pw = 2.2064e5 * np.exp(1.0 / (1.0 - theta) * ( c1 * theta +
                                                       c2 * theta ** 1.5 +
                                                       c3 * theta ** 3 +
                                                       c4 * theta ** 3.5 +
                                                       c5 * theta ** 4 +
                                                       c6 * theta ** 7.5 ))
    elif T > 173.15:
        theta = T / 273.16
        Pw = 6.11657 * np.exp(a0 * (1.0 - theta ** -1.5) +
                              a1 * (1.0 - theta ** -1.25))
    else:
        Pw = 0.0

    return Pw


if __name__ == "__main__":
    pressure = 796.22906
    temperature = 270.40
    humidity = 40.0
    angle = 40.8
    co2 = 368.5
    o3 = 0.039
    ch4 = 4.0
    co = 0.15
    o2 = 2.2e5

    lowwave = 620
    highwave = 650
    lowfreq = 1e7 / highwave
    highfreq = 1e7 / lowwave
    modeler = Modeler(debug=False)
    modeler.MakeModel(pressure=pressure, temperature=temperature, humidity=humidity, lowfreq=lowfreq, highfreq=highfreq,
                      angle=angle, o2=o2, alt=2.1, save=True)
          

