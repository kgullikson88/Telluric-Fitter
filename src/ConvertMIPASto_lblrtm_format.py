"""
    This file will convert the MIPAS atmosphere profile (which exists in each
    run directory) into the format that LBLRTM expects for an atmosphere. You
    should never need to call this directly!

    
    
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


def Main(lines, jcharP="A", jcharT="A", jlong="L", jchar="", num_molecules=12):
    # ignore first 22 lines, which are just comments
    num_levels = int(lines[23].split()[0])
    data = {}  #set up a dictionary

    """
      Loop over the lines. The type of profile is denoted by lines
        that start with a '*'
    """
    for i in range(24, len(lines)):
        line = lines[i]
        if "*" in line and 'END' not in line:
            #We have found a new profile! Read it in.
            key = line.split('*')[1].split()[0]
            j = i + 1
            temp = []
            while "*" not in lines[j]:
                numbers = lines[j].split()
                for num in numbers:
                    temp.append(float(num))
                j = j + 1
            data[key] = temp


    #Now that we have read the atmosphere in, prepare the output lines
    outputlines = []
    for j in range(num_molecules):
        jchar = jchar + "A"
    for i in range(num_levels):
        altitude = "%.3e" % data['HGT'][i]
        pressure = "%.3e" % data["PRE"][i]
        temp = "%.3e" % data["TEM"][i]
        record35 = altitude.rjust(10) + pressure.rjust(10) + temp.rjust(
            10) + "     " + jcharP + jcharT + " " + jlong + " " + jchar
        record36 = ""
        for j in range(num_molecules):
            number = "%.8e" % data[mol[j]][i]
            record36 = record36 + number.rjust(15)
        outputlines.append(record35)
        outputlines.append(record36)

    return outputlines


# Some constant data:
mol = []
mol.append("H2O")
mol.append("CO2")
mol.append("O3")
mol.append("N2O")
mol.append("CO")
mol.append("CH4")
mol.append("O2")
mol.append("NO")
mol.append("SO2")
mol.append("NO2")
mol.append("NH3")
mol.append("HNO3")
mol.append("OH")
mol.append("HF")
mol.append("HCL")
mol.append("HBR")
mol.append("HI")

if __name__ == '__main__':
    Main("MIPAS_atmosphere_profile")
