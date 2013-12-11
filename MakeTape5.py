"""
    This file will read in the file 'ParameterFile', which is a human-readable
      version of the input for LBLRTM (called TAPE5). It will convert that into
      TAPE5 for use by LBLRTM. You should never need to call this directly!




    This file is part of the TelluricFitter program.

    TelluricFitter is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TelluricFitter is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TelluricFitter.  If not, see <http://www.gnu.org/licenses/>.



"""

import sys
from fortranformat import FortranRecordWriter as writer
from fortranformat import FortranRecordReader as reader
from numpy import pi, cos, tan, arctan


def Convert(parfile="ParameterFile", output="TAPE5"):
  #Read parameter file
  infile = open(parfile)
  lines = infile.readlines()
  infile.close()
  
  #open output file
  outfile = open(output, "w")
  
  #Record 1.1:
  description = lines[2][24:]
  outfile.write(writer(("A%i" %min(80, len(description)+1))).write(("$%s" %description, )))
  
  #Record 1.2
  ihirac = int(lines[4][24:])
  ilblf4 = int(lines[5][24:])
  icntnm = int(lines[6][24:])
  iaersl = int(lines[7][24:])
  iemit = int(lines[8][24:])
  iscan = int(lines[9][24:])
  ifilter = int(lines[10][24:])
  iplot = int(lines[11][24:])
  itest = int(lines[12][24:])
  iatm = int(lines[13][24:])
  imrg = int(lines[14][24:])
  ilas = int(lines[15][24:])
  iod = int(lines[16][24:])
  ixsec = int(lines[17][24:])
  mpts = int(lines[18][24:])
  npts = int(lines[19][24:])
  outfile.write(writer("(10(4x,i1),3x,i2, 3(4x,i1),2(1x,i4))").write((ihirac, ilblf4, icntnm, iaersl, iemit, iscan, ifilter, iplot, itest, iatm, imrg, ilas, iod, ixsec, mpts, npts)) + "\n")
  
  if icntnm == 6:
    print "icntnm = 6"
  if iemit == 2:
    print "iemit = 2"
    
  #Record 1.3
  v1 = float(lines[23][24:])
  v2 = float(lines[24][24:])
  sample = float(lines[25][24:])
  dvset = float(lines[26][24:])
  alfal0 = float(lines[27][24:])
  avmass = float(lines[28][24:])
  dptmin = float(lines[29][24:])
  dptflg = float(lines[30][24:])
  ilnflg = int(lines[31][24:])
  dvout = float(lines[32][24:])
  nmol_scal = int(lines[33][24:])
  outfile.write(writer("(8(F10.3),4x,I1,5x,F10.3,3x,I2)").write((v1, v2, sample, dvset, alfal0, avmass, dptmin, dptflg, ilnflg, dvout, nmol_scal)) + "\n")
 
  #Record 1.3.a and 1.3.b (TODO!)
  if nmol_scal > 0:
    print nmol_scal
    
  #Record 1.4
  if iemit == 1 or (iemit == 2 and iotflg == 2):
    tbound = float(lines[37][24:])
    sremis1 = float(lines[38][24:])
    sremis2 = float(lines[39][24:])
    sremis3 = float(lines[40][24:])
    srrefl1 = float(lines[41][24:])
    srrefl2 = float(lines[42][24:])
    srrefl3 = float(lines[43][24:])
    surf_refl = lines[44][24:].strip()
    outfile.write(writer("(7(F10.3),4x,A1)").write((tbound, sremis1, sremis2, sremis3, srrefl1, srrefl2, srrefl3, surf_refl)) + "\n")
    
  #Record 3.1
  if iatm == 1:
    model = int(lines[49][24:])
    itype = int(lines[50][24:])
    ibmax = int(lines[51][24:])
    zero = int(lines[52][24:])
    noprnt = int(lines[53][24:])
    nmol = int(lines[54][24:])
    ipunch = int(lines[55][24:])
    ifxtyp = int(lines[56][24:])
    munits = int(lines[57][24:])
    re = float(lines[58][24:])
    hspace = float(lines[59][24:])
    vbar = float(lines[60][24:])
    if vbar < 0:
      vbar = (v1 + v2) / 2.0
    ref_lat = float(lines[61][24:])
    outfile.write(writer("(7I5,I2,1X,I2,3F10.3,10x,f10.3)").write((model, itype, ibmax, zero, noprnt, nmol, ipunch, ifxtyp, munits, re, hspace, vbar, ref_lat)) + "\n")
    
  #Record 3.2
  if itype == 2 or itype == 3:
    h1 = float(lines[63][24:])
    h2 = float(lines[64][24:])
    angle = float(lines[65][24:])
    length = float(lines[66][24:])
    beta = float(lines[67][24:])
    path = float(lines[68][24:])
    hobs = float(lines[69][24:])
    if h2 == 0.0 and itype == 3:
      length = hspace/cos(angle*pi/180.)
      beta = arctan(hspace*tan(angle*pi/180.)/(hspace+re))
    outfile.write(writer("(5F10.3,I5,5x,F10.3)").write((h1, h2, angle, length, beta, path, hobs)) + "\n")
    
  #Record 3.3a
  if ibmax == 0:
    avtrat = float(lines[71][24:])
    tdiff1 = float(lines[72][24:])
    tdiff2 = float(lines[73][24:])
    altd1 = float(lines[74][24:])
    altd2 = float(lines[75][24:])
    outfile.write(writer("(5F10.3)").write((avtrat, tdiff1, tdiff2, altd1, altd2)) + "\n")
    
  #Record 3.4
  if model == 0:
    immax = int(lines[78][24:])
    hmod = lines[79][24:].strip()
    outfile.write(writer("(I5,A24)").write((immax, hmod)) + "\n")
    
  #Record 3.5
  #NOTE: parameter file should already have this in the correct format!
  if model == 0:
    for i in range(immax):
      zm, pm, tm, jcharp, jchart, jlong, jchars = reader("(3G10.4, 5x, 2A1, 1x, A1, 1x, A39)").read(lines[81+2*i])
      #print zm, pm, tm, jcharp, jchart, jlong, jchars
      fmt = writer("(3E10.3, 5x, 2A1, 1x, A1, 1x, A39)")
      #print writer("(3E10.4)").write((zm, pm, tm))
      #print fmt.write((zm, pm, tm, jcharp, jchart, jlong, jchars))
      outfile.write(writer("(3E10.4, 5x, 2A1, 1x, A1, 1x, A39)").write((zm, pm, tm, jcharp, jchart, jlong, jchars)) + "\n")
      vmol = []
      line = lines[81 + 2*i + 1]
      if jlong == "L":
        vmol = reader("(%iE15.8)" %nmol).read(lines[81+2*i+1])
        outfile.write(writer("(%iE15.8)" %nmol).write(vmol) + '\n')
      else:
        vmol = reader("(%iE10.3)" %nmol).read(lines[81+2*i+1])
        outfile.write(writer("(%iE10.3)" %nmol).write(vmol) + '\n')
      
  #Finish of the file
  stopflg = -1.0
  outfile.write(writer("(F4.1)").write((stopflg,)) + "\n")
  outfile.write(writer("(F4.1)").write((stopflg,)) + "\n")
  outfile.write("%")
    
  

def ReadParFile(parameterfile="ParameterFile"):
  infile = open(parameterfile)
  lines = infile.readlines()
  infile.close()

  pars = []
  for line in lines:
    if not line.startswith("!") and not line.startswith("#"):
      pars.append(line[24:].strip())

  return pars


def WriteTape5(pars, atmosphere=None, output="TAPE5"):
  outfile = open(output, "w")

  description = pars[0]
  outfile.write(writer(("A%i" %min(80, len(description)+1))).write(("$%s" %description, )) + "\n")

  #Record 1.2
  ihirac = int(pars[1])
  ilblf4 = int(pars[2])
  icntnm = int(pars[3])
  iaersl = int(pars[4])
  iemit = int(pars[5])
  iscan = int(pars[6])
  ifilter = int(pars[7])
  iplot = int(pars[8])
  itest = int(pars[9])
  iatm = int(pars[10])
  imrg = int(pars[11])
  ilas = int(pars[12])
  iod = int(pars[13])
  ixsec = int(pars[14])
  mpts = int(pars[15])
  npts = int(pars[16])
  outfile.write(writer("(10(4x,i1),3x,i2, 3(4x,i1),2(1x,i4))").write((ihirac, ilblf4, icntnm, iaersl, iemit, iscan, ifilter, iplot, itest, iatm, imrg, ilas, iod, ixsec, mpts, npts)) + "\n")

  if icntnm == 6:
    sys.exit("Error! icntnm = 6 not supported!")
  if iemit == 2:
    sys.exit("Error! iemit = 2 not supported!")
  
  #Record 1.3
  v1 = float(pars[17])
  v2 = float(pars[18])
  sample = float(pars[19])
  dvset = float(pars[20])
  alfal0 = float(pars[21])
  avmass = float(pars[22])
  dptmin = float(pars[23])
  dptflg = float(pars[24])
  ilnflg = int(pars[25])
  dvout = float(pars[26])
  nmol_scal = int(pars[27])
  outfile.write(writer("(8(F10.3),4x,I1,5x,F10.3,3x,I2)").write((v1, v2, sample, dvset, alfal0, avmass, dptmin, dptflg, ilnflg, dvout, nmol_scal)) + "\n")

  #Record 1.3.a and 1.3.b (TODO!)
  if nmol_scal > 0:
    sys.exit("Error! nmol_scal > 0 not supported!")

  #Record 1.4
  if iemit == 1 or (iemit == 2 and iotflg == 2):
    tbound = float(pars[28])
    sremis1 = float(pars[29])
    sremis2 = float(pars[30])
    sremis3 = float(pars[31])
    srrefl1 = float(pars[32])
    srrefl2 = float(pars[33])
    srrefl3 = float(pars[34])
    surf_refl = pars[35]
    outfile.write(writer("(7(F10.3),4x,A1)").write((tbound, sremis1, sremis2, sremis3, srrefl1, srrefl2, srrefl3, surf_refl)) + "\n")

  #Record 3.1
  if iatm == 1:
    model = int(pars[36])
    itype = int(pars[37])
    ibmax = int(pars[38])
    zero = int(pars[39])
    noprnt = int(pars[40])
    nmol = int(pars[41])
    ipunch = int(pars[42])
    ifxtyp = int(pars[43])
    munits = int(pars[44])
    re = float(pars[45])
    hspace = float(pars[46])
    vbar = float(pars[47])
    if vbar < 0:
      vbar = (v1 + v2)/2.0
    ref_lat = float(pars[48])
    outfile.write(writer("(7I5,I2,1X,I2,3F10.3,10x,f10.3)").write((model, itype, ibmax, zero, noprnt, nmol, ipunch, ifxtyp, munits, re, hspace, vbar, ref_lat)) + "\n")

  #Record 3.2
  if itype ==2 or itype == 3:
    h1 = float(pars[49])
    h2 = float(pars[50])
    angle = float(pars[51])
    length = float(pars[52])
    beta = float(pars[53])
    path = float(pars[54])
    hobs = float(pars[55])
    if h2 == 0.0 and itype == 3:
      length = hspace / cos(angle*pi/180.)
      beta = arctan(hspace*tan(angle*pi/180.) / (hspace + re))
    outfile.write(writer("(5F10.3,I5,5x,F10.3)").write((h1, h2, angle, length, beta, path, hobs)) + "\n")

  #Record 3.3a
  if ibmax == 0:
    avtrat = float(pars[56])
    tdiff1 = float(pars[57])
    tdiff2 = float(pars[58])
    altd1 = float(pars[59])
    altd2 = float(pars[60])
    outfile.write(writer("(5F10.3)").write((avtrat, tdiff1, tdiff2, altd1, altd2)) + "\n")

  #Record 3.4
  if model == 0:
    immax = int(pars[61])
    hmod = pars[62]
    outfile.write(writer("(I5,A24)").write((immax, hmod)) + "\n")
    if atmosphere == None:
      sys.exit("Error! User atmosphere requested but no atmosphere supplied!")

    #Record 3.5
    layers = sorted(atmosphere.keys())
    jcharp = "A"
    jchart = "A"
    jlong = "L"
    jchars = "A" * nmol
    for i in range(immax):
      zm = layers[i]
      pm, tm = atmosphere[zm][0], atmosphere[zm][1]
      abundances = atmosphere[zm][2]
      if len(abundances) != nmol:
        sys.exit("Error! len(abundances) != nmol in atmosphere. Exiting")
      outfile.write(writer("(3E10.4, 5x, 2A1, 1x, A1, 1x, A%i)" %min(39, nmol)).write((zm, pm, tm, jcharp, jchart, jlong, jchars)) + "\n")
      #outfile.write(writer("(%iE15.8)" %nmol).write(abundances) + '\n')
      start = 0
      end = min(start+8, nmol)
      done = False
      while not done:
        outfile.write(writer("(%iE15.8)" %(end-start)).write(abundances[start:end]) + '\n')
	start = end
	end = min(start+8, nmol)
	if start == end:
	  done = True


  #Finish the file
  stopflg = -1.0
  outfile.write(writer("(F4.1)").write((stopflg,)) + "\n")
  outfile.write(writer("(F4.1)").write((stopflg,)) + "\n")
  outfile.write("%")

  
  
if __name__ == "__main__":
  Convert(output="test")
  pars = ReadParFile()
  
  #Read in atmosphere
  infile = open("atmosphere")
  lines = infile.readlines()
  infile.close()

  Atmosphere = {}
  nmol = 12
  for i in range(0, len(lines), 2):
    zm, pm, tm, jcharp, jchart, jlong, jchars = reader("(3G10.4, 5x, 2A1, 1x, A1, 1x, A39)").read(lines[i])
    abundances = reader("(%iE15.8)" %nmol).read(lines[i+1])
    Atmosphere[zm] = (pm, tm, abundances)

  
  WriteTape5(pars, atmosphere = Atmosphere, output="test2")
