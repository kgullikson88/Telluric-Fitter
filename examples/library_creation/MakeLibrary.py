import MakeModel


humidity = range(10, 90, 10)   #Let humidity go from 10% to 90%, in steps of 10%
co2 = range(350, 400, 5)  #Let the CO2 surface abundance range from 350-400 ppmv, in steps of 5
wavestart = 1500.0
waveend = 1550.0

#Make an instance of the Modeler class
modeler = MakeModel.Modeler()
for h20 in humidity:
  for co2val in co2:
    modeler.MakeModel(h2o=h2o, co2=co2val, save=True, libfile="Library_Files.dat", lowfreq=1e7/waveend, highfreq=1e7/wavestart)

#The filenames of your telluric library are given in the file 'Library_Files.dat' in the current directory.
