import MakeModel


humidity = range(10, 90, 10)   #Let humidity go from 10% to 90%, in steps of 10%
o2 = [1.5e5, 2.1e5, 2.5e5]  #Let the O2 surface abundance range vary too

# Set the start and end wavelength, in nm!
wavestart = 500
waveend = 900

#Make an instance of the Modeler class
modeler = MakeModel.Modeler()
for h2o in humidity:
  for o2val in o2:
    modeler.MakeModel(humidity=h2o, o2=o2val, save=True, libfile="Library_Files.dat", lowfreq=1e7/waveend, highfreq=1e7/wavestart)

print "The filenames of your telluric library are given in the file 'Library_Files.dat' in the current directory."
print "You can move them all here with the script move_library_here.sh"
