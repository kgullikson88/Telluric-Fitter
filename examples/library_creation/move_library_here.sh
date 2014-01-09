#After running the python script to make the library, you can move it here by executing this

for fname in `cat Library_Files.dat`
do
  mv $fname .
done
