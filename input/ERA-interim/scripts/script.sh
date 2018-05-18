#Small script to select a single point from grb data and convert it into ascii
#See pythonscript.py for how to prepare the data to be read into SAMSIM
for x in *.grb
do
  echo $x
  cdo selindexbox,1,1,1,1 $x $x.cut
  cdo output $x.cut > $x.txt
done

