#!/bin/csh -f
# Run this script to get the U. Heidelberg CH4C13 data
rm d13Cdata*.txt
wget http://www.iup.uni-heidelberg.de/institut/forschung/groups/kk/Data/d13Cdata_Alert.txt
wget http://www.iup.uni-heidelberg.de/institut/forschung/groups/kk/Data/d13Cdata_Izana.txt
wget http://www.iup.uni-heidelberg.de/institut/forschung/groups/kk/Data/d13Cdata_Neumayer.txt
exit(0)
