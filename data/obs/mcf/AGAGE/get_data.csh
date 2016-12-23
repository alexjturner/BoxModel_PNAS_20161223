#!/bin/csh -f
# Run this script to get the AGAGE MCF data
rm ./*-gcmd.mon
wget -r --no-parent -A '*-gcmd.mon' http://agage.eas.gatech.edu/data_archive/agage/gc-md/monthly/
mv agage.eas.gatech.edu/data_archive/agage/gc-md/monthly/*-gcmd.mon ./.
rm -rf agage.eas.gatech.edu
exit(0)
