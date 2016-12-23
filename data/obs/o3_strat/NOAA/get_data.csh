#!/bin/csh -f
# Run this script to get the NOAA stratospheric ozone data
rm ./*TOTo3.txt
wget -r --no-parent -A '*TOTo3.txt' ftp://aftp.cmdl.noaa.gov/data/ozwv/Dobson/
mv aftp.cmdl.noaa.gov/data/ozwv/Dobson/*TOTo3.txt ./.
rm -rf aftp.cmdl.noaa.gov
exit(0)
