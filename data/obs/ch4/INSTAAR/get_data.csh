#!/bin/csh -f
# Run this script to get the NOAA CH4 data
rm -rf ./month ; mkdir ./month
wget -r --no-parent -A 'ch4_*_month.txt' ftp://aftp.cmdl.noaa.gov/data/trace_gases/ch4/flask/surface/
mv aftp.cmdl.noaa.gov/data/trace_gases/ch4/flask/surface/ch4_*_month.txt ./month/.
rm -rf aftp.cmdl.noaa.gov
exit(0)
