#!/bin/csh -f
# Run this script to get the NOAA CH4C13 data
rm -rf ./month ; mkdir month
wget -r --no-parent -A 'ch4c13_*_month.txt' ftp://aftp.cmdl.noaa.gov/data/trace_gases/ch4c13/flask/surface/
mv aftp.cmdl.noaa.gov/data/trace_gases/ch4c13/flask/surface/ch4c13_*_month.txt ./month/.
rm -rf aftp.cmdl.noaa.gov
exit(0)
