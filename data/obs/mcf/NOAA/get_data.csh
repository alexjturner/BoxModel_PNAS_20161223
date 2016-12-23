#!/bin/csh -f
# Run this script to get the NOAA MCF data
rm -rf ./combined ; mkdir ./combined
wget "ftp://aftp.cmdl.noaa.gov/data/hats/solvents/CH3CCl3/combined/HATS_global_MC.txt"
mv HATS_global_MC.txt ./combined/.
exit(0)
