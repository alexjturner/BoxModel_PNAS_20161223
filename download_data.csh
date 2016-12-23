#!/bin/csh -f
# Run this script to download all the datasets

### Click the WDCGG link so it starts archiving the VOC data
open "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/archiver.cgi?tool=gzip&archive=VOCS_EVENT"
set bDir = `pwd`

### Get the datasets
# CH4
cd ${bDir}/data/obs/ch4/INSTAAR/        ; ./get_data.csh
# CH4C13
cd ${bDir}/data/obs/ch4c13/INSTAAR/     ; ./get_data.csh
cd ${bDir}/data/obs/ch4c13/Heidelberg/  ; ./get_data.csh
# MCF
cd ${bDir}/data/obs/mcf/NOAA/           ; ./get_data.csh
cd ${bDir}/data/obs/mcf/AGAGE/          ; ./get_data.csh
cd ${bDir}/data/obs/mcf/GAGE/           ; ./get_data.csh
# Stratospheric ozone
cd ${bDir}/data/obs/o3_strat/NOAA/      ; ./get_data.csh
# C2H6
cd ${bDir}/data/obs/c2h6/WDCGG/         ; ./get_data.csh

### Done
exit(0)
