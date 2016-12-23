#!/bin/csh -f
# This script should download the VOC data from the WDCGG page.
# However, getting the VOC data requires submitting a query to the WDCGG page.
# The script *should* submit the query but sometimes it's a bit slow to archive the data.
# You can manually download the tarball from the WDCGG page here: "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/download.cgi?para=VOCs"
# Choose the "event" data here: "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/archiver.cgi?tool=gzip&archive=VOCS_EVENT"
# This should give you a VOCS_EVENT.tgz tarball.
# Untar the tarball and manually run the commands from this script.

rm -rf ./event VOCS_EVENT.tgz ; mkdir ./event
open "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/archiver.cgi?tool=gzip&archive=VOCS_EVENT"
wget "http://ds.data.jma.go.jp/gmd/wdcgg/tmp/ftptmp/VOCS_EVENT.tgz"
tar -xvf VOCS_EVENT.tgz
mv vocs/ethane/event/* ./event/.
rm ./event/bsl999900.niwa.am.fl.ethane.nl.ev.dat # Remove site w/ different format
rm -rf vocs VOCS_EVENT.tgz
exit(0)
