# BoxModel
README file for the 2-box model from Turner et al. (2017)

Alex Turner

December 23, 2016


# Info
 * The published paper can be found here: "https://doi.org/10.1073/pnas.1616020114".
 * The MATLAB code can be run using the "DriverScript.m" file.
 * The inversions use publicly available datasets.  As such, we have not included those datasets in our tarball.
 * We have included scripts to download all of the datasets though.  This can be done by running the "download_data.csh" script.
 * Users could also manually download the data by navigating the websites (see below), contacting the PIs for those datasets, contacting Alex Turner (aturner@fas.harvard.edu), or Christian Frankenberg (cfranken@caltech.edu).
 * The code will still run without the original datasets and will just use the hemispheric averages that were bootstrapped from Turner et al.
 * The code in the tarball is consistent with the code on the github (as of December 23, 2016).


# Public datasets:
 * CH4 from NOAA/ESRL
     * PI: Ed Dlugokencky
     * URL: "ftp://aftp.cmdl.noaa.gov/data/trace_gases/ch4/flask/"
 * MCF from NOAA/ESRL
     * PI: Steve Montzka
     * URL: "ftp://aftp.cmdl.noaa.gov/data/hats/solvents/CH3CCl3/"
 * MCF from GAGE/AGAGE
     * PI: Ron Prinn
     * URL: "http://agage.mit.edu/"
 * CH4C13 from NOAA/ESRL
     *  PI: James White
     *  URL: "ftp://aftp.cmdl.noaa.gov/data/trace_gases/ch4c13/flask/"
 * CH4C13 from U. Heidelberg
     *  PI: Ingeborg Levin
     *  URL: "http://www.iup.uni-heidelberg.de/institut/forschung/groups/kk/Data_html"
 * CH4C13 from U. Washington
     *  Data is included in the tarball
 * CH4C13 from UC Irvine
     *  Data is included in the tarball
 * C2H6 data (wrote an observation operator but did not use in the manuscript)
     *  URL: "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/download.cgi?para=VOCs"
