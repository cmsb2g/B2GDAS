B2GDAS
======


CMS Data Analysis School (CMSDAS) exercise for the
Beyond Two Generations Physics Analysis Group (B2G PAG)

To run :


`cmsrel CMSSW_7_6_4`

`cd CMSSW_7_6_4/src`

`git clone https://github.com/cmsb2g/B2GDAS.git Analysis/B2GDAS`

`cd Analysis/B2GDAS`

`scram b -j 10`

`cd test`

`voms-proxy-init`

`python b2gdas_fwlite.py --files=inputfiles/rsgluon_ttbar_3TeV.txt --outname=rsgluon_ttbar_3TeV.root`
