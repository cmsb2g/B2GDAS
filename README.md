B2GDAS
======


CMS Data Analysis School (CMSDAS) exercise for the
Beyond Two Generations Physics Analysis Group (B2G PAG)

To run :


<verbatim>
cmsrel CMSSW\_7\_3\_0
cd CMSSW\_7\_3\_0/src
mkdir Analysis
cd Analysis
git clone https://github.com/cmsb2g/B2GDAS.git
scram b -j 10
cd test
voms-proxy-init (to get a grid proxy for xrootd)
python b2gdas\_fwlite.py --files=rsgluon\_ttbar\_3TeV.txt --outname=rsgluon\_ttbar\_3TeV.root
</verbatim>
