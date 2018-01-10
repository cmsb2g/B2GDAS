#! /usr/bin/env python


## _________                _____.__                            __  .__               
## \_   ___ \  ____   _____/ ____\__| ____  __ ______________ _/  |_|__| ____   ____  
## /    \  \/ /  _ \ /    \   __\|  |/ ___\|  |  \_  __ \__  \\   __\  |/  _ \ /    \ 
## \     \___(  <_> )   |  \  |  |  / /_/  >  |  /|  | \// __ \|  | |  (  <_> )   |  \
##  \______  /\____/|___|  /__|  |__\___  /|____/ |__|  (____  /__| |__|\____/|___|  /
##         \/            \/        /_____/                   \/                    \/ 
import sys
import array as array
from plot_mttbar import plot_mttbar
import subprocess

import os

path = '/eos/uscms/store/user/cmsdas/2018/long_exercises/B2GTTbar/'

# Dictionaries
filenames = {
	'QCD' : [],
	'singleMuon' : [],
	'SingleElectron' : [],
	'WJets' : [],
	'rsg' : [],
	'ttbar' : [],
	'singletop' : []
}

outnames = {
	'QCD' : [],
	'singleMuon' : [],
	'SingleElectron' : [],
	'WJets' : [],
	'rsg' : [],
	'ttbar' : [],
	'singletop' : []
}


# Extract file names
for name in filenames.keys():
	files, outfiles = [], []
	batcmd="xrdfs root://cmseos.fnal.gov ls -u /store/user/cmsdas/2018/long_exercises/B2GTTbar/"
	temps = subprocess.check_output(batcmd, shell=True)
	for file in temps.split("\n"):
	    #print file.split("/")
	    if file.split("/")[-1].startswith(name) :
	        filenames[name].append(file)
	        outnames[name].append(file.split("/")[-1][0:-5])

# Compile function inputs
ins = []
for leptype in ['mu', 'ele']:
	for typ in filenames.keys(): 
		for i, n in enumerate(filenames[typ]):
			in_file = filenames[typ][i]
			out_file = outnames[typ][i]+"_plots_"+leptype+".root"
			ins.append(["--file_in", in_file, "--file_out", out_file, "--lepton", leptype])  # can include --jer up/down or --jec up/down

# Run in parallel
from multiprocessing import Pool
p = Pool(10)
print(p.map(plot_mttbar, ins))
