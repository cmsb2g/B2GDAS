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


ttjetsStr = ["--file_in", "ttjets.root", "--file_out", "ttjets_plots.root"]
rsgluon3TeVStr = ["--file_in", "rsgluon_ttbar_3TeV.root", "--file_out", "rsluon3TeV_plots.root"]
singlemuStr = ["--file_in", "singlemuon.root", "--file_out", "singlemuon_plots.root", "--isData"]

plot_mttbar( ttjetsStr )
plot_mttbar( rsgluon3TeVStr )
plot_mttbar( singlemuStr )


