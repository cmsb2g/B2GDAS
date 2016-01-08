#! /usr/bin/env python


## _________                _____.__                            __  .__               
## \_   ___ \  ____   _____/ ____\__| ____  __ ______________ _/  |_|__| ____   ____  
## /    \  \/ /  _ \ /    \   __\|  |/ ___\|  |  \_  __ \__  \\   __\  |/  _ \ /    \ 
## \     \___(  <_> )   |  \  |  |  / /_/  >  |  /|  | \// __ \|  | |  (  <_> )   |  \
##  \______  /\____/|___|  /__|  |__\___  /|____/ |__|  (____  /__| |__|\____/|___|  /
##         \/            \/        /_____/                   \/                    \/ 
import sys
from array import array
from optparse import OptionParser


def makepuhist(argv) : 
    parser = OptionParser()

    parser.add_option('--file_data', type='string', action='store',
                      dest='file_data',
                      help='Input file for data')

    parser.add_option('--file_mc', type='string', action='store',
                      dest='file_mc',
                      help='Input file for MC')

    parser.add_option('--file_out', type='string', action='store',
                      dest='file_out',
                      help='Output file')
    
    (options, args) = parser.parse_args(argv)
    argv = []

    print '===== Command line options ====='
    print options
    print '================================'

    import ROOT

    f1 = ROOT.TFile(options.file_data)
    f2 = ROOT.TFile(options.file_mc)
    f3 = ROOT.TFile(options.file_out, "RECREATE")

    h1 = f1.Get("pileup")
    h2 = f2.Get("pileup")

    h1.Sumw2()
    h1.Scale(1.0 / h1.Integral() )
    h2.Sumw2()
    h2.Scale(1.0 / h2.Integral() )

    h1.Divide(h2)

    f3.cd()
    h1.Write()
    f3.Close()



if __name__ == "__main__" :
    makepuhist(sys.argv)


