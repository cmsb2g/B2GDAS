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


def makepu_fwlite(argv) : 
    parser = OptionParser()

    parser.add_option('--files', type='string', action='store',
                      dest='files',
                      help='Input files')

    parser.add_option('--outname', type='string', action='store',
                      default='pumc.root',
                      dest='outname',
                      help='Name of output file')

    parser.add_option('--verbose', action='store_true',
                      default=False,
                      dest='verbose',
                      help='Print debugging info')

    parser.add_option('--maxevents', type='int', action='store',
                      default=-1,
                      dest='maxevents',
                      help='Number of events to run. -1 is all events')

    (options, args) = parser.parse_args(argv)
    argv = []

    print '===== Command line options ====='
    print options
    print '================================'


    ## _____________      __.____    .__  __             _________ __          _____  _____ 
    ## \_   _____/  \    /  \    |   |__|/  |_  ____    /   _____//  |_ __ ___/ ____\/ ____\
    ##  |    __) \   \/\/   /    |   |  \   __\/ __ \   \_____  \\   __\  |  \   __\\   __\ 
    ##  |     \   \        /|    |___|  ||  | \  ___/   /        \|  | |  |  /|  |   |  |   
    ##  \___  /    \__/\  / |_______ \__||__|  \___  > /_______  /|__| |____/ |__|   |__|   
    ##      \/          \/          \/             \/          \/                           

    import ROOT
    import sys
    from DataFormats.FWLite import Events, Handle
    ROOT.gROOT.Macro("rootlogon.C")

    vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
    
    ##   ___ ___ .__          __                                             
    ##  /   |   \|__| _______/  |_  ____   ________________    _____   ______
    ## /    ~    \  |/  ___/\   __\/  _ \ / ___\_  __ \__  \  /     \ /  ___/
    ## \    Y    /  |\___ \  |  | (  <_> ) /_/  >  | \// __ \|  Y Y  \\___ \ 
    ##  \___|_  /|__/____  > |__|  \____/\___  /|__|  (____  /__|_|  /____  >
    ##        \/         \/             /_____/            \/      \/     \/

    f = ROOT.TFile(options.outname, "RECREATE")
    f.cd()
    # and also make a few 1-d histograms
    pileup = ROOT.TH1F("pileup", "Pileup", 50, 0, 50)
        
    ## ___________                    __    .____                         
    ## \_   _____/__  __ ____   _____/  |_  |    |    ____   ____ ______  
    ##  |    __)_\  \/ // __ \ /    \   __\ |    |   /  _ \ /  _ \\____ \ 
    ##  |        \\   /\  ___/|   |  \  |   |    |__(  <_> |  <_> )  |_> >
    ## /_______  / \_/  \___  >___|  /__|   |_______ \____/ \____/|   __/ 
    ##         \/           \/     \/               \/            |__|    


    # IMPORTANT : Run one FWLite instance per file. Otherwise,
    # FWLite aggregates ALL of the information immediately, which
    # can take a long time to parse. 
    filelist = file( options.files )
    filesraw = filelist.readlines()
    files = []
    nevents = 0
    for ifile in filesraw :
        if len( ifile ) > 2 : 
            s = 'root://cmsxrootd-site.fnal.gov/' + ifile.rstrip()
            files.append( s )
            print 'Added ' + s


    # loop over files
    for ifile in files :
        print 'Processing file ' + ifile
        events = Events (ifile)
        if options.maxevents > 0 and nevents > options.maxevents :
            break

        # loop over events in this file
        i = 0
        for iev,event in enumerate(events):
            evWeight = 1.0
            
            if options.maxevents > 0 and nevents > options.maxevents :
                break
            i += 1
            nevents += 1

            if nevents % 1000 == 0 : 
                print '    ---> Event ' + str(nevents)
            if options.verbose :
                print '==============================================='
                print '    ---> Event ' + str(nevents)


            ## ____   ____             __                    _________      .__                 __  .__               
            ## \   \ /   /____________/  |_  ____ ___  ___  /   _____/ ____ |  |   ____   _____/  |_|__| ____   ____  
            ##  \   Y   // __ \_  __ \   __\/ __ \\  \/  /  \_____  \_/ __ \|  | _/ __ \_/ ___\   __\  |/  _ \ /    \ 
            ##   \     /\  ___/|  | \/|  | \  ___/ >    <   /        \  ___/|  |_\  ___/\  \___|  | |  (  <_> )   |  \
            ##    \___/  \___  >__|   |__|  \___  >__/\_ \ /_______  /\___  >____/\___  >\___  >__| |__|\____/|___|  /
            ##               \/                 \/      \/         \/     \/          \/     \/                    \/ 


            event.getByLabel(vertexLabel, vertices)
            # Vertices
            NPV = len(vertices.product())
            if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
                if options.verbose : 
                    print "Event has no good primary vertex."
                continue
            else:
                PV = vertices.product()[0]
                if options.verbose : 
                    print "PV at x,y,z = %+5.3f, %+5.3f, %+6.3f (ndof %.1f)" % (PV.x(), PV.y(), PV.z(), PV.ndof())
            pileup.Fill( NPV )

    ## _________ .__                                     
    ## \_   ___ \|  |   ____ _____    ____  __ ________  
    ## /    \  \/|  | _/ __ \\__  \  /    \|  |  \____ \ 
    ## \     \___|  |_\  ___/ / __ \|   |  \  |  /  |_> >
    ##  \______  /____/\___  >____  /___|  /____/|   __/ 
    ##         \/          \/     \/     \/      |__|    
    f.cd()
    f.Write()
    f.Close()


if __name__ == "__main__" :
    makepu_fwlite(sys.argv)


