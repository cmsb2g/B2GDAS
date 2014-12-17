#! /usr/bin/env python


## _________                _____.__                            __  .__               
## \_   ___ \  ____   _____/ ____\__| ____  __ ______________ _/  |_|__| ____   ____  
## /    \  \/ /  _ \ /    \   __\|  |/ ___\|  |  \_  __ \__  \\   __\  |/  _ \ /    \ 
## \     \___(  <_> )   |  \  |  |  / /_/  >  |  /|  | \// __ \|  | |  (  <_> )   |  \
##  \______  /\____/|___|  /__|  |__\___  /|____/ |__|  (____  /__| |__|\____/|___|  /
##         \/            \/        /_____/                   \/                    \/ 

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--files', type='string', action='store',
                  dest='files',
                  help='Input files')

parser.add_option('--outname', type='string', action='store',
                  default='outplots.root',
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

parser.add_option('--maxjets', type='int', action='store',
                  default=999,
                  dest='maxjets',
                  help='Number of jets to plot. To plot all jets, set to a big number like 999')


parser.add_option('--bDiscMin', type='float', action='store',
                  default=0.200,
                  dest='bDiscMin',
                  help='Minimum b discriminator')

parser.add_option('--minAK4Pt', type='float', action='store',
                  default=30.,
                  dest='minAK4Pt',
                  help='Minimum PT for AK4 jets')

parser.add_option('--maxAK4Rapidity', type='float', action='store',
                  default=2.4,
                  dest='maxAK4Rapidity',
                  help='Maximum AK4 rapidity')

parser.add_option('--minAK8Pt', type='float', action='store',
                  default=400.,
                  dest='minAK8Pt',
                  help='Minimum PT for AK8 jets')

parser.add_option('--maxAK8Rapidity', type='float', action='store',
                  default=2.4,
                  dest='maxAK8Rapidity',
                  help='Maximum AK8 rapidity')

(options, args) = parser.parse_args()
argv = []


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
from leptonic_nu_z_component import solve_nu_tmass, solve_nu

muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
electrons, electronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"
photons, photonLabel = Handle("std::vector<pat::Photon>"), "slimmedPhotons"
taus, tauLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
jets, jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
ak8jets, ak8jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"
mets, metLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"


##   ___ ___ .__          __                                             
##  /   |   \|__| _______/  |_  ____   ________________    _____   ______
## /    ~    \  |/  ___/\   __\/  _ \ / ___\_  __ \__  \  /     \ /  ___/
## \    Y    /  |\___ \  |  | (  <_> ) /_/  >  | \// __ \|  Y Y  \\___ \ 
##  \___|_  /|__/____  > |__|  \____/\___  /|__|  (____  /__|_|  /____  >
##        \/         \/             /_____/            \/      \/     \/

f = ROOT.TFile(options.outname, "RECREATE")
f.cd()


h_mttbar = ROOT.TH1F("h_mttbar", ";m_{t#bar{t}} (GeV)", 600, 0, 6000)

h_ptAK8 = ROOT.TH1F("h_ptAK8", "AK8 Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
h_etaAK8 = ROOT.TH1F("h_etaAK8", "AK8 Jet #eta;#eta", 120, -6, 6)
h_yAK8 = ROOT.TH1F("h_yAK8", "AK8 Jet Rapidity;y", 120, -6, 6)
h_phiAK8 = ROOT.TH1F("h_phiAK8", "AK8 Jet #phi;#phi (radians)",100,-ROOT.Math.Pi(),ROOT.Math.Pi())
h_mAK8 = ROOT.TH1F("h_mAK8", "AK8 Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mprunedAK8 = ROOT.TH1F("h_mprunedAK8", "AK8 Pruned Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mfilteredAK8 = ROOT.TH1F("h_mfilteredAK8", "AK8 Filtered Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mtrimmedAK8 = ROOT.TH1F("h_mtrimmedAK8", "AK8 Trimmed Jet Mass;Mass (GeV)", 100, 0, 1000)
h_minmassAK8 = ROOT.TH1F("h_minmassAK8", "AK8 CMS Top Tagger Min Mass Paring;m_{min} (GeV)", 100, 0, 1000)
h_nsjAK8 = ROOT.TH1F("h_nsjAK8", "AK8 CMS Top Tagger N_{subjets};N_{subjets}", 5, 0, 5)
h_tau21AK8 = ROOT.TH1F("h_tau21AK8", "AK8 Jet #tau_{2} / #tau_{1};Mass#tau_{21}", 100, 0, 1.0)
h_tau32AK8 = ROOT.TH1F("h_tau32AK8", "AK8 Jet #tau_{3} / #tau_{2};Mass#tau_{32}", 100, 0, 1.0)


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
        s = 'root://cmsxrootd.fnal.gov/' + ifile.rstrip()
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
    for event in events:
        if options.maxevents > 0 and nevents > options.maxevents :
            break
        i += 1
        nevents += 1

        if i % 1000 == 0 :
            print '    ---> Event ' + str(i)


        ## ____   ____             __                    _________      .__                 __  .__               
        ## \   \ /   /____________/  |_  ____ ___  ___  /   _____/ ____ |  |   ____   _____/  |_|__| ____   ____  
        ##  \   Y   // __ \_  __ \   __\/ __ \\  \/  /  \_____  \_/ __ \|  | _/ __ \_/ ___\   __\  |/  _ \ /    \ 
        ##   \     /\  ___/|  | \/|  | \  ___/ >    <   /        \  ___/|  |_\  ___/\  \___|  | |  (  <_> )   |  \
        ##    \___/  \___  >__|   |__|  \___  >__/\_ \ /_______  /\___  >____/\___  >\___  >__| |__|\____/|___|  /
        ##               \/                 \/      \/         \/     \/          \/     \/                    \/ 


        event.getByLabel(vertexLabel, vertices)

        # Vertices
        if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
            print "Event has no good primary vertex."
            continue
        else:
            PV = vertices.product()[0]
            if options.verbose : 
                print "PV at x,y,z = %+5.3f, %+5.3f, %+6.3f (ndof %.1f)" % (PV.x(), PV.y(), PV.z(), PV.ndof())

            
            
        ## .____                  __                    _________      .__                 __  .__               
        ## |    |    ____ _______/  |_  ____   ____    /   _____/ ____ |  |   ____   _____/  |_|__| ____   ____  
        ## |    |  _/ __ \\____ \   __\/  _ \ /    \   \_____  \_/ __ \|  | _/ __ \_/ ___\   __\  |/  _ \ /    \ 
        ## |    |__\  ___/|  |_> >  | (  <_> )   |  \  /        \  ___/|  |_\  ___/\  \___|  | |  (  <_> )   |  \
        ## |_______ \___  >   __/|__|  \____/|___|  / /_______  /\___  >____/\___  >\___  >__| |__|\____/|___|  /
        ##         \/   \/|__|                    \/          \/     \/          \/     \/                    \/ 



        event.getByLabel( muonLabel, muons )
        event.getByLabel( electronLabel, electrons )




        # Select tight good muons
        goodmuons = []
        if len(muons.product()) > 0 :
            for i,muon in enumerate( muons.product() ) :
                if muon.pt() > 45.0 and abs(muon.eta()) < 2.1 and muon.muonBestTrack().dz(PV.position()) < 5.0 and muon.isTightMuon(PV) :
                    goodmuons.append( muon )
                if options.verbose :
                    print "muon %2d: pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d." % (
                        i, muon.pt(), muon.muonBestTrack().dz(PV.position()), muon.isLooseMuon(), muon.isTightMuon(PV))

        # Select tight good electrons
        goodelectrons = []
        if len(electrons.product()) > 0 :
            for i,electron in enumerate( electrons.product() ) :

                if electron.pt() < 45.0 :
                    continue

                dEtaIn = electron.deltaEtaSuperClusterTrackAtVtx()
                dPhiIn = electron.deltaPhiSuperClusterTrackAtVtx()
                hOverE = electron.hcalOverEcal()
                sigmaIetaIeta = electron.sigmaIetaIeta()
                full5x5_sigmaIetaIeta = electron.full5x5_sigmaIetaIeta()
                if abs(electron.ecalEnergy()) > 0.0 : 
                    ooEmooP = abs(1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy() )
                else :
                    ooEmooP = 1.0e30

                d0 = (-1.0) * electron.gsfTrack().dxy(PV.position())
                dz = electron.gsfTrack().dz(PV.position())
                pfIso = electron.pfIsolationVariables()
                absIso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
                relIso = absIso/electron.pt()
                #expectedMissingInnerHits = electron.gsfTrack().trackerExpectedHitsInner().numberOfLostHits()
                passConversionVeto = electron.passConversionVeto

                goodElectron = False
                # Barrel ECAL cuts
                if abs(electron.superCluster().eta()) < 1.479 :                    
                    goodElectron = \
                      abs( dEtaIn ) < 0.0091 and \
                      abs(dPhiIn) < 0.031 and \
                      full5x5_sigmaIetaIeta < 0.0106 and \
                      hOverE < 0.0532 and \
                      abs(d0) < 0.0126 and \
                      abs(dz) < 0.0116 and \
                      abs( ooEmooP ) < 0.0609 and \
                      relIso < 0.1649 and \
                      passConversionVeto
                # Endcap ECAL cuts
                elif abs(electron.superCluster().eta()) < 2.5 and abs(electron.superCluster().eta()) > 1.479 :
                    goodElectron = \
                      abs( dEtaIn ) < 0.0106 and \
                      abs(dPhiIn) < 0.0359 and \
                      full5x5_sigmaIetaIeta < 0.0305 and \
                      hOverE < 0.0835 and \
                      abs(d0) < 0.0163 and \
                      abs(dz) < 0.5999 and \
                      abs( ooEmooP ) < 0.1126 and \
                      relIso < 0.2075 and \
                      passConversionVeto
                
                if goodElectron == True :
                    goodelectrons.append( electron )
                if options.verbose :
                    print "elec %2d: pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d" % (
                        i, electron.pt(), electron.superCluster().eta(), electron.sigmaIetaIeta(), electron.full5x5_sigmaIetaIeta(), electron.gsfTrack().trackerExpectedHitsInner().numberOfLostHits(), electron.passConversionVeto())
        


        # Veto on dilepton events
        if len(goodmuons) + len(goodelectrons) != 1 :
            continue
        elif len(goodmuons) > 0 :
            theLepton = ROOT.TLorentzVector( goodmuons[0].px(),
                                             goodmuons[0].py(),
                                             goodmuons[0].pz(),
                                             goodmuons[0].energy() )
        else :
            theLepton = ROOT.TLorentzVector( goodelectrons[0].px(),
                                             goodelectrons[0].py(),
                                             goodelectrons[0].pz(),
                                             goodelectrons[0].energy() )



        ##      ____.       __      _________      .__                 __  .__               
        ##     |    | _____/  |_   /   _____/ ____ |  |   ____   _____/  |_|__| ____   ____  
        ##     |    |/ __ \   __\  \_____  \_/ __ \|  | _/ __ \_/ ___\   __\  |/  _ \ /    \ 
        ## /\__|    \  ___/|  |    /        \  ___/|  |_\  ___/\  \___|  | |  (  <_> )   |  \
        ## \________|\___  >__|   /_______  /\___  >____/\___  >\___  >__| |__|\____/|___|  /
        ##               \/               \/     \/          \/     \/                    \/ 

        #
        #
        #
        # Here, we have TWO jet collections, AK4 and AK8. The
        # AK4 jets are used for b-tagging, while the AK8 jets are used
        # for top-tagging.
        # In the future, the AK8 jets will contain "subjet b-tagging" in
        # miniAOD but as of now (Dec 2014) this is not ready so we
        # need to adjust.
        #
        #
        
        # use getByLabel, just like in cmsRun
        event.getByLabel (jetLabel, jets)          # For b-tagging
        event.getByLabel (ak8jetLabel, ak8jets)    # For top-tagging
        
        # loop over jets and fill hists
        ijet = 0

        # These will hold all of the jets we need for the selection
        ak4JetsGood = []
        ak8JetsGood = []

        # For selecting leptons, look at 2-d cut of dRMin, ptRel of
        # lepton and nearest jet that has pt > 30 GeV
        dRMin = 9999.0
        inearestJet = -1    # Index of nearest jet
        nearestJet = None   # Nearest jet

        ############################################
        # First get the AK4 jet nearest the lepton :
        ############################################ 
        for i,jet in enumerate(jets.product()) :

            if jet.pt() < options.minAK4Pt or abs(jet.rapidity()) > options.maxAK4Rapidity :
                continue
            # perform jet ID with UNCORRECTED jet energy
            jetP4 = ROOT.TLorentzVector( jet.px(), jet.py(), jet.pz(), jet.energy() )
            p4Raw = jetP4 * jet.jecFactor(0)

            nhf = jet.neutralHadronEnergy() / p4Raw.E()
            nef = jet.neutralEmEnergy() / p4Raw.E()
            chf = jet.chargedHadronEnergy() / p4Raw.E()
            cef = jet.chargedEmEnergy() / p4Raw.E()
            nconstituents = jet.numberOfDaughters()
            nch = jet.chargedMultiplicity()
            goodJet = \
              nhf < 0.99 and \
              nef < 0.99 and \
              chf > 0.00 and \
              cef < 0.99 and \
              nconstituents > 1 and \
              nch > 0

            if goodJet:
                dR = jetP4.DeltaR(theLepton )
                ak4JetsGood.append(jet)
                if dR < dRMin :
                    inearestJet = ijet
                    nearestJet = jet
                    dRMin = dR
                    
                    
        # Require at least one leptonic-side jet, and 2d isolation cut
        if nearestJet == None :
            continue
        theLepJet = ROOT.TLorentzVector( nearestJet.px(), nearestJet.py(), nearestJet.pz(), nearestJet.energy() )
        ptRel = theLepJet.Perp( theLepton.Vect() )
        pass2D = ptRel > 20.0 or dRMin > 0.5
        if pass2D == False :
            continue

        ############################################
        # Second, get the AK8 jet away from the lepton
        ############################################
        for i,jet in enumerate(ak8jets.product()) :

            if jet.pt() < options.minAK8Pt or abs(jet.rapidity()) > options.maxAK8Rapidity :
                continue
            # perform jet ID with UNCORRECTED jet energy
            jetP4 = ROOT.TLorentzVector( jet.px(), jet.py(), jet.pz(), jet.energy() )
            p4Raw = jetP4 * jet.jecFactor(0)

            nhf = jet.neutralHadronEnergy() / p4Raw.E()
            nef = jet.neutralEmEnergy() / p4Raw.E()
            chf = jet.chargedHadronEnergy() / p4Raw.E()
            cef = jet.chargedEmEnergy() / p4Raw.E()
            nconstituents = jet.numberOfDaughters()
            nch = jet.chargedMultiplicity()
            goodJet = \
              nhf < 0.99 and \
              nef < 0.99 and \
              chf > 0.00 and \
              cef < 0.99 and \
              nconstituents > 1 and \
              nch > 0

            if goodJet:
                dR = jetP4.DeltaR(theLepton )
                if dR > ROOT.TMath.Pi()/2.0 :
                    ak8JetsGood.append(jet)
        


        ############################################
        # Investigate the b-tagging and t-tagging
        ############################################
        if len(ak4JetsGood) < 1 or len(ak8JetsGood) < 1 :
            continue

        nbtags = 0
        bJetIndices = []
        for ijet,jet in enumerate(ak4JetsGood) :
            bdisc = jet.bDiscriminator('')
            if bdisc > options.bDiscMin :
                nbtags += 1
                bJetIndices.append( ijet )
            

        nttags = 0
        tJets = []
        for ijet,jet in enumerate(ak8JetsGood) : 
            if jet.pt() < options.minAK8Pt :
                continue

            mAK8Pruned = jet.userFloat('ak8PFJetsCHSPrunedLinks')
            mAK8Filtered = jet.userFloat('ak8PFJetsCHSFilteredLinks')
            mAK8Trimmed = jet.userFloat('ak8PFJetsCHSTrimmedLinks')
            # Make sure there are top tags if we want to plot them
            minMass = -1.0
            nsubjets = -1                
            tagInfoLabels = jet.tagInfoLabels()
            hasTopTagInfo = 'caTop' in tagInfoLabels 
            if hasTopTagInfo :
                minMass = jet.tagInfo('caTop').properties().minMass
                nsubjets = jet.tagInfo('caTop').properties().nSubJets
            # Get n-subjettiness "tau" variables
            tau1 = jet.userFloat('NjettinessAK8:tau1')
            tau2 = jet.userFloat('NjettinessAK8:tau2')
            tau3 = jet.userFloat('NjettinessAK8:tau3')
            if tau1 > 0.0001 :
                tau21 = tau2 / tau1
                h_tau21AK8.Fill( tau21 )
            else :
                h_tau21AK8.Fill( -1.0 )
            if tau2 > 0.0001 :
                tau32 = tau3 / tau2
                h_tau32AK8.Fill( tau32 )
            else :
                h_tau32AK8.Fill( -1.0 )

            # Perform CMS top tagging with trimmed jet mass
            print 'minMass = {0:6.2f}, trimmed mass = {1:6.2f}, tau32 = {2:6.2f}'.format(
                minMass, mAK8Trimmed, tau32
                ), 
            if minMass > 50.0 and mAK8Trimmed > 100. and mAK8Trimmed < 250. and tau32 > 0.4 :
                nttags += 1
                tJets.append( jet )
                print ' --->Tagged jet!'
            else :
                print ''


        # Now we do our kinematic calculation based on the categories of the
        # number of top and bottom tags

        if nttags > 0 :
            # Finally get the METs
            event.getByLabel( metLabel, mets )
            met = mets.product()[0]
            
            hadTopCandP4 = ROOT.TLorentzVector( tJets[0].px(), tJets[0].py(), tJets[0].pz(), tJets[0].energy() )
            lepTopCandP4 = None
            
            # Check if the nearest jet to the lepton is b-tagged
            if inearestJet not in bJetIndices :
                print 'closest jet not tagged'
            else  :


                # Get the z-component of the lepton from the W mass constraint
                bJetCandP4 = ROOT.TLorentzVector( ak4JetsGood[inearestJet].px(),
                                                  ak4JetsGood[inearestJet].py(),
                                                  ak4JetsGood[inearestJet].pz(),
                                                  ak4JetsGood[inearestJet].energy())
                nuCandP4 = ROOT.TLorentzVector( met.px(), met.py(), 0, met.energy() )

                nuz1 = None
                nuz2 = None
                solution = solve_nu( vlep=theLepton,
                                     vnu=nuCandP4,
                                     nuz1=nuz1,
                                     nuz2=nuz2
                                     )
                # If there is at least one real solution, pick it up
                if solution :
                    print '--- Have a solution --- '
                    nuCandP4.SetPz( nuz1 )
                    lepTopCandP4 = nuCandP4 + theLepton + bJetCandP4

                    ttbarCand = hadTopCandP4 + lepTopCandP4
                    h_mttbar.Fill( ttbarCand.M() )
                else :
                    print '--- No solution for neutrino z ---'            
                    

        

## _________ .__                                     
## \_   ___ \|  |   ____ _____    ____  __ ________  
## /    \  \/|  | _/ __ \\__  \  /    \|  |  \____ \ 
## \     \___|  |_\  ___/ / __ \|   |  \  |  /  |_> >
##  \______  /____/\___  >____  /___|  /____/|   __/ 
##         \/          \/     \/     \/      |__|    
f.cd()
f.Write()
f.Close()
