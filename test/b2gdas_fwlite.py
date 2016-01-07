#! /usr/bin/env python



##   ___ ___         .__                        ___________                   __  .__                      
##  /   |   \   ____ |  | ______   ___________  \_   _____/_ __  ____   _____/  |_|__| ____   ____   ______
## /    ~    \_/ __ \|  | \____ \_/ __ \_  __ \  |    __)|  |  \/    \_/ ___\   __\  |/  _ \ /    \ /  ___/
## \    Y    /\  ___/|  |_|  |_> >  ___/|  | \/  |     \ |  |  /   |  \  \___|  | |  (  <_> )   |  \\___ \ 
##  \___|_  /  \___  >____/   __/ \___  >__|     \___  / |____/|___|  /\___  >__| |__|\____/|___|  /____  >
##        \/       \/     |__|        \/             \/             \/     \/                    \/     \/ 

def getJER(jetEta, sysType) :
    """
    Here, jetEta should be the jet pseudorapidity, and sysType is :
        nominal : 0
        down    : -1
        up      : +1
    """

    jerSF = 1.0

    if ( (sysType==0 or sysType==-1 or sysType==1) == False):
        print "ERROR: Can't get JER! use type=0 (nom), -1 (down), +1 (up)"
        return float(jerSF)

    # Values from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    etamin = [0.0,0.8,1.3,1.9,2.5,3.0,3.2]
    etamax = [0.8,1.3,1.9,2.5,3.0,3.2,5.0]
    scale_nom =    [1.061,1.088,1.106,1.126,1.343,1.303,1.320]
    scale_uncert = [0.023,0.029,0.030,0.094,0.123,0.111,0.286]
 
    # old 8 TeV
    # etamin = [0.0,0.5,1.1,1.7,2.3,2.8,3.2]
    # etamax = [0.5,1.1,1.7,2.3,2.8,3.2,5.0]
    # scale_nom = [1.079,1.099,1.121,1.208,1.254,1.395,1.056]
    # scale_dn  = [1.053,1.071,1.092,1.162,1.192,1.332,0.865]
    # scale_up  = [1.105,1.127,1.150,1.254,1.316,1.458,1.247]

    for iSF in range(0,len(scale_nom)) :
        if abs(jetEta) >= etamin[iSF] and abs(jetEta) < etamax[iSF] :
            if sysType < 0 :
                jerSF = scale_nom[iSF] - scale_uncert[iSF]
            elif sysType > 0 :
                jerSF = scale_nom[iSF] + scale_uncert[iSF]
            else :
                jerSF = scale_nom[iSF]
            break

    return float(jerSF)


## _________                _____.__                            __  .__               
## \_   ___ \  ____   _____/ ____\__| ____  __ ______________ _/  |_|__| ____   ____  
## /    \  \/ /  _ \ /    \   __\|  |/ ___\|  |  \_  __ \__  \\   __\  |/  _ \ /    \ 
## \     \___(  <_> )   |  \  |  |  / /_/  >  |  /|  | \// __ \|  | |  (  <_> )   |  \
##  \______  /\____/|___|  /__|  |__\___  /|____/ |__|  (____  /__| |__|\____/|___|  /
##         \/            \/        /_____/                   \/                    \/ 
import sys
from array import array
from optparse import OptionParser


def b2gdas_fwlite(argv) : 
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

    parser.add_option('--writeTree', action='store_true',
                      default=False,
                      dest='writeTree',
                      help='Write TTree?')    

    parser.add_option('--isData', action='store_true',
                      default=False,
                      dest='isData',
                      help='Is this data?')

    parser.add_option('--maxevents', type='int', action='store',
                      default=-1,
                      dest='maxevents',
                      help='Number of events to run. -1 is all events')

    parser.add_option('--maxjets', type='int', action='store',
                      default=999,
                      dest='maxjets',
                      help='Number of jets to plot. To plot all jets, set to a big number like 999')

    parser.add_option('--syst', type='string', action='store',
                      default='NOM',
                      dest='syst',
                      help='Jet systematic to run, options are NOM,JECUP,JECDN,JERUP,JERDN')

    
    parser.add_option('--bdisc', type='string', action='store',
                      default='pfCombinedInclusiveSecondaryVertexV2BJetTags',
                      dest='bdisc',
                      help='Name of output file')


    parser.add_option('--bdiscMin', type='float', action='store',
                      default=0.679,
                      dest='bDiscMin',
                      help='Minimum b discriminator')

    parser.add_option('--minMuonPt', type='float', action='store',
                      default=30.,
                      dest='minMuonPt',
                      help='Minimum PT for muons')

    parser.add_option('--maxMuonEta', type='float', action='store',
                      default=2.1,
                      dest='maxMuonEta',
                      help='Maximum muon pseudorapidity')

    parser.add_option('--minElectronPt', type='float', action='store',
                      default=30.,
                      dest='minElectronPt',
                      help='Minimum PT for electrons')

    parser.add_option('--maxElectronEta', type='float', action='store',
                      default=2.5,
                      dest='maxElectronEta',
                      help='Maximum electron pseudorapidity')


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

    (options, args) = parser.parse_args(argv)
    argv = []



    # Keep a quick dictionary for systematics
    NOM = 0
    JECUP = 1
    JECDN = 2
    JERUP = 3
    JERDN = 4

    systematics = {
        'NOM':NOM,
        'JECUP':JECUP,
        'JECDN':JECDN,
        'JERUP':JERUP,
        'JERDN':JERDN        
        }

    try : 
        syst = systematics[options.syst]
    except :
        print 'Invalid systematic option.'
    
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
    import copy

    muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
    electrons, electronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"
    photons, photonLabel = Handle("std::vector<pat::Photon>"), "slimmedPhotons"
    taus, tauLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
    jets, jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
    ak8jets, ak8jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"
    mets, metLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
    vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
    rhos, rhoLabel = Handle("double"), "fixedGridRhoAll"
    gens, genLabel = Handle("std::vector<reco::GenParticle>"), "prunedGenParticles"

    ##   ___ ___ .__          __                                             
    ##  /   |   \|__| _______/  |_  ____   ________________    _____   ______
    ## /    ~    \  |/  ___/\   __\/  _ \ / ___\_  __ \__  \  /     \ /  ___/
    ## \    Y    /  |\___ \  |  | (  <_> ) /_/  >  | \// __ \|  Y Y  \\___ \ 
    ##  \___|_  /|__/____  > |__|  \____/\___  /|__|  (____  /__|_|  /____  >
    ##        \/         \/             /_____/            \/      \/     \/

    f = ROOT.TFile(options.outname, "RECREATE")
    f.cd()

    # Actually to make life easy, we're going to make "N-dimensional histograms" aka Ntuples
    if options.writeTree : 
        TreeSemiLept = ROOT.TTree("TreeSemiLept", "TreeSemiLept")
        SemiLeptTrig        = array('i', [0]  )
        SemiLeptWeight      = array('f', [0.] )
        FatJetPt            = array('f', [-1.])
        FatJetEta           = array('f', [-1.])
        FatJetPhi           = array('f', [-1.])
        FatJetRap           = array('f', [-1.])
        FatJetEnergy        = array('f', [-1.])
        FatJetBDisc         = array('f', [-1.])
        FatJetMass          = array('f', [-1.])
        FatJetMassSoftDrop  = array('f', [-1.])
        FatJetTau32         = array('f', [-1.])
        FatJetTau21         = array('f', [-1.]) 
        FatJetSDbdiscW      = array('f', [-1.])
        FatJetSDbdiscB      = array('f', [-1.])
        FatJetSDsubjetWpt   = array('f', [-1.])
        FatJetSDsubjetWmass = array('f', [-1.])
        FatJetSDsubjetBpt   = array('f', [-1.])
        FatJetSDsubjetBmass = array('f', [-1.])
        LeptonType          = array('i', [-1])
        LeptonPt            = array('f', [-1.])
        LeptonEta           = array('f', [-1.])
        LeptonPhi           = array('f', [-1.])
        LeptonEnergy        = array('f', [-1.])
        LeptonIso           = array('f', [-1.])
        LeptonPtRel         = array('f', [-1.])
        LeptonDRMin         = array('f', [-1.])
        SemiLepMETpt        = array('f', [-1.])
        SemiLepMETphi       = array('f', [-1.])
        SemiLepNvtx         = array('f', [-1.])
        SemiLepTTmass       = array('f', [-1.])
        DeltaPhiLepFat      = array('f', [-1.]) 
        AK4bDisc            = array('f', [-1.])
        NearestAK4JetPt     = array('f', [-1.])
        NearestAK4JetEta    = array('f', [-1.])
        NearestAK4JetPhi    = array('f', [-1.])
        NearestAK4JetMass   = array('f', [-1.])
        SemiLeptRunNum        = array('f', [-1.])   
        SemiLeptLumiBlock     = array('f', [-1.])   
        SemiLeptEventNum      = array('f', [-1.])   


        TreeSemiLept.Branch('SemiLeptTrig'        , SemiLeptTrig        ,  'SemiLeptTrig/I'        )
        TreeSemiLept.Branch('SemiLeptWeight'      , SemiLeptWeight      ,  'SemiLeptWeight/F'      )    
        TreeSemiLept.Branch('FatJetPt'            , FatJetPt            ,  'FatJetPt/F'            )
        TreeSemiLept.Branch('FatJetEta'           , FatJetEta           ,  'FatJetEta/F'           )
        TreeSemiLept.Branch('FatJetPhi'           , FatJetPhi           ,  'FatJetPhi/F'           )
        TreeSemiLept.Branch('FatJetRap'           , FatJetRap           ,  'FatJetRap/F'           )
        TreeSemiLept.Branch('FatJetEnergy'        , FatJetEnergy        ,  'FatJetEnergy/F'        )
        TreeSemiLept.Branch('FatJetBDisc'         , FatJetBDisc         ,  'FatJetBDisc/F'         )
        TreeSemiLept.Branch('FatJetMass'          , FatJetMass          ,  'FatJetMass/F'          ) 
        TreeSemiLept.Branch('FatJetMassSoftDrop'  , FatJetMassSoftDrop  ,  'FatJetMassSoftDrop/F'  )
        TreeSemiLept.Branch('FatJetTau32'         , FatJetTau32         ,  'FatJetTau32/F'         )
        TreeSemiLept.Branch('FatJetTau21'         , FatJetTau21         ,  'FatJetTau21/F'         )
        TreeSemiLept.Branch('FatJetSDbdiscW'      , FatJetSDbdiscW      ,  'FatJetSDbdiscW/F'      )
        TreeSemiLept.Branch('FatJetSDbdiscB'      , FatJetSDbdiscB      ,  'FatJetSDbdiscB/F'      )        
        TreeSemiLept.Branch('FatJetSDsubjetWpt'   , FatJetSDsubjetWpt   ,  'FatJetSDsubjetWpt/F'   )
        TreeSemiLept.Branch('FatJetSDsubjetWmass' , FatJetSDsubjetWmass ,  'FatJetSDsubjetWmass/F' )
        TreeSemiLept.Branch('FatJetSDsubjetBpt'   , FatJetSDsubjetBpt   ,  'FatJetSDsubjetBpt/F'   )
        TreeSemiLept.Branch('FatJetSDsubjetBmass' , FatJetSDsubjetBmass ,  'FatJetSDsubjetBmass/F' )
        TreeSemiLept.Branch('LeptonType'          , LeptonType          ,  'LeptonType/I'          )
        TreeSemiLept.Branch('LeptonPt'            , LeptonPt            ,  'LeptonPt/F'            )
        TreeSemiLept.Branch('LeptonEta'           , LeptonEta           ,  'LeptonEta/F'           )
        TreeSemiLept.Branch('LeptonPhi'           , LeptonPhi           ,  'LeptonPhi/F'           )
        TreeSemiLept.Branch('LeptonEnergy'        , LeptonEnergy        ,  'LeptonEnergy/F'        )
        TreeSemiLept.Branch('LeptonIso'           , LeptonIso           ,  'LeptonIso/F'           )
        TreeSemiLept.Branch('LeptonPtRel'         , LeptonPtRel         ,  'LeptonPtRel/F'         )
        TreeSemiLept.Branch('LeptonDRMin'         , LeptonDRMin         ,  'LeptonDRMin/F'         )        
        TreeSemiLept.Branch('SemiLepMETpt'        , SemiLepMETpt        ,  'SemiLepMETpt/F'        )
        TreeSemiLept.Branch('SemiLepMETphi'       , SemiLepMETphi       ,  'SemiLepMETphi/F'       )
        TreeSemiLept.Branch('SemiLepNvtx'         , SemiLepNvtx         ,  'SemiLepNvtx/F'         )
        TreeSemiLept.Branch('SemiLepTTmass'       , SemiLepTTmass       ,  'SemiLepTTmass/F'      )
        TreeSemiLept.Branch('DeltaPhiLepFat'      , DeltaPhiLepFat      ,  'DeltaPhiLepFat/F'      )
        TreeSemiLept.Branch('AK4bDisc'            ,AK4bDisc             ,  'AK4bDisc/F'            )
        TreeSemiLept.Branch('NearestAK4JetPt'     ,NearestAK4JetPt      ,  'NearestAK4JetPt/F'     )
        TreeSemiLept.Branch('NearestAK4JetEta'    ,NearestAK4JetEta     ,  'NearestAK4JetEta/F'    )
        TreeSemiLept.Branch('NearestAK4JetPhi'    ,NearestAK4JetPhi     ,  'NearestAK4JetPhi/F'    )
        TreeSemiLept.Branch('NearestAK4JetMass'   ,NearestAK4JetMass    ,  'NearestAK4JetMass/F'   )
        TreeSemiLept.Branch('SemiLeptRunNum'         ,  SemiLeptRunNum       ,  'SemiLeptRunNum/F'          )
        TreeSemiLept.Branch('SemiLeptLumiBlock'      ,  SemiLeptLumiBlock    ,  'SemiLeptLumiBlock/F'       )
        TreeSemiLept.Branch('SemiLeptEventNum'       ,  SemiLeptEventNum     ,  'SemiLeptEventNum/F'        )

    

    # and also make a few 1-d histograms
    h_mttbar = ROOT.TH1F("h_mttbar", ";m_{t#bar{t}} (GeV)", 200, 0, 6000)
    h_mttbar_true = ROOT.TH1F("h_mttbar_true", "True m_{t#bar{t}};m_{t#bar{t}} (GeV)", 200, 0, 6000)

    h_ptLep = ROOT.TH1F("h_ptLep", "Lepton p_{T};p_{T} (GeV)", 100, 0, 1000)
    h_etaLep = ROOT.TH1F("h_etaLep", "Lepton #eta;p_{T} (GeV)#eta", 100, 0, ROOT.TMath.TwoPi() )
    h_met = ROOT.TH1F("h_met", "Missing p_{T};p_{T} (GeV)", 100, 0, 1000)
    h_ptRel = ROOT.TH1F("h_ptRel", "p_{T}^{REL};p_{T}^{REL} (GeV)", 100, 0, 100)
    h_dRMin = ROOT.TH1F("h_dRMin", "#Delta R_{MIN};#Delta R_{MIN}", 100, 0, 5.0)
    h_2DCut = ROOT.TH2F("h_2DCut", "2D Cut;#Delta R;p_{T}^{REL}", 20, 0, 5.0, 20, 0, 100 )

    h_ptAK4 = ROOT.TH1F("h_ptAK4", "AK4 Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
    h_etaAK4 = ROOT.TH1F("h_etaAK4", "AK4 Jet #eta;#eta", 120, -6, 6)
    h_yAK4 = ROOT.TH1F("h_yAK4", "AK4 Jet Rapidity;y", 120, -6, 6)
    h_phiAK4 = ROOT.TH1F("h_phiAK4", "AK4 Jet #phi;#phi (radians)",100,-ROOT.Math.Pi(),ROOT.Math.Pi())
    h_mAK4 = ROOT.TH1F("h_mAK4", "AK4 Jet Mass;Mass (GeV)", 100, 0, 1000)
    h_bdiscAK4 = ROOT.TH1F("h_bdiscAK4", "AK4 b discriminator;b discriminator", 100, 0, 1.0)

    h_ptAK8 = ROOT.TH1F("h_ptAK8", "AK8 Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
    h_etaAK8 = ROOT.TH1F("h_etaAK8", "AK8 Jet #eta;#eta", 120, -6, 6)
    h_yAK8 = ROOT.TH1F("h_yAK8", "AK8 Jet Rapidity;y", 120, -6, 6)
    h_phiAK8 = ROOT.TH1F("h_phiAK8", "AK8 Jet #phi;#phi (radians)",100,-ROOT.Math.Pi(),ROOT.Math.Pi())
    h_mAK8 = ROOT.TH1F("h_mAK8", "AK8 Jet Mass;Mass (GeV)", 100, 0, 1000)
    h_msoftdropAK8 = ROOT.TH1F("h_msoftdropAK8", "AK8 Softdrop Jet Mass;Mass (GeV)", 100, 0, 1000)
    h_mprunedAK8 = ROOT.TH1F("h_mprunedAK8", "AK8 Pruned Jet Mass;Mass (GeV)", 100, 0, 1000)
    h_mfilteredAK8 = ROOT.TH1F("h_mfilteredAK8", "AK8 Filtered Jet Mass;Mass (GeV)", 100, 0, 1000)
    h_mtrimmedAK8 = ROOT.TH1F("h_mtrimmedAK8", "AK8 Trimmed Jet Mass;Mass (GeV)", 100, 0, 1000)
    h_minmassAK8 = ROOT.TH1F("h_minmassAK8", "AK8 CMS Top Tagger Min Mass Paring;m_{min} (GeV)", 100, 0, 1000)
    h_nsjAK8 = ROOT.TH1F("h_nsjAK8", "AK8 CMS Top Tagger N_{subjets};N_{subjets}", 5, 0, 5)
    h_tau21AK8 = ROOT.TH1F("h_tau21AK8", "AK8 Jet #tau_{2} / #tau_{1};Mass#tau_{21}", 100, 0, 1.0)
    h_tau32AK8 = ROOT.TH1F("h_tau32AK8", "AK8 Jet #tau_{3} / #tau_{2};Mass#tau_{32}", 100, 0, 1.0)



    ##      ____.       __    _________                                     __  .__                      
    ##     |    | _____/  |_  \_   ___ \  __________________   ____   _____/  |_|__| ____   ____   ______
    ##     |    |/ __ \   __\ /    \  \/ /  _ \_  __ \_  __ \_/ __ \_/ ___\   __\  |/  _ \ /    \ /  ___/
    ## /\__|    \  ___/|  |   \     \___(  <_> )  | \/|  | \/\  ___/\  \___|  | |  (  <_> )   |  \\___ \ 
    ## \________|\___  >__|    \______  /\____/|__|   |__|    \___  >\___  >__| |__|\____/|___|  /____  >
    ##               \/               \/                          \/     \/                    \/     \/ 
    ROOT.gSystem.Load('libCondFormatsJetMETObjects')


    if options.isData :
        print 'Getting L2L3 for AK4'
        L2L3JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L2L3Residual_AK4PFchs.txt");
        print 'Getting L3 for AK4'
        L3JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L3Absolute_AK4PFchs.txt");
        print 'Getting L2 for AK4'
        L2JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L2Relative_AK4PFchs.txt");
        print 'Getting L1 for AK4'
        L1JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L1FastJet_AK4PFchs.txt");
        # for data only :
        #ResJetParAK4 = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L2L3Residual_AK4PFchs.txt");

        print 'Getting L2L3 for AK8'
        L2L3JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L2L3Residual_AK8PFchs.txt");    
        print 'Getting L3 for AK8'
        L3JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L3Absolute_AK8PFchs.txt");
        print 'Getting L2 for AK8'
        L2JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L2Relative_AK8PFchs.txt");
        print 'Getting L1 for AK8'
        L1JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L1FastJet_AK8PFchs.txt");
        # for data only :
        #ResJetParAK8 = ROOT.JetCorrectorParameters("JECs/74X_dataRun2_v5_L2L3Residual_AK8PFchs.txt");
        #  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
        vParJecAK4 = ROOT.vector('JetCorrectorParameters')()
        vParJecAK4.push_back(L1JetParAK4)
        vParJecAK4.push_back(L2JetParAK4)
        vParJecAK4.push_back(L3JetParAK4)
        vParJecAK4.push_back(L2L3JetParAK4)

        ak4JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK4)

        vParJecAK8 = ROOT.vector('JetCorrectorParameters')()
        vParJecAK8.push_back(L1JetParAK8)
        vParJecAK8.push_back(L2JetParAK8)
        vParJecAK8.push_back(L3JetParAK8)
        vParJecAK8.push_back(L2L3JetParAK8)

        ak8JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK8)


        jecParUncStrAK4 = ROOT.std.string('JECs/74X_mcRun2_asymptotic_v4_Uncertainty_AK4PFchs.txt')
        jecUncAK4 = ROOT.JetCorrectionUncertainty( jecParUncStrAK4 )
        jecParUncStrAK8 = ROOT.std.string('JECs/74X_mcRun2_asymptotic_v4_Uncertainty_AK8PFchs.txt')
        jecUncAK8 = ROOT.JetCorrectionUncertainty( jecParUncStrAK8 )    

        
    else :
        print 'Getting L3 for AK4'
        L3JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L3Absolute_AK4PFchs.txt");
        print 'Getting L2 for AK4'
        L2JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L2Relative_AK4PFchs.txt");
        print 'Getting L1 for AK4'
        L1JetParAK4  = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L1FastJet_AK4PFchs.txt");
        # for data only :
        #ResJetParAK4 = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L2L3Residual_AK4PFchs.txt");

        print 'Getting L3 for AK8'
        L3JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L3Absolute_AK8PFchs.txt");
        print 'Getting L2 for AK8'
        L2JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L2Relative_AK8PFchs.txt");
        print 'Getting L1 for AK8'
        L1JetParAK8  = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L1FastJet_AK8PFchs.txt");
        # for data only :
        #ResJetParAK8 = ROOT.JetCorrectorParameters("JECs/74X_mcRun2_asymptotic_v4_L2L3Residual_AK8PFchs.txt"); 

        #  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
        vParJecAK4 = ROOT.vector('JetCorrectorParameters')()
        vParJecAK4.push_back(L1JetParAK4)
        vParJecAK4.push_back(L2JetParAK4)
        vParJecAK4.push_back(L3JetParAK4)
        # for data only :
        #vParJecAK4.push_back(ResJetPar)

        ak4JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK4)

        vParJecAK8 = ROOT.vector('JetCorrectorParameters')()
        vParJecAK8.push_back(L1JetParAK8)
        vParJecAK8.push_back(L2JetParAK8)
        vParJecAK8.push_back(L3JetParAK8)
        # for data only :
        #vParJecAK8.push_back(ResJetPar)

        ak8JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK8)

        jecParUncStrAK4 = ROOT.std.string('JECs/74X_mcRun2_asymptotic_v4_Uncertainty_AK4PFchs.txt')
        jecUncAK4 = ROOT.JetCorrectionUncertainty( jecParUncStrAK4 )
        jecParUncStrAK8 = ROOT.std.string('JECs/74X_mcRun2_asymptotic_v4_Uncertainty_AK8PFchs.txt')
        jecUncAK8 = ROOT.JetCorrectionUncertainty( jecParUncStrAK8 )    

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
        for event in events:
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


            ##   ________                __________.__          __          
            ##  /  _____/  ____   ____   \______   \  |   _____/  |_  ______
            ## /   \  ____/ __ \ /    \   |     ___/  |  /  _ \   __\/  ___/
            ## \    \_\  \  ___/|   |  \  |    |   |  |_(  <_> )  |  \___ \ 
            ##  \______  /\___  >___|  /  |____|   |____/\____/|__| /____  >
            ##         \/     \/     \/                                  \/
            if not options.isData: 
                haveGenSolution = False
                isGenPresent = event.getByLabel( genLabel, gens )
                if isGenPresent : 
                    topQuark = None
                    antitopQuark = None
                    for igen,gen in enumerate( gens.product() ) :
                        #if options.verbose :
                        #    print 'GEN id=%.1f, pt=%+5.3f' % ( gen.pdgId(), gen.pt() )
                        if gen.pdgId() == 6 :
                            topQuark = gen
                        elif gen.pdgId() == -6 :
                            antitopQuark = gen

                    if topQuark != None and antitopQuark != None : 
                        ttbarCandP4 = topQuark.p4() + antitopQuark.p4()
                        h_mttbar_true.Fill( ttbarCandP4.mass() )
                        haveGenSolution = True
                    else :
                        if options.verbose :
                            print 'No top quarks, not filling mttbar'
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


            ## __________.__             ____   ____      .__                 
            ## \______   \  |__   ____   \   \ /   /____  |  |  __ __   ____  
            ##  |       _/  |  \ /  _ \   \   Y   /\__  \ |  | |  |  \_/ __ \ 
            ##  |    |   \   Y  (  <_> )   \     /  / __ \|  |_|  |  /\  ___/ 
            ##  |____|_  /___|  /\____/     \___/  (____  /____/____/  \___  >
            ##         \/     \/                        \/                 \/ 
            event.getByLabel(rhoLabel, rhos)
            # Rhos
            if len(rhos.product()) == 0 :
                print "Event has no rho values."
                continue
            else:
                rho = rhos.product()[0]
                if options.verbose : 
                    print 'rho = {0:6.2f}'.format( rho )


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
                    if muon.pt() > options.minMuonPt and abs(muon.eta()) < options.maxMuonEta and muon.muonBestTrack().dz(PV.position()) < 5.0 and muon.isTightMuon(PV) :
                        goodmuons.append( muon )
                        if options.verbose :
                            print "muon %2d: pt %4.1f, eta %+5.3f phi %+5.3f dz(PV) %+5.3f, POG loose id %d, tight id %d." % (
                                i, muon.pt(), muon.eta(), muon.phi(), muon.muonBestTrack().dz(PV.position()), muon.isLooseMuon(), muon.isTightMuon(PV))

            # Select tight good electrons
            goodelectrons = []
            if len(electrons.product()) > 0 :
                for i,electron in enumerate( electrons.product() ) :

                    if electron.pt() < electron.pt() and abs(electron.eta()) < options.maxElectronEta and electron.electronID("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight") == 5 :
                        goodelectrons.append( electron )
                        if options.verbose :
                            print "elec %2d: pt %4.1f, supercluster eta %+5.3f, phi %+5.3f sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), pass conv veto %d" % \
                            ( i, electron.pt(), electron.superCluster().eta(), electron.phi(), electron.sigmaIetaIeta(), electron.full5x5_sigmaIetaIeta(), electron.passConversionVeto())



            # Veto on dilepton events
            # Also keep track of the PF index of the lepton
            # for lepton-jet cleaning (see below)
            theLeptonObjKey = -1

            if len(goodmuons) + len(goodelectrons) != 1 :
                continue
            elif len(goodmuons) > 0 :
                theLeptonCand = goodmuons[0]
                theLepton = ROOT.TLorentzVector( goodmuons[0].px(),
                                                 goodmuons[0].py(),
                                                 goodmuons[0].pz(),
                                                 goodmuons[0].energy() )
                theLeptonObjKey = goodmuons[0].originalObjectRef().key()
                leptonType = 13
            else :
                theLeptonCand = goodelectrons[0]
                theLepton = ROOT.TLorentzVector( goodelectrons[0].px(),
                                                 goodelectrons[0].py(),
                                                 goodelectrons[0].pz(),
                                                 goodelectrons[0].energy() )
                theLeptonObjKey = goodelectrons[0].originalObjectRef().key()
                leptonType = 11

            # Get the "footprint" of the lepton. That is, all of the candidates making up the lepton.

            # now get a list of the PF candidates used to build this lepton, so to exclude them
            footprint = set()
            for i in xrange(theLeptonCand.numberOfSourceCandidatePtrs()):
                footprint.add(theLeptonCand.sourceCandidatePtr(i).key()) # the key is the index in the pf collection


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
            # In addition, we must perform "lepton-jet" cleaning.
            # This is because the PF leptons are actually counted in the
            # list of particles sent to the jet clustering.
            # Therefore, we need to loop over the jet constituents and
            # remove the lepton. 

            # use getByLabel, just like in cmsRun
            event.getByLabel (jetLabel, jets)          # For b-tagging
            event.getByLabel (ak8jetLabel, ak8jets)    # For top-tagging

            # loop over jets and fill hists
            ijet = 0

            # These will hold all of the jets we need for the selection
            ak4JetsGood = []
            ak8JetsGood = []
            ak4JetsGoodP4 = []
            ak8JetsGoodP4 = []


            # For selecting leptons, look at 2-d cut of dRMin, ptRel of
            # lepton and nearest jet that has pt > 30 GeV
            dRMin = 9999.0
            inearestJet = -1    # Index of nearest jet
            nearestJet = None   # Nearest jet


            ############################################
            # Get the AK4 jet nearest the lepton :
            ############################################ 
            for i,jet in enumerate(jets.product()) :
                # Get the jet p4
                jetP4Raw = ROOT.TLorentzVector( jet.px(), jet.py(), jet.pz(), jet.energy() )
                # Get the correction that was applied at RECO level for MINIADO
                jetJECFromMiniAOD = jet.jecFactor(0)
                # Remove the old JEC's to get raw energy
                jetP4Raw *= jetJECFromMiniAOD
                # Apply jet ID
                nhf = jet.neutralHadronEnergy() / jetP4Raw.E()
                nef = jet.neutralEmEnergy() / jetP4Raw.E()
                chf = jet.chargedHadronEnergy() / jetP4Raw.E()
                cef = jet.chargedEmEnergy() / jetP4Raw.E()
                nconstituents = jet.numberOfDaughters()
                nch = jet.chargedMultiplicity()
                goodJet = \
                  nhf < 0.99 and \
                  nef < 0.99 and \
                  chf > 0.00 and \
                  cef < 0.99 and \
                  nconstituents > 1 and \
                  nch > 0

                if not goodJet :
                    if options.verbose : 
                        print 'bad jet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}, bdisc = {4:6.2f}'.format (
                            jetP4Raw.Perp(), jetP4Raw.Rapidity(), jetP4Raw.Phi(), jetP4Raw.M(), jet.bDiscriminator( options.bdisc )
                            )
                    continue



                if options.verbose :
                    print 'raw jet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}, bdisc = {4:6.2f}'.format (
                        jetP4Raw.Perp(), jetP4Raw.Rapidity(), jetP4Raw.Phi(), jetP4Raw.M(), jet.bDiscriminator( options.bdisc )
                        )


                # Remove the lepton from the list of constituents for lepton/jet cleaning
                # Speed up computation, only do this for DR < 0.6
                cleaned = False
                if theLepton.DeltaR(jetP4Raw) < 0.6:
                    # Check all daughters of jets close to the lepton
                    pfcands = jet.daughterPtrVector()
                    for ipf,pf in enumerate( pfcands ) :

                        # If any of the jet daughters matches the good lepton, remove the lepton p4 from the jet p4
                        if pf.key() in footprint :
                            if options.verbose :
                                print 'REMOVING LEPTON, pt/eta/phi = {0:6.2f},{1:6.2f},{2:6.2f}'.format(
                                    theLepton.Perp(), theLepton.Eta(), theLepton.Phi()
                                    )
                            if jetP4Raw.Energy() > theLepton.Energy() : 
                                jetP4Raw -= theLepton
                            else :
                                jetP4Raw -= theLepton
                                jetP4Raw.SetEnergy(0.0)
                            cleaned = True
                            break



                # Apply new JEC's
                ak4JetCorrector.setJetEta( jetP4Raw.Eta() )
                ak4JetCorrector.setJetPt ( jetP4Raw.Perp() )
                ak4JetCorrector.setJetE  ( jetP4Raw.E() )
                ak4JetCorrector.setJetA  ( jet.jetArea() )
                ak4JetCorrector.setRho   ( rho )
                ak4JetCorrector.setNPV   ( NPV )            
                newJEC = ak4JetCorrector.getCorrection()
                

                # Get uncertainties
                jecUncAK4.setJetPhi(  jetP4Raw.Phi()  )
                jecUncAK4.setJetEta(  jetP4Raw.Eta()  )
                jecUncAK4.setJetPt(   jetP4Raw.Perp() * newJEC    )
                corrDn = newJEC - jecUncAK4.getUncertainty(0)
                jecUncAK4.setJetPhi(  jetP4Raw.Phi()  )
                jecUncAK4.setJetEta(  jetP4Raw.Eta()  )
                jecUncAK4.setJetPt(   jetP4Raw.Perp() * newJEC    )
                corrUp = newJEC + jecUncAK4.getUncertainty(1)



                # If MC, get jet energy resolution
                ptsmear   = 1.0
                ptsmearUp = 1.0
                ptsmearDn = 1.0
                if not options.isData:
                    # ---------------------------------------
                    # JER
                    # ---------------------------------------
                    
                    smear     = getJER( AK8P4Corr.Eta(),  0) 
                    smearUp   = getJER( AK8P4Corr.Eta(),  1) 
                    smearDn   = getJER( AK8P4Corr.Eta(), -1) 
                    recopt    = jetP4Raw.Perp() * newJEC
                    genpt     = jet.genJet().pt()
                    deltapt   = (recopt-genpt)*(smear-1.0)
                    deltaptUp = (recopt-genpt)*(smearUp-1.0)
                    deltaptDn = (recopt-genpt)*(smearDn-1.0)
                    ptsmear   = max(0.0, (recopt+deltapt)/recopt)
                    ptsmearUp = max(0.0, (recopt+deltaptUp)/recopt)
                    ptsmearDn = max(0.0, (recopt+deltaptDn)/recopt)



                if syst == NOM :
                    jetP4 = jetP4Raw * newJEC * ptsmear   # Nominal JEC, nominal JER
                elif syst == JECUP :
                    jetP4 = jetP4Raw * corrUp * ptsmear   # JEC up, nominal JER
                elif syst == JECDN :
                    jetP4 = jetP4Raw * corrDn * ptsmear   # JEC dn, nominal JER
                elif syst == JERUP :
                    jetP4 = jetP4Raw * newJEC * ptsmearUp # Nominal JEC, JER up
                elif syst == JERDN :
                    jetP4 = jetP4Raw * newJEC * ptsmearDn # Nominal JEC, JER dn
                else :
                    print "Invalid systematic option, using nominal. You should reconfigure or else I will talk a lot."
                    jetP4 = jetP4Raw * newJEC * ptsmear   # Nominal JEC, nominal JER
                    
                # Now perform jet kinematic cuts
                if jetP4.Perp() < options.minAK4Pt or abs(jetP4.Rapidity()) > options.maxAK4Rapidity :
                    continue

                # Get the jet nearest the lepton
                dR = jetP4.DeltaR(theLepton )
                ak4JetsGood.append(jet)
                ak4JetsGoodP4.append( jetP4 )
                if options.verbose :
                    print 'corrjet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}, bdisc = {4:6.2f}'.format (
                        jetP4.Perp(), jetP4.Rapidity(), jetP4.Phi(), jetP4.M(), jet.bDiscriminator( options.bdisc )
                        )

                if dR < dRMin :
                    inearestJet = ijet
                    nearestJet = jet
                    nearestJetP4 = jetP4
                    dRMin = dR


            ############################################
            # Require at least one leptonic-side jet, and 2d isolation cut
            ############################################ 
            if nearestJet == None :
                continue

            # Finally get the METs
            event.getByLabel( metLabel, mets )
            met = mets.product()[0]

            theLepJet = nearestJetP4
            theLepJetBDisc = nearestJet.bDiscriminator( options.bdisc )

            # Fill some plots related to the jets
            h_ptAK4.Fill( theLepJet.Perp() )
            h_etaAK4.Fill( theLepJet.Eta() )
            h_yAK4.Fill( theLepJet.Rapidity() )
            h_mAK4.Fill( theLepJet.M() )
            h_bdiscAK4.Fill( theLepJetBDisc )
            # Fill some plots related to the lepton, the MET, and the 2-d cut
            ptRel = theLepJet.Perp( theLepton.Vect() )
            h_ptLep.Fill(theLepton.Perp())
            h_etaLep.Fill(theLepton.Eta())
            h_met.Fill(met.pt())
            h_ptRel.Fill( ptRel )
            h_dRMin.Fill( dRMin )
            h_2DCut.Fill( dRMin, ptRel )
            pass2D = ptRel > 20.0 or dRMin > 0.4
            if options.verbose : 
                print '2d cut : dRMin = {0:6.2f}, ptRel = {1:6.2f}, pass = {2:6d}'.format( dRMin, ptRel, pass2D )
            if pass2D == False :
                continue

            ############################################
            # Get the AK8 jet away from the lepton
            ############################################
            for i,jet in enumerate(ak8jets.product()) :


                # perform jet ID with UNCORRECTED jet energy
                jetP4 = ROOT.TLorentzVector( jet.px(), jet.py(), jet.pz(), jet.energy() )
                jetP4Raw = copy.copy(jetP4)
                jetP4Raw *= jet.jecFactor(0)

                nhf = jet.neutralHadronEnergy() / jetP4Raw.E()
                nef = jet.neutralEmEnergy() / jetP4Raw.E()
                chf = jet.chargedHadronEnergy() / jetP4Raw.E()
                cef = jet.chargedEmEnergy() / jetP4Raw.E()
                nconstituents = jet.numberOfDaughters()
                nch = jet.chargedMultiplicity()
                goodJet = \
                  nhf < 0.99 and \
                  nef < 0.99 and \
                  chf > 0.00 and \
                  cef < 0.99 and \
                  nconstituents > 1 and \
                  nch > 0

                if not goodJet :
                    continue


                # Apply new JEC's
                ak8JetCorrector.setJetEta( jetP4Raw.Eta() )
                ak8JetCorrector.setJetPt ( jetP4Raw.Perp() )
                ak8JetCorrector.setJetE  ( jetP4Raw.E() )
                ak8JetCorrector.setJetA  ( jet.jetArea() )
                ak8JetCorrector.setRho   ( rho )
                ak8JetCorrector.setNPV   ( NPV )            
                newJEC = ak8JetCorrector.getCorrection()

                # Get uncertainties
                jecUncAK8.setJetPhi(  jetP4Raw.Phi()  )
                jecUncAK8.setJetEta(  jetP4Raw.Eta()  )
                jecUncAK8.setJetPt(   jetP4Raw.Perp() * newJEC   )
                corrDn = newJEC - jecUncAK8.getUncertainty(0)
                jecUncAK8.setJetPhi(  jetP4Raw.Phi()  )
                jecUncAK8.setJetEta(  jetP4Raw.Eta()  )
                jecUncAK8.setJetPt(   jetP4Raw.Perp() * newJEC   )
                corrUp = newJEC + jecUncAK8.getUncertainty(1)



                # If MC, get jet energy resolution
                ptsmear   = 1.0
                ptsmearUp = 1.0
                ptsmearDn = 1.0
                if not options.isData:
                    # ---------------------------------------
                    # JER
                    # ---------------------------------------
                    
                    smear     = getJER( AK8P4Corr.Eta(),  0) 
                    smearUp   = getJER( AK8P4Corr.Eta(),  1) 
                    smearDn   = getJER( AK8P4Corr.Eta(), -1) 
                    recopt    = jetP4Raw.Perp() * newJEC
                    genpt     = jet.genJet().pt()
                    deltapt   = (recopt-genpt)*(smear-1.0)
                    deltaptUp = (recopt-genpt)*(smearUp-1.0)
                    deltaptDn = (recopt-genpt)*(smearDn-1.0)
                    ptsmear   = max(0.0, (recopt+deltapt)/recopt)
                    ptsmearUp = max(0.0, (recopt+deltaptUp)/recopt)
                    ptsmearDn = max(0.0, (recopt+deltaptDn)/recopt)



                if syst == NOM :
                    jetP4 = jetP4Raw * newJEC * ptsmear   # Nominal JEC, nominal JER
                elif syst == JECUP :
                    jetP4 = jetP4Raw * corrUp * ptsmear   # JEC up, nominal JER
                elif syst == JECDN :
                    jetP4 = jetP4Raw * corrDn * ptsmear   # JEC dn, nominal JER
                elif syst == JERUP :
                    jetP4 = jetP4Raw * newJEC * ptsmearUp # Nominal JEC, JER up
                elif syst == JERDN :
                    jetP4 = jetP4Raw * newJEC * ptsmearDn # Nominal JEC, JER dn
                else :
                    print "Invalid systematic option, using nominal. You should reconfigure or else I will talk a lot."
                    jetP4 = jetP4Raw * newJEC * ptsmear   # Nominal JEC, nominal JER

                
                # Now perform jet kinematic cuts
                if jetP4.Perp() < options.minAK8Pt or abs(jetP4.Rapidity()) > options.maxAK8Rapidity :
                    continue

                # Only keep AK8 jets "away" from the lepton, so we do not need
                # lepton-jet cleaning here. There's no double counting. 
                dR = jetP4.DeltaR(theLepton )
                if dR > ROOT.TMath.Pi()/2.0 :
                    ak8JetsGood.append(jet)
                    ak8JetsGoodP4.append( jetP4 )

            ## ___________                     .__                
            ## \__    ___/____     ____   ____ |__| ____    ____  
            ##   |    |  \__  \   / ___\ / ___\|  |/    \  / ___\ 
            ##   |    |   / __ \_/ /_/  > /_/  >  |   |  \/ /_/  >
            ##   |____|  (____  /\___  /\___  /|__|___|  /\___  / 
            ##                \//_____//_____/         \//_____/  

            ############################################
            # Investigate the b-tagging and t-tagging
            ############################################
            if len(ak4JetsGoodP4) < 1 or len(ak8JetsGoodP4) < 1 :
                continue


            tJets = []
            for ijet,jet in enumerate(ak8JetsGood) : 
                if jet.pt() < options.minAK8Pt :
                    continue

                mAK8Softdrop = jet.userFloat('ak8PFJetsCHSSoftDropMass')
                mAK8Pruned = jet.userFloat('ak8PFJetsCHSPrunedMass')
                mAK8Filtered = jet.userFloat('ak8PFJetsCHSFilteredMass')
                mAK8Trimmed = jet.userFloat('ak8PFJetsCHSTrimmedMass')


                h_ptAK8.Fill( jet.pt() )
                h_etaAK8.Fill( jet.eta() )
                h_yAK8.Fill( jet.rapidity() )
                h_mAK8.Fill( jet.mass() )
                h_msoftdropAK8.Fill( mAK8Softdrop )
                h_mprunedAK8.Fill( mAK8Pruned )
                h_mfilteredAK8.Fill( mAK8Filtered )
                h_mtrimmedAK8.Fill( mAK8Trimmed )



                tJets.append( jet )



            ##  ____  __.__                              __  .__         __________                     
            ## |    |/ _|__| ____   ____   _____ _____ _/  |_|__| ____   \______   \ ____   ____  ____  
            ## |      < |  |/    \_/ __ \ /     \\__  \\   __\  |/ ___\   |       _// __ \_/ ___\/  _ \ 
            ## |    |  \|  |   |  \  ___/|  Y Y  \/ __ \|  | |  \  \___   |    |   \  ___/\  \__(  <_> )
            ## |____|__ \__|___|  /\___  >__|_|  (____  /__| |__|\___  >  |____|_  /\___  >\___  >____/ 
            ##         \/       \/     \/      \/     \/             \/          \/     \/     \/       

            # Now we do our kinematic calculation based on the categories of the
            # number of top and bottom tags
            mttbar = -1.0
            if len(tJets) == 0 :
                if options.verbose : 
                    print 'No top tags'
            else :


                hadTopCandP4 = ROOT.TLorentzVector( tJets[0].px(), tJets[0].py(), tJets[0].pz(), tJets[0].energy() )
                lepTopCandP4 = None
                # Get the z-component of the lepton from the W mass constraint
                bJetCandP4 = ak4JetsGoodP4[inearestJet]
                nuCandP4 = ROOT.TLorentzVector( met.px(), met.py(), 0, met.energy() )


                solution, nuz1, nuz2 = solve_nu( vlep=theLepton, vnu=nuCandP4 )
                # If there is at least one real solution, pick it up
                if solution :
                    if options.verbose : 
                        print '--- Have a solution --- '
                    nuCandP4.SetPz( nuz1 )
                else :
                    if options.verbose : 
                        print '--- No solution for neutrino z ---'
                    nuCandP4.SetPz( nuz1.real )

                lepTopCandP4 = nuCandP4 + theLepton + bJetCandP4

                ttbarCand = hadTopCandP4 + lepTopCandP4
                mttbar = ttbarCand.M()
                h_mttbar.Fill( mttbar )
                if not options.isData and  haveGenSolution == False :
                    print 'Very strange. No gen solution, but it is a perfectly good event. mttbar = ' + str(ttbarCand.M() )


            ## ___________.__.__  .__    ___________                      
            ## \_   _____/|__|  | |  |   \__    ___/______   ____   ____  
            ##  |    __)  |  |  | |  |     |    |  \_  __ \_/ __ \_/ __ \ 
            ##  |     \   |  |  |_|  |__   |    |   |  | \/\  ___/\  ___/ 
            ##  \___  /   |__|____/____/   |____|   |__|    \___  >\___  >
            ##      \/                                          \/     \/ 
            if options.writeTree :
                candToPlot = 0

                # Make sure there are top tags if we want to plot them               
                tagInfoLabels = ak8JetsGood[candToPlot].tagInfoLabels()
                # Get n-subjettiness "tau" variables
                tau1 = ak8JetsGood[candToPlot].userFloat('NjettinessAK8:tau1')
                tau2 = ak8JetsGood[candToPlot].userFloat('NjettinessAK8:tau2')
                tau3 = ak8JetsGood[candToPlot].userFloat('NjettinessAK8:tau3')
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

                # Get the subjets from the modified mass drop algorithm
                # aka softdrop with beta=0.
                # The heaviest should correspond to the W and the lightest
                # should correspond to the b. 
                subjets = ak8JetsGood[candToPlot].subjets('SoftDrop')
                subjetW = None
                subjetB = None
                if len(subjets) >= 2 : 
                    if subjets[0].mass() > subjets[1].mass() :
                        subjetW = subjets[0]
                        subjetB = subjets[1]                    
                    else :
                        subjetB = subjets[0]
                        subjetW = subjets[1]
                else :
                    continue 
                
                
                SemiLeptTrig        [0] = 0
                SemiLeptWeight      [0] = evWeight
                FatJetPt            [0] = ak8JetsGoodP4[candToPlot].Perp()
                FatJetEta           [0] = ak8JetsGoodP4[candToPlot].Eta()
                FatJetPhi           [0] = ak8JetsGoodP4[candToPlot].Phi()
                FatJetRap           [0] = ak8JetsGoodP4[candToPlot].Rapidity()
                FatJetEnergy        [0] = ak8JetsGoodP4[candToPlot].Energy()
                FatJetMass          [0] = ak8JetsGoodP4[candToPlot].M()
                FatJetBDisc         [0] = ak8JetsGood[candToPlot].bDiscriminator(options.bdisc)
                FatJetMassSoftDrop  [0] = ak8JetsGood[candToPlot].userFloat('ak8PFJetsCHSSoftDropMass')
                FatJetTau32         [0] = tau32
                FatJetTau21         [0] = tau21
                if subjetW != None : 
                    FatJetSDbdiscW      [0] = subjetW.bDiscriminator(options.bdisc)
                    FatJetSDbdiscB      [0] = subjetB.bDiscriminator(options.bdisc)
                    FatJetSDsubjetWpt   [0] = subjetW.pt()
                    FatJetSDsubjetWmass [0] = subjetW.mass()
                    FatJetSDsubjetBpt   [0] = subjetB.pt()
                    FatJetSDsubjetBmass [0] = subjetB.mass()
                LeptonType          [0] = leptonType
                LeptonPt            [0] = theLepton.Perp()  
                LeptonEta           [0] = theLepton.Eta()
                LeptonPhi           [0] = theLepton.Phi()
                LeptonEnergy        [0] = theLepton.E()
                LeptonPtRel         [0] = nearestJetP4.Perp(theLepton.Vect())
                LeptonDRMin         [0] = nearestJetP4.DeltaR(theLepton)
                SemiLepMETpt        [0] = met.pt()
                SemiLepMETphi       [0] = met.phi()   
                SemiLepNvtx         [0] = NPV
                SemiLepTTmass       [0] = mttbar
                DeltaPhiLepFat      [0] = ak8JetsGoodP4[candToPlot].DeltaR(theLepton)
                AK4bDisc            [0] = nearestJet.bDiscriminator(options.bdisc)
                NearestAK4JetPt     [0] = nearestJetP4.Perp()
                NearestAK4JetEta    [0] = nearestJetP4.Eta()
                NearestAK4JetPhi    [0] = nearestJetP4.Phi()
                NearestAK4JetMass   [0] = nearestJetP4.M()
                SemiLeptRunNum      [0] = event.object().id().run()
                SemiLeptLumiBlock   [0] = event.object().luminosityBlock()
                SemiLeptEventNum    [0] = event.object().id().event()

                TreeSemiLept.Fill()                        

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
    b2gdas_fwlite(sys.argv)


