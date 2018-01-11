#! /usr/bin/env python


## _________                _____.__                            __  .__               
## \_   ___ \  ____   _____/ ____\__| ____  __ ______________ _/  |_|__| ____   ____  
## /    \  \/ /  _ \ /    \   __\|  |/ ___\|  |  \_  __ \__  \\   __\  |/  _ \ /    \ 
## \     \___(  <_> )   |  \  |  |  / /_/  >  |  /|  | \// __ \|  | |  (  <_> )   |  \
##  \______  /\____/|___|  /__|  |__\___  /|____/ |__|  (____  /__| |__|\____/|___|  /
##         \/            \/        /_____/                   \/                    \/ 
import sys
import array as array
from optparse import OptionParser
import numpy as np


def plot_mttbar(argv) :
    parser = OptionParser()

    parser.add_option('--file_in', type='string', action='store',
                      dest='file_in',
                      help='Input file')

    parser.add_option('--file_out', type='string', action='store',
                      dest='file_out',
                      help='Output file')

    parser.add_option('--jer', type='string', action='store',
                      dest='jer',
                      default = None,
                      help='Choice of up and down for Jet Energy Resolution (only set jer or jec')

    parser.add_option('--jec', type='string', action='store',
                      dest='jec',
                      default = None,
                      help='Choice of up and down for Jet Energy Correction (only set jer or jec')
    
    parser.add_option('--lepton', type='string', action='store',
	                  dest='lepton',
					  default='mu',
                      help='Choice of lepton (mu or ele)')
					  
    parser.add_option('--isData', action='store_true',
                      dest='isData',
                      default = False,
                      help='Is this Data?')
        
    (options, args) = parser.parse_args(argv)
    argv = []

    #write to temp file
    #fh = open("num.txt", "a")

    #print '===== Command line options ====='
    #print options
    #print '================================'

    import ROOT

    from leptonic_nu_z_component import solve_nu_tmass, solve_nu  #load Z_momentum with these functions

    histogramSuffix = ''
    #Adding name change to histograms based on --jec options and --jer options
    if options.lepton == 'ele':
        histogramSuffix = '_el'
    if options.lepton == 'mu':
        histogramSuffix = '_mu'

    if options.jec is not None and options.jer is not None:
        print 'You are trying to do two systematics at once! Please fix your input options'

    if options.jec == 'up' :
        histogramSuffix += '_jec_Up'

    if options.jec == 'down' :
	    histogramSuffix += '_jec_Down'

    if options.jer == 'up' :
        histogramSuffix += '_jer_Up'

    if options.jer == 'down' :
	    histogramSuffix += '_jer_Down'

    fpileup = ROOT.TFile.Open('purw.root', 'read')
    h_pileupWeight = fpileup.Get('pileup') 

    fout= ROOT.TFile(options.file_out, "RECREATE")
    fout.cd()
    h_cuts = ROOT.TH1F("Cut_flow", "", 4,0,4)

    h_mttbar = ROOT.TH1F("h_mttbar"+histogramSuffix, ";m_{t#bar{t}} (GeV);Number", 100, 0, 5000)#invariant ttbar mass
    h_mttbar_control = ROOT.TH1F("h_mttbar_control"+histogramSuffix, ";m_{t#bar{t}} (GeV);Number", 100, 0, 5000)#invariant ttbar mass
    h_mtopHad = ROOT.TH1F("h_mtopHad"+histogramSuffix, ";m_{jet} (GeV);Number", 100, 0, 400)
    h_mtopHadGroomed = ROOT.TH1F("h_mtopHadGroomed"+histogramSuffix, ";Groomed m_{jet} (GeV);Number", 100, 0, 400)

    #ele plots
    h_lepPt  = ROOT.TH1F("h_lepPt"+histogramSuffix, "; lep_{pt}(GeV);Number", 100, 0, 1200)
    h_lepEta = ROOT.TH1F("h_lepEta"+histogramSuffix, ";lep_{#eta};Number", 100, -2.5, 2.5)
    h_lepPhi = ROOT.TH1F("h_lepPhi"+histogramSuffix, ";lep_{#phi};Number", 100, -3.5, 3.5)

    #AK8

    h_AK8Pt = ROOT.TH1F("h_AK8Pt"+histogramSuffix  , ";AK8_{pt} (GeV);Number", 100, 400, 2500)
    h_AK8Eta = ROOT.TH1F("h_AK8Eta"+histogramSuffix, ";AK8_{#eta} ;Number", 100, -2.5, 2.5)
    h_AK8Phi = ROOT.TH1F("h_AK8Phi"+histogramSuffix, ";AK8_{#phi} ;Number", 100, -3.5, 3.5)
    h_AK8Tau32 = ROOT.TH1F("h_AK8Tau32"+histogramSuffix,";AK8_{#tau_{32}};Number", 50, 0.0, 1.0)
    h_AK8Tau32PreSel = ROOT.TH1F("h_AK8Tau32PreSel"+histogramSuffix,";AK8_{#tau_{32}};Number", 50, 0.0, 1.0)
    h_AK8Tau21 = ROOT.TH1F("h_AK8Tau21"+histogramSuffix,";AK8_{#tau_{21}};Number", 50, 0.0, 1.0)
    h_AK8Tau21PreSel = ROOT.TH1F("h_AK8Tau21PreSel"+histogramSuffix,";AK8_{#tau_{21}};Number", 50, 0.0, 1.0)
    
	#AK8
    h_AK4Pt    = ROOT.TH1F("h_AK4Pt"+histogramSuffix,";ak4jet_{pT} (GeV);Number", 100, 0, 1500)
    h_AK4Eta   = ROOT.TH1F("h_AK4Eta"+histogramSuffix,";ak4jet_{#eta};Number", 100, -2.5, 2.5)
    h_AK4Phi   = ROOT.TH1F("h_AK4Phi"+histogramSuffix,";ak4jet_{#phi};Number", 100, -3.5, 3.5)
    h_AK4M     = ROOT.TH1F("h_AK4M"+histogramSuffix,";ak4jet_{mass};Number", 100, 0, 400)
    h_AK4Bdisc = ROOT.TH1F("h_AK4Bdisc"+histogramSuffix,";ak4jet_{bdisc};Number", 100, 0, 1.0)
    h_AK4BdiscPreSel = ROOT.TH1F("h_AK4BdiscPreSel"+histogramSuffix,";ak4jet_{bdisc};Number", 100, 0, 1.0)


    h_drAK4AK8    = ROOT.TH1F("h_drAK4AK8"+histogramSuffix,";#DeltaR_{AK4, AK8} ;Number", 100, 0, 5)
#    h_drLepAK8    = ROOT.TH1F("h_drLepAK8",";{#delta r}_{lep, AK8} ;Number", 100, 0, 1500)
    h_drLepAK4    = ROOT.TH1F("h_drLepAK4"+histogramSuffix,";#DeltaR_{lep, AK4} ;Number", 100, 0, 5)
#    h_dPhiLepAK8 = ROOT.TH1F("h_dPhiLepAK8"+histogramSuffix,";#Delta#phi_{l,AK8};Number", 100, 0.0, 1.0) #Not actually filled in any of the ntuples

    #Following lines is to make sure that the statistical errors are kept and stored
    h_mttbar.Sumw2()
    h_mtopHad.Sumw2()
    h_mtopHadGroomed.Sumw2()
    h_lepPt.Sumw2()
    h_lepEta.Sumw2()
    h_lepPhi.Sumw2()
    h_AK8Pt.Sumw2()
    h_AK8Eta.Sumw2()
    h_AK8Phi.Sumw2()
    h_AK8Tau32.Sumw2()
    h_AK8Tau21.Sumw2()
    h_AK4Pt.Sumw2()
    h_AK4Eta.Sumw2()
    h_AK4Phi.Sumw2()
    h_AK4M.Sumw2()
    h_AK4Bdisc.Sumw2()
    h_drAK4AK8.Sumw2()
    h_drLepAK4.Sumw2()
    
    fin = ROOT.TFile.Open(options.file_in)

    trees = [ fin.Get("TreeSemiLept") ]

    tot_entries, count = 0, 0
    cut1, cut2, cut3, cut4 = 0 ,0 ,0 ,0

    for itree,t in enumerate(trees) :

        #if options.isData : 
        SemiLeptTrig        =  ROOT.vector('int')()
        SemiLeptWeight      = array.array('f', [0.] )
        PUWeight            = array.array('f', [0.] )
        GenWeight           = array.array('f', [0.] )
        FatJetPt            = array.array('f', [-1.])
        FatJetEta           = array.array('f', [-1.])
        FatJetPhi           = array.array('f', [-1.])
        FatJetRap           = array.array('f', [-1.])
        FatJetEnergy        = array.array('f', [-1.])
        FatJetBDisc         = array.array('f', [-1.])
        FatJetMass          = array.array('f', [-1.])
        FatJetMassSoftDrop  = array.array('f', [-1.])
        FatJetTau32         = array.array('f', [-1.])
        FatJetTau21         = array.array('f', [-1.]) 
        FatJetSDBDiscW      = array.array('f', [-1.])
        FatJetSDBDiscB      = array.array('f', [-1.])
        FatJetSDsubjetWpt   = array.array('f', [-1.])
        FatJetSDsubjetWmass = array.array('f', [-1.])
        FatJetSDsubjetBpt   = array.array('f', [-1.])
        FatJetSDsubjetBmass = array.array('f', [-1.])
        FatJetJECUpSys      = array.array('f', [-1.])
        FatJetJECDnSys      = array.array('f', [-1.])
        FatJetJERUpSys      = array.array('f', [-1.])
        FatJetJERDnSys      = array.array('f', [-1.])
        LeptonType          = array.array('i', [-1])
        LeptonPt            = array.array('f', [-1.])
        LeptonEta           = array.array('f', [-1.])
        LeptonPhi           = array.array('f', [-1.])
        LeptonEnergy        = array.array('f', [-1.])
        LeptonIso           = array.array('f', [-1.])
        LeptonPtRel         = array.array('f', [-1.])
        LeptonDRMin         = array.array('f', [-1.])
        SemiLepMETpt        = array.array('f', [-1.])
        SemiLepMETphi       = array.array('f', [-1.])
        SemiLepNvtx         = array.array('i', [0])
        FatJetDeltaPhiLep   = array.array('f', [-1.]) 
        NearestAK4JetBDisc  = array.array('f', [-1.])
        NearestAK4JetPt     = array.array('f', [-1.])
        NearestAK4JetEta    = array.array('f', [-1.])
        NearestAK4JetPhi    = array.array('f', [-1.])
        NearestAK4JetMass   = array.array('f', [-1.])
        NearestAK4JetJECUpSys = array.array('f', [-1.])
        NearestAK4JetJECDnSys = array.array('f', [-1.])
        NearestAK4JetJERUpSys = array.array('f', [-1.])
        NearestAK4JetJERDnSys = array.array('f', [-1.])
        SemiLeptRunNum        = array.array('f', [-1.])   
        SemiLeptLumiNum     = array.array('f', [-1.])   
        SemiLeptEventNum      = array.array('f', [-1.])   


        #if options.isData : 
        t.SetBranchAddress('SemiLeptTrig'        , SemiLeptTrig )
        t.SetBranchAddress('SemiLeptWeight'      , SemiLeptWeight      ) #Combined weight of all scale factors (lepton, PU, generator) relevant for the smeileptonic event selection
        t.SetBranchAddress('PUWeight'            , PUWeight            )
        t.SetBranchAddress('GenWeight'           , GenWeight               )
        t.SetBranchAddress('FatJetPt'            , FatJetPt            )
        t.SetBranchAddress('FatJetEta'           , FatJetEta           )
        t.SetBranchAddress('FatJetPhi'           , FatJetPhi           )
        t.SetBranchAddress('FatJetRap'           , FatJetRap           )
        t.SetBranchAddress('FatJetEnergy'        , FatJetEnergy        )
        t.SetBranchAddress('FatJetBDisc'         , FatJetBDisc         )
        t.SetBranchAddress('FatJetMass'          , FatJetMass           )
        t.SetBranchAddress('FatJetMassSoftDrop'  , FatJetMassSoftDrop  )
        t.SetBranchAddress('FatJetTau32'         , FatJetTau32         )
        t.SetBranchAddress('FatJetTau21'         , FatJetTau21         )
        t.SetBranchAddress('FatJetSDBDiscW'      , FatJetSDBDiscW      )
        t.SetBranchAddress('FatJetSDBDiscB'      , FatJetSDBDiscB              )
        t.SetBranchAddress('FatJetSDsubjetWpt'   , FatJetSDsubjetWpt   )
        t.SetBranchAddress('FatJetSDsubjetWmass' , FatJetSDsubjetWmass )
        t.SetBranchAddress('FatJetSDsubjetBpt'   , FatJetSDsubjetBpt   )
        t.SetBranchAddress('FatJetSDsubjetBmass' , FatJetSDsubjetBmass )
        t.SetBranchAddress('FatJetJECUpSys'      , FatJetJECUpSys      )
        t.SetBranchAddress('FatJetJECDnSys'      , FatJetJECDnSys      )
        t.SetBranchAddress('FatJetJERUpSys'      , FatJetJERUpSys      )
        t.SetBranchAddress('FatJetJERDnSys'      , FatJetJERDnSys      )
        t.SetBranchAddress('LeptonType'          , LeptonType          )
        t.SetBranchAddress('LeptonPt'            , LeptonPt            )
        t.SetBranchAddress('LeptonEta'           , LeptonEta           )
        t.SetBranchAddress('LeptonPhi'           , LeptonPhi           )
        t.SetBranchAddress('LeptonEnergy'        , LeptonEnergy        )
        t.SetBranchAddress('LeptonIso'           , LeptonIso           )
        t.SetBranchAddress('LeptonPtRel'         , LeptonPtRel         )
        t.SetBranchAddress('LeptonDRMin'         , LeptonDRMin         )
        t.SetBranchAddress('SemiLepMETpt'        , SemiLepMETpt        )
        t.SetBranchAddress('SemiLepMETphi'       , SemiLepMETphi       )
        t.SetBranchAddress('SemiLepNvtx'         , SemiLepNvtx         )
        t.SetBranchAddress('FatJetDeltaPhiLep'      , FatJetDeltaPhiLep      )
        t.SetBranchAddress('NearestAK4JetBDisc'            ,NearestAK4JetBDisc             )
        t.SetBranchAddress('NearestAK4JetPt'     ,NearestAK4JetPt      )
        t.SetBranchAddress('NearestAK4JetEta'    ,NearestAK4JetEta     )
        t.SetBranchAddress('NearestAK4JetPhi'    ,NearestAK4JetPhi     )
        t.SetBranchAddress('NearestAK4JetMass'   ,NearestAK4JetMass    )
        t.SetBranchAddress('NearestAK4JetJECUpSys'      , NearestAK4JetJECUpSys)
        t.SetBranchAddress('NearestAK4JetJECDnSys'      , NearestAK4JetJECDnSys)
        t.SetBranchAddress('NearestAK4JetJERUpSys'      , NearestAK4JetJERUpSys)
        t.SetBranchAddress('NearestAK4JetJERDnSys'      , NearestAK4JetJERDnSys)
        t.SetBranchAddress('SemiLeptRunNum'         ,  SemiLeptRunNum       )
        t.SetBranchAddress('SemiLeptLumiNum'      ,  SemiLeptLumiNum    )
        t.SetBranchAddress('SemiLeptEventNum'       ,  SemiLeptEventNum     )
        


        t.SetBranchStatus ('*', 0)
        t.SetBranchStatus ('SemiLeptWeight', 1)
        t.SetBranchStatus ('PUWeight', 1)
        t.SetBranchStatus ('GenWeight', 1)
        t.SetBranchStatus ('FatJetPt', 1)
        t.SetBranchStatus ('FatJetEta', 1)
        t.SetBranchStatus ('FatJetPhi', 1)
        t.SetBranchStatus ('FatJetMass', 1)
        t.SetBranchStatus ('FatJetMassSoftDrop', 1)
        t.SetBranchStatus ('FatJetTau32', 1)
        t.SetBranchStatus ('FatJetTau21', 1)
        t.SetBranchStatus ('SemiLeptTrig', 1)
        t.SetBranchStatus ('NearestAK4JetBDisc', 1)
        t.SetBranchStatus ('NearestAK4JetPt'   ,1 )
        t.SetBranchStatus ('NearestAK4JetEta'  ,1 )
        t.SetBranchStatus ('NearestAK4JetPhi'  ,1 )
        t.SetBranchStatus ('NearestAK4JetMass' ,1 )
        t.SetBranchStatus ('SemiLepMETpt' , 1 )
        t.SetBranchStatus ('SemiLepMETphi' , 1 )
        t.SetBranchStatus ('LeptonType'          , 1 )
        t.SetBranchStatus ('LeptonPt'            , 1)
        t.SetBranchStatus ('LeptonEta'           , 1)
        t.SetBranchStatus ('LeptonPhi'           , 1)
        t.SetBranchStatus ('LeptonEnergy'        , 1)
        t.SetBranchStatus ('LeptonIso'           , 1)
        t.SetBranchStatus ('LeptonPtRel'         , 1)
        t.SetBranchStatus ('LeptonDRMin'         , 1)
        t.SetBranchStatus ('FatJetJECUpSys'      , 1)
        t.SetBranchStatus ('FatJetJECDnSys'      , 1)
        t.SetBranchStatus ('FatJetJERUpSys'      , 1)
        t.SetBranchStatus ('FatJetJERDnSys'      , 1)
        t.SetBranchStatus ('NearestAK4JetJECUpSys', 1)
        t.SetBranchStatus ('NearestAK4JetJECDnSys' , 1)
        t.SetBranchStatus ('NearestAK4JetJERUpSys' , 1)
        t.SetBranchStatus ('NearestAK4JetJERDnSys' , 1)
        t.SetBranchStatus ('SemiLepNvtx' , 1)
        
        entries = t.GetEntriesFast()
        tot_entries +=entries

        #print 'Processing tree ' + str(itree)

        eventsToRun = entries
        for jentry in xrange( eventsToRun ):
            #if jentry % 100000 == 0 :
            #    print 'processing ' + str(jentry)
            # get the next tree in the chain and verify
            ientry = t.GetEntry( jentry )
            if ientry < 0:
                break
			#Triggers that are available: 
			#   0: "HLT_Mu50"
			#   1: "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"
			#   2: "HLT_Ele115_CaloIdVT_GsfTrkIdT"
			#   3: "HLT_PFHT800"

            if options.lepton=='mu' :

                if LeptonType[0] != 13 :
                    continue

                # Muon triggers only for now (use HLT_Mu50 with index 0)

                if SemiLeptTrig[0] == 0  :
                    continue

            if options.lepton=='ele' :
    
                if LeptonType[0] != 11 :
                    continue


                # Muon triggers only for now (use HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 with index 1 and HLT_Ele115_CaloIdVT_GsfTrkIdT with index 2)
                if SemiLeptTrig[1] == 0 and SemiLeptTrig[2] == 0 :
                    continue

            # Hadronic top
            hadTopCandP4 = ROOT.TLorentzVector()
            hadTopCandP4.SetPtEtaPhiM( FatJetPt[0], FatJetEta[0], FatJetPhi[0], FatJetMass[0])#set up with lead Ak8 jet in event
            bJetCandP4 = ROOT.TLorentzVector()
            bJetCandP4.SetPtEtaPhiM( NearestAK4JetPt[0], NearestAK4JetEta[0], NearestAK4JetPhi[0], NearestAK4JetMass[0])
            # MET
            nuCandP4 = ROOT.TLorentzVector( )
            nuCandP4.SetPtEtaPhiM( SemiLepMETpt[0], 0, SemiLepMETphi[0], SemiLepMETpt[0] )
            # Leptoon
            theLepton = ROOT.TLorentzVector()
            theLepton.SetPtEtaPhiE( LeptonPt[0], LeptonEta[0], LeptonPhi[0], LeptonEnergy[0] ) # Assume massless

            tau32 = FatJetTau32[0]
            mass_sd = FatJetMassSoftDrop[0]
            bdisc = NearestAK4JetBDisc[0]
            pileupWeight=  h_pileupWeight.GetBinContent(SemiLepNvtx[0]+1)
            #Weights
            weight = pileupWeight
            if options.jec =='up':
                weight = 1*NearestAK4JetJECUpSys[0]*FatJetJECUpSys[0]
            if options.jec =='down':
                weight = 1*NearestAK4JetJECDnSys[0]*FatJetJECDnSys[0]
            if options.jer =='up':
                weight = 1*NearestAK4JetJERUpSys[0]*FatJetJECUpSys[0]
            if options.jer =='down':
                weight = 1*NearestAK4JetJERDnSys[0]*FatJetJECDnSys[0]

            #print weight
            
            #preselection histos            
            h_AK4BdiscPreSel.Fill( NearestAK4JetBDisc[0], weight )
            h_AK8Tau32PreSel.Fill(FatJetTau32[0], weight )
            h_AK8Tau21PreSel.Fill(FatJetTau21[0], weight )

            passKin = hadTopCandP4.Perp() > 100.
            passTopTag = tau32 < 0.8 and mass_sd > 110. and mass_sd < 250.
            pass2DCut = LeptonPtRel[0] > 20. or LeptonDRMin[0] > 0.4
            passBtag = bdisc > 0.7

            # Invariant mass calculation
            def calculate_m():
                lepTopCandP4 = None
                # Get the z-component of the lepton from the W mass constraint
                solution, nuz1, nuz2 = solve_nu( vlep=theLepton, vnu=nuCandP4 )
                # If there is at least one real solution, pick it up
                if solution :
                    nuCandP4.SetPz( nuz1 )
                else :
                    nuCandP4.SetPz( nuz1.real )

                lepTopCandP4 = nuCandP4 + theLepton + bJetCandP4

                ttbarCand = hadTopCandP4 + lepTopCandP4
                mttbar = ttbarCand.M()
                return mttbar

            # Applying and counting cuts
            if not passKin: 
                continue
            else:
                cut1 +=1
            if not pass2DCut: 
                continue
            else:
                cut2 +=1
            if not passBtag: 
                # Fill control region
                if not passTopTag: 
                    continue 
                else:
                    mttbar = calculate_m()
                    h_mttbar_control.Fill( mttbar, weight )
                continue                
            else:
                cut3 +=1
            if not passTopTag: 
                continue                
            else:
                cut4 +=1


            ##  ____  __.__                              __  .__         __________                     
            ## |    |/ _|__| ____   ____   _____ _____ _/  |_|__| ____   \______   \ ____   ____  ____  
            ## |      < |  |/    \_/ __ \ /     \\__  \\   __\  |/ ___\   |       _// __ \_/ ___\/  _ \ 
            ## |    |  \|  |   |  \  ___/|  Y Y  \/ __ \|  | |  \  \___   |    |   \  ___/\  \__(  <_> )
            ## |____|__ \__|___|  /\___  >__|_|  (____  /__| |__|\___  >  |____|_  /\___  >\___  >____/ 
            ##         \/       \/     \/      \/     \/             \/          \/     \/     \/       

            # Now we do our kinematic calculation based on the categories of the
            # number of top and bottom tags
            mttbar = -1.0

            mttbar = calculate_m()
            # Filling plots
            count +=1
            h_mttbar.Fill( mttbar, weight )
            h_mtopHadGroomed.Fill( mass_sd, weight )
            h_mtopHad.Fill( hadTopCandP4.M(), weight )
            
            #fill lepton histos
            h_lepPt.Fill(LeptonPt[0], weight )
            h_lepEta.Fill(LeptonEta[0], weight )
            h_lepPhi.Fill(LeptonPhi[0], weight )

            #fill Jet histos
            h_AK8Pt.Fill(FatJetPt[0], weight )
            h_AK8Eta.Fill(FatJetEta[0], weight )
            h_AK8Phi.Fill(FatJetPhi[0], weight )
    
            h_AK4Pt.Fill( NearestAK4JetPt[0], weight )
            h_AK4Eta.Fill( NearestAK4JetEta[0], weight )
            h_AK4Phi.Fill( NearestAK4JetPhi[0], weight )
            h_AK4M.Fill( NearestAK4JetMass[0], weight )
            h_AK4Bdisc.Fill( NearestAK4JetBDisc[0], weight )
      
            #dr's
            h_drAK4AK8.Fill(bJetCandP4.DeltaR(hadTopCandP4), weight)
            h_drLepAK4.Fill(theLepton.DeltaR(bJetCandP4), weight)
            h_AK8Tau32.Fill(FatJetTau32[0], weight )
            h_AK8Tau21.Fill(FatJetTau21[0], weight )

            #h_dPhiLepAK8.Fill(FatJetDeltaPhiLep[0], weight )

    # Fill cut-flow
    h_cuts.SetBinContent(1, cut1)
    h_cuts.SetBinContent(2, cut2)
    h_cuts.SetBinContent(3, cut3)
    h_cuts.SetBinContent(4, cut4)
    h_cuts.GetXaxis().SetBinLabel(1, "passKin")
    h_cuts.GetXaxis().SetBinLabel(2, "pass2DCut" )
    h_cuts.GetXaxis().SetBinLabel(3, "passBtag")
    h_cuts.GetXaxis().SetBinLabel(4, "passTopTag")

    print options.file_out, " : ", count, "/", tot_entries, ", Percentage:", round(float(count)/(float(tot_entries+1))*100,3), "%", \
     "Cut_flow: [", cut1, cut2, cut3, cut4, "]"

    #fh.write(options.file_in)
    
   
    #fh.close

    fout.cd()
    fout.Write()
    fout.Close()

if __name__ == "__main__" :
    plot_mttbar(sys.argv)
