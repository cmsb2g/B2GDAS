from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'B2GDAS_RSKKGluon3TeV'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'

config.Data.inputDataset = '/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Site.storageSite = 'T2_DE_DESY'

config.JobType.scriptExe = 'execute_for_crab.sh'

config.JobType.outputFiles = ['outplots.root']
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'b2gdas_fwlite.py', 'leptonic_nu_z_component.py', 'JECs', 'purw.root', 'egammaEffi.txt_SF2D.root', 'MuonID_Z_2016runB_2p6fb.root', 'general_tracks_and_early_general_tracks_corr_ratio.root' ]
