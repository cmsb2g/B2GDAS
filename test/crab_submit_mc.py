from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'B2GDAS_RSGluonToTT_M-3000'
config.General.workArea = 'B2GDAS'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'

config.Data.inputDataset = '/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/CMSDAS' % (getUsernameFromSiteDB())
config.Data.publication = False

config.Site.storageSite = 'T3_US_FNALLPC'

config.JobType.scriptExe = 'execute_for_crab.sh'
config.JobType.outputFiles = ['output.root']
config.JobType.sendExternalFolder = True
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'b2gdas_fwlite.py', 'leptonic_nu_z_component.py', 'JECs', 'purw.root', 'egammaEffi.txt_SF2D.root', 'MuonID_Z_RunBCD_prompt80X_7p65.root', 'general_tracks_and_early_general_tracks_corr_ratio.root' ]
