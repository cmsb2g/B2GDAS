from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'B2GDAS_SingleMuon'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'

config.Data.inputDataset = '/SingleMuon/jpilot-crab_CMSDAS_skim_mu-4f50d602a6e6a83d24076c7043579a32/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Site.storageSite = 'T3_US_FNALLPC'

config.JobType.scriptExe = 'execute_for_crab_data.sh'

config.JobType.outputFiles = ['outplots.root']
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'b2gdas_fwlite.py', 'leptonic_nu_z_component.py', 'JECs', 'purw.root' ]
