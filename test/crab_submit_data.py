from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'B2GDAS_SingleMuon'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'

config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.runRange = '273403-273404'

config.Site.storageSite = 'T2_DE_DESY'

config.JobType.scriptExe = 'execute_for_crab_data.sh'

config.JobType.outputFiles = ['output.root']
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab_data.py', 'b2gdas_fwlite.py', 'leptonic_nu_z_component.py', 'JECs' ]
