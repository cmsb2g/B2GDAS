import FWCore.ParameterSet.Config as cms

process = cms.Process('NoSplit')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring([
            '/store/mc/RunIISpring16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/B8AEB5CC-9638-E611-BD77-D4AE526A048B.root',
            '/store/mc/RunIISpring16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/BA012C20-9638-E611-9FB4-782BCB161FC2.root',
            '/store/mc/RunIISpring16MiniAODv2/RSGluonToTT_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/F058BE0A-9638-E611-8115-0025905C42FE.root'
            ]))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
