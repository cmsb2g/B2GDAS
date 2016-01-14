#! /usr/bin/env python
import sys
import array as array
from plot_mttbar import plot_mttbar


#ttjetsStr = ["--file_in", "ttjets.root", "--file_out", "ttjets_plots.root"]
#rsgluon3TeVStr = ["--file_in", "rskkgluon3TeV.root", "--file_out", "rskkgluon3TeV_plots.root"]
singlemuStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleMuon/crab_B2GDAS_SingleMuon_PromptRecov4/160114_144645/0000/MuData_PromptRecov4.root", "--file_out", "MuData_PromptRecov4.root", "--isData"]
#singlemu_elStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleMuon/crab_B2GDAS_SingleMuon_PromptRecov4/160114_144645/0000/MuData_PromptRecov4.root", "--file_out", "MuData_PromptRecov4_Electron.root", "--isData","--isElectron"]
singlemu_OctStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleMuon/crab_B2GDAS_SingleMuon/160113_235122/0000/MuData_Oct5.root", "--file_out", "MuData_Oct5_plots.root", "--isData"]
#singlemu_Oct_elStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleMuon/crab_B2GDAS_SingleMuon/160113_235122/0000/MuData_Oct5.root", "--file_out", "MuData_Oct5_plots_Electron.root", "--isData","--isElectron"]


singleelStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleElectron/crab_B2GDAS_SingleElectron_PromptRecov4/160114_144229/0000/EleData_PromptRecov4.root", "--file_out", "EleData_PromptRecov4.root", "--isData","--isElectron"]
#singleel_muStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleElectron/crab_B2GDAS_SingleElectron_PromptRecov4/160114_144229/0000/EleData_PromptRecov4.root", "--file_out", "EleData_PromptRecov4_Muon.root", "--isData"]
singleel_OctStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleElectron/crab_B2GDAS_SingleElectron/160113_234934/0000/EleData_Oct5.root", "--file_out", "EleData_Oct5_plots.root", "--isData","--isElectron"]
#singleel_mu_OctStr = ["--file_in", "root://cmsxrootd.fnal.gov//store/user/pilot/SingleElectron/crab_B2GDAS_SingleElectron/160113_234934/0000/EleData_Oct5.root", "--file_out", "EleData_Oct5_plots_Muon.root", "--isData"]

plot_mttbar( singlemuStr )
#plot_mttbar( singlemu_elStr  )
plot_mttbar( singlemu_OctStr )
#plot_mttbar( singlemu_Oct_elStr )
plot_mttbar( singleelStr )
#plot_mttbar( singleel_muStr )
plot_mttbar( singleel_OctStr )
#plot_mttbar( singleel_mu_OctStr )
