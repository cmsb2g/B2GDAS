import os

files = os.listdir("root_files")
names = []
for file in files:
	if file.endswith("_plots_ele_jec_down.root"): names.append(file[:-(len("_plots_ele_jec_down.root"))])

for name in names:
	os.system("hadd root_files/hadd_" + name+".root root_files/" + name +"_plots*")