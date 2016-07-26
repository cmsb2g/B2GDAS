import sys
import PSet
files = []
outfile = file( 'filesToProcess.txt', 'w')
for ifile in PSet.process.source.fileNames :    
    outfile.write( ifile + '\n' )


outfile.close()

sys.argv.append('--input')
sys.argv.append('filesToProcess.txt')
sys.argv.append('--isCrabRun')
sys.argv.append('--trigProc')
sys.argv.append('HLT2')

print sys.argv

from b2gdas_fwlite import *

b2gdas_fwlite(sys.argv)

