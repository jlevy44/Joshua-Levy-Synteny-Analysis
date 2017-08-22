import os, subprocess, sys
"""python com1_2 sample refText CDS version"""
sample = sys.argv[1]
refText = sys.argv[2] # ['','','','']
CDS = sys.argv[3]
version = sys.argv[4]
references = refText.strip(']').strip('[').replace("'",'').split(',')
from pipelineFunctions import *
os.chdir(version+'/'+sample)
if all([os.path.isfile('%s.%s.lifted.anchors' %(sample, ref)) and os.stat('%s.%s.lifted.anchors' %(sample, ref)).st_size > 0 for ref in references]) == 0:
    subprocess.call('sh constructv1_1.sh',shell=True)
    sampleCount = 0
    for ref in references:
        try:
            sampleCount = replaceGeneNames(sample, ref, sampleCount)
        except:
            print sample, ref, sampleCount
try:
    tiling2bed('%snuc.tiling'%CDS, CDS, sample, sample+'_%ssyn'%CDS+'.bed')
    replaceGeneNames(sample, CDS, 0, 1)
except:
    print 'Unable to finish nucmify'
try:
    BB2bed('BBmapped.bed',CDS, sample,'../../referenceGenomes/%s/centromere.bed'%CDS)
    replaceGeneNames(sample, CDS, 0, 0, 1)
except:
    print 'Unable to finish BBMapping'
try:
    subprocess.call('sh constructv1_2.sh',shell=True)
except:
    print 'Unable to finish sh constructv1_2.sh for ' + sample
