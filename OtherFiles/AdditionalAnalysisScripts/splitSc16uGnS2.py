from pyfaidx import Fasta
from statistics import mean, stdev
import numpy as np


gffFiles = ['546.gff3','547.gff3']

for file in gffFiles:
    inputFile = open(file,'r')
    key = file.split('.')[0]
    minVal=[]
    maxVal=[]
    for line in inputFile:
        if line and line.split()[0] == 'Sc16uGnS2' and 'longest=1' in line and 'mRNA' in line:
            minVal.append(int(line.split()[3]))
            maxVal.append(int(line.split()[4]))

    minValArr = np.array(minVal)
    maxValArr = np.array(maxVal)

    print '%s\tmin=%d,median=%d,15th=%d,82.5th=%d\tmax=%d,median=%d,15th=%d,82.5th=%d\t\n'%(key,np.min(minValArr),
        np.median(minValArr),np.percentile(minValArr,15,interpolation='nearest'),np.percentile(minValArr,82.5,interpolation='nearest'),np.max(maxValArr),np.median(maxValArr),np.percentile(maxValArr,15,interpolation='nearest'),np.percentile(maxValArr,82.5,interpolation='nearest'),)