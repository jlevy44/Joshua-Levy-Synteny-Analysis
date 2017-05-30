import os, numpy as np, itertools

from collections import defaultdict
directories = ['./20Species','./27Species']
for directory in directories:
    speciesCNS = defaultdict(list)
    for file in [file for file in os.listdir(directory) if file.endswith('AllCNSElements.bed')]:
        with open(directory+'/'+file,'r') as f:
            speciesCNS[file.split('_')[0]] = set(np.vectorize(lambda line: line[line.rfind('\t')+1:].split(';')[-4].replace('SegmentID=',''))(f.readlines()))
    CNSList = set(itertools.chain.from_iterable(speciesCNS.values()))
    lenCNSList = float(len(CNSList))
    with open(directory+'/'+'CNSLoss.txt','w') as f:
        for species in speciesCNS.keys():
            f.write(species+'\nCNSRetentionRatio=%f Percent\n\n'%(float(len(speciesCNS[species]))/lenCNSList*100.))

