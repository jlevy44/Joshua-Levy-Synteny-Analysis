import numpy as np
import matplotlib.pyplot as plt
import os,sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import cPickle as pickle
pickleLoad = 0
"""
try:
    if 'savefile.p' in os.listdir('.'):
        with open('savefile.p','rb') as save:
            (x,y,z) = pickle.load(save)
            pickleLoad = 1
except:
    pass
"""



def CNS(directory):
    print directory
    MASegDict = defaultdict(list)
    seqCount = Counter()
    numFeatures = defaultdict(list)
    speciesDistributionMaster = defaultdict(list)
    for species in [file for file in os.listdir(directory) if file.endswith('.bed')]:
        try:
            print directory+species
            seqCount[species] = 0
            speciesDistribution = Counter()
            with open(directory+species,'r') as f:
                lines = f.readlines()
                numFeatures[species] = [len(lines)]
                if species.endswith('ConservedElements.bed'):
                    for line in lines:
                        if line:
                            lineList = line.split('\t')
                            lineList2 = lineList[-1].split(';')
                            lineList3 = lineList2[1].split(',')
                            tempDict = {word.split(':')[0]:int(word.split(':')[1] != '0') for word in lineList3}
                            MASegDict[lineList2[2].replace('SegmentID=','')] = sum(tempDict.values())
                            seqCount[species] += int(lineList[2])-int(lineList[1])
                            for species2 in tempDict.keys():
                                if species2 not in speciesDistribution.keys():
                                    speciesDistribution[species2] = 0
                                else:
                                    speciesDistribution[species2] += tempDict[species2]
                else:
                    for line in lines:
                        if line:
                            lineList = line.split('\t')
                            lineList2 = lineList[-1].split(';')
                            lineList3 = lineList2[1].split(',')
                            tempDict = {word.split(':')[0]:int(word.split(':')[1] != '0') for word in lineList3}
                            seqCount[species] += int(lineList[2])-int(lineList[1])
                            for species2 in tempDict.keys():
                                if species2 not in speciesDistribution.keys():
                                    speciesDistribution[species2] = 0
                                else:
                                    speciesDistribution[species2] += tempDict[species2]
                speciesDistributionMaster[species] = speciesDistribution
                #print speciesDistributionMaster
                #print numFeatures
                #print ','.join('%s:%d'%(key,speciesDistributionMaster[species][key]) for key in speciesDistributionMaster[species].keys())
        except:
            print 'Error with ' + species
    with open(directory+'CNSStatistics.txt','w') as f:
        for species in sorted(numFeatures.keys()):
            if species:
                try:
                    f.write(species+'\nTotalSequenceAmount=%dbps\nNumberOfElements=%d\n%s\n\n'%(seqCount[species],numFeatures[species][0],'SpeciesDistribution='+','.join('%s:%d'%(key,speciesDistributionMaster[species][key]) for key in speciesDistributionMaster[species].keys())))#FIXME Add species number and graph
                except:
                    print 'Error writing ' + species
    plt.figure()
    plt.hist(MASegDict.values(),bins=np.arange(0,int(np.max(MASegDict.values()))) + 0.5)
    plt.title('Distribution of Number of Species for Conserved Segments')
    plt.ylabel('Count')
    plt.xlabel('Number of species in Conserved Segment')
    plt.savefig(directory+'SpeciesNumberDistribution.png')

if __name__ == '__main__':
    p = Pool(5)
    p.map(CNS, ['./20Species/','./27Species/'])