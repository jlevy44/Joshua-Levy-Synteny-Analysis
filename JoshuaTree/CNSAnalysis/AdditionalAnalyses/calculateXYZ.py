import sys,os
from collections import defaultdict, Counter
from pyfaidx import Faidx, Fasta
from pybedtools import BedTool
import subprocess
import numpy as np
import shutil
import operator
from ete2 import Tree
import matplotlib.pyplot as plt
import cPickle as pickle
pickleLoad = 0

try:
    if 'savefile.p' in os.listdir('.'):
        with open('savefile.p','rb') as save:
            (x,y,z) = pickle.load(save)
            pickleLoad = 1
except:
    pass

findlen = lambda x: sum([int(line.split('\t')[2])-int(line.split('\t')[1]) for line in str(x).split('\n') if line])
def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of relevant filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list of files or list information
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    configFile.seek(0)
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or info
            configFile.seek(0)
            return line.split()[-1].strip('\n')
    configFile.seek(0)
with open('configCNSAnalysis.txt','r') as f:
    inputSpecies = parseConfigFindList('masterListSpecies', f)

protId = defaultdict(list)
for line in inputSpecies:
    if line:
        print line
        protId[line.split('_')[0]] = line.split('_')[1]


inputList = protId.keys()
analysisRun = os.getcwd().split('/')[-2] + ' Run: '
specId = {k:v for v,k in protId.iteritems()}
fastaFolders = ['../CactusRun/output/'+folder+'/' for folder in os.listdir('../CactusRun/output') if folder.startswith('FastaOut')]
if pickleLoad == 0:
    x = Counter()  # amount of sequence in genome
    y = Counter()  # amount of sequence passed into cactus
    z = defaultdict(list)  # amount of sequence passed out of Cactus
    for species in inputList:
        print species
        x[species] = 0
        y[species] = 0
        z[species] = []
        for file in os.listdir('.'):
            if file.endswith('.fai') and protId[species] in file:
                with open(file,'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if line:
                            x[species] += int(line.split('\t')[1])#abs(int(line.split('\t')[3]) - int(line.split('\t')[2]))#max(map(int,line.split('\t')[2:4]))-min(map(int,line.split('\t')[2:4]))
        for folder in [folder2 for folder2 in fastaFolders if species+'.fa' in os.listdir(folder2)]:
            try:
                fa = Fasta(folder+species+'.fa')
                #bedText = '\n'.join('\t'.join(['_'.join(line.split('_')[:-2])] + line.split('_')[-2:]) for line in fa.keys())
                y[species] += sum([len(fa[key][:].seq) for key in fa.keys()])#findlen(BedTool(bedText, from_string=True))
            except:
                print 'Error for ' + folder+species+'.fa'
            print y[species]

    """
    with open('finalSyntenyMultipleSpecies.bed','r') as f:
        print 'Bed Open...'
        lineOut = []
        for line in f.readlines():
            lineOut.append('-'.join(line.split('\t')[0:4])+'|'+line[line.rfind('\t')+1:])
        for line in lineOut:
            for seq in line.split('|'):
                y[specId[seq.split('-')[0]]] += abs(int(seq.split('-')[3]) - int(seq.split('-')[2]))
    """

    for mafFile in [file for file in os.listdir('.') if file.endswith('.maf')]:
        print 'maf Start'
        with open(mafFile,'r') as f:
            for line in f.readlines():
                if line and line.startswith('s') and 'Anc' not in line.split('\t')[1].split('.')[0]:
                    z[line.split('\t')[1].split('.')[0]].append(int(line.split('\t')[3]))

    with open('savefile.p','wb') as save:
        pickle.dump((x,y,z),save)

with open('sequenceTest.txt','w') as f:
    print 'write output...'
    for species in inputList:
        f.write('%s\nGenomeSequence(x)=%dbps\nSequenceInputToCactus(y)=%dbps\nSequenceOutputFromCactus(z)=%dbps\nAverageConservedSequenceLength=%fbps\ny/x=%f\nz/y=%f\n\n'%(species,x[species],y[species],sum(z[species]),np.mean(z[species]),float(y[species])/float(x[species]),float(sum(z[species]))/float(y[species])))
        plt.hist(z[species], bins=np.arange(0, 100, 5))
        plt.axvline(x=15, color='r', linewidth=5)
        plt.title(analysisRun+'Histogram of ' + species+ ' Conserved Sequence Lengths')
        plt.xlabel('CS Length (bps)')
        plt.ylabel('Count')
        plt.savefig('%s_%s_LengthThreshold.png'%(analysisRun.replace(' ','').strip(':'),species), bbox_inches='tight')

numberSpecies = len(inputList)
X = tuple([x[species] for species in sorted(z.keys())])
Y = tuple([y[species] for species in sorted(z.keys())])
Z = tuple([sum(z[species]) for species in sorted(z.keys())])

fig, ax = plt.subplots()
index = np.arange(numberSpecies)
bar_width = 0.3
opacity = 0.8
#FIXME Ordering is wrong!!!
rects1 = plt.bar(index - bar_width, X, bar_width,
                 alpha=opacity,
                 color='b',
                 label='Total Sequence in Genome')

rects2 = plt.bar(index, Y, bar_width,
                 alpha=opacity,
                 color='g',
                 label='Total Sequence Input to Cactus')

rects3 = plt.bar(index + bar_width, Z, bar_width,
                 alpha=opacity,
                 color='r',
                 label='Total Sequence Output from Cactus')

plt.xlabel('Species')
plt.ylabel('Total Sequence Length')
plt.title('Amount of Sequence Passing Through Analysis')
plt.xticks(index + bar_width, tuple([protId[species] for species in sorted(z.keys())]), rotation = 'vertical')
plt.legend()
plt.savefig(analysisRun.replace(' ','').strip(':') + '_sequenceFlow.png',bbox_inches='tight')
plt.tight_layout()

plt.show()
