import re
import sys
from re import sub
for path in sys.path:
    if path and 'anaconda' in path:
        sys.path.remove(path)

import numpy as np
from pybedtools import *
from pyfaidx import Fasta
import subprocess, os, shutil
from collections import *
import time
import dill as pickle
#from multiprocessing import Pool


from difflib import SequenceMatcher

def softmask(a):
    return

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
# help from online sources ^^^^, helps verify backtracking from MAFs to the original fastas

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or info
            configFile.seek(0)
            return line.split()[-1].strip('\n')
    configFile.seek(0)
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
start = time.clock()
pickleSkip = 0
print 'Loading CNS configuration file...','time=',time.clock()-start
configFile = open('configCNSAnalysis.txt','r')
rootfolder = parseConfigFindPath('root_folder',configFile)
pathPython = parseConfigFindPath('pathPython',configFile)
# species names and IDs, FIXME please can add this info to a config file
masterListSpecies = parseConfigFindList('masterListSpecies',configFile)
checkValidity = parseConfigFindPath('checkValidity',configFile)
intragenus = parseConfigFindList('intragenus',configFile)
intergenus = parseConfigFindList('intergenus',configFile)
subgenome = parseConfigFindList('subgenome',configFile)
conservedFastaPath = parseConfigFindPath('conservedFastaPath',configFile)
pickleSkip = parseConfigFindPath('pickleSkip',configFile)
pickleName = parseConfigFindPath('pickleName',configFile)
fasta2phylip = parseConfigFindPath('fasta2phylip',configFile)
phyML = parseConfigFindPath('PhyML',configFile)
bootstrap = parseConfigFindPath('bootstrap',configFile)
treeFile = parseConfigFindPath('treeFile',configFile)
treeOut = parseConfigFindPath('treeOut',configFile)
ratioCopy = parseConfigFindPath('ratioCopy',configFile)
outputTreeImages = parseConfigFindPath('outputTreeImages',configFile)
configFile.close()

if phyML == '1':
    phyML = 1
else:
    phyML = 0

if outputTreeImages == '1':
    outputTreeImages = 1
else:
    outputTreeImages = 0

if ratioCopy == '1':
    ratioCopy = 1
else:
    ratioCopy = 0

if fasta2phylip == '1':
    fasta2phylip = 1
else:
    fasta2phylip = 0

if treeOut == '1':
    treeOut = 1
else:
    treeOut = 0


if pickleSkip == '1':
    pickleSkip = 1
else:
    pickleSkip = 0

if checkValidity == '0':
    checkValidity = 0

sys.path.append(pathPython) # add python path



class speciesClass(): # add information about species that stores name and protyome ID, genome .fa file Fasta object, conserved Bed Element files
    #and you can generate bed files for genes of species and CDS of species
    def __init__(self,speciesNumber,genomeFileList,gffFileList,speciesName,speciesShortName):
        self.speciesNumber = speciesNumber
        for file in genomeFileList:
            if self.speciesNumber in file:
                self.genome = Fasta(file)
        for file in gffFileList:
            if self.speciesNumber in file and 'PAC' in file:
                self.gffFile = file
        self.speciesName = speciesName
        self.speciesShortName = speciesShortName
        self.conservedElementsBed = '%s_ConservedElements.bed'%self.speciesName

        #self.conservedElementsBedFile = open(self.conservedElementsBed, 'w')




speciesInfo = {}

conditionalDictionary = defaultdict(list)



# list all files in analysis directory
listALLFiles = str(subprocess.Popen('ls', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()).split('\n')

# FIXME ['ls', '%s' % '']


print 'Generating File List','time=',time.clock()-start
# generate list of MAF, GFF and .fa files
listMAFfiles = []
listGFFFiles = []
listGenomeFiles = []
for file in listALLFiles:
    if file.endswith('.maf'):
        listMAFfiles.append(file.strip('\n'))
    if file.endswith('.gff') or file.endswith('.gff3'):
        listGFFFiles.append(file.strip('\n'))
    if file.endswith('.fa') or file.endswith('.fasta'):
        listGenomeFiles.append(file.strip('\n'))

print 'Initializing instances of species class...','time=',time.clock()-start
# generate speciesClass objects with relevant info seen above, for each species on masterListSpecies
for species in masterListSpecies:
    speciesInfo[species] = speciesClass(species.split('_')[1], listGenomeFiles, listGFFFiles, species.split('_')[0], species.split('_')[2])
"""
def turnSixteen(x):
    if x == -1:
        return 16
    else:
        return x
"""
print 'Generating list of intragenus, intergenus species and complete species list','time=',time.clock()-start
listIntraGenus = []
listInterGenus = []
listSubgenomes = []
for shortName in intragenus:
    listIntraGenus.append([species.split('_')[0] for species in masterListSpecies if shortName == species.split('_')[-1].strip('\n')][0])
    if shortName in subgenome:
        listSubgenomes.append(
            [species.split('_')[0] for species in masterListSpecies if shortName == species.split('_')[-1].strip('\n')][
                0])
for shortName in intergenus:
    listInterGenus.append([species.split('_')[0] for species in masterListSpecies if shortName == species.split('_')[-1].strip('\n')][0])
listIncludedSpecies = listIntraGenus+listInterGenus

def findBadCharPosition(strSeq):
    """for each MAF sequence, output maximum number of valid characters in a row, exclude duplicates/lowercase/N/softmask <- invalid
    only accept sequence in analysis if at least 15 valid characters in a row"""
    #minVal = np.min(np.vectorize(lambda y: turnSixteen(y))(np.vectorize(lambda x: strSeq.find(x))(np.array(['a','t','c','g','N']))))
    if 'a' in strSeq or 'c' in strSeq or 'N' in strSeq or 'g' in strSeq or 't' in strSeq:
        return np.max(np.vectorize(lambda x: len(x))(np.array(strSeq.replace('a','N').replace('c','N').replace('t','N').replace('g','N').strip('-').split('N'))))
    else:
        return 16

# if checking whether original fasta sequences are valid, MAF info can backtrack to find bed information/fasta DNA sequence start stop info
if checkValidity:
    open('CheckValidity.txt','w').close()
    checkValidFile = open('CheckValidity.txt','w')

segmentCount = 0
# original conditionals, two-copy ratios between species, NKH 111 is PvirN has 1 sequence, K has 1, Hallii has 1
# can design something similar for Bhybridum analysis... eg. BdDSBs1011 D and S are subgenomes

# check out each MAF
mafAnalysisStructure = defaultdict(list)
writeConservedBedFiles = dict.fromkeys(masterListSpecies,)
validityCount = Counter()
cacOutCount = Counter()
for species in listIncludedSpecies:
    validityCount[species] = 0

for file in listMAFfiles:
    print file
    inputMAF = open(file,'r')

    # for each  segment in MAF file
    for segment in inputMAF.read().split('\n\n'):
        # dont skip analyzing segment or writing lines
        skipSegment = 0

        speciesList = []
        outputInfo = []
        if '#' in segment:
            skipSegment = 1
        for line in segment.split('\n'):
            if line and skipSegment == 0:
                if line[0] == 's' and 'Anc' not in line.split()[1]: # if this is a sequence from an actual species
                    # if length of sequence is >= 20 and over 15 valid characters in row
                    if int(line.split()[3]) >= 20 and findBadCharPosition(line.split()[-1]) >= 15:
                        lineList = line.split()
                        lineList2 = lineList[1].split('.')
                        speciesName = lineList2[0]
                        if len(lineList2) > 3:
                            lineList2 = lineList2[0:1] + ['.'.join(lineList2[2:])]
                        lineList3 = lineList2[-1].split('_')
                        if lineList[4] == '-': # negative orientation of sequence, start and end coordinates found from right end of chromosome
                            startCoord,endCoord = int(lineList3[-1])-int(lineList[2])-int(lineList[3]),int(lineList3[-1])-int(lineList[2])
                        else: # positive orientation of sequence, start end coords found from left end chromosomes
                            startCoord, endCoord = int(lineList3[-2])+int(lineList[2]),int(lineList3[-2])+int(lineList[2])+int(lineList[3])

                        # using MAF header info, see if stripped MAF sequence is same as sequence grabbed from original Fasta to make sure backtracking works
                        for species in masterListSpecies:
                            if speciesName in species:
                                MAFseq = lineList[-1].replace('-','')
                                if lineList[4] == '-': # reverse complement find for negative orientation
                                    conservedSeq = str(speciesInfo[species].genome[lineList2[-1][:lineList2[-1].find(lineList3[-2])-1]][startCoord:endCoord].reverse.complement)
                                else:
                                    conservedSeq = str(speciesInfo[species].genome[lineList2[-1][:lineList2[-1].find(lineList3[-2])-1]][startCoord:endCoord])
                                # should have 100% validity
                                cacOutCount[speciesName]+=len(MAFseq)
                                validityCount[speciesName]+=similar(MAFseq,sub("[a-z]",'N',conservedSeq))*100.*len(MAFseq)

with open('backtrackCheck.txt','w') as f:
    for speciesName in listIncludedSpecies:
        print speciesName
        #print validityCount[speciesName]
        #print cacOutCount[speciesName]
        f.write('%s\nCactus Output First Filtering = %dbps\nAverage backtrack percentage = %f Percent\n\n'%(speciesName,cacOutCount[speciesName],float(validityCount[speciesName])/float(cacOutCount[speciesName])))