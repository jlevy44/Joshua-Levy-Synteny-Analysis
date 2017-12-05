import re
import sys

import numpy as np
from pybedtools import *
from pyfaidx import Fasta
import subprocess, os, shutil
from collections import *
import time
import dill as pickle
from ete3 import Tree,TreeStyle,NodeStyle
#from multiprocessing import Pool


from difflib import SequenceMatcher



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
shortFastaOutputName = parseConfigFindPath('short_fasta_output_name',configFile)
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

if shortFastaOutputName == '1':
    shortFastaOutputName = 1
else:
    shortFastaOutputName = 0

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

    def gff2bedobjs(self):
        outputfilename = self.gffFile.replace('.gff', '_Genes.bed').replace('.gff3', '_Genes.bed')
        outputfilename2 = self.gffFile.replace('.gff', '_CDS.bed').replace('.gff3', '_CDS.bed')
        open(outputfilename, 'w').close()
        open(outputfilename2, 'w').close()
        outputFile = open(outputfilename, 'w')
        outputFile2 = open(outputfilename2, 'w')
        CDS_read = 0
        for line in open(self.gffFile, 'r'):
            if line:
                if 'mRNA' in line and 'longest=1' in line and int(line.split('\t')[3] >= 0):
                    geneName = line.split()[-1].split(';')[1].replace('Name=','')
                    outputFile.write('%s\t%d\t%s\t%s\n'%(line.split()[0],int(line.split()[3]) - 1,line.split()[4],geneName))
                    CDS_read = 1
                elif 'mRNA' in line and 'longest=1' not in line:
                    CDS_read = 0

                if 'CDS' in line:
                    if CDS_read and not geneName:
                        geneName = 'NoCNSName'
                    outputFile2.write('%s\t%d\t%s\t%s\n' % (line.split()[0], int(line.split()[3]) -1, line.split()[4], geneName))
                #outputFile.write('%s\t%d\t%s\n' % (line.split()[1],int(line.split()[2])-1,line.split()[3]))
        outputFile.close()
        outputFile2.close()
        self.bedGenes = BedTool(outputfilename).sort().merge(c=4,o='distinct',delim='|')
        self.bedCDS = BedTool(outputfilename2).sort().merge(c=4,o='distinct',delim='|')


def count2Conditional(countSeq,speciesList):
    global speciesInfo
    return ''.join('%s%d'%(speciesInfo[[speciesName for speciesName in speciesInfo.keys() if species in speciesName][0]].speciesShortName,countSeq[species]) for species in speciesList)

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
    if file.endswith('.maf') and file.startswith('FastaOut'):
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

if pickleSkip == 0:
    print 'Reading MAF Segments and inputting them into MAF analysis structure','time=',time.clock()-start
    for species in masterListSpecies:
        open(speciesInfo[species].conservedElementsBed, 'w').close()
        writeConservedBedFiles[species] = open(speciesInfo[species].conservedElementsBed,'w')
    for file in listMAFfiles:
        inputMAF = open(file,'r')

        # for each  segment in MAF file
        for segment in inputMAF.read().split('\n\n'):
            writeTheseLines = defaultdict(list)
            countSeq = Counter()
            # dont skip analyzing segment or writing lines
            skipSegment = 0
            # set count of # seq for each species in segment to 0, can have multiple copies throughout genome
            for species in masterListSpecies:
                countSeq[species.split('_')[0]] = 0# FIXME change _ to - the newest MAF files will likely use - instead of _
                # , may have to very code to reflect this
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
                            countSeq[speciesName] += 1
                            if len(lineList2) > 3:
                                lineList2 = lineList2[0:1] + ['.'.join(lineList2[2:])]
                            lineList3 = lineList2[-1].split('_')
                            if lineList[4] == '-': # negative orientation of sequence, start and end coordinates found from right end of chromosome
                                startCoord,endCoord = int(lineList3[-1])-int(lineList[2])-int(lineList[3]),int(lineList3[-1])-int(lineList[2])
                            else: # positive orientation of sequence, start end coords found from left end chromosomes
                                startCoord, endCoord = int(lineList3[-2])+int(lineList[2]),int(lineList3[-2])+int(lineList[2])+int(lineList[3])
                            writeTheseLines[speciesName].append('%s\t%d\t%d\t%s;'
                                        %(lineList2[-1][:lineList2[-1].find(lineList3[-2])-1],startCoord,endCoord,lineList[4]))
                            # writing to a bed file for a particular species, generate list to write for each species (chr,xi,xf,orientation)!!!!
                        else: # if dont meet reqs, skip analysis of this segment
                            skipSegment = 1
                        # using MAF header info, see if stripped MAF sequence is same as sequence grabbed from original Fasta to make sure backtracking works
                        if checkValidity:
                            for species in masterListSpecies:
                                if speciesName in species:
                                    MAFseq = lineList[-1].replace('-','')
                                    if lineList[4] == '-': # reverse complement find for negative orientation
                                        conservedSeq = str(speciesInfo[species].genome[lineList2[-1][:lineList2[-1].find(lineList3[-2])-1]][startCoord:endCoord].reverse.complement)
                                    else:
                                        conservedSeq = str(speciesInfo[species].genome[lineList2[-1][:lineList2[-1].find(lineList3[-2])-1]][startCoord:endCoord])
                                    # should have 100% validity
                                    checkValidFile.write('%s\tMAF=%s\tConservedSequence=%s\tOrientation=%s\tPercentSimilar=%d%%\n' %(speciesName,MAFseq,conservedSeq,lineList[4],similar(MAFseq,conservedSeq)*100.))
            if len([key for key in countSeq.keys() if countSeq[key] >= 1]) < 2:
                skipSegment = 1
            # FIXME need INCLUDE AT LEAST TWO SPECIES ELSE SKIP SEGMENT ^^^ CAN FIX ABOVE
            if skipSegment == 0: # if staying with segment for all species, write the lines to output bed and add respective
                #  number of sequences per species for each segment, also output MASegment ID number
                for speciesName in countSeq.keys():
                    for line2 in writeTheseLines[speciesName]:
                        writeConservedBedFiles[[species for species in masterListSpecies if speciesName in species][0]].write(line2+''.join('%s:%d'%(speciesName2,countSeq[speciesName2]) + ',' for speciesName2 in countSeq.keys())[:-1]+';SegmentID=MASeg%d'%segmentCount+'\n')
                # FIXME Start here build new dictionary structure, it will be like this for now
                # if NKH111, output NKH MAF sequences to fasta file in NKH111 folder
                if all(countSeq[species] < 2 for species in listIncludedSpecies) and sum([countSeq[key] for key in listIncludedSpecies]) > 1:
                    conditional = count2Conditional(countSeq,listIncludedSpecies)
                    if conditionalDictionary[conditional] == []:
                        conditionalDictionary[conditional] = countSeq
                        for species in conditionalDictionary[conditional].keys():
                            if species not in listIncludedSpecies:
                                del conditionalDictionary[conditional][species]
                    if mafAnalysisStructure[conditional] == None or mafAnalysisStructure[conditional] == []: #FIXME
                        mafAnalysisStructure[conditional] = {'CS': defaultdict(list),'ProteinCoding_CS': [],'Mixed_CS': [],'NonProteinCoding_CS': []}
                    mafAnalysisStructure[conditional]['CS']['MASeg%d'%segmentCount] = ''.join('>%s_%s_%s_%s\n%s\n'%(tuple(line.split()[1:5])+(line.split()[-1],)) for species in [species for species in listIncludedSpecies if species in listIncludedSpecies and countSeq[species] == 1] for line in segment.split('\n') if line and species in line)
                elif all(countSeq[species] < 2 for species in listIntraGenus) and sum([countSeq[key] for key in listIntraGenus]) > 1:
                    conditional = count2Conditional(countSeq,listIntraGenus)
                    if conditionalDictionary[conditional] == []:
                        #FIXME and fix this!!!
                        conditionalDictionary[conditional] = countSeq
                        for species in conditionalDictionary[conditional].keys():
                            if species not in listIntraGenus:
                                del conditionalDictionary[conditional][species]
                    if mafAnalysisStructure[conditional] == None or mafAnalysisStructure[conditional] == []: #FIXME
                        mafAnalysisStructure[conditional] = {'CS': defaultdict(list),'ProteinCoding_CS': [],'Mixed_CS': [],'NonProteinCoding_CS': []}
                    mafAnalysisStructure[conditional]['CS']['MASeg%d'%segmentCount] = ''.join('>%s_%s_%s_%s\n%s\n'%(tuple(line.split()[1:5])+(line.split()[-1],)) for species in [species for species in listIntraGenus if countSeq[species] == 1] for line in segment.split('\n') if line and species in line)
                    a=1

            #FIXME _ to -


            segmentCount+=1 # add count to number of segments to generate multiple alignment segment IDs


                                    #, lineList[4],''.join(speciesName2 + ',' for speciesName2 in speciesList)[:-1]))
                            #FIXME '_' or '-' ^^^^^
            #for infoTuple in outputInfo:
            #write -->    (infoTuple+(''.join(speciesName + ',' for speciesName in speciesList)[:-1],)))#str(Pvir[infoTuple[0]][infoTuple[1]:infoTuple[2]]))))
        a=1
        inputMAF.close()
    for species in writeConservedBedFiles.keys():
        writeConservedBedFiles[species].close()


    a=1



    if checkValidity:
        checkValidFile.close()

    bedCNSIntergenic = {}
    bedCNSIntronic = {}
    bedConservedCDS = {}
    bedCNS = {}
    print 'Sorting sequences into conserved elements, CDS, CNS and intronic/intergenic elements...','time=',time.clock()-start
    for species in masterListSpecies:

        conservedElementsBed = open(speciesInfo[species].conservedElementsBed,'r')

        # create bed object from conserved elements for each species, conserved elements found from MAF files done above
        speciesInfo[species].conservedElementsBed = BedTool(conservedElementsBed).sort().merge(c=4,o='distinct',delim='|')
        conservedElementsBed.close()
        # generate CDS and gene bed objects for each species
        speciesInfo[species].gff2bedobjs()

        # This is both CNS Intronic + intergenic, [CNS total] = CS-CDS = CS - CDS
        bedCNS[species] = speciesInfo[species].conservedElementsBed.subtract(speciesInfo[species].bedCDS)\
            .sort().merge(c=4,o='distinct',delim='|')

        #bedCNS = {}
        # find CNS Intronic genes = CS + (genes-CDS)  + is intersection, - is subtraction of sequences
        bedCNSIntronic[species] = bedCNS[species].intersect(speciesInfo[species].bedGenes).sort().merge(c=4,o='distinct',delim='|').closest(speciesInfo[species].bedGenes,nonamecheck = True)
        # CNS intergenic = CS - genes
        bedCNSIntergenic[species] = speciesInfo[species].conservedElementsBed.subtract(speciesInfo[species].bedGenes).sort() \
            .merge(c=4, o='distinct', delim='|').closest(speciesInfo[species].bedGenes,k=1,d=True,nonamecheck = True)

        # Conserved CDS = CS +CDS
        bedConservedCDS[species] = speciesInfo[species].conservedElementsBed.intersect(speciesInfo[species].bedCDS)\
            .sort().merge(c=4,o='distinct',delim='|').closest(speciesInfo[species].bedGenes,nonamecheck = True)

        bedCNS[species] = bedCNS[species].closest(speciesInfo[species].bedGenes,k=1,d=True,nonamecheck = True)



    # .merge(o='distinct',c=4,delim = '|')

    #bedCNSFinal = bedCNS[bedCNS.keys()[0]].sort()

        #,str(Pvir[infoTuple[0]][infoTuple[1]:infoTuple[2]]

        #for key in bedCNS.keys()[1:]:
        #    bedCNSFinal = bedCNSFinal.cat(bedCNS[key],postmerge = False).sort()
        #FIXME start here!!!
        #bedCNSFinal = bedCNSFinal.merge(o='distinct',c=4,delim = '|').sort()
        CNSOutFiles = [('%s_CNSElements_Intronic.bed'%speciesInfo[species].speciesName,bedCNSIntronic[species].filter(lambda x: len(x) > 1)),
                       ('%s_CNSElements_Intergenic.bed'%speciesInfo[species].speciesName,bedCNSIntergenic[species].filter(lambda x: len(x) > 1)),
                       ('%s_Conserved_CDS.bed' % speciesInfo[species].speciesName,bedConservedCDS[species].filter(lambda x: len(x) > 1)),
                       ('%s_AllCNSElements.bed' % speciesInfo[species].speciesName,bedCNS[species].filter(lambda x: len(x) > 1))]
        for bedOut in CNSOutFiles: # for each of the above described files, append the sequences from the original fastas and include closest gene/ distance information
            open(bedOut[0],'w').close()
            if 0:
                bedOutFile = open(bedOut[0],'w')

                for line in str(bedOut[1]).split('\n'):
                    if line:
                        lineList = line.split()
                        # need to have a long enough sequence to get outputted...
                        if len(outputSequence) >= 15 and findBadCharPosition(outputSequence) >= 15:
                            if 'CNSElements_Intronic' in bedOut[0] or 'Conserved_CDS' in bedOut[0]:
                                bedOutFile.write('%s\t%s\t%s\t%s;geneID=%s\n'%tuple(lineList[0:4]+[lineList[-1]]))
                            elif 'CNSElements_Intergenic' in bedOut[0] or 'AllCNSElements' in bedOut[0]:
                                bedOutFile.write('%s\t%s\t%s\t%s;closestGene=%s;distance=%s\n'%tuple(lineList[0:4]+[lineList[-2:],lineList[-1].strip('\n')]))

            if 1: #FIXME removed for now, please add later
                bedOutFile = open(bedOut[0],'w')

                for line in str(bedOut[1]).split('\n'):
                    if line:
                        lineList = line.split()
                        try:
                            outputSequence = str(speciesInfo[species].genome[lineList[0]][int(lineList[1]):int(lineList[2])])
                            # need to have a long enough sequence to get outputted...
                            if len(outputSequence) >= 15 and findBadCharPosition(outputSequence) >= 15:
                                if 'CNSElements_Intronic' in bedOut[0] or 'Conserved_CDS' in bedOut[0]:
                                    bedOutFile.write('%s\t%s\t%s\t%s;geneID=%s;%s\n'%(tuple(lineList[0:4])+(lineList[-1],outputSequence)))
                                elif 'CNSElements_Intergenic' in bedOut[0] or 'AllCNSElements' in bedOut[0]:
                                    bedOutFile.write('%s\t%s\t%s\t%s;closestGene=%s;distance=%s;%s\n'%(tuple(lineList[0:4])+(lineList[-2],lineList[-1].strip('\n'),outputSequence)))
                        except:
                            print species,lineList


        #bedOut.write(str(bedCNSFinal))

            bedOutFile.close()

        conservedElementsBed.close()



    # run secondary analysis that checks out conservation of two-copy genes etc... can generate trees from this secondary analysis
    """Second Analysis, used to be manip files test..."""
    def bed2dict(bedfname): # create a dictionary that has MA Segments as keys
        bedDict={}
        for line in open(bedfname,'r'):
            if line:
                MASeg = line.split()[-1].split(';')[2].replace('SegmentID=','')
                if bedDict.has_key(MASeg) and tuple(line.split()[0:3]) not in bedDict[MASeg]:
                    bedDict[MASeg].append(tuple(line.split()[0:3]))
                else:
                    bedDict[MASeg] = [tuple(line.split()[0:3])]
        return bedDict

    ultimateBedDicts = {'Conserved_CDS.bed':{},'AllCNSElements.bed':{}} # split into ~CDS or ~CNS
    print 'Generating ultimate CNS, CDS structure...','time=',time.clock()-start
    for file in os.listdir('.'): # FIXME for some reason had problem saving
        if file.endswith('Conserved_CDS.bed') or file.endswith('AllCNSElements.bed'):
            ultimateBedDicts[file[file.find('_')+1:]][file[:file.find('_')]] = bed2dict(file) # a dictionary with keys [species] and inside a key of MASEG
        # hierarchy is ultimate[CDS or CNS][species][MASeg]
    print 'Writing data to %s'%pickleName,'time=',time.clock()-start
    pickle.dump(tuple([ultimateBedDicts,mafAnalysisStructure, speciesInfo, conditionalDictionary]), open(pickleName, 'wb'),protocol=2)


if pickleSkip == 1:
    print 'Skipping MAF file reading and loading bed file structures, species class instances, conditional dictionary, and maf analysis structure','time=',time.clock()-start
    (ultimateBedDicts,mafAnalysisStructure, speciesInfo, conditionalDictionary) = pickle.load(open(pickleName,'rb'),protocol=2)

def conditional2SpeciesList(conditionalDict): # code obsolete
    global masterListSpecies
    # get dictionary keys that have 1, convert those short names to species names
    return np.vectorize(lambda shortName: [species.split('_')[0] for species in masterListSpecies if shortName ==
                                           species.split('_')[-1].strip('\n')][0])\
        ([shortname for shortname in conditionalDict.keys() if conditionalDict[shortname]])

ratio2CopyStats = dict.fromkeys(mafAnalysisStructure.keys(),{'CS':0,'ProteinCoding_CS':0,
                     'Mixed_CS':0,'NonProteinCoding_CS':0})
print 'Generating final fasta files...','time=',time.clock()-start
# generate final fastas
#print conditionalDictionary
for conditional in mafAnalysisStructure.keys():
    if type(conditionalDictionary[conditional] )==type(list):
        a=1
    speciesListCDSTest = [species for species in conditionalDictionary[conditional].keys() if conditionalDictionary[conditional][species] == 1]
    conservedSeqs = {'CS':dict.fromkeys(speciesListCDSTest,''),'ProteinCoding_CS':dict.fromkeys(speciesListCDSTest,''),
                     'Mixed_CS':dict.fromkeys(speciesListCDSTest,''),'NonProteinCoding_CS':dict.fromkeys(speciesListCDSTest,'')}
    for MASeg in mafAnalysisStructure[conditional]['CS']:
        conservedSeq = dict.fromkeys(speciesListCDSTest,'')
        if type(mafAnalysisStructure[conditional]['CS'][MASeg]) == type(list):
            a=1
        for sequence in mafAnalysisStructure[conditional]['CS'][MASeg].split('>'):#FIXME error here
            if sequence:
                speciesName = sequence.split('.')[0]
                conservedSeq[speciesName] = sequence.split('\n')[1].strip('\n')
        for species in speciesListCDSTest:
            conservedSeqs['CS'][species] += conservedSeq[species]
        CDS_conditional = np.vectorize(
            lambda speciesName: ultimateBedDicts['Conserved_CDS.bed'][speciesName].has_key(MASeg))(speciesListCDSTest)
        if np.all(CDS_conditional):
            # Fixme change mafAnalysisStructure, probably do not need to include other analysis types, but may need to keep this to calculate two-copy ratios
            mafAnalysisStructure[conditional]['ProteinCoding_CS'].append(MASeg)
            for species in speciesListCDSTest:
                conservedSeqs['ProteinCoding_CS'][species] += conservedSeq[species]
        elif np.any(CDS_conditional) and not np.all(CDS_conditional):
            mafAnalysisStructure[conditional]['Mixed_CS'].append(MASeg)
            for species in speciesListCDSTest:
                conservedSeqs['Mixed_CS'][species] += conservedSeq[species]
        elif not np.any(CDS_conditional):
            mafAnalysisStructure[conditional]['NonProteinCoding_CS'].append(MASeg)
            for species in speciesListCDSTest:
                conservedSeqs['NonProteinCoding_CS'][species] += conservedSeq[species]
    for analysisType in ['CS','ProteinCoding_CS','Mixed_CS','NonProteinCoding_CS']:
        try:
            if shortFastaOutputName:
                #print conditional, conditionalDictionary[conditional]
                conditional = ''.join(map(str,conditionalDictionary[conditional].values()))
                #print conditional
            finalFastaFileName = conditional+'_'+analysisType+'.fasta'
            open(conservedFastaPath+finalFastaFileName,'w').close()
            outputFinalFastaFile = open(conservedFastaPath+finalFastaFileName,'w')
            outputFinalFastaFile.write(''.join('>%s\n%s\n'%(species,conservedSeqs[analysisType][species]) for species in conservedSeqs[analysisType].keys() if conservedSeqs[analysisType][species]))
            outputFinalFastaFile.close()
        except:
            print conditional
        """if analysisType == 'CS': #FIXME removed for now
            ratio2CopyStats[conditional][analysisType]= len(mafAnalysisStructure[conditional][analysisType].keys())
        else:
            ratio2CopyStats[conditional][analysisType]= len(mafAnalysisStructure[conditional][analysisType])"""
#FIXME everything below needs to be fixed...
# output two copy statistics


if ratioCopy: #FIXME defunct for now...
    print 'Performing 2 copy ratio statistics...','time=',time.clock()-start
    #species Lists
    listConditionals = conditionalDictionary.keys()
    setIntraNotSubgenomesSpecies = set(listIntraGenus) - set(listSubgenomes) # intragenus, but not subgenomes
    setOfSpecies = set(listIncludedSpecies)#set([conditionalDictionary[conditional].keys() for conditional in listConditionals if len(conditionalDictionary[conditional].keys()) == len(masterListSpecies)][0]) #list of species
    setOnlyExtragenus = set(listInterGenus)#setOfSpecies - set(listIntraGenus) # only extra/inter genus species


    #Conditional Lists

    listIntraNotSubgenomeConditionals = set([conditional for conditional in listConditionals if
                                       all(conditionalDictionary[conditional][species] == 1 for species in
                                           setIntraNotSubgenomesSpecies)]) # conditionals that are make sure intragenus species that are not subgenomes are all one copy, contains subset and total list of species conditionals, can be intersected
    allIntraGenusLowSpeciesConditionals = set([conditional for conditional in listConditionals if
                                       any(conditionalDictionary[conditional][species] == 1 for species in
                                           listIntraGenus) and set(conditionalDictionary[conditional].keys()) == set(listIntraGenus)]) #FIXME conditionals of only a subset of the total species, only of the intragenus
    allIntraInterGenusSpeciesConditionals = set(listConditionals)-allIntraGenusLowSpeciesConditionals # conditionals of all of the species, excludes the subset above
    extraGenusAllOnesConditionals = set([conditional for conditional in listConditionals if all(conditionalDictionary[conditional][species] == 1 for species in setOnlyExtragenus)])
    extraGenusMixOnesConditionals = set([conditional for conditional in listConditionals if any(conditionalDictionary[conditional][species] == 1 for species in setOnlyExtragenus)])-extraGenusAllOnesConditionals
    extraGenusAllZerosConditionals = set([conditional for conditional in listConditionals if all(conditionalDictionary[conditional][species] == 0 for species in setOnlyExtragenus)])
    twoCopyList = (set([conditional for conditional in listIntraNotSubgenomeConditionals if
                   all(conditionalDictionary[conditional][species] == 1 for species in listSubgenomes)])).intersection(listIntraNotSubgenomeConditionals)
    oneCopyList = (set([conditional for conditional in listIntraNotSubgenomeConditionals if
                         any(conditionalDictionary[conditional][species] == 1 for species in listSubgenomes)])- twoCopyList).intersection(listIntraNotSubgenomeConditionals)
    twoCopy = {'Intra': twoCopyList.intersection(allIntraInterGenusSpeciesConditionals.intersection(extraGenusAllZerosConditionals)),
               'Inter': twoCopyList.intersection(allIntraInterGenusSpeciesConditionals.intersection(extraGenusAllOnesConditionals)),
               'Any': twoCopyList.intersection(allIntraInterGenusSpeciesConditionals.intersection(extraGenusMixOnesConditionals)),'Mixed Intra':twoCopyList.intersection(allIntraGenusLowSpeciesConditionals)}
    oneCopy = {'Intra': oneCopyList.intersection(allIntraInterGenusSpeciesConditionals.intersection(extraGenusAllZerosConditionals)),
               'Inter': oneCopyList.intersection(allIntraInterGenusSpeciesConditionals.intersection(extraGenusAllOnesConditionals)),
               'Any': oneCopyList.intersection(allIntraInterGenusSpeciesConditionals.intersection(extraGenusMixOnesConditionals)),'Mixed Intra':oneCopyList.intersection(allIntraGenusLowSpeciesConditionals)}

    if subgenome and subgenome != [] and len(subgenome) == 2: #FIXME if one entry is missing in onecopy two copy, then dont output that analysis unless ratio is truly one
        for analysisType in ['CS','ProteinCoding_CS','Mixed_CS','NonProteinCoding_CS']:
            #FIXME start HERE
            a=1
            twoCopyNumber = dict.fromkeys(twoCopy.keys(),[])
            oneCopyNumber = dict.fromkeys(oneCopy.keys(), [])
            for key in twoCopy.keys():
                RatioCopyFile = analysisType+'_ConservationRatios.txt'
                twoCopyNumber[key] = [ratio2CopyStats[conditional][analysisType] for conditional in twoCopy[key]]
                oneCopyNumber[key] = [ratio2CopyStats[conditional][analysisType] for conditional in oneCopy[key]]
            open(conservedFastaPath + RatioCopyFile, 'w').close()
            RatioCopy = open(conservedFastaPath + RatioCopyFile, 'w')
            RatioCopy.write('\n'.join('%s Genus Conservation Ratio: %d'%(key,sum(twoCopyNumber[key])/(sum(twoCopyNumber[key])+sum(oneCopyNumber[key])) ) for key in twoCopy.keys()))
            RatioCopy.write('\n' + '\n'.join('%s Genus Conditionals Included in Calculation:(%s)/(%s)'%(key,
                            '+'.join(conditional for conditional in twoCopy[key]),'+'.join(conditional for conditional
                                    in twoCopy[key])+'+'.join(conditional for conditional in oneCopy[key])) for key in twoCopy.keys()))
            RatioCopy.close()
    # FIXME add weighted ratio copy analysis, and analyzing how much DNA goes into Synteny and comes out of Cactus or my analysis
# FIXME why getting same conservation ratios???

conservedFastaPathFiles = os.listdir(conservedFastaPath)#str(subprocess.Popen(['ls', conservedFastaPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #.stdout.read()).split('\n')
conservedFastaFiles = [fasta for fasta in conservedFastaPathFiles if fasta.endswith('.fasta')]

#FIXME need to be able to run PhyML and tree rendering
if fasta2phylip:
    print 'Fasta to Phylip conversion...','time=',time.clock()-start
    for fasta in conservedFastaFiles:
        with open(conservedFastaPath + fasta, 'r') as f:
            lengthCheck = len(set([len(line) for line in f.read().split('\n')  if line and not line.startswith('>')  ])) == 1
        if lengthCheck:
            subprocess.call(['perl', 'Fasta2Phylip.pl', conservedFastaPath+fasta, conservedFastaPath+fasta.replace('fasta','phylip')])
if phyML:
    print 'Now running PhyML on working outputs... Producing ancestral trees...','time=',time.clock()-start
    phymlLine = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin/phyml'
    phylipFiles = [phylip for phylip in os.listdir(conservedFastaPath) if phylip.endswith('.phylip')]
    #if treeOut == 0:
    for phylip in phylipFiles:
        if phylip:
            subprocess.call([phymlLine, '-i', conservedFastaPath+phylip, '-s', 'BEST', '-q', '-b', bootstrap, '-m', 'GTR']) #FIXME works for now...
            tree = phylip+'_phyml_tree.txt'
            try:
                with open(conservedFastaPath + tree, 'r') as f:  # FIXME boot_trees verus phyml_tree
                    t = Tree(f.read())
                    ts = TreeStyle()
                    ns = NodeStyle()
                    ns['size']=0
                    ts.show_leaf_name = True
                    ts.show_branch_length = False
                    ts.show_branch_support = True
                    for n in t.traverse():
                        n.set_style(ns)
                    #t.show(tree_style=ts)
                    t.render( conservedFastaPath+'/'+tree.replace('_phyml_tree.txt', '.png'),tree_style = ts)
            except:
                pass
    #else:
    #    for phylip in phylipFiles:
    #        if phylip:
    #            subprocess.call(['PhyML', '-i', conservedFastaPath+phylip, '-s', 'BEST', '-q', '-b', bootstrap, '-m', 'GTR','-u',treeFile])
    if 0:#outputTreeImages: #FIXME need to install SIP PyQt4
        print "Generating images for produced trees...",'time=',time.clock()-start
        with open('runTree.sh','w') as f:
            f.write('export PATH=~/anaconda_ete/bin:$PATH\npython treeImage.py')
        with open('treeImage.py','w') as f:
            f.write("""from ete3 import Tree,TreeStyle,NodeStyle
import os
conservedFastaPath = '%s'
for tree in [file for file in os.listdir(conservedFastaPath) if file and file.endswith('_phyml_tree.txt')]:
    try:
        with open(conservedFastaPath + tree, 'r') as f:  # FIXME boot_trees verus phyml_tree
            t = Tree(open(tree,'r').read())
            ts = TreeStyle()
            ns = NodeStyle()
            ns['size']=0
            ts.show_leaf_name = True
            ts.show_branch_length = False
            ts.show_branch_support = True
            for n in t.traverse():
                n.set_style(ns)
            #t.show(tree_style=ts)
            t.render( conservedFastaPath+'/'+tree.replace('_phyml_tree.txt', '.png'),tree_style = ts)
    except:
        pass"""%(conservedFastaPath))
        subprocess.call(['sh', 'runTree.sh'])


# FIXME add user specified output location
# FIXME add tree figure generation PDF!!


print 'End..','time=',time.clock()-start