import re
import numpy as np
from pybedtools import *
from pyfaidx import Fasta
import subprocess, os, sys, shutil
from collections import *
import time
import dill as pickle


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

pickleSkip = 0

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
ratioCopy = parseConfigFindPath('ratioCopy',configFile)
configFile.close()

if phyML == '1':
    phyML = 1
else:
    phyML = 0

if ratioCopy == '1':
    ratioCopy = 1
else:
    ratioCopy = 0

if fasta2phylip == '1':
    fasta2phylip = 1
else:
    fasta2phylip = 0

if treeFile != '0' or treeFile != '':
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
            if self.speciesNumber in file:
                self.gffFile = file
        self.speciesName = speciesName
        self.speciesShortName = speciesShortName
        self.conservedElementsBed = '%s_ConservedElements.bed'%self.speciesName
        open(self.conservedElementsBed, 'w').close()
        #self.conservedElementsBedFile = open(self.conservedElementsBed, 'w')

    def gff2bedobjs(self):
        outputfilename = self.gffFile.replace('.gff', '_Genes.bed').replace('.gff3', '_Genes.bed')
        outputfilename2 = self.gffFile.replace('.gff', '_CDS.bed').replace('.gff3', '_CDS.bed')
        open(outputfilename, 'w').close()
        open(outputfilename2, 'w').close()
        outputFile = open(outputfilename, 'w')
        outputFile2 = open(outputfilename2, 'w')
        for line in open(self.gffFile, 'r'):
            if line:
                if 'mRNA' in line and 'longest=1' in line:
                    geneName = line.split()[-1].split(';')[1].replace('Name=','')
                    outputFile.write('%s\t%d\t%s\t%s\n'%(line.split()[0],int(line.split()[3]) - 1,line.split()[4],geneName))

                if 'CDS' in line:
                    if not geneName:
                        geneName = 'NoCNSName'
                    outputFile2.write('%s\t%d\t%s\t%s\n' % (line.split()[0], int(line.split()[3]) -1, line.split()[4], geneName))
                #outputFile.write('%s\t%d\t%s\n' % (line.split()[1],int(line.split()[2])-1,line.split()[3]))
        outputFile.close()
        outputFile2.close()
        self.bedGenes = BedTool(outputfilename).sort()
        self.bedCDS = BedTool(outputfilename2).sort()


def count2Conditional(countSeq,speciesList):
    global speciesInfo
    return ''.join('%s%d'%(speciesInfo[[speciesName for speciesName in speciesInfo.keys() if species in speciesName][0]].speciesShortName,countSeq[species]) for species in speciesList)

speciesInfo = {}

conditionalDictionary = defaultdict(list)



# list all files in analysis directory
listALLFiles = str(subprocess.Popen('ls', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()).split('\n')

# FIXME ['ls', '%s' % '']


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
    """
    if minVal == -1:
        return 16
    else:
        return minVal
    """
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
paths = ['%s/CS/NKH111/'%rootfolder,
         '%s/CS/NKH011/'%rootfolder,
         '%s/CS/NKH101/'%rootfolder,
         '%s/CS/HDS111/'%rootfolder]
#FIXME get rid of this!!! ^^^
# check out each MAF
mafAnalysisStructure = defaultdict(list)
start = time.clock()
writeConservedBedFiles = dict.fromkeys(masterListSpecies,)
if pickleSkip == 0:
    for species in masterListSpecies:
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
                    if mafAnalysisStructure[conditional] == None or mafAnalysisStructure[conditional] == []: #FIXME
                        mafAnalysisStructure[conditional] = {'CS': defaultdict(list),'ProteinCoding_CS': [],'Mixed_CS': [],'NonProteinCoding_CS': []}
                    mafAnalysisStructure[conditional]['CS']['MASeg%d'%segmentCount] = ''.join('>%s_%s_%s_%s\n%s\n'%(tuple(line.split()[1:5])+(line.split()[-1],)) for species in [species for species in countSeq.keys() if countSeq[species] == 1] for line in segment.split('\n') if line and species in line)
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
                    """
                    for line in segment.split('\n'):
                        for species in [species for species in listIntraGenus if countSeq[species] == 1]:
                            if line and species in line:
                    """
                    mafAnalysisStructure[conditional]['CS']['MASeg%d'%segmentCount] = ''.join('>%s_%s_%s_%s\n%s\n'%(tuple(line.split()[1:5])+(line.split()[-1],)) for species in [species for species in listIntraGenus if countSeq[species] == 1] for line in segment.split('\n') if line and species in line)
                    a=1
                """
                if countSeq['PvirgatumN'] == countSeq['PvirgatumK'] == countSeq['Phallii'] == 1:
                    open(paths[0]+'MASeg%d.fasta'%segmentCount,'w').close()
                    MAsegFile = open(paths[0]+'MASeg%d.fasta'%segmentCount,'w')
                    for line in segment.split('\n'):
                        if line and ('Pvirgatum' in line or 'Phallii' in line):
                            MAsegFile.write('>%s_%s_%s_%s\n%s\n'%(tuple(line.split()[1:5])+(line.split()[-1],)))
                    MAsegFile.close()
                # similar for NKH011
                elif countSeq['PvirgatumK'] == countSeq['Phallii'] == 1 and countSeq['PvirgatumN'] == 0:
                    open(paths[1] + 'MASeg%d.fasta' % segmentCount, 'w').close()
                    MAsegFile = open(paths[1] + 'MASeg%d.fasta' % segmentCount, 'w')
                    for line in segment.split('\n'):
                        if line and ('Pvirgatum' in line or 'Phallii' in line):
                            MAsegFile.write('>%s_%s_%s_%s\n%s\n' % (tuple(line.split()[1:5]) + (line.split()[-1],)))
                    MAsegFile.close()
                # NKH101
                elif countSeq['PvirgatumN'] == countSeq['Phallii'] == 1 and countSeq['PvirgatumK'] == 0:
                    open(paths[2] + 'MASeg%d.fasta' % segmentCount, 'w').close()
                    MAsegFile = open(paths[2] + 'MASeg%d.fasta' % segmentCount, 'w')
                    for line in segment.split('\n'):
                        if line and ('Pvirgatum' in line or 'Phallii' in line):
                            MAsegFile.write('>%s_%s_%s_%s\n%s\n' % (tuple(line.split()[1:5]) + (line.split()[-1],)))
                    MAsegFile.close()
                # HDS111 intragenus
                if countSeq['Bdistachyon'] == countSeq['Bstacei'] == countSeq['Phallii'] == 1:
                    open(paths[3] + 'MASeg%d.fasta' % segmentCount, 'w').close()
                    MAsegFile = open(paths[3] + 'MASeg%d.fasta' % segmentCount, 'w')
                    for line in segment.split('\n'):
                        if line and ('Bdistachyon' in line or 'Bstacei' in line or 'Phallii' in line):
                            MAsegFile.write('>%s_%s_%s_%s\n%s\n' % (tuple(line.split()[1:5]) + (line.split()[-1],)))
                    MAsegFile.close()
                    """

            #FIXME _ to -


            segmentCount+=1 # add count to number of segments to generate multiple alignment segment IDs


                                    #, lineList[4],''.join(speciesName2 + ',' for speciesName2 in speciesList)[:-1]))
                            #FIXME '_' or '-' ^^^^^
            #for infoTuple in outputInfo:
            #write -->    (infoTuple+(''.join(speciesName + ',' for speciesName in speciesList)[:-1],)))#str(Pvir[infoTuple[0]][infoTuple[1]:infoTuple[2]]))))



            """a
            s	Anc0.Anc0refChr0	115	8	+	10715	TCAGCTTA
            s	Anc1.Anc1refChr0	123	8	+	57143	TCAGCTTA
            s	Anc3.Anc3refChr0	123	8	+	22297	TCAGCTTA
            s	Anc4.Anc4refChr0	123	8	+	55673	TCAGCTTA
            s	Sbicolor.Sbicolor.Chr09_952912_1185299	152	8	+	232387	TCAGCTTA
            s	Anc6.Anc6refChr2	143	8	-	27358	ACAGCTTA
            s	Anc7.Anc7refChr0	167	8	+	30340	TCAGCTTA
            s	Sitalica.Sitalica.scaffold_3_4088912_4264959	147	8	-	176047	ACAGCTTA
            s	Pvirgatum.Pvirgatum.Chr01K_14565140_14836634	217	8	+	271494	TCAGCTTA
            s	Phallii.Phallii.Chr03_5677168_5837565	166	8	+	160397	TCAGCTTA
            s	Anc2.Anc2refChr0	147	8	+	38180	TCAGCTTA
            s	Osativa.Osativa.Chr5_516915_603321	172	8	+	86406	TCAGCTTA
            s	Anc5.Anc5refChr7	69	8	-	22721	TCAGCTTA
            s	Bdistachyon.Bdistachyon.Bd2_39074157_39126388	133	8	-	52231	TCAGCTTA
            s	Bstacei.Bstacei.Chr08_20211117_20262991	83	8	-	51874	TCAGCTTA"""
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
    bedCSMinusCDS = {}

    for species in masterListSpecies:

        conservedElementsBed = open(speciesInfo[species].conservedElementsBed,'r')

        # create bed object from conserved elements for each species, conserved elements found from MAF files done above
        speciesInfo[species].conservedElementsBed = BedTool(conservedElementsBed).sort()
        conservedElementsBed.close()
        # generate CDS and gene bed objects for each species
        speciesInfo[species].gff2bedobjs()

        #bedCNS = {}
        # find CNS Intronic genes = CS + (genes-CDS)  + is intersection, - is subtraction of sequences
        bedCNSIntronic[species] = speciesInfo[species].conservedElementsBed.intersect(speciesInfo[species].bedGenes.
                subtract(speciesInfo[species].bedCDS).sort(),u=True).sort().closest(speciesInfo[species].bedGenes,nonamecheck = True)
        # CNS intergenic = CS - genes
        bedCNSIntergenic[species] = speciesInfo[species].conservedElementsBed.subtract(speciesInfo[species].bedGenes).sort()\
                                                .sort().closest(speciesInfo[species].bedGenes,k=1,d=True,nonamecheck = True)

        # Conserved CDS = CS +CDS
        bedConservedCDS[species] = speciesInfo[species].conservedElementsBed.intersect(speciesInfo[species].bedCDS,
                                                u=True).sort().closest(speciesInfo[species].bedGenes,nonamecheck = True)

        # This is both CNS Intronic + intergenic, [CNS total] = CS-CDS = CS - CDS
        bedCSMinusCDS[species] = speciesInfo[species].conservedElementsBed.subtract(speciesInfo[species].bedCDS,
                                                u=True).sort().closest(speciesInfo[species].bedGenes,k=1,d=True,nonamecheck = True)


    # .merge(o='distinct',c=4,delim = '|')

    #bedCNSFinal = bedCNS[bedCNS.keys()[0]].sort()

        #,str(Pvir[infoTuple[0]][infoTuple[1]:infoTuple[2]]

        #for key in bedCNS.keys()[1:]:
        #    bedCNSFinal = bedCNSFinal.cat(bedCNS[key],postmerge = False).sort()
        #FIXME start here!!!
        #bedCNSFinal = bedCNSFinal.merge(o='distinct',c=4,delim = '|').sort()
        CNSOutFiles = [('%s_CNSElements_Intronic.bed'%speciesInfo[species].speciesName,bedCNSIntronic[species]),
                       ('%s_CNSElements_Intergenic.bed'%speciesInfo[species].speciesName,bedCNSIntergenic[species]),
                       ('%s_Conserved_CDS.bed' % speciesInfo[species].speciesName,bedConservedCDS[species]),
                       ('%s_CS_Minus_CDS.bed' % speciesInfo[species].speciesName,bedCSMinusCDS[species])]
        for bedOut in CNSOutFiles: # for each of the above described files, append the sequences from the original fastas and include closest gene/ distance information
            open(bedOut[0],'w').close()

            bedOutFile = open(bedOut[0],'w')

            for line in str(bedOut[1]).split('\n'):
                if line:
                    lineList = line.split()
                    outputSequence = str(speciesInfo[species].genome[lineList[0]][int(lineList[1]):int(lineList[2])])
                    # need to have a long enough sequence to get outputted...
                    if len(outputSequence) >= 15 and findBadCharPosition(outputSequence) >= 15:
                        if 'CNSElements_Intronic' in bedOut[0] or 'Conserved_CDS' in bedOut[0]:
                            bedOutFile.write('%s\t%s\t%s\t%s;geneID=%s;%s\n'%(tuple(lineList[0:4])+(lineList[-1],outputSequence)))
                        elif 'CNSElements_Intergenic' in bedOut[0] or 'CS_Minus_CDS' in bedOut[0]:
                            bedOutFile.write('%s\t%s\t%s\t%s;closestGene=%s;distance=%s;%s\n'%(tuple(lineList[0:4])+(lineList[-2],lineList[-1].strip('\n'),outputSequence)))



        #bedOut.write(str(bedCNSFinal))

            bedOutFile.close()

        conservedElementsBed.close()

        """
        open('%s_Conserved_CDS.bed' % speciesInfo[species].speciesName, 'w').close()
        conservedCDSFile = open('%s_Conserved_CDS.bed'%speciesInfo[species].speciesName,'w')
        for line in str(bedConservedCDS[species]).split('\n'):
            if line:
                conservedCDSFile.write('%s\t%s\t%s\t%s;geneID=%s\n'%(tuple(line.split()[0:4]) + (line.split()[-1],)))

        conservedCDSFile.close()
        """

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

    def grabPartMAFSeq(MAFSeq,x_coords): # irrelevant script now...
        coordtransform = np.array(range(len(MAFSeq)))
        coordtransform = np.delete(coordtransform,[val.start(0) for val in re.finditer('-',MAFSeq)])
        print coordtransform, coordtransform[x_coords[1]-1], x_coords[1]
        return MAFSeq[coordtransform[x_coords[0]]:coordtransform[x_coords[1]-1]]

    ultimateBedDicts = {'Conserved_CDS.bed':{},'CS_Minus_CDS.bed':{}} # split into ~CDS or ~CNS

    for file in listALLFiles:
        if file.endswith('Conserved_CDS.bed') or file.endswith('CS_Minus_CDS.bed'):
            ultimateBedDicts[file[file.find('_')+1:]][file[:file.find('_')]] = bed2dict(file) # a dictionary with keys [species] and inside a key of MASEG
        # hierarchy is ultimate[CDS or CNS][species][MASeg]
    pickle.dump(tuple([ultimateBedDicts,mafAnalysisStructure, speciesInfo, conditionalDictionary]), open(pickleName, 'wb'))


if pickleSkip == 1:
    (ultimateBedDicts,mafAnalysisStructure, speciesInfo, conditionalDictionary) = pickle.load(open(pickleName,'rb'))

print 'time=',time.clock()-start,'size=',sys.getsizeof(mafAnalysisStructure)
def conditional2SpeciesList(conditionalDict): # code obsolete
    global masterListSpecies
    # get dictionary keys that have 1, convert those short names to species names
    return np.vectorize(lambda shortName: [species.split('_')[0] for species in masterListSpecies if shortName ==
                                           species.split('_')[-1].strip('\n')][0])\
        ([shortname for shortname in conditionalDict.keys() if conditionalDict[shortname]])

ratio2CopyStats = dict.fromkeys(mafAnalysisStructure.keys(),{'CS':0,'ProteinCoding_CS':0,
                     'Mixed_CS':0,'NonProteinCoding_CS':0})
# generate final fastas
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
        finalFastaFileName = conditional+'_'+analysisType+'.fasta'
        open(conservedFastaPath+finalFastaFileName,'w').close()
        outputFinalFastaFile = open(conservedFastaPath+finalFastaFileName,'w')
        outputFinalFastaFile.write(''.join('>%s\n%s\n'%(species,conservedSeqs[analysisType][species]) for species in conservedSeqs[analysisType].keys() if conservedSeqs[analysisType][species]))
        outputFinalFastaFile.close()
        if analysisType == 'CS':
            ratio2CopyStats[conditional][analysisType]= len(mafAnalysisStructure[conditional][analysisType].keys())
        else:
            ratio2CopyStats[conditional][analysisType]= len(mafAnalysisStructure[conditional][analysisType])
#FIXME everything below needs to be fixed...
# output two copy statistics


if ratioCopy:
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

conservedFastaPathFiles = str(subprocess.Popen(['ls', conservedFastaPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        .stdout.read()).split('\n')
conservedFastaFiles = [fasta for fasta in conservedFastaPathFiles if fasta.endswith('.fasta')]

if fasta2phylip:
    for fasta in conservedFastaFiles:
        subprocess.call(['perl', 'Fasta2Phylip.pl', conservedFastaPath+fasta, conservedFastaPath+fasta.replace('fasta','phylip')])
if phyML:
    phylipFiles = [fasta.replace('fasta','phylip') for fasta in conservedFastaFiles]
    if treeOut == 0:
        for phylip in phylipFiles:
            subprocess.call(['PhyML', '-i', conservedFastaFiles+phylip, '-s', 'BEST', '-q', '-b', bootstrap, '-m', 'GTR'])
    else:
        for phylip in phylipFiles:
            subprocess.call(['PhyML', '-i', conservedFastaFiles+phylip, '-s', 'BEST', '-q', '-b', bootstrap, '-m', 'GTR','-u',treeFile])
# FIXME add user specified output location


"""
        for species in conservedSeqs.keys():
            outputFinalFastaFile.write('>%s\n%s\n'%(species,''.join(conservedSeq for conservedSeq in conservedSeqs[species])))

"""
#FIXME start here
"""
#below is obselete
def concatFastas(fname,HDSPath,NKHPath,FinalPath): # if intergenus conservation, concatenate fastas together and add it to relevant path
    NKHFileStringList = open(NKHPath+fname,'r').read().split('>')
    HDSFileStringList = open(HDSPath+fname,'r').read().split('>')
    outputStringList = list(set(NKHFileStringList+HDSFileStringList))
    outputFile = open(FinalPath+fname,'w')
    outputFile.write('>'.join(seqs for seqs in outputStringList))
    outputFile.close()

#below is obsolete
def grabConservedSeq(PathAndfname,conservedSeqs): # build a list of conserved sequences that are either Conserved sequence,
    # protein coding, nonprotein, or mixed and will join list together
    inputFastaReadSegments = open(PathAndfname,'r').read().split('>')
    for sequence in inputFastaReadSegments:
        if sequence:
            speciesName = sequence.split('.')[0]
            if 'K' == sequence.split('_')[0][-1]: #FIXME cahnge _ to -
                speciesName += 'K'
            if 'N' == sequence.split('_')[0][-1]:
                speciesName += 'N'
            if conservedSeqs.has_key(speciesName):
                conservedSeqs[speciesName].append(sequence.split('\n')[1].strip('\n')) # add sequence to main list to generate main fasta file output
            else:
                conservedSeqs[speciesName] = [sequence.split('\n')[1].strip('\n')]
    return conservedSeqs

for analysisType in ['CS','ProteinCoding_CS','Mixed_CS','NonProteinCodingN_CS']: # run analysis on these types and sort them into relevant conditionals and output fastas for PhyML analysis and conservation ratios
    # see more details about ratio conservation in various documents...
    outFname = ['NKHDS11111.txt','NKHDS01111.txt','NKHDS10111.txt','NKHDS11100.txt','NKHDS01100.txt','NKHDS10100.txt','NKH111.txt','NKH011.txt','NKH101.txt','HDS111.txt']
    if analysisType == 'CS':
        analysisType2 = ''
        arrayFiles = []
        for path in paths:
            arrayFiles.append(str(subprocess.Popen(['ls', '%s' % path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    .stdout.read()).split('\n'))#.replace('CS',analysisType)
        #print analysisPath
        for i in range(3):
            np.savetxt(analysisPath+outFname[i].replace('.txt','_%s.txt'%analysisType),np.intersect1d(arrayFiles[3],arrayFiles[i]),fmt='%s',comments='')
            np.savetxt(analysisPath+outFname[i+3].replace('.txt','_%s.txt'%analysisType),np.setdiff1d(arrayFiles[i],arrayFiles[3]),fmt='%s',comments='')
            np.savetxt(analysisPath+outFname[i+6].replace('.txt','_%s.txt'%analysisType),arrayFiles[i],fmt='%s',comments='')

        np.savetxt(analysisPath+outFname[9].replace('.txt','_%s.txt'%analysisType),arrayFiles[3],fmt='%s',comments='')
        # the above script finds the intersection and differences
    else:
        analysisType2 = analysisType
        for fname in outFname[0:6]:
            outFile = open(analysisPath+fname.replace('.txt','_%s.txt'%analysisType),'w')
            outFile.write(str(subprocess.Popen('ls %s/%s/%s | grep "%s"'%(rootfolder,'CS',fname.replace('.txt',''),
                    analysisType),shell=True , stdout=subprocess.PIPE).stdout.read()))
            outFile.close() # the CS analysis separates NHK files into various NKHDS folders, now they are named according to type of analysis

    if analysisType == 'CS':
        # if in CS, sort NKH files into various NKHDS folders according to conditional statements and which files fall under each conditional
        for fname in outFname[0:6]:
            inputFile = open(analysisPath+fname.replace('.txt','_%s.txt'%analysisType),'r')
            conditionalState = fname.split('.')[0].strip('NKHDS')
            if len(conditionalState) == 5 and conditionalState.endswith('11'):
                HDSpath = '%s/%s/HDS111/'%(rootfolder,'CS')
                NKHpath = '%s/%s/NKH%s/'%(rootfolder,'CS',conditionalState[0:3])
                finalPath = '%s/%s/NKHDS%s/'%(rootfolder,'CS',conditionalState)
                for line in inputFile:
                    if line and line != '\n' and ':' not in line and '/' not in line:
                        inputNKHDSfname = line.strip('\n')
                        concatFastas(inputNKHDSfname,HDSpath,NKHpath,finalPath)
            elif len(conditionalState) == 5 and conditionalState.endswith('00'):
                NKHpath = '%s/%s/NKH%s/'%(rootfolder,'CS',conditionalState[0:3])
                finalPath = '%s/%s/NKHDS%s/'%(rootfolder,'CS',conditionalState)
                for line in inputFile:
                    if line and line != '\n' and ':' not in line and '/' not in line:
                        inputNKHDSfname = line.strip('\n')
                        #print line, fname.replace('.txt','_%s.txt'%analysisType)
                        #print NKHpath, inputNKHDSfname, finalPath
                        shutil.copy(NKHpath+inputNKHDSfname,finalPath)
            inputFile.seek(0)
            inputFile.close()
    for fname in outFname[0:6]: # compute twocopy statistics (ratio) between NKHDS, see output text
        # output final fasta for NKHDS conditional type and analysis type seen above, which generate trees for each
        fileArray2 = open(analysisPath+fname.replace('.txt','_%s.txt'%analysisType),'r').read().splitlines()
        speciesRatios = fname.replace('.txt','')
        TwoCopyStatNumbers[speciesRatios] = float(len(fileArray2))
        conservedSeqs = {}
        for file in fileArray2:
            if file and ':' not in file and '/' not in file:
                try:
                    conservedSeqs = grabConservedSeq(analysisPath+speciesRatios+'/'+file,conservedSeqs)
                except:
                    print analysisType
        outputFinalFasta = fname.replace('.txt','_%s.fasta'%analysisType)
        open(analysisPath+outputFinalFasta,'w').close()
        outputFinalFastaFile = open(analysisPath+outputFinalFasta, 'w')
        for species in conservedSeqs.keys():
            outputFinalFastaFile.write('>%s\n%s\n'%(species,''.join(conservedSeq for conservedSeq in conservedSeqs[species])))
        outputFinalFastaFile.close()

    ratio2copyStat = [TwoCopyStatNumbers['NKHDS11111']/(TwoCopyStatNumbers['NKHDS10111']+TwoCopyStatNumbers['NKHDS01111']+TwoCopyStatNumbers['NKHDS11111']),
                      TwoCopyStatNumbers['NKHDS11100'] / (TwoCopyStatNumbers['NKHDS10100'] + TwoCopyStatNumbers['NKHDS01100']+TwoCopyStatNumbers['NKHDS11100'])]
    open(analysisPath+'Twocopyratios_%s.txt'%analysisType,'w').close()
    twoCopyFile = open(analysisPath+'Twocopyratios_%s.txt'%analysisType,'w')
    twoCopyFile.write('Intergenus Conservation [NKHDS11111/(NKHDS11111+NKHDS10111+NKHDS01111)]: %f\n'
                      'Intragenus Conservation [NKHDS11100/(NKHDS11100+NKHDS10100+NKHDS01100)]: %f\n'%(tuple(ratio2copyStat)))
    twoCopyFile.close()
    if analysisType == 'CS':
        for path in ['%s/%s/%s/'%(rootfolder,'CS',fname.replace('.txt','')) for fname in outFname[0:6]]:
            for file in str(subprocess.Popen(['ls', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                    .stdout.read()).split('\n'):
                if file and file.endswith('.fasta'):
                    """"""">Pvirgatum.Pvirgatum.Chr01K_14565140_14836634_328_23_+
                    AGACAGCTAACAGCTTACTTTCT
                    >Phallii.Phallii.Chr03_5677168_5837565_290_23_+
                    AGACAGCTAACAGCTTACTTTCT
                    """"""
                    inputFasta = open(path + file, 'r')  # 011
                    MASeg = file.split('.')[0]
                    speciesList = []
                    for line in inputFasta:
                        if '>' in line:
                            lineList = line[1:].split('.')
                            speciesName = lineList[0]
                            lineList2 = lineList[-1].split('_')
                            scaffold_name = lineList[-1][:lineList[-1].find(lineList2[-5]) - 1]
                            # orientation = lineList2[-1].strip('\n')
                            if 'K' == scaffold_name[-1]:
                                speciesName += 'K'
                            if 'N' == scaffold_name[-1]:
                                speciesName += 'N'
                            """"""
                            if orientation == '+':
                                x = (int(lineList2[-5])+int(lineList2[-3]),int(lineList2[-5])+int(lineList2[-3])+int(lineList2[-2]))
                                print (scaffold_name,),x
                            else:
                                x = (int(lineList2[-4])-int(lineList2[-3])-int(lineList2[-2]),int(lineList2[-4])-int(lineList2[-3]))
                                print (scaffold_name,),x
                            """"""
                        else:
                            speciesList.append(speciesName)
                            # for key in ultimateBedDicts.keys():
                            #    if ultimateBedDicts[key][speciesName].has_key(MASeg):
                            #        for CNSCDS in ultimateBedDicts[key][speciesName][MASeg]:
                            # xRel = tuple(np.vectorize(lambda y: int(y))(CNSCDS[1:])-x[0])
                            # print ultimateBedDicts[key][speciesName][MASeg],conservedSeg, len(conservedSeg), xRel,grabPartMAFSeq(conservedSeg,xRel)
                    CDS_conditional = np.vectorize(
                        lambda speciesName: ultimateBedDicts['Conserved_CDS.bed'][speciesName].has_key(MASeg))(speciesList)
                    # in each folder for the NKHDS conditionals, change file names to match relevant analysis, eg. whether
                    # sequence intersects a CDS, is not, or some do and some don't between different species
                    if np.all(CDS_conditional):
                        os.rename(path + file, path + file.replace('.fasta', '_ProteinCoding_CS.fasta'))
                    elif np.any(CDS_conditional) and not np.all(CDS_conditional):
                        os.rename(path + file, path + file.replace('.fasta', '_Mixed_CS.fasta'))
                    elif not np.any(CDS_conditional):
                        os.rename(path + file, path + file.replace('.fasta', '_NonProteinCodingN_CS.fasta'))
                    inputFasta.close()
"""
print 'time=',time.clock()-start