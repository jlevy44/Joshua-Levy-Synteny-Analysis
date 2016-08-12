#import re
import numpy as np
from pybedtools import *
from pyfaidx import Fasta
import subprocess, os, sys, shutil

from difflib import SequenceMatcher

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
# help from online sources ^^^^

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname
            configFile.seek(0)
            return line.split()[-1].strip('\n')
configFile = open('configCNSAnalysis.txt','r')
rootfolder = parseConfigFindPath('root_folder',configFile)
pathPython = parseConfigFindPath('pathPython',configFile)
configFile.close()

sys.path.append(pathPython)


class speciesClass():
    def __init__(self,speciesNumber,genomeFileList,gffFileList,speciesName):
        self.speciesNumber = speciesNumber
        for file in genomeFileList:
            if self.speciesNumber in file:
                self.genome = Fasta(file)
        for file in gffFileList:
            if self.speciesNumber in file:
                self.gffFile = file
        self.speciesName = speciesName
        self.conservedElementsBed = '%s_ConservedElements.bed'%self.speciesName
        open(self.conservedElementsBed, 'w').close()
        self.conservedElementsBedFile = open(self.conservedElementsBed, 'w')

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


speciesInfo = {}

masterListSpecies = ['Sbicolor_313','Sitalica_312','PvirgatumK_383','PvirgatumN_383','Phallii_308','Osativa_323','Bdistachyon_314','Bstacei_316']

checkValidity = 0

listALLFiles = str(subprocess.Popen('ls', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()).split('\n')

# FIXME ['ls', '%s' % '']


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

for species in masterListSpecies:
    speciesInfo[species] = speciesClass(species.split('_')[1], listGenomeFiles, listGFFFiles, species.split('_')[0])
"""
def turnSixteen(x):
    if x == -1:
        return 16
    else:
        return x
"""
def findBadCharPosition(strSeq):
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
if checkValidity:
    open('CheckValidity.txt','w').close()
    checkValidFile = open('CheckValidity.txt','w')

segmentCount = 0
paths = ['%s/CS/NKH111/'%rootfolder,
         '%s/CS/NKH011/'%rootfolder,
         '%s/CS/NKH101/'%rootfolder,
         '%s/CS/HDS111/'%rootfolder]
for file in listMAFfiles:
    inputMAF = open(file,'r')

    for segment in inputMAF.read().split('\n\n'):
        writeTheseLines = {}
        countSeq = {}
        skipSegment = 0
        for species in masterListSpecies:
            countSeq[species.split('_')[0]] = 0# FIXME change
        speciesList = []
        outputInfo = []
        for line in segment.split('\n'):
            if line and skipSegment == 0:
                if line[0] == 's' and 'Anc' not in line.split()[1]:
                    if int(line.split()[3]) >= 20 and findBadCharPosition(line.split()[-1]) >= 15:
                        lineList = line.split()
                        lineList2 = lineList[1].split('.')
                        speciesName = lineList2[0]
                        if line.split('_')[0].endswith('K'):  # FIXME change _ to -
                            speciesName += 'K'
                        elif line.split('_')[0].endswith('N'):
                            speciesName += 'N'
                        if speciesName in speciesList:
                            countSeq[speciesName] += 1
                        else:
                            speciesList.append(speciesName)
                            writeTheseLines[speciesName] = []
                            countSeq[speciesName] = 1
                        lineList3 = lineList2[-1].split('_')
                        if lineList[4] == '-':
                            startCoord,endCoord = int(lineList3[-1])-int(lineList[2])-int(lineList[3]),int(lineList3[-1])-int(lineList[2])
                        else:
                            startCoord, endCoord = int(lineList3[-2])+int(lineList[2]),int(lineList3[-2])+int(lineList[2])+int(lineList[3])
                        writeTheseLines[speciesName].append('%s\t%d\t%d\t%s;'
                                    %(lineList2[-1][:lineList2[-1].find(lineList3[-2])-1],startCoord,endCoord,lineList[4]))
                    else:
                        skipSegment = 1
                    if checkValidity:
                        for species in masterListSpecies:
                            if speciesName in species:
                                MAFseq = lineList[-1].replace('-','')
                                if lineList[4] == '-':
                                    conservedSeq = str(speciesInfo[species].genome[lineList2[-1][:lineList2[-1].find(lineList3[-2])-1]][startCoord:endCoord].reverse.complement)
                                else:
                                    conservedSeq = str(speciesInfo[species].genome[lineList2[-1][:lineList2[-1].find(lineList3[-2])-1]][startCoord:endCoord])
                                checkValidFile.write('%s\tMAF=%s\tConservedSequence=%s\tOrientation=%s\tPercentSimilar=%d%%\n' %(speciesName,MAFseq,conservedSeq,lineList[4],similar(MAFseq,conservedSeq)*100.))
        if skipSegment == 0:
            for speciesName in speciesList:
                for species in masterListSpecies:
                    if speciesName in species:
                        for line in writeTheseLines[speciesName]:
                            speciesInfo[species].conservedElementsBedFile.write(line+''.join('%s:%d'%(speciesName2,countSeq[speciesName2]) + ',' for speciesName2 in speciesList)[:-1]+';SegmentID=MASeg%d'%segmentCount+'\n')

            if countSeq['PvirgatumN'] == countSeq['PvirgatumK'] == countSeq['Phallii'] == 1:
                open(paths[0]+'MASeg%d.fasta'%segmentCount,'w').close()
                MAsegFile = open(paths[0]+'MASeg%d.fasta'%segmentCount,'w')
                for line in segment.split('\n'):
                    if line and ('Pvirgatum' in line or 'Phallii' in line):
                        MAsegFile.write('>%s_%s_%s_%s\n%s\n'%(tuple(line.split()[1:5])+(line.split()[-1],)))
                MAsegFile.close()
            elif countSeq['PvirgatumK'] == countSeq['Phallii'] == 1 and countSeq['PvirgatumN'] == 0:
                open(paths[1] + 'MASeg%d.fasta' % segmentCount, 'w').close()
                MAsegFile = open(paths[1] + 'MASeg%d.fasta' % segmentCount, 'w')
                for line in segment.split('\n'):
                    if line and ('Pvirgatum' in line or 'Phallii' in line):
                        MAsegFile.write('>%s_%s_%s_%s\n%s\n' % (tuple(line.split()[1:5]) + (line.split()[-1],)))
                MAsegFile.close()
            elif countSeq['PvirgatumN'] == countSeq['Phallii'] == 1 and countSeq['PvirgatumK'] == 0:
                open(paths[2] + 'MASeg%d.fasta' % segmentCount, 'w').close()
                MAsegFile = open(paths[2] + 'MASeg%d.fasta' % segmentCount, 'w')
                for line in segment.split('\n'):
                    if line and ('Pvirgatum' in line or 'Phallii' in line):
                        MAsegFile.write('>%s_%s_%s_%s\n%s\n' % (tuple(line.split()[1:5]) + (line.split()[-1],)))
                MAsegFile.close()

            if countSeq['Bdistachyon'] == countSeq['Bstacei'] == countSeq['Phallii'] == 1:
                open(paths[3] + 'MASeg%d.fasta' % segmentCount, 'w').close()
                MAsegFile = open(paths[3] + 'MASeg%d.fasta' % segmentCount, 'w')
                for line in segment.split('\n'):
                    if line and ('Bdistachyon' in line or 'Bstacei' in line or 'Phallii' in line):
                        MAsegFile.write('>%s_%s_%s_%s\n%s\n' % (tuple(line.split()[1:5]) + (line.split()[-1],)))
                MAsegFile.close()
        #FIXME _ to -


        segmentCount+=1


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



    inputMAF.close()

if checkValidity:
    checkValidFile.close()

bedCNSIntergenic = {}
bedCNSIntronic = {}
bedConservedCDS = {}
bedCSMinusCDS = {}

for species in masterListSpecies:

    speciesInfo[species].conservedElementsBedFile.close()

    speciesInfo[species].conservedElementsBed = BedTool(speciesInfo[species].conservedElementsBed).sort()

    speciesInfo[species].gff2bedobjs()

    #bedCNS = {}
    bedCNSIntronic[species] = speciesInfo[species].conservedElementsBed.intersect(speciesInfo[species].bedGenes.
            subtract(speciesInfo[species].bedCDS).sort(),u=True).sort().closest(speciesInfo[species].bedGenes,nonamecheck = True)
    bedCNSIntergenic[species] = speciesInfo[species].conservedElementsBed.subtract(speciesInfo[species].bedGenes).sort()\
                                            .sort().closest(speciesInfo[species].bedGenes,k=1,d=True,nonamecheck = True)
    bedConservedCDS[species] = speciesInfo[species].conservedElementsBed.intersect(speciesInfo[species].bedCDS,
                                            u=True).sort().closest(speciesInfo[species].bedGenes,nonamecheck = True)
    bedCSMinusCDS[species] = speciesInfo[species].conservedElementsBed.subtract(speciesInfo[species].bedCDS,
                                            u=True).sort().closest(speciesInfo[species].bedGenes,k=1,d=True,nonamecheck = True)


    # .merge(o='distinct',c=4,delim = '|')

    #bedCNSFinal = bedCNS[bedCNS.keys()[0]].sort()

    #,str(Pvir[infoTuple[0]][infoTuple[1]:infoTuple[2]]

    #for key in bedCNS.keys()[1:]:
    #    bedCNSFinal = bedCNSFinal.cat(bedCNS[key],postmerge = False).sort()

    #bedCNSFinal = bedCNSFinal.merge(o='distinct',c=4,delim = '|').sort()
    CNSOutFiles = [('%s_CNSElements_Intronic.bed'%speciesInfo[species].speciesName,bedCNSIntronic[species]),
                   ('%s_CNSElements_Intergenic.bed'%speciesInfo[species].speciesName,bedCNSIntergenic[species]),
                   ('%s_Conserved_CDS.bed' % speciesInfo[species].speciesName,bedConservedCDS[species]),
                   ('%s_CS_Minus_CDS.bed' % speciesInfo[species].speciesName,bedCSMinusCDS[species])]
    for bedOut in CNSOutFiles:
        open(bedOut[0],'w').close()

        bedOutFile = open(bedOut[0],'w')

        for line in str(bedOut[1]).split('\n'):
            if line:
                lineList = line.split()
                outputSequence = str(speciesInfo[species].genome[lineList[0]][int(lineList[1]):int(lineList[2])])
                if len(outputSequence) >= 15 and findBadCharPosition(outputSequence) >= 15:
                    if 'CNSElements_Intronic' in bedOut[0] or 'Conserved_CDS' in bedOut[0]:
                        bedOutFile.write('%s\t%s\t%s\t%s;geneID=%s;%s\n'%(tuple(lineList[0:4])+(lineList[-1],outputSequence)))
                    elif 'CNSElements_Intergenic' in bedOut[0] or 'CS_Minus_CDS' in bedOut[0]:
                        bedOutFile.write('%s\t%s\t%s\t%s;closestGene=%s;distance=%s;%s\n'%(tuple(lineList[0:4])+(lineList[-2],lineList[-1].strip('\n'),outputSequence)))



    #bedOut.write(str(bedCNSFinal))

        bedOutFile.close()

    """
    open('%s_Conserved_CDS.bed' % speciesInfo[species].speciesName, 'w').close()
    conservedCDSFile = open('%s_Conserved_CDS.bed'%speciesInfo[species].speciesName,'w')
    for line in str(bedConservedCDS[species]).split('\n'):
        if line:
            conservedCDSFile.write('%s\t%s\t%s\t%s;geneID=%s\n'%(tuple(line.split()[0:4]) + (line.split()[-1],)))

    conservedCDSFile.close()
    """


execfile('manipulateFilesTest.py')




