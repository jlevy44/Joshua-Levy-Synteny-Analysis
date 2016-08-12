import subprocess, os
from pybedtools import *

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname
            configFile.seek(0)
            return line.split()[-1].strip('\n')

def printFirst10Lines(bed):
    count = 0
    for first10lines in str(bed).split('\n'):
        if first10lines:
            if count >= 10:
                break
            count += 1
            print first10lines, count

configFile = open('ghostConfig.txt','r')

findInfoList = ['pathBed','ghostNonRefPath','gffNonRefPath','ghostAnalysisOutputPathAll','ghostAnalysisOutputPathExon',
                'mainSpecies']
            # list of query strings into config path finder

subprocess.call('module load bedtools/2.25.0',shell=True)


for i in range(len(findInfoList)):  # find the paths of above queries
    findInfoList[i] = parseConfigFindPath(findInfoList[i], configFile)
(pathBed,ghostNonRefPath,gffNonRefPath,ghostAnalysisOutputPathAll,ghostAnalysisOutputPathExon,mainSpecies) = \
    tuple(findInfoList)


def findFilesonPath(path,extensionName):
    listRelevantFiles = []
    listFiles = str(subprocess.Popen(['ls', '%s' % path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    .stdout.read()).split('\n')
    for file in listFiles:
        if file.endswith(extensionName):
            listRelevantFiles.append(file)
    return listRelevantFiles


filesPairwiseBed = tuple(findFilesonPath(pathBed,'.bed'))
filesNonrefGFF = tuple(findFilesonPath(gffNonRefPath,'.gff'))


def gffNonRef2bed(gfffile,gffNonrefPath,nonrefPath):
    inputFile = open(gffNonrefPath+gfffile,'r')
    open(nonrefPath + gfffile[:gfffile.rfind('.')] + '.bed', 'w').close()
    outputFile = open(nonrefPath+gfffile[:gfffile.rfind('.')]+'.bed','w')
# psuedomolecule_1	32775	33556	ID=PAC4GC:2552982;Name=Brdisv1Bd2-31000009m;longest=1;Parent=Brdisv1Bd2-31000009m.g

    for line in inputFile:
        if line:
            lineList = line.split()

            outputTuple = (gfffile.split('.')[1],lineList[0],lineList[1],lineList[2],
                           lineList[-1].split(';')[1].replace('Name=',''))
            outputFile.write('%s-%s\t%s\t%s\t%s\n'%outputTuple)

    outputFile.close()
    inputFile.close()

mainSpeciesFileList = []

for file in filesNonrefGFF:
    if mainSpecies in file:
        if 'all' in file:
            # Bd4	JGI	mRNA	6597379	6598760	.	-	.	ID=alift968089m.265;Name=Brdisv1ABR21024638m;Parent=alift968089g.264;peptide_similarity=0.966183574879227;longest=1;pepEntropy=4.12
            outputFile = open(ghostNonRefPath + file[:file.rfind('.')] + '.bed', 'w')
            for line in open(gffNonRefPath + file, 'r'):
                if line:
                    lineList = line.split()
                    outputTuple = ('283', lineList[0], lineList[3], lineList[4],
                                   lineList[8].split(';')[1].replace('Name=',''))
                    outputFile.write('%s-%s\t%s\t%s\t%s\n' % outputTuple)
            outputFile.close()
            mainSpeciesFileList.append(file[:file.rfind('.')] + '.bed')

        if 'exonerate' in file:
            # Bd3	exonerate	match	31061733	31062363	697	+	.	Match 120308_Brdisv1ABR51039317m.faa;start 1;stop 211;coverage 1.00;subjectLen 211;locusId Bd3:1:31061672:31063417
            # get rid of first underscore and .faa
            outputFile = open(ghostNonRefPath + file[:file.rfind('.')] + '.bed', 'w')
            for line in open(gffNonRefPath+file,'r'):
                if line:
                    lineList = line.split()
                    outputTuple = ('283',lineList[0],lineList[3],lineList[4],
                                   lineList[9][lineList[9].find('_')+1:lineList[9].find(';')].replace('.faa',''))
                    outputFile.write('%s-%s\t%s\t%s\t%s\n' % outputTuple)
            outputFile.close()
            mainSpeciesFileList.append(file[:file.rfind('.')] + '.bed')

    else:
        gffNonRef2bed(file, gffNonRefPath, ghostNonRefPath)


filesGhostNonrefBed = tuple(findFilesonPath(ghostNonRefPath,'.bed'))



def reverseBed(bedFilename,bedpath,mainSpecies,targetSpecies):
    #283-Bd1	35904179	35973723	167-ChrC-97998-98793
    reverseBedFilename = bedFilename.replace(mainSpecies,'$').replace(targetSpecies,mainSpecies).replace('$',targetSpecies)
    open(bedpath+reverseBedFilename,'w').close()
    inputBed = open(bedpath+bedFilename,'r')
    global reverseBedFile
    reverseBedFile = open(bedpath+reverseBedFilename,'w')

    for line in inputBed:
        if line:
            lineList = line.split()

            lineList2 = lineList[-1].split('-')
            outputTuple = ('%s-%s'%(lineList2[0],lineList2[1]),lineList2[2],lineList2[3],lineList[0]+
                           '-%s-%s'%(lineList[1],lineList[2]))

            reverseBedFile.write('%s\t%s\t%s\t%s\n'%outputTuple)
    reverseBedFile.close()
    inputBed.close()

    #print open(bedpath+reverseBedFilename,'r').read()

    return reverseBedFilename

listGhostAnalysis = []

# PAC4GC.523-PAC2_0.312_5.bed

for file in filesPairwiseBed:
    if file.endswith('283_10.bed') == 0:
        targetSpecies = file[file.find('-')+1:file.rfind('_')].replace('PAC2_0.','').replace('PAC4GC.','')
        for genesFile in filesGhostNonrefBed:
            if targetSpecies in genesFile:
                targetSpeciesFile = genesFile



        #######
        reverseBedFilename = file.replace(mainSpecies, '$').replace(targetSpecies, mainSpecies).replace('$',targetSpecies)
        open(pathBed + reverseBedFilename, 'w').close()
        inputBed = open(pathBed + file, 'r')
        reverseBedFile = open(pathBed + reverseBedFilename, 'w')

        for line in inputBed:
            if line:
                lineList = line.split()

                lineList2 = lineList[-1].split('-')
                outputTuple = ('%s-%s' % (lineList2[0], lineList2[1]), lineList2[2], lineList2[3], lineList[0] +
                               '-%s-%s' % (lineList[1], lineList[2]))

                reverseBedFile.write('%s\t%s\t%s\t%s\n' % outputTuple)
        reverseBedFile.close()
        inputBed.close()


    #print open(pathBed + reverseBedFilename, 'r').read()

    #####

    listGhostAnalysis.append((targetSpeciesFile,file,reverseBedFilename,
                              targetSpecies))



    #printFirst10Lines(bedPossibleGhostGenes)

    #printFirst10Lines(bedPossibleGhostGenes)
bedPossibleGhostGenes = {}
ghostAnalysisOutputPath = {}
for mainSpeciesFile in mainSpeciesFileList:

    if 'all' in mainSpeciesFile:
        ghostAnalysisOutputPath[mainSpeciesFile] = ghostAnalysisOutputPathAll
    elif 'exonerate' in mainSpeciesFile:
        ghostAnalysisOutputPath[mainSpeciesFile] = ghostAnalysisOutputPathExon
    else:
        os.error('Cannot find proper output Path')
        exit()

    bedPossibleGhostGenes[mainSpeciesFile] = BedTool(ghostNonRefPath + mainSpeciesFile).sort()

for eachGhostAnalysis in listGhostAnalysis:

    bedNonReferencedGenesFile = open(ghostNonRefPath+eachGhostAnalysis[0],'r')
    #print str(bedNonReferencedGenesFile)
    bedPairwiseSynteny = BedTool(pathBed+eachGhostAnalysis[1]).sort()
    bedPairwiseSyntenyReverse = BedTool(pathBed+eachGhostAnalysis[2]).sort()

    #printFirst10Lines(bedPairwiseSynteny)
    #printFirst10Lines(bedPairwiseSyntenyReverse)
    ghostAnalysisOutputFile = {}
    for mainSpeciesFile in mainSpeciesFileList:
        open(ghostAnalysisOutputPath[mainSpeciesFile] + '%s.ghostOutput.txt' % eachGhostAnalysis[-1], 'w').close()
        ghostAnalysisOutputFile[mainSpeciesFile] = open(ghostAnalysisOutputPath[mainSpeciesFile] + '%s.ghostOutput.txt' % eachGhostAnalysis[-1], 'w')

    for line in bedNonReferencedGenesFile:
        if line:
            bedNonRefGene = BedTool(line.strip('\n'),from_string=True)

            #print bedNonRefGene, 'This is nonref'
            #printFirst10Lines(bedNonRefGene)

            #print bedPairwiseSyntenyReverse, 'Reverse pairwise'

            nonrefIntersectReversal = bedPairwiseSyntenyReverse.intersect(bedNonRefGene,wa = True)
            if nonrefIntersectReversal.count():
                #if nonrefIntersectReversal.count() > 1:
                #    print nonrefIntersectReversal.count()
                #printFirst10Lines(nonrefIntersectReversal)
                #print 'Nonref intersect with reversal'

                reverseReversalIntersectionBedText = ''
                for lineRev in str(nonrefIntersectReversal).split('\n'):
                    if lineRev:
                        lineList = lineRev.split()
                        lineList2 = lineList[-1].split('-')
                        outputTuple = ('%s-%s' % (lineList2[0], lineList2[1]), lineList2[2], lineList2[3], lineList[0] +
                                       '-%s-%s' % (lineList[1], lineList[2]))
                        reverseReversalIntersectionBedText += '%s\t%s\t%s\t%s\n' % outputTuple

                reverseReversalIntersectionBed = BedTool(reverseReversalIntersectionBedText,from_string=True)

                #printFirst10Lines(reverseReversalIntersectionBed)

                #intersectWithOriginalBed = bedPairwiseSynteny.intersect(reverseReversalIntersectionBed)
                # is above necessary?

                #printFirst10Lines(intersectWithOriginalBed)

                #printFirst10Lines(bedPossibleGhostGenes.intersect(reverseReversalIntersectionBed,wa=True))


                for mainSpeciesFile in mainSpeciesFileList:
                    intersectionWithGhosts = bedPossibleGhostGenes[mainSpeciesFile].intersect(reverseReversalIntersectionBed,wa=True)
                    listGhostGenes = []


                    #printFirst10Lines(intersectionWithGhosts)
                    if str(intersectionWithGhosts).split('\n')[0]:
                        for lineFinal in str(intersectionWithGhosts).split('\n'):
                            if lineFinal:
                                #print lineFinal
                                listGhostGenes.append(lineFinal.split()[-1])
                        #print line
                        listGhostGenes = list(set(listGhostGenes))
                        ghostAnalysisOutputFile[mainSpeciesFile].write(line.strip('\n') + ';%d;NonRefGhostGeneMatch=%s;%s;%s\n'
                                %(len(listGhostGenes),str(line.split()[-1].strip('\n') in listGhostGenes),
                                ''.join('%s|%s-%s-%s,'%(synSeq.split()[-1],synSeq.split()[0],synSeq.split()[1],
                                synSeq.split()[2]) if synSeq else '' for synSeq in str(reverseReversalIntersectionBedText).split('\n'))[:-1],
                                ''.join(ghostGene+',' for ghostGene in listGhostGenes)[:-1]) )
                    else:
                        ghostAnalysisOutputFile[mainSpeciesFile].write(line.strip('\n') + ';%d\n' % 0)


            else:
                for mainSpeciesFile in mainSpeciesFileList:
                    ghostAnalysisOutputFile[mainSpeciesFile].write(line.strip('\n') + ';%d\n' %0)
    for mainSpeciesFile in mainSpeciesFileList:
        ghostAnalysisOutputFile[mainSpeciesFile].close()
    bedNonReferencedGenesFile.close()

"""


#FIXME need to remove these files after every run
raw_input()
os.chdir(pathBed)
subprocess.call(['rm','*283_10.bed'])
"""