from pyfaidx import Fasta
import pandas as pd
import numpy as np
import os, sys
from pybedtools import *

"""Initialize the syntenic structure; return is the initialized structure..."""

def generateFastaObjectStructure(listOfGenomeFiles,genomePath):
    """Generates dictionary of different species that contains a list of chromosomes for that species and a fasta object
    for that species"""
    fastaObjectStructure={}

    for genomeFile in listOfGenomeFiles:
        #generate a fasta object for given genome
        fastaGenome=Fasta(genomePath+genomeFile)
        speciesName = genomeFile[genomeFile.find('_') + 1:genomeFile.rfind('_')]
        # access list of chromosomes and write them to dictionary listing for species...
        genomeFile2=genomeFile+'.fai'
        #faiOpenFile=open(genomePath+genomeFile2,'r')
        genomeNameFind=genomeFile2.split('_')
        fastaObjectStructure[speciesName] = fastaGenome
        # add list of different chromosomes for species
        #for line in faiOpenFile:
        #    lineList=line.split()
        #    fastaObjectStructure[genomeNameFind[0]][1]+=[lineList[0]]
    return fastaObjectStructure

def isbetween(start_coord_query,end_coord_query,start_coord_higherlvl,end_coord_higherlvl):
    """Compare function for if the query start or end coordinate is between the higher level start/end coordinates"""
    if start_coord_higherlvl <= start_coord_query <= end_coord_higherlvl or \
                            start_coord_higherlvl <= end_coord_query <= end_coord_higherlvl or \
                            start_coord_query <= start_coord_higherlvl <= end_coord_query:
        return True
    else:
        return False


def pairComparisonSynteny(syntenicInputFiles,pathUnOut,pathSort):
    """Creates a list of syntenies for a specific pair wise comparison from the .unout files and .sort files.
    Syntenies are of form [(species1chromosome,start_coord,end_coord),(species2chromosome,start_coord,end_coord)]"""

    # NOTE: syntenicInputTuple=(SyntenicFile,species1SortFile,species2SortFile)
    Syntenies=[]
    # now lets try to parse the .unout...
    syntenyFile = open(pathUnOut+syntenicInputFiles[0], 'r')

    # Find species Names (numbers)
    speciesName = [syntenicInputFiles[0][syntenicInputFiles[0].find('.')+1:syntenicInputFiles[0].find('-')],
                   syntenicInputFiles[0][syntenicInputFiles[0].find('PAC2'):syntenicInputFiles[0].find('_5')][7:]]

    # read the line that has the proper header to be parsed
    reading = 1
    """
    lines=[]
    for line in syntenyFile:
        if '[' in line:
            lines+=[line]
    """
    # the following will grab start and end genes for each species for each synteny
    chunkList = []
    newChunk=[0,0,0,0]
    grabLine = ['$', '$','$','$']
    # find lines that match description of syntenic sequences
    for line in syntenyFile:
        for i in range(4):
            if grabLine[i] in line:
                newChunk[i] = line
        if 0 not in newChunk:
            chunkList += [newChunk]
            newChunk=[0,0,0,0]
            grabLine = ['$', '$','$','$']
        if '[' in line:
            lineList = line.split()
            if int(lineList[-6]) >= 4: # if number loci greater than or equal to 4
                grabLine = [lineList[2], lineList[4],lineList[-11],lineList[-9]]
                #if 'Minus' in line:
                #    grabLine[2],grabLine[3] = grabLine[3],grabLine[2]
                #if 'Plus' in line:
                #    tag = ' Plus'
                #elif 'Minus' in line:
                #    tag = ' Minus'
                #else: tag=''

    syntenyFile.seek(0)
    global syntenicSequences
    syntenicSequences = []
    # Bd2,510683,Bradi2g41430.1,Chr1,2751500,LOC_Os01g38960.1
    # list of all of the genes in all syntenies for pair; used to speed up processing time when adding coord for gene
    geneListSpecies = [[],[]] #species 1 and 2 gene lists respecitvely
    # adding start and end genes here to particular synteny (will add coords later)
    for chunk in chunkList:
        #lineList=range(4)
        lineImportantInfo=range(4)
        for i in range(4):
            #lineList[i] = chunk[i].split()
            lineImportantInfo[i] = chunk[i].split()[0].split(',')
        #if 0 in speciesName:
        #    speciesName[0], speciesName[1] = line.split()



        startEndGenes =[[speciesName[0]+'_'+lineImportantInfo[0][0], lineImportantInfo[0][2], lineImportantInfo[1][2]],
                       [speciesName[1]+'_'+lineImportantInfo[0][3], lineImportantInfo[2][5], lineImportantInfo[3][5]]]
        #for i in range(2):
            #if startEndGenes[i][1] > startEndGenes[i][2]: #if gene1 is ordered after gene2, switch genes....
             #   startEndGenes[i][1],startEndGenes[i][2] = startEndGenes[i][2],startEndGenes[i][1]

        # now switch order of genes according to plus or minus rules
        syntenicSequences += [startEndGenes]
        geneListSpecies[0] += [lineImportantInfo[0][2], lineImportantInfo[1][2]] #species 1 specific genes....
        geneListSpecies[1] += [lineImportantInfo[2][5], lineImportantInfo[3][5]] # species 2

    #print lines

    syntenyFile.close()

    def removePlusMinus(speciesChromosomeName):
        # GET RID OF PLUS OR MINUS... FOR NOW keeping it!!!
        if 'Plus' in speciesChromosomeName:
            speciesChromosomeName = speciesChromosomeName[:speciesChromosomeName.find('Plus')-1]
        if 'Minus' in speciesChromosomeName:
            speciesChromosomeName = speciesChromosomeName[:speciesChromosomeName.find('Minus')-1]
        return speciesChromosomeName

    # import the sorted file to find start and end locations of each gene
    sortFiles = [syntenicInputFiles[1],syntenicInputFiles[2]]
    # open species file to find positions in the genome where syntenic sequence is located in between
    sortFileOpen = [open(pathSort+sortFiles[0], 'r'), open(pathSort+sortFiles[1], 'r')]
    """
    #test code
    testOut= open('testout.txt','w')
    for line in lines:
        testOut.write(line+'\n')

    exit()
    """
    # recreate structure here!!!!! START HERE!!!
    # more parsing, generate syntenies from each line

    def checkForGene(geneLine,speciesNumber):
        """Checks particular line of sort file and searches through syntenic sequence list structure and outputs start
        or end coordinate of gene"""
        global syntenicSequences
        geneCoordInfos = []
        for i in range(len(syntenicSequences)):
            for k in [1,2]:
                #if findGenes[i][speciesNumber][k] == 'Pavir.9KG291000.1' and geneLine.split()[0] == 'Pavir.9KG291000.1': #REMOVE
                 #   a=1
                if str(syntenicSequences[i][speciesNumber][k]) == geneLine.split()[0]:
                    syntenicSequences[i][speciesNumber][k] = int(geneLine.split()[k+1])
                """   if type(syntenicSequences[i][speciesNumber][1]) == \
                            type(syntenicSequences[i][speciesNumber][2]) and \
                                    syntenicSequences[i][speciesNumber][2] < \
                                    syntenicSequences[i][speciesNumber][1]:
                        # check if end coordinate less than start coord (find out why??)
                """

                    #geneCoordInfos += [(i, k, geneLine.split()[k+1])]
        #if geneCoordInfos == []:
        #    return [('NaN', 'NaN', 'NaN')]
        #else:
        #    return geneCoordInfos


# species 1 and 2
    for speciesNumber in range(2):
        # check each line in sort file and..
        for line in sortFileOpen[speciesNumber]:
            # see if gene is a gene in a list of all syntenic genes
            if line.split()[0] in geneListSpecies[speciesNumber]:
                # if so, replace gene with its corresponding start/stop coordinate in the genome
                # note: output multiple gene coords!!!!
                checkForGene(line,speciesNumber)
                #geneCoordInfos=checkForGene(line,syntenicSequences,speciesNumber)
                """
                for geneCoord in geneCoordInfos:
                    if geneCoord[0] != 'NaN':
                        syntenicSequences[geneCoord[0]][speciesNumber][geneCoord[1]]=int(geneCoord[2].strip('\n')) # use int for now
                        if type(syntenicSequences[geneCoord[0]][speciesNumber][1]) == \
                            type(syntenicSequences[geneCoord[0]][speciesNumber][2]) and \
                                syntenicSequences[geneCoord[0]][speciesNumber][2] < \
                                        syntenicSequences[geneCoord[0]][speciesNumber][1]:
                            errorFile.write(line.split()[0]+' '+str(syntenicSequences[geneCoord[0]])+('\n')) #check if end coordinate less than start coord (find out why??)
#Seita.4G005900.1	scaffold_4	406794	407432	213	1	U	1509211
        #comment and read into greater detail!!!!!
                """


    #for syntenySequence in syntenicSequences:
        """
        string1list = line[line.find('[') + 1:line.find(']')].split(' ')
        string2list = line[line.rfind('[') + 1:line.rfind(']')].split(' ')

        # grabs the beginning and end genes of the syntenic sequence for each species
        speciesgenes = [[string1list[1], string1list[3]], [string2list[1], string2list[3]]]

        # grabs species+chromosome names
        speciesChr1 = line[0:line.find('[') - 1]
        speciesChr2 = line[line.find(':') + 2:line.rfind('[') - 1]
        speciesChr1 = removePlusMinus(speciesChr1)
        speciesChr2 = removePlusMinus(speciesChr2)
        """
        #PLEASE CHANGE THIS... WILL TAKE TOO LONG TO READ THROUGH FILES

        """
        j = 0
        genePositionList = [[], []]
        for i in range(2):
            # for each line in the sorted file
            for line2 in sortFileOpen[i]:
                lineList = line2.split()
                # determine if there is reverse sequencing in the syntenic sequence
                if int(lineList[2]) >= int(lineList[3]):
                    print 'Reverse sequencing found @:', line2
                # locate the start/end gene within the sorted list of genes, add start and end coordinate info
                if speciesgenes[i][0] == lineList[7] or speciesgenes[i][1] == lineList[7]:
                    if speciesChr1 not in lineList[1] and speciesChr2 not in lineList[1]:
                        print line2, speciesChr1, speciesChr2, speciesgenes #DEBUGGING!!!
                    genePositionList[i] += [int(lineList[2]), int(lineList[3])]
                if len(genePositionList[i])==4:
                    # keep only start coord of first gene and end coord of last gene to capture entire sequence
                    del genePositionList[i][1]
                    del genePositionList[i][1]
                    break
            """
                #if i == 0:
                    ########$^#$%^#$%&#^ AND HERE!!!

            # read again from the top (FOR NOW)
            #sortFileOpen[i].seek(0)

        #START HERE!!!!
        """
        try:
            Syntenies+=[[(speciesChr1,genePositionList[0][0],genePositionList[0][1]),
                        (speciesChr2,genePositionList[1][0],genePositionList[1][1])]]
        except:
            print (genePositionList[0][0],genePositionList[0][1],genePositionList,sortFiles,
                    speciesgenes[1][1])
            print speciesChr1, speciesChr2, speciesgenes
            exit()
        """

    # close sort files...
    for i in range(2):
        sortFileOpen[i].seek(0)
        sortFileOpen[i].close()

    ###########

    # create bedtool files for AB and AC pairwise analysis to test intersections...
    open('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')]),'w').close()
    bedoutfile = open('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')]),'w')
    for syntenicSequence in syntenicSequences:
        if int(syntenicSequence[1][2])<int(syntenicSequence[1][1]): #FIXME
            print syntenicSequence[1]
        bedoutfile.write('%s\t%s\t%s\t%s_%s_%s \n' %tuple(syntenicSequence[0]+syntenicSequence[1]))
    bedoutfile.close()
    bedAnalyze = BedTool('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')])).sort()


    stop=1 #pause here and drag bed file to new location

    ###########
    return bedAnalyze



def syntenicStructure(subgenome): #input N or K for subgenome

    # NOTE: FASTA structures must match syntenies!!! ~~~~~~~
    # generate the fasta data structure
    listOfGenomeFiles=[]
    fastaFindFile = open('syntenicTuplesList.txt','r')
    readFasta = 0 # so far do not input
    for line in fastaFindFile:
        if 'genomePath' in line:
            genomePath = line.split()[-1]
        if readFasta == 1:
            if 'Stop' in line:
                break
            listOfGenomeFiles += [line.strip('\n')]
        if 'listOfGenomeFiles:' in line:
            readFasta = 1
    fastaFindFile.seek(0)
    fastaFindFile.close()
    fastaObjectStructure = generateFastaObjectStructure(listOfGenomeFiles,genomePath)



    # generate all of the syntenies from the following input tuples 'PAC2_0.283-PAC2_0.323_5.unout'
    # generate list of SyntenicTuples from unout files (in correct order for analysis)... will change!
    # do analysis on N subgenome for now...
    listOfSyntenicTuples = []

    #let's create syntenic tuples list from file
    synTupFile = open('syntenicTuplesList.txt','r')
    read = 0 # so far, not reading lines to create tuples
    for line in synTupFile:
        if 'pathUnOut' in line:
            pathUnOut = line.split()[-1]
        if 'pathSort' in line:
            pathSort = line.split()[-1]
        if read == 1:
            if 'Stop' in line:
                break
            syntenicList = line.split('   ')
            syntenicList[2] = syntenicList[2][:-1]
            listOfSyntenicTuples += [tuple(syntenicList)]
        if '%s Test'%subgenome in line: # if performing N or K Test analysis
            read = 1
            pathFastaOutput = line.split()[-1].strip('\n')

    #for file in os.listdir(pathUnOut):

    listOfPairedComparisonSyntenies=[]
    #^^^ add more above! Species A and B first
    # analyze pairwise comparison of syntenies
    # each of these pairwise comparison syntenies are bedtools
    for syntenicInputTuple in listOfSyntenicTuples:
        listOfPairedComparisonSyntenies+=[pairComparisonSynteny(syntenicInputTuple,pathUnOut,pathSort)]

    finalSyntenyStructureBed = listOfPairedComparisonSyntenies[0]

    for pairwiseComparison in listOfPairedComparisonSyntenies[1:]:
        finalSyntenyStructureBed = finalSyntenyStructureBed.cat(pairwiseComparison,postmerge=False)

    finalSyntenyStructureBed = finalSyntenyStructureBed.sort().merge(o='distinct',c=4,delim='|',d=100000)

    print finalSyntenyStructureBed

    count = 1
    open('ERRTEST.txt', 'w').close()
    errorFile = open('ERRTEST.txt', 'w')
    # OUTPUTTING TO FASTA FILES!!!
    for line in str(finalSyntenyStructureBed).split('\n'):
        if line:
            open(pathFastaOutput+'FastaOut%d.fasta'%count,'w').close()
            fastaOutFile = open(pathFastaOutput+'FastaOut%d.fasta'%count,'w')
            syntenicSequenceParseList = line.split('\t')
            speciesAOut = syntenicSequenceParseList[0:3]
            if speciesAOut[0].split('_')[0] in ['523','524']:
                speciesAName = '383'
            else:
                speciesAName = speciesAOut[0].split('_')[0]
            #Target and query A species have tuples that describe syntenic sequence (Species, chromosome, xi, xf)
            speciesAOutTuple = (speciesAName,speciesAOut[0].split('_')[1],speciesAOut[1],speciesAOut[2])
            fastaOutputTuples = [speciesAOutTuple]
            if '|' in syntenicSequenceParseList[3]:
                listOfParseTargetSpecies = syntenicSequenceParseList[3].split('|')
            else:
                listOfParseTargetSpecies = [syntenicSequenceParseList[3]]
            for targetSpecies in listOfParseTargetSpecies:
                outputTargetList = targetSpecies.split('_')
                fastaOutputTuples += [(outputTargetList[0],
                               targetSpecies[targetSpecies.find('_')+1:targetSpecies.find(outputTargetList[-2])-1],
                               outputTargetList[-2], outputTargetList[-1])]
                if targetSpecies[targetSpecies.find('_')+1:targetSpecies.find(outputTargetList[-2])-1]=='scaffold_':
                    a=1
            for fastaOutputTuple in fastaOutputTuples:
                try:
                    if int(fastaOutputTuple[2]) > int(fastaOutputTuple[3]): #FIXME!!!
                        print fastaOutputTuple
                        fastaOutputTuple[2],fastaOutputTuple[3] = fastaOutputTuple[3],fastaOutputTuple[2]
                    fastaOutFile.write('> %s %s %s %s\n' %fastaOutputTuple)
                    fastaOutFile.write(str(fastaObjectStructure[fastaOutputTuple[0]][fastaOutputTuple[1]][int(fastaOutputTuple[2])-1:int(fastaOutputTuple[3])])+'\n')
                except:
                    errorFile.write(str(fastaOutputTuple) + ('\n'))
            fastaOutFile.close()
        count += 1
    errorFile.close()

    """
    # first generate list of syntenic dictionaries from species A&B from the first synteny
    # create initial synteny dictionary structure from species A & B
    finalSyntenyStructure=[]


    def isSyn2dict(syntenyList,eachSyntenyDictionary):
        #Determine whether to convert syntenic sequence represented by the second tuple of a list of tuples
        #into a dictionary object.
        queryName=syntenyList[0][0] #Always Species A

        resultName=syntenyList[1][0] # Species B, C, D, E... etc

        if eachSyntenyDictionary.has_key(queryName):
            for item in eachSyntenyDictionary[queryName]:
                if isbetween(int(syntenyList[0][1]),int(syntenyList[0][2]), int(item[0]),int(item[1])):
                    return True
            return False

        else:
            return False


    # initialize structure for final syntenies, this will combine similar syntenies across pair comparisons
    for eachPairedComparisonSyntenies in listOfPairedComparisonSyntenies[0]:
        # initialize the structure with the paired comparison between species A & B, paired comparison converted to dict
        finalSyntenyStructure+=[{eachPairedComparisonSyntenies[0][0]:[(eachPairedComparisonSyntenies[0][1],
                                                                     eachPairedComparisonSyntenies[0][2])],
                                 eachPairedComparisonSyntenies[1][0]:[(eachPairedComparisonSyntenies[1][1],
                                                                     eachPairedComparisonSyntenies[1][2])]}]

    # let's add the other paired comparisons to the dictionary structure under the premise that we'll combine similar
    # syntenies
    if len(listOfPairedComparisonSyntenies) > 1: # if we are performing analysis on more than 2 species
        for i in range(1,len(listOfPairedComparisonSyntenies)): # already initialized structure using 0th entry
            # goes through all syntenies in particular paired comparison
            for eachSynteny in listOfPairedComparisonSyntenies[i]:
                # finds out if particular synteny [(species1,xi,xf),(species2,yi,yf)] can be added to dictionary for
                # particular synteny type that combines all species into dictionary, synt dict incorporate all species
                if isSyn2dict(eachSynteny,finalSyntenyStructure[i]):
                    # either add onto existing dictionary/synteny type or create new key...
                    if finalSyntenyStructure[i].has_key(eachSynteny[1][0]):
                        finalSyntenyStructure[i][eachSynteny[1][0]].append((eachSynteny[1][1],eachSynteny[1][2]))
                    else:
                            finalSyntenyStructure[i][eachSynteny[1][0]]=[(eachSynteny[1][1], eachSynteny[1][2])]

    # Test item here to see if dictionaries work....????
    for item in finalSyntenyStructure:
        if len(item)>2:
            print item
    a=1
    """

syntenicStructure('K')
syntenicStructure('N')
# maybe input inputSyntenicFilesTuples and genome .fa file list.... not sure

# repair process that looks for start and end coordinates

# output syntenic files to new folder!







#fastaObjectStructure['Bdistachyon'][0]['Bd2'][0:10] indexes like python... Take this into account!!!!
#listOfGenomeFiles=['Bdistachyon_314_v3.0.softmasked.fa','Osativa_323_v7.0.softmasked.fa']
#generateFastaObjectStructure(listOfGenomeFiles)
# note: accessing fasta in form str(speciesAccess['Chr'][0]['Chr1'][1:10]) grabs 10 basepairs of chromosome 1 of Chr
#syntenicInputTuple = ('PAC2_0.283-PAC2_0.323_5.unout','q.PAC2_0.283.sort.txt','t.PAC2_0.323.sort.txt')
#pairComparisonSynteny(syntenicInputTuple)
#^^^ need to generate list of tuples above