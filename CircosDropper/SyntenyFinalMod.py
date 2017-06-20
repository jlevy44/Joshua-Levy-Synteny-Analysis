import os, sys
from pyfaidx import Fasta
from pybedtools import *
from collections import defaultdict

"""Creates Fasta object structure of genomes and multiple bedtools files/objects for each pairwise comparison
(another structure created for this. The bedtool object structure of pairwise comparison is merged to stitch together
all of these comparisons to create syntenies between more than two species. The resulting BedTool object is parsed and
referenced to access the Fasta structure to pull information/sequences from the genomes themselves and export to
fasta files to be used by a MSA such as Cactus."""

def generateFastaObjectStructure(listOfGenomeFiles,genomePath):
    """Generates dictionary of different species that contains a fasta object (genome) for that species
    Inputs list of genome files in .fa format and the path to find them. Outputs fasta object structure, dictionary of
    Fasta objects that will be referenced in final analysis"""
    fastaObjectStructure={}

    for genomeFile in listOfGenomeFiles:
        #generate a fasta object for given genome
        fastaGenome=Fasta(genomePath+genomeFile)
        # find species name (species number)
        #speciesName = genomeFile[genomeFile.find('_') + 1:genomeFile.rfind('_')]
        speciesName = genomeFile.split('_')[1]
        try:
            int(speciesName)
        except:
            speciesName = genomeFile.split('.')[1] # if in PAC2_0.275.fa format
        # access add fasta genome object to fasta object structure, searchable under species name
        fastaObjectStructure[speciesName] = fastaGenome

        # removed code
        #genomeFile2=genomeFile+'.fai'
        #faiOpenFile=open(genomePath+genomeFile2,'r')
        #genomeNameFind=genomeFile2.split('_')
        # add list of different chromosomes for species
        #for line in faiOpenFile:
        #    lineList=line.split()
        #    fastaObjectStructure[genomeNameFind[0]][1]+=[lineList[0]]
    return fastaObjectStructure

# irrelevant code that was originally used to try to stitch together syntenic sequences
def isbetween(start_coord_query,end_coord_query,start_coord_higherlvl,end_coord_higherlvl):
    """Compare function for if the query start or end coordinate is between the higher level start/end coordinates"""
    if start_coord_higherlvl <= start_coord_query <= end_coord_higherlvl or \
                            start_coord_higherlvl <= end_coord_query <= end_coord_higherlvl or \
                            start_coord_query <= start_coord_higherlvl <= end_coord_query:
        return True
    else:
        return False


def unout2bed(syntenicInputFiles,pathUnOut='./',pathSort='./',Loci_Threshold=4):
    """Creates a list of syntenies for a specific pair wise comparison from the .unout files and .sort2 files and turns
    this list into a bedfile.
    Inputs paired comparison files in the form of a syntenic tuple (.unout file, species 1 sort2 file, species2 sort2)
    and the unout and sort path locations.
    Outputs a BedTool object created from the paired comparison that lists all of the syntenies for the paired
    comparison with formatting for each line: Chr is chromosome
    'Species1Name_Species1Chr startCoordSyntenicSeq endCoordSynSeq species2Name_Chr_startCoord_endCoord'
    Bed files that are converted into BedTool objects were created from
    syntenies of form [[species1chromosome,start_coord,end_coord],[species2chromosome,start_coord,end_coord]]"""

    # NOTE: syntenicInputTuple=(SyntenicFileUnOut,species1SortFile,species2SortFile)
    Syntenies=[]
    # now lets try to parse the .unout...
    syntenyFile = open(pathUnOut+syntenicInputFiles[0], 'r')

    # Find species Names (numbers) #FIXME
    speciesName = [syntenicInputFiles[1].split('.')[-2],syntenicInputFiles[2].split('.')[-2]]

    # unused code below
    # read the line that has the proper header to be parsed
    #reading = 1
    """
    lines=[]
    for line in syntenyFile:
        if '[' in line:
            lines+=[line]
    """

    # the following will grab start and end genes for each species for each syntenic block in the unout
    # Unout file is the syntenyFile object listed below
    chunkList = []
    newChunk=[0,0,0,0]
    lineCheck = {0:1,1:1,2:4,3:4}
    grabLine = ['$', '$','$','$']
    # find lines that match description of syntenic sequences in the unout file
    # a Chunk is a list of lines read from the the unout file
    # newChunk contains list of [line containing species 1 start gene, line containing species 1 end gene, line
    #containing species 2 start gene, line contain species 2 end ]
    # first it parses the header line containing a '[' to find which lines to grab...
    for line in syntenyFile:
        for i in range(4):
            if ',' in line.split()[0] and grabLine[i] == line.split()[0].split(',')[lineCheck[i]]:
                newChunk[i] = line
        if 0 not in newChunk:
            chunkList += [newChunk]
            newChunk=[0,0,0,0]
            grabLine = ['$', '$','$','$']
        if '[' in line:
            lineList = line.split()
            if int(lineList[-6]) >= Loci_Threshold: # if number loci greater than or equal to 4
                # parsing to find the start/end genes in the genes listed under each header in unout file
                grabLine = [lineList[2], lineList[4],lineList[-11],lineList[-9]]

                # removed code that attempted to deal with gene ordering
                #if 'Minus' in line:
                #    grabLine[2],grabLine[3] = grabLine[3],grabLine[2]
                #if 'Plus' in line:
                #    tag = ' Plus'
                #elif 'Minus' in line:
                #    tag = ' Minus'
                #else: tag=''

    syntenyFile.seek(0)


    global syntenicSequences # creating syntenic sequences as aformentioned
    # [[species1chromosome,start_gene,end_gene],[species2chromosome,start_gene,end_gene]]
    # the genes will later be turned into their corresponding coordinates as the mutable list of syntenic sequence
    # is accessed through the sort2 file

    syntenicSequences = []

    # example line from a Chunk
    # Bd2,510683,Bradi2g41430.1,Chr1,2751500,LOC_Os01g38960.1   other stuff over here

    # list of all of the genes in all syntenies for pair; used to speed up processing time when adding coord for gene
    geneListSpecies = [[],[]] #species 1 and 2 gene lists respecitvely
    # the above is just a list of genes found in the unout files so when we open the sort2 files, we can know what we
    # are looking for and can speed up processing

    #FIXME
    def endsWithM(gene):
        if gene.endswith('m'):
            return gene[:-1]

    # adding start and end genes here to particular synteny (will add coords later)
    for chunk in chunkList:

        # Bd2,510683,Bradi2g41430.1,Chr1,2751500,LOC_Os01g38960.1
        # each lineImportant info will contain above info, gets rid of other stuff from a line in chunk and splits
        # above info into a list
        lineImportantInfo=range(4)
        for i in range(4):
            lineImportantInfo[i] = chunk[i].split()[0].split(',')

        # unused code
        #if 0 in speciesName:
        #    speciesName[0], speciesName[1] = line.split()

        # this is a syntenic sequence here, derived from parsing through the line of important info and retaining
        # species names
        #if 'Brsylv' in chunk[0] and 'Brsylv' in chunk[1]: # fix output for bad naming for this new species FIXME future
        #    for i in range(2):
        #        lineImportantInfo[i][2] = lineImportantInfo[i][2][:-1]
        startEndGenes =[[speciesName[0]+'-'+lineImportantInfo[0][0].replace('-','~'), lineImportantInfo[0][2], lineImportantInfo[1][2]],
                    [speciesName[1]+'-'+lineImportantInfo[0][3].replace('-','~'), lineImportantInfo[2][5], lineImportantInfo[3][5]]]

        # Removed code that tried to deal with gene ordering
        #for i in range(2):
            #if startEndGenes[i][1] > startEndGenes[i][2]: #if gene1 is ordered after gene2, switch genes....
             #   startEndGenes[i][1],startEndGenes[i][2] = startEndGenes[i][2],startEndGenes[i][1]

        # added a syntenic sequence to the list of syntenic sequences, keeping gene names/IDs for now but will find
        # coords for the gene IDs later
        syntenicSequences += [startEndGenes]
        # adding to list of genes for a species to speed up processing later on when finding start/end coords for gene
        geneListSpecies[0] += [lineImportantInfo[0][2], lineImportantInfo[1][2]] #species 1 specific genes....
        geneListSpecies[1] += [lineImportantInfo[2][5], lineImportantInfo[3][5]] # species 2

    #print lines

    syntenyFile.close()

    # import the sorted file to find start and end locations/coordinates of each gene in syntenic sequence
    sortFiles = [syntenicInputFiles[1],syntenicInputFiles[2]]
    sortDictionary = dict.fromkeys(sortFiles, defaultdict(list))
    for file in sortFiles:
        openFile = open(pathSort+file,'r')
        for line in openFile.read().split('\n'):
            if line:
                sortDictionary[file][line.split()[0]] = line.split()[1:]
        openFile.seek(0)
        openFile.close()

    def checkForGene(geneLine,speciesNumber):
        """Checks particular line of sort file and searches through syntenic sequence list structure and outputs start
        or end coordinate of gene. Will turn:
        [[species1chromosome,start_gene,end_gene],[species2chromosome,start_gene,end_gene]]
        into: (start_coord is start coordinate of the starting gene/ first gene of syntenic sequence)
        [[species1chromosome,start_coord,end_coord],[species2chromosome,start_coord,end_coord]]

        Inputs a specified line from the sort2 and speciesNumber (first/second species represented by [0] or [1] index
        and replaces in the global variable syntenicSequences the start/end gene with its corresponding start/end coord
        in the sort2 file. Searches the syntenic sequences structure to find start/end coord referenced in the sort2."""
        global syntenicSequences

        # searching each syntenic sequence (ith) in the structure, only for speciesNumber (species 1 or 2) specified
        # searches 'k', the species start (1) or end (2) gene and replaces it with start/end coordinate
        for i in range(len(syntenicSequences)):
            for k in [1,2]:
                # unused code testing purposes
                #if findGenes[i][speciesNumber][k] == 'Pavir.9KG291000.1' and geneLine.split()[0] == 'Pavir.9KG291000.1': #REMOVE
                 #   a=1

                # if the gene found matches gene in sort2 file
                if str(syntenicSequences[i][speciesNumber][k]) == geneLine.split()[0]:
                    # replace start/end gene with corresponding start/end coordinate
                    syntenicSequences[i][speciesNumber][k] = int(geneLine.split()[k+1])

    for j in range(len(syntenicSequences)):
        for i in range(2):
            for k in range(2):
                syntenicSequences[j][i][k+1]  = sortDictionary[sortFiles[i]][syntenicSequences[j][i][k+1]][k+1]
                a=1


    ###########
    open(pathUnOut[:-1][:pathUnOut[:-1].rfind('/')+1]+'ErrorFiles/ErrFile%s.txt'%speciesName[1],'w').close()
    errfile = open(pathUnOut[:-1][:pathUnOut[:-1].rfind('/') + 1] + 'ErrorFiles/ErrFile%s.txt' % speciesName[1], 'w')
    # create bedtool file from the syntenicSequences structure
    open('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')]),'w').close()
    bedoutfile = open('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')]),'w')
    for syntenicSequence in syntenicSequences:
        # once again, code commented out below was used to check gene/gene coordinate ordering and make sure start coord
        # is less than end coord for syntenic sequence
        #if int(syntenicSequence[1][2])<int(syntenicSequence[1][1]):
        #    print syntenicSequence[1]
        try:
            anyWrong = 0
            for i in range(2):
                syntenicSequence[i][1] = str(int(syntenicSequence[i][1])-1)
                syntenicSequence[i][2] = str(syntenicSequence[i][2])
                if int(syntenicSequence[i][1]) > int(syntenicSequence[i][2]):
                    errfile.write(str(syntenicSequence[i]) + '\n')
                    anyWrong += 1

            # write each sequence to bed file, A/query species will be accessed by BedTools, B/C/D/etc target species is
            # added on as a feature to be preserved when the bedtool objects are merged
            if anyWrong == 0:
                bedoutfile.write('%s\t%s\t%s\t%s-%s-%s \n' %tuple(syntenicSequence[0]+syntenicSequence[1]))
        except:
            errfile.write(str(tuple(syntenicSequence[0]+syntenicSequence[1]))+'\n')
    bedoutfile.close()
    errfile.close()
    # convert the pairwise comparison bed file into BedTool object to be manipulated when BedTool obj structure created
    bedAnalyze = BedTool('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')])).sort().saveas('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')]))

    # old debug method...
    #stop=1 #pause here and drag bed file to new location

    ###########
    return '%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')])



def syntenicStructure(): #input N or K for subgenome, will generalize later
    """Input 'N' or 'K' for now in oder to perform a synteny analysis on either the N or K subgenome. This will compare
    the subgenome to its closest species relatives and builds a Fasta object structure and a BedTool Object structure to
    access the Fasta structure, which accesses the genomes themselves.
    Fasta files are outputted from each analysis that contain information about each merged syntenic sequence and the
    actual basepairs written to the fasta files from the genome .fa files. Synthesizes multiple pairwise comparisons
    between species and creates syntenies between more than 2 species..."""


    # NOTE: FASTA structures must match syntenies!!! ~~~~~~~
    # generate the fasta data structure
    listOfGenomeFiles=[]
    # open config file to grab list of genomes to access
    fastaFindFile = open('syntenicConfig.txt','r') # keep config in running directory
    readFasta = 0 # so far do not input
    for line in fastaFindFile:
        if line:
            # Reads the type of Analysis to be performed; modify config file. this name must match name used for output path
            if 'NameAnalysis' in line:
                nameAnalysis = line.split()[-1].strip('\n')
            if 'writeFastaOut' in line:
                writeFastaOut = line.split()[-1].strip('\n')
            if 'pathPythonModules' in line:
                sys.path.append(line.split()[-1].strip('\n'))
            if 'Loci' in line: # establish loci threshold
                Loci_Threshold = int(line.split()[-1].strip('\n'))
            # in order to run on the server, list path of installed modules^^^
            if 'genomePath' in line:
                # grabs path of genome files locations
                genomePath = line.split()[-1]
            if readFasta == 1:
                # if reading lines in config, pull the genome file names out to be used in setting up fasta structure
                # will generate .fai files in the genome path if not already created...
                if 'Stop' in line:
                    break # stop reading file
                listOfGenomeFiles += [line.strip('\n')]
            if 'listOfGenomeFiles:' in line:
                # begin reading lines to pull genome file names
                readFasta = 1
    fastaFindFile.seek(0) # reset where reading config file to beginning
    fastaFindFile.close()

    try:
        int(writeFastaOut)
    except:
        writeFastaOut = 1

    # generate the Fasta Object structure, a dictionary of Fasta objects referenced by species names
    fastaObjectStructure = generateFastaObjectStructure(listOfGenomeFiles,genomePath)



    # generate all of the syntenies from the following input tuples (unout file, species1 sort2, species2 sort2)
    # generate list of SyntenicTuples from config file 'SyntenicInputTuples.txt'
    # do analysis on N subgenome for now...
    listOfSyntenicTuples = []

    #let's create syntenic tuples list from file, read config file
    synTupFile = open('syntenicConfig.txt','r')
    read = 0 # so far, not reading lines to create tuples
    for line in synTupFile:
        if line:
            if 'pathUnOut' in line: # pull path of unout synteny comparison file
                pathUnOut = line.split()[-1]
            if 'pathSort' in line: # pull path of sort2 files
                pathSort = line.split()[-1]
            if 'BPsMergeDist' in line:
                BPsMergeDist = line.split()[-1]
                print line
            if 'softMasked' in line:
                softmask = line.split()[-1]
            if read == 1:
                if 'Stop' in line:
                    break # read lines until stop
                # PAC4GC.524-PAC2_0.383_5.unout   q.PAC4GC.524.sort2   t.PAC2_0.383.sort2
                syntenicList = line.split('   ') # not tab delimited, may change, but read syntenic tuples in this manner
                syntenicList[2] = syntenicList[2].strip('\n') # remove end line for last entry
                listOfSyntenicTuples += [tuple(syntenicList)]
            if '%s Test'%nameAnalysis in line: # if performing N or K Test analysis, find in line and read lines until stop
                read = 1
                pathFastaOutput = line.split()[-1].strip('\n') # add output path for fasta files

    # if a loci threshold does not exist, set the loci threshold = 4
    try:
        Loci_Threshold
    except:
        Loci_Threshold = 4

    try:
        BPsMergeDist = int(BPsMergeDist)
    except:
        BPsMergeDist = 100000
    try:
        softmask = int(softmask)
    except:
        softmask = 0
    print BPsMergeDist
    synTupFile.close()

    #for file in os.listdir(pathUnOut):
    listOfPairedComparisonSyntenies=[] # this is out BEDTool object structure which contains bed objects that are pair-
    #wise comparisons
    # each of these pairwise comparison syntenies are bedtools
    # generate structure (list of bedtool objects) by going through each paired syntenic comparison
    for syntenicInputTuple in listOfSyntenicTuples:
        listOfPairedComparisonSyntenies+=[pairComparisonSynteny(syntenicInputTuple,pathUnOut,pathSort,Loci_Threshold)]
    if int(writeFastaOut):
        # initialize the final synteny structure (just a bed object that will be concatenated and merged with other beds)
        # first paired comparison bed is initial bed object
        finalSyntenyStructureBed = listOfPairedComparisonSyntenies[0]

        # concatenate remaining paired comparison bed objects together
        for pairwiseComparison in listOfPairedComparisonSyntenies[1:]:
            finalSyntenyStructureBed = finalSyntenyStructureBed.cat(pairwiseComparison,postmerge=False)

        # merge these concatenated bed objects and format, query species do not have to overlap, can be 100000 BPs apart...
        # one line should look like:
        # 523-Chr09N	118260675	118632318	308-Chr07-33574614-33759677|308-Chr09-1325484-1436481| ...etc
        # SpeciesA-Chrom xi xf     SpeciesB-Chr-xi-xf|SpeciesC-Chr-xi-xf|SpeciesC#2-Chr-xi-xf...etc
        finalSyntenyStructureBed = finalSyntenyStructureBed.sort().merge(o='distinct',c=4,delim='|',d=BPsMergeDist)

        #output bed object to file
        open(pathFastaOutput+'finalSyntenyMultipleSpecies.bed','w').close()
        finalSyntenyBedFile = open(pathFastaOutput+'finalSyntenyMultipleSpecies.bed','w')
        finalSyntenyBedFile.write(str(finalSyntenyStructureBed))
        finalSyntenyBedFile.close()

        # hide for now
        #print finalSyntenyStructureBed

        # there exists an error file to make sure program runs to end, output errors onto this file, and search the sort2
        # files for faulty genes
        # error file will contain species-chr, xi, xf for all bad outputs, can look them up in sort and match to unout
        count = 1
        open('ERRTEST.txt', 'w').close()
        errorFile = open('ERRTEST.txt', 'w')
        # OUTPUTTING TO FASTA FILES!!!
        # for each merged syntenic sequence (according to final structure bed object), each line corresponds to one sequence

        def softmaskedOutput(synSeqStr): # format softmasked assembly by turning lower case letters to 'N' for Cactus align
            for ch in ['c','t','a','g']:
                if ch in synSeqStr:
                    synSeqStr=synSeqStr.replace(ch,'N')
            return synSeqStr

        for line in str(finalSyntenyStructureBed).split('\n'):
            # if line exists
            if line:
                open(pathFastaOutput+'FastaOut%d.fasta'%count,'w').close() # create new fasta output file

                # open for writing
                fastaOutFile = open(pathFastaOutput+'FastaOut%d.fasta'%count,'w')
                # parse sequence line bed object eg.
                # 523-Chr09N	118260675	118632318	308-Chr07-33574614-33759677|308-Chr09-1325484-1436481|
                # split into [523-Chr09N, 118260675, 118632318,308-Chr07-33574614-33759677|308-Chr09-1325484-1436481|..]
                syntenicSequenceParseList = line.split('\t')
                #speciesAOut would look like [523-Chr09N, 118260675, 118632318]
                speciesAOut = syntenicSequenceParseList[0:3]
                # if the speciesName (523 in above example) is either 523 or 524 species Name, change name to 383 for naming
                if speciesAOut[0].split('-')[0] in ['523','524']:
                    speciesAName = '383' # P.virgatum if query species
                else:
                    speciesAName = speciesAOut[0].split('-')[0] # hold old naming convention
                #Target and query A species have tuples that describe syntenic sequence (Species, chromosome, xi, xf)
                speciesAOutTuple =(speciesAName,speciesAOut[0].split('-')[1].replace('~','-'),speciesAOut[1],speciesAOut[2])
                # initialize final output structure for fasta file, which uses list of output tuples
                fastaOutputTuples = [speciesAOutTuple]
                # split up the last entry of the syntenicSequenceParseList by the delimiter, these create outputs for
                # bed features/ more distant species relatives (target species)
                if '|' in syntenicSequenceParseList[3]:
                    listOfParseTargetSpecies = syntenicSequenceParseList[3].split('|')
                else:
                    listOfParseTargetSpecies = [syntenicSequenceParseList[3]]
                # for each of the target species (not species A) found from splitting delimiter into list, create similar
                # output tuples as those used for species A output and add all of them to the fastaOutputTuples
                for targetSpecies in listOfParseTargetSpecies:
                    outputTargetList = targetSpecies.split('-')
                    if outputTargetList[0] in ['523','524']: #SPECIFIC TO N OR K ANALYSIS... MAY NEED TO CHANGE
                        speciesBName = '383'
                    else:
                        speciesBName = outputTargetList[0]
                    fastaOutputTuples += [(speciesBName, targetSpecies[targetSpecies.find('-')+1:
                                            targetSpecies.find(outputTargetList[-2])-1].replace('~','-'),
                                            outputTargetList[-2], outputTargetList[-1])]
                    # removed code for testing
                    #if targetSpecies[targetSpecies.find('-')+1:targetSpecies.find(outputTargetList[-2])-1]=='scaffold_':
                    #    a=1
                # for each output tuple for one Synteny across the merged comparison final bed structure, output info/
                # actual DNA sequence to fasta under appropriate header
                if softmask:
                    for fastaOutputTuple in fastaOutputTuples:
                        try:
                            # eliminate reordering of gene coordinates for now... unless drawing error file and need debug
                            #if int(fastaOutputTuple[2]) > int(fastaOutputTuple[3]): #FIXME!!!
                            #    print fastaOutputTuple
                            #    fastaOutputTuple[2],fastaOutputTuple[3] = fastaOutputTuple[3],fastaOutputTuple[2]
                            # header to outputted sequence > Species Chromosome start coord (xi) end coord xf
                            fastaOutFile.write('> %s %s %s %s\n' %fastaOutputTuple)
                            # writing the actual sequence
                            fastaOutFile.write(softmaskedOutput(str(fastaObjectStructure[fastaOutputTuple[0]]
                                                [fastaOutputTuple[1]][int(fastaOutputTuple[2]):int(fastaOutputTuple[3])]))+'\n')
                        except: # if not able to write
                            errorFile.write(str(fastaOutputTuple) + ('\n')) # search through sort2 file to debug error
                else:
                    for fastaOutputTuple in fastaOutputTuples:
                        try:
                            # eliminate reordering of gene coordinates for now... unless drawing error file and need debug
                            #if int(fastaOutputTuple[2]) > int(fastaOutputTuple[3]): #FIXME!!!
                            #    print fastaOutputTuple
                            #    fastaOutputTuple[2],fastaOutputTuple[3] = fastaOutputTuple[3],fastaOutputTuple[2]
                            # header to outputted sequence > Species Chromosome start coord (xi) end coord xf
                            fastaOutFile.write('> %s %s %s %s\n' %fastaOutputTuple)
                            # writing the actual sequence
                            fastaOutFile.write(str(fastaObjectStructure[fastaOutputTuple[0]]
                                                [fastaOutputTuple[1]][int(fastaOutputTuple[2]):int(fastaOutputTuple[3])])+'\n')
                        except: # if not able to write
                            errorFile.write(str(fastaOutputTuple) + ('\n')) # search through sort2 file to debug error
                fastaOutFile.close()
            count += 1 # change cound to change title of fasta out file
        errorFile.close()


    # old code, removed, does not work!! prior attempt without bedtools
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

# run K and N subgenome analysis
#syntenicStructure()








#MISC notes/concerns that may or may not have been addressed
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