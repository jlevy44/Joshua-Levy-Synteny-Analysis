from pybedtools import BedTool
from collections import Counter
from collections import defaultdict

def pairComparisonSynteny(syntenicInputFiles,pathUnOut,pathSort,Loci_Threshold):
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
    bedAnalyze = BedTool('%s.bed'%(syntenicInputFiles[0][:syntenicInputFiles[0].rfind('.')])).sort()

    # old debug method...
    #stop=1 #pause here and drag bed file to new location

    ###########
    return bedAnalyze # this is our BedTool object

bedSynt = pairComparisonSynteny(('PAC4GC.001-PAC4GC.003_5.unout','q.PAC4GC.001.sort2','q.PAC4GC.003.sort2'),'/Users/jlevy/Desktop/Projects/Synteny/OtherProjects/FixBhydTest/','/Users/jlevy/Desktop/Projects/Synteny/OtherProjects/FixBhydTest/',1)

a=1