from collections import defaultdict
import os, subprocess
from pybedtools import BedTool
def pairComparisonSynteny(unoutFile,Loci_Threshold):
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
    syntenyFile = open(unoutFile, 'r')
    # Find species Names (numbers) #FIXME Change to parse through unout
    speciesName = unoutFile.replace('PAC2_0.','').replace('PAC4GC.','').split('-')
    speciesName = [speciesName[0],speciesName[1][:speciesName[1].find('_')]]

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
            newChunk += [score,orientation]
            chunkList += [newChunk]
            newChunk=[0,0,0,0]
            grabLine = ['$', '$','$','$']
        if '[' in line:
            lineList = line.split()
            if int(lineList[-6]) >= Loci_Threshold: # if number loci greater than or equal to 4
                # parsing to find the start/end genes in the genes listed under each header in unout file
                grabLine = [lineList[2], lineList[4],lineList[-11],lineList[-9]]
            score = lineList[-3]
            if 'Plus' in line:
                orientation = '+'
            elif 'Minus' in line:
                orientation = '-'
            else:
                orientation = '+'
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

    with open('%s-%s.simple'%(speciesName[0],speciesName[1]),'w') as f:
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
            startEndGenes =[ lineImportantInfo[0][2], lineImportantInfo[1][2],
                         lineImportantInfo[2][5], lineImportantInfo[3][5]] + chunk[4:]
            f.write('\t'.join(gene for gene in startEndGenes)+'\n')
            # Removed code that tried to deal with gene ordering
            #for i in range(2):
                #if startEndGenes[i][1] > startEndGenes[i][2]: #if gene1 is ordered after gene2, switch genes....
                 #   startEndGenes[i][1],startEndGenes[i][2] = startEndGenes[i][2],startEndGenes[i][1]

        #print lines

        syntenyFile.close()
        f.close()
    return '%s-%s.simple'%(speciesName[0],speciesName[1])

speciesChromosomes = defaultdict(list)

def reverseChromosomes(bedFile,chromosomes,speciesChromosomes):
    with open(bedFile,'r') as f:
        bedlines = f.readlines()
    bedlines2 = []
    for line in bedlines:
        if line:
            lineList = line.split('\t')
            for chromosome in chromosomes:
                if chromosome in lineList[0]:
                    #print list(map(str,[speciesChromosomes[chromosome]-int(lineList[2]),speciesChromosomes[chromosome]-int(lineList[1])]))
                    #print '\t'.join([lineList[0]] + list(map(str,[speciesChromosomes[chromosome]-int(lineList[2]),speciesChromosomes[chromosome]-int(lineList[1])])) + lineList[3:])+'\n'
                    try:
                        bedlines2.append('\t'.join([lineList[0]] + list(map(str,[speciesChromosomes[chromosome]-int(lineList[2]),speciesChromosomes[chromosome]-int(lineList[1])])) + lineList[3:])+'\n')
                        break
                    except:
                        print line
                if all(chromosome not in line.split('\t')[0] for chromosome in chromosomes):
                    bedlines2.append(line+'\n')

    with open(bedFile+'2','w') as f:
        f.writelines(bedlines2)
    BedTool(bedFile+'2').sort().saveas(bedFile+'2')
    return 0



layoutTxt = defaultdict(list)
simpleFiles = []
for file in os.listdir('.'):
    if file.endswith('.gff3'):
        if file.replace('.gff3','.bed') not in os.listdir('.'):
            subprocess.call('python -m jcvi.formats.gff bed --type=mRNA --key=Name %s -o %s'%(file,file.replace('.gff3','.bed')),shell=True)
    if file.endswith('.fai'):
        speciesname = file.split('_')[1]
        layoutTxt[speciesname] = ','.join([line.split()[0] for line in open(file,'r').readlines()])
        speciesChromosomes[speciesname] = {line.split('\t')[0]:int(line.split('\t')[1]) for line in open(file,'r').readlines() if line}
    if file.endswith('.unout'):
        simpleFiles.append(pairComparisonSynteny(file,4))
seqid = '\n'.join(layoutTxt[species] for species in ['323','316','314','001','002'])
print speciesChromosomes
layout = """# y, xstart, xend, rotation, color, label, va,  bed
.8,     .35,    .65,      0,      , Rice, top, t.PAC2_0.323.bed2
.6,     .35,    .65,      0,      , Bstacei, middle, t.PAC2_0.316.bed2
.4,     .1,    .4,      -15,      , Bdistachyon, top, t.PAC2_0.314.bed
.2,     .1,    .4,      0,      , ABR113D, bottom, q.PAC4GC.001.bed
.3,     .6,    .9,      15,      , ABR113S, bottom, q.PAC4GC.002.bed2
# edges
e, 1, 0, 316-323.simple
e, 1, 2, 316-314.simple
e, 1, 4, 316-002.simple
e, 2, 3, 314-001.simple
#"""
with open('layout','w') as f:
    f.write(layout)
with open('seqid','w') as f:
    f.write("""Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,ChrUn,ChrSy
Chr01,Chr04,Chr02,Chr09,Chr08,Chr07,Chr06,Chr03,Chr05,Chr10
Bd2,Bd3,Bd1,Bd5,Bd4
BhD2,BhD3,BhD1,BhD5,BhD4
BhS1,BhS4,BhS2,BhS9,BhS8,BhS7,BhS6,BhS3,BhS5,BhS10""")
reverseChromosomes('t.PAC2_0.323.bed',['Chr5','Chr7'],speciesChromosomes['323'])
reverseChromosomes('t.PAC2_0.316.bed',['Chr01','Chr04'],speciesChromosomes['316'])
reverseChromosomes('q.PAC4GC.002.bed',['BhS1','BhS4'],speciesChromosomes['002'])


    #f.write(seqid)
subprocess.call('python -m jcvi.graphics.karyotype seqid layout',shell = True)