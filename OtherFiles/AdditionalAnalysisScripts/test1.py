import sys
import commands
from pybedtools import *
from Bio import SeqIO
from pyfaidx import Fasta

# lets try g

# test with bedtools, can try this with pybedtools
#commands.getoutput("bedtools getfasta -fi test.fasta -bed test.gff -fo testout.txt")
#commands.getoutput("cat")

# test with biopython; will learn libraries more
from Bio.Seq import Seq

#create a sequence object
my_seq = Seq('CATGTAGACTAG')

#print out some details about it
print 'seq %s is %i bases long' % (my_seq, len(my_seq))
print 'reverse complement is %s' % my_seq.reverse_complement()
print 'protein translation is %s' % my_seq.translate()

# now lets try to parse the test.unout...

testFile=open('test.unout','r')


#read the line that has the proper header to be parsed
reading=1
while reading==1:
    line=testFile.readline()
    if '[' in line:
        reading=0

testFile.close()

print line

# more parsing
string1list=line[line.find('[')+1:line.find(']')].split(' ')
string2list=line[line.rfind('[')+1:line.rfind(']')].split(' ')

# grabs the beginning and end genes of the syntenic sequence for each species
speciesgenes=[[string1list[1],string1list[3]],[string2list[1],string2list[3]]]


# grabs species names
species1=line[0:line.find('[')-1]
species2=line[line.find(':')+2:line.rfind('[')-1]
# GET RID OF PLUS OR MINUS... FOR NOW!!!
species2=species2[:species2.find('Plus')]
species2=species2[:species2.find('Minus')]


print 'Syntenic sequence 1 is located between the following genes:', species1, speciesgenes[0],species2,speciesgenes[1]

# import the sorted file to find start and end locations of each gene
sortFiles=['q.PAC2_0.283.sort.txt','t.PAC2_0.323.sort.txt']

sortFileOpen=[0,0]
j=0
genePositionList=[[],[]]
for i in range(2):
    # open species file to find positions in the genome where syntenic sequence is located in between
    sortFileOpen[i]=open(sortFiles[i],'r')
    # for each line in the sorted file
    for line in sortFileOpen[i]:
        lineList=line.split()
        # determine if there is reverse sequencing in the syntenic sequence
        if int(lineList[2])>=int(lineList[3]):
            print 'Reverse sequencing found @:', line
        # locate the start/end gene within the sorted list of genes, add start and end coordinate info
        if speciesgenes[i][0] == lineList[7] or speciesgenes[i][1] == lineList[7]:
            genePositionList[i]+=[int(lineList[2]),int(lineList[3])]
    # keep only start coord of first gene and end coord of last gene to capture entire sequence
    del genePositionList[i][1]
    del genePositionList[i][2]
    sortFileOpen[i].close()

print genePositionList

# now let's use bedtools to verify above data (?) and fetch these syntenic sequences and output to fasta
# create a bedtool
#bedpool1=BedTool('q.PAC2_0.283.sort.txt')

#print bedpool1

# let's try to read the .fa file into small chuncks using biopython
#record = SeqIO.read("Bdistachyon_314_v3.0.softmasked.fa", "fasta")
genome1=Fasta('Bdistachyon_314_v3.0.softmasked.fa')
genome2=Fasta('Osativa_323_v7.0.softmasked.fa')

genome1Final=str(genome1[species1][genePositionList[0][0]:genePositionList[0][1]])
genome2Final=str(genome2[species2][genePositionList[1][0]:genePositionList[1][1]])

testOutputFile=open('testOutJosh.fasta','w')

#format the species and chromosome name a little bit!..... CHANGE VARIABLE NAMES AND GENERALIZE!!!
testOutputFile.write('\n\n>'+species1[:-1]+'_'+species1[-1]+' '+str(genePositionList[0][0])+' '+str(genePositionList[0][1])
                    +'\n'+genome1Final+'\n\n')
testOutputFile.write('>'+species2[:-1]+'_'+species2[-1]+' '+str(genePositionList[1][0])+' '+str(genePositionList[1][1])
                    +'\n'+genome2Final+'\n\n')
testOutputFile.close()