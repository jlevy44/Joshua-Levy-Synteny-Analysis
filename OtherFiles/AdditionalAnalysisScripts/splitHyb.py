from pyfaidx import Fasta
import numpy as np

faiFile = open('hybridum_550_v1.fasta.fai','r')

chrDict = {}
chrDictNum = {}

for line in faiFile:
    if line:
        chrDict[line.split()[0]]={'546':[],'547':[]}
        chrDictNum[line.split()[0]] = {'546': 0, '547': 0}
faiFile.close()

gffFiles = ['546.gff3','547.gff3']

for file in gffFiles:
    inputFile = open(file,'r')
    key = file.split('.')[0]

    for line in inputFile:
        if line:
            if 'mRNA' in line and 'longest=1' in line:
                chrDict[line.split()[0]][key].append(int(line.split()[4])-int(line.split()[3]))

    inputFile.close()

open('hybridum_546_v1.fasta','w').close()
open('hybridum_547_v1.fasta','w').close()

file546 = open('hybridum_546_v1.fasta','w')
file547 = open('hybridum_547_v1.fasta','w')

bhybridWhole = Fasta('hybridum_550_v1.fasta')

for key in chrDict.keys():
    for key2 in chrDict[key].keys():
        if not chrDict[key][key2]:
            chrDict[key][key2].append(0)
    if key != 'Sc16uGnS2':
        if sum(chrDict[key]['546']) > sum(chrDict[key]['547']):
            file546.write('>%s\n%s\n'%(key,str(bhybridWhole[key][:])))
        else:
            file547.write('>%s\n%s\n' % (key, str(bhybridWhole[key][:])))

cutoff = 15481775
key = 'Sc16uGnS2'
file546.write('>%sd\n%s\n'%(key,str(bhybridWhole[key][:cutoff])))
file547.write('>%ss\n%s\n' % (key, str(bhybridWhole[key][cutoff:])))

file546.close()
file547.close()

Fasta('hybridum_546_v1.fasta')
Fasta('hybridum_547_v1.fasta')