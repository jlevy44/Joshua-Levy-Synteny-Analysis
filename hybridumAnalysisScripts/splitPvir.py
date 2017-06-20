from pyfaidx import Fasta
import numpy as np

faiFile = open('Pvirgatum_450_v4.0.softmasked.fa.fai','r')

PvirK=[]
PvirN=[]

for line in faiFile:
    if line:
        if line.split('\t')[0].endswith('K'):
            PvirK.append(line.split('\t')[0])
        elif line.split('\t')[0].endswith('N'):
            PvirN.append(line.split('\t')[0])
faiFile.close()

gffFiles = {'N':'t.PAC2_0.451.gff3','K':'t.PAC2_0.452.gff3'}

with open('Pvirgatum_450_v4.1.gene.gff3','r') as f:
    GFFlines = f.readlines()

with open(gffFiles['N'],'w') as fN:
    with open(gffFiles['K'],'w') as fK:
        for line in GFFlines:
            if line:
                if line.split('\t')[0].endswith('K'):
                    fK.write(line)
                if line.split('\t')[0].endswith('N'):
                    fN.write(line )


PvirG = Fasta('Pvirgatum_450_v4.0.softmasked.fa')

with open('PvirgatumN_451_v1.fa','w') as f:
    for key in PvirN:
        f.write('>%s\n%s\n' % (key, PvirG[key][:].seq))
with open('PvirgatumK_452_v1.fa','w') as f:
    for key in PvirK:
        f.write('>%s\n%s\n' % (key, PvirG[key][:].seq))


Fasta('PvirgatumN_451_v1.fa')
Fasta('PvirgatumK_452_v1.fa')