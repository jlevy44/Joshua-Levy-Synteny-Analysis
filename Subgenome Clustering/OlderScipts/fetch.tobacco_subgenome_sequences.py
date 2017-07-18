from pybedtools import *
from pyfaidx import Fasta
import numpy as np
import subprocess, os, shutil, sys
import string

__author__ = 'sgordon007'


seqListA = open('subgenome.A.chrom.txt','r')
seqListB = open('subgenome.B.chrom.txt','r')

open('Nitab_v4.5_genome_Chr_S.fasta','w').close()
file_S = open('Nitab_v4.5_genome_Chr_S.fasta','w')

wholeGenome = Fasta('Nitab_v4.5_genome_Chr_Edwards2017.fasta')

for line in seqListA:
    # line.strip()
    scaffold = line.strip().replace('>','')
    print scaffold
    file_S.write('>%s\n%s\n' % (scaffold, str(wholeGenome[scaffold][:])))

seqListA.close()
file_S.close()

open('Nitab_v4.5_genome_Chr_T.fasta','w').close()
file_T = open('Nitab_v4.5_genome_Chr_T.fasta','w')

for line in seqListB:
    # line.strip()
    scaffold = line.strip().replace('>','')
    print scaffold
    file_T.write('>%s\n%s\n' % (scaffold, str(wholeGenome[scaffold][:])))

seqListB.close()
file_T.close()

