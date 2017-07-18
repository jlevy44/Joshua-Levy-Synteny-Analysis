from pybedtools import *
from pyfaidx import Fasta
import numpy as np
import subprocess, os, shutil, sys
import string

__author__ = 'sgordon007'

seqListA = open('subgenome.A.chrom.txt','r')
seqListB = open('subgenome.B.chrom.txt','r')

open('Pvirgatum_450_K.fasta','w').close()
file_K = open('Pvirgatum_450_K.fasta','w')

wholeGenome = Fasta('Pvirgatum_450_v4.0.fa')

for line in seqListA:
    # line.strip()
    scaffold = line.strip().replace('>','')
    print scaffold
    file_K.write('>%s\n%s\n' % (scaffold, str(wholeGenome[scaffold][:])))

seqListA.close()
file_K.close()

open('Pvirgatum_450_N.fasta','w').close()
file_N = open('Pvirgatum_450_N.fasta','w')

for line in seqListB:
    # line.strip()
    scaffold = line.strip().replace('>','')
    print scaffold
    file_N.write('>%s\n%s\n' % (scaffold, str(wholeGenome[scaffold][:])))

seqListB.close()
file_N.close()



