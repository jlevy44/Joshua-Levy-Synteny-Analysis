import sys, os
import numpy as np


refSamp = sys.argv[1] # r or s

fname = sys.argv[2]+'.cds'

version = sys.argv[3] # v0, v1, ...


if refSamp == 'r':
    folder = './referenceGenomes/'
else:
    folder = './'+version+'/'+sys.argv[2]+'/'

os.chdir(folder)

with open(fname,'r') as f:
    writeLines = np.vectorize(lambda line: line[:line.find('.')].replace('.','').replace('-','').strip('\n') + '\n')(f.readlines())
with open(fname,'w') as f2:
    f2.writelines(writeLines)
