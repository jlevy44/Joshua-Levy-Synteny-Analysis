import numpy as np
import multiprocessing as mp
import subprocess
import sys
from pyfaidx import Fasta
from pybedtools import BedTool
fastaFile = Fasta(sys.argv[1])

global inputStr
global key

def grabLine(positions):
    global inputStr
    global key
    if positions[-1] != 'end':
        #return inputStr[positions[0]:positions[1]]
        return '%s\t%s\t%s\t%s\n'%tuple([key]+map(str,positions[0:2])+['_'.join([key]+map(str,positions[0:2]))])
    else:
        return '%s\t%s\t%s\t%s\n' % (key, str(positions[0]), str(len(inputStr)), '_'.join([key, str(positions[0]), str(len(inputStr))]))
        #print inputStr[positions[0]:]
        #return inputStr[positions[0]:]
def split(inputStr='',width=60):
    if inputStr:
        positions = np.arange(0,len(inputStr),width)
        posFinal = []
        if len(inputStr) > width:
            for i in range(len(positions[:-1])):
                posFinal += [(positions[i],positions[i+1])]
        posFinal += [(positions[-1],'end')]
        try:
            print posFinal[0:4]
            print key
        except:
            print 'length scaffold very small'
        if __name__ == '__main__':
            p=mp.Pool(processes=8)
            splitLines = p.map(grabLine,posFinal)
            p.close()
            p.join()
            #print wrappedLines[0:10]
        return splitLines
    else:
        return ''
bedText = []
for key in fastaFile.keys():
    inputStr = fastaFile[key][:].seq
    bedText += split(inputStr,75000)
with open('correspondence.bed','w') as f:
    f.write('\n'.join(bedText))
BedTool('correspondence.bed').sort().saveas('correspondence.bed')
subprocess.call('bedtools getfasta -fi %s -fo %s -bed %s -name'%(sys.argv[1],sys.argv[1].split('.fa')[0]+'_split.fa','correspondence.bed'),shell=True)
Fasta(sys.argv[1].split('.fa')[0]+'_split.fa')