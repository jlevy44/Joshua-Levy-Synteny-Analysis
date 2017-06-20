from pybedtools import BedTool
from pyfaidx import Fasta
import pandas as pd
import numpy as np
from textwrap import fill
import multiprocessing as mp
from collections import defaultdict
import matplotlib.pyplot as plt
import subprocess
from pybedtools import BedTool
import time
mays = pd.read_csv('maize_subgenome_filled.csv')
"""Make file from above
M1_1  chr xi xf
M2_1

etc

"""

"""create gff from xi data and split into two gff"""

"""extract fastas"""

def grabLine(positions):
    global inputStr
    if positions[-1] != 'end':
        return inputStr[positions[0]:positions[1]]
    else:
        print inputStr[positions[0]:]
        return inputStr[positions[0]:]
def wraplines(width=60):
    positions = np.arange(0,len(inputStr),width)
    posFinal = []
    for i in range(len(positions[:-1])):
        posFinal += [(positions[i],positions[i+1])]
    posFinal += [(positions[-1],'end')]
    print posFinal[0:10]
    if __name__ == '__main__':
        p=mp.Pool(processes=8)
        wrappedLines = p.map(grabLine,posFinal)
        p.close()
        p.join()
        #print wrappedLines[0:10]
    return '\n'.join(wrappedLines)

def df2bed(df):
    npDf = np.array(df)
    return BedTool('\n'.join('%s\t%d\t%d'%(row[1].strip(' '),row[2]-1,row[3]) for row in npDf),from_string=True).sort()

mays = mays[mays['subgenome'].isnull() == False]
maysNpChr = [np.array(mays[mays['chr'] == str(i+1)]) for i in range(10)]
maysNpChrSplit = []
for chrom in maysNpChr:
    chromM = {'maize1':[],'maize2':[]}
    chromNew = []
    lineControl = chrom[0][-1]
    #print chrom[0]
    chromNew.append(chrom[0])
    for line in chrom[1:]:
        if line[-1] == 'maize1' and lineControl == 'maize1':
            chromNew.append(line)
        if line[-1] == 'maize1' and lineControl == 'maize2':
            #print line
            lineControl = 'maize1'
            chromM['maize2'].append(np.array(chromNew))
            chromNew=[]
            chromNew.append(line)
        if line[-1] == 'maize2' and lineControl == 'maize2':
            chromNew.append(line)
        if line[-1] == 'maize2' and lineControl == 'maize1':
            lineControl = 'maize2'
            chromM['maize1'].append(np.array(chromNew))
            chromNew = []
            chromNew.append(line)
    chromM[lineControl].append(np.array(chromNew))
    maysNpChrSplit.append(chromM)
#print maysNpChrSplit

mays1 = mays[mays['subgenome'] == 'maize1']
mays2 = mays[mays['subgenome'] == 'maize2']

#a = BedTool('284.bed')

#Mays1bed = df2bed(mays1)
#Mays2bed = df2bed(mays2)

#Mays1bed.saveas('m1.bed')
#Mays2bed.saveas('m2.bed')

#mays1Int = a.intersect(Mays1bed,wa=True)
#mays2Int = a.intersect(Mays2bed,wa=True)


mays1Chr = [mays1[mays1['chr'] == str(i+1)] for i in range(10)]
mays2Chr = [mays2[mays2['chr'] == str(i+1)] for i in range(10)]
#from Tkinter import Tk
#tk = Tk()
maysFasta = Fasta('Zmays_284_AGPv3.softmasked.fa')
if 0:
    for i in range(10):
        plt.figure()
        if mays1Chr[i].empty == 0 and mays2Chr[i].empty == 1:
            df = mays1Chr[i]
            df['xi'].plot.hist(color = 'b',bins=50,range=[0, maysFasta[str(i+1)][:].end])
            plt.title('M1_'+str(i+1))
        elif mays2Chr[i].empty == 0 and mays1Chr[i].empty == 0:
            plt.subplot(211)
            df = mays1Chr[i]
            df['xi'].plot.hist(color='b',bins=50,range=[0, maysFasta[str(i+1)][:].end])
            plt.title('M1_' + str(i+1))
            plt.subplot(212)
            df = mays2Chr[i]
            df['xi'].plot.hist(color='r',bins=50,range=[0, maysFasta[str(i+1)][:].end])
            plt.title('M2_' + str(i+1))
        elif mays2Chr[i].empty == 0 and mays1Chr[i].empty == 1:
            df = mays2Chr[i]
            df['xi'].plot.hist(color='r',bins=50,range=[0, maysFasta[str(i+1)][:].end])
            plt.title('M2_' + str(i+1))
        plt.savefig(str(i+1)+'.png')

#npMays1Chr = [np.array(mays1Chr[i]) for i in range(len(mays1Chr)) if mays1Chr[i].empty == 0]
#npMays2Chr = [np.array(mays2Chr[i]) for i in range(len(mays2Chr)) if mays2Chr[i].empty == 0]

#mays1ChrXiXfBed = [('M1_%s'%(npMays1Chr[i][0,1]),npMays1Chr[i][0,1],npMays1Chr[i][0,2]-1,npMays1Chr[i][len(npMays1Chr[i])-1,3]) for i in range(len(npMays1Chr))]
#mays2ChrXiXfBed = [('M2_%s'%(npMays2Chr[i][0,1]),npMays2Chr[i][0,1],npMays2Chr[i][0,2]-1,npMays2Chr[i][len(npMays2Chr[i])-1,3]) for i in range(len(npMays2Chr))]


mays1ChrXiXfBed = []
mays2ChrXiXfBed = []

for chrom in maysNpChrSplit:
    if 'maize1' in chrom.keys():
        mays1ChrXiXfBed += [('_'.join([npChr[0,1],str(npChr[0,2]),str(npChr[len(npChr)-1,3])]),npChr[0,1],npChr[0,2]-1,npChr[len(npChr)-1,3]) for npChr in chrom['maize1']]
    if 'maize2' in chrom.keys():
        mays2ChrXiXfBed += [('_'.join([npChr[0,1],str(npChr[0,2]),str(npChr[len(npChr)-1,3])]),npChr[0,1],npChr[0,2]-1,npChr[len(npChr)-1,3]) for npChr in chrom['maize2']]


with open('Maize1_SplitInfo.txt','w') as f:
    f.write('\t'.join(['SubgenomeChr','MaizeOriginalChr','xi','xf'])+'\n')
    f.write('\n'.join('\t'.join([str(x) for x in bedTuple]) for bedTuple in mays1ChrXiXfBed) + '\n')

with open('Maize2_SplitInfo.txt','w') as f:
    f.write('\t'.join(['SubgenomeChr','MaizeOriginalChr','xi','xf'])+'\n')
    f.write('\n'.join('\t'.join([str(x) for x in bedTuple]) for bedTuple in mays2ChrXiXfBed) + '\n')


if 0:
    global inputStr
    with open('ZmaysM1_301_AGPv3.fa','w') as f:
        for chrom in mays1ChrXiXfBed:
            inputStr = maysFasta[chrom[1]][chrom[2]:chrom[3]].seq
            f.write('>%s\n%s\n'%(chrom[0],wraplines()))#fill(,60, break_on_hyphens = False)))
    with open('ZmaysM2_302_AGPv3.fa', 'w') as f:
        for chrom in mays2ChrXiXfBed:
            inputStr = maysFasta[chrom[1]][chrom[2]:chrom[3]].seq
            f.write('>%s\n%s\n'%(chrom[0],wraplines()))#fill(maysFasta[chrom[1]][chrom[2]:chrom[3]].seq,60, break_on_hyphens = False)))

    mays1Fasta= Fasta('ZmaysM1_301_AGPv3.fa')

    mays2Fasta = Fasta('ZmaysM2_302_AGPv3.fa')

    #print mays1Fasta['M1_1'][:].seq == maysFasta['1'][mays1ChrXiXfBed[0][2]:mays1ChrXiXfBed[0][3]].seq

#b = BedTool('284.bed')

b = BedTool('q.PAC2_0.284.gff3')

def changeGFFCoordinates(gffstring,xi,chrName):
    gffOut = []
    start = 0
    for line in gffstring.split('\n'):
        if line and '##' not in line:
            if start == 1:
                lineList = line.split('\t')
                gffOut.append(chrName + '\t' + '\t'.join(lineList[1:3]) + '\t%d\t%d\t'%tuple(map(lambda x: int(x)-xi,lineList[3:5])) + '\t'.join(lineList[5:]) + '\n')
            if line.split('\t')[2] == 'gene' and start == 0:
                start = 1
                print line
                lineList = line.split('\t')
                gffOut.append(chrName + '\t' + '\t'.join(lineList[1:3]) + '\t%d\t%d\t' % tuple(map(lambda x: int(x) - xi, lineList[3:5])) + '\t'.join(lineList[5:]) + '\n')
        else:
            if '##' in line:
                #print line
                gffOut.append(line + '\n')
    return gffOut
M1gffLines = []
M2gffLines = []

gffOutLines = []
for chrom in mays1ChrXiXfBed:
    gffstr = str(b.intersect(BedTool('\t'.join([chrom[1],str(chrom[2]+1),str(chrom[3])]),from_string=True),wa=True,f=0.95))
    gffOutLines += changeGFFCoordinates(gffstr,chrom[2],chrom[0])
with open('301_presort.gff3','w') as f:
    f.writelines(gffOutLines)

#BedTool('q.PAC2_0.301.gff3').sort().saveas('q.PAC2_0.301.gff3')

gffOutLines = []
for chrom in mays2ChrXiXfBed:
    gffstr = str(b.intersect(BedTool('\t'.join([chrom[1],str(chrom[2]+1),str(chrom[3])]),from_string=True),wa=True,f=0.95))
    gffOutLines += changeGFFCoordinates(gffstr,chrom[2],chrom[0])
with open('302_presort.gff3','w') as f:
    f.writelines(gffOutLines)

runCommand = lambda x: subprocess.call(x, shell=True)
for command in ['python -m jcvi.formats.gff sort %s -o %s'%(gff[0],gff[1]) for gff in [('301_presort.gff3','q.PAC2_0.301.gff3'),('302_presort.gff3','q.PAC2_0.302.gff3')]]:
    runCommand(command)
    #BedTool('q.PAC2_0.302.gff3').sort().saveas('q.PAC2_0.302.gff3')