from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
sortList = open('q.PAC4GC.001.sort2','r').read().split('\n')
sortDict = defaultdict(list)
for line in sortList:
    if line:
        lineList = line.split()
        sortDict[lineList[0]] = [lineList[1]] + [int(item) for item in lineList[2:3]]
geneList = []
geneCompareFiles = open('BhD001003Analysis.txt','r').read().split('\n')
genelinesList = [line.split()[0].split(',')[2] for line in geneCompareFiles if len(line.split())>2 and ':' in line.split()[1] and line.split()[0].split(',')[0] == 'BhD1']

geneList = set(genelinesList)

totalGeneList = set(sortDict.keys())

finalGeneList = list(geneList & totalGeneList)

coordinategrab = []
for gene in finalGeneList:
    coordinategrab += sortDict[gene][1:2]
hist,bins = np.histogram(coordinategrab,bins=80)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)
plt.title('Chromosome 1 on ABR113 for 001-003 Mapping UnOut')
plt.xlabel('Position on Chromosome')
plt.ylabel('Count')
plt.show()
