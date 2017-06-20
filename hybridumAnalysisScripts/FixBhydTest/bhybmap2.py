from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
sortList = open('hyb1MAP2.txt','r').read().split('\n')
sortDict = defaultdict(list)
coordinategrab = []
for line in sortList:
    if line:
        lineList = line.split()
        coordinategrab += list(np.arange(int(lineList[1]),int(lineList[2]),1000))


hist,bins = np.histogram(coordinategrab,bins=80)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)
plt.title('Chromosome 1 on ABR113 for 001-003 Mapping-Bed File')
plt.xlabel('Position on Chromosome')
plt.ylabel('Count')
plt.show()
