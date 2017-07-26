import os
from pybedtools import BedTool
import numpy as np
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
"""python generate"""

# search for fai file and transposonDensity file

faiFile = [file for file in os.listdir('.') if file.endswith('.fai')][0]
transposonGFF = 'transposonDensity.gff'

histInterval = defaultdict(list)
with open(faiFile, 'r') as f:
    chromCount = 0
    for line in f:
        chromCount += 1
        # bedread = '\n'.join('\t'.join('%s\t%d\t%d'%tuple([line.split('\t')[0]]+sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3])))))
        interval = sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3]))
        histInterval[line.split('\t')[0]] = list(np.arange(0., interval[-1], 250000.)) + [interval[-1]]
        if chromCount > 22:
            break
    bedHist = BedTool('\n'.join('\n'.join(
        '\t'.join([key] + [str(int(x)) for x in [histInterval[key][i], histInterval[key][i + 1]]]) for i in
        range(len(histInterval[key]) - 1)) for key in histInterval.keys()), from_string=True)

bedTrans = bedHist.coverage(BedTool(transposonGFF)).sort()

transposonDensity = []
for line in str(bedTrans).splitlines():
    if line:
        lineList = line.split('\t')
        transposonDensity.append(lineList[0:3]+[float(lineList[-1])])
transposonDensityArr = np.array(transposonDensity)

for chrom in histInterval:
    arraySubset = transposonDensityArr[transposonDensityArr[:,0] == chrom]
    filtered_density = savgol_filter(arraySubset[:,-1])
    fig = plt.figure()
    plt.plot(arraySubset[:,-1],)
    plt.plot()
