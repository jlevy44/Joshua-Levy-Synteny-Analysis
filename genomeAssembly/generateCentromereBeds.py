import os
from pybedtools import BedTool
import numpy as np
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, hilbert
from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import argrelextrema

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
        #interval = sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3]))
        interval = line.split('\t')[1]
        histInterval[line.split('\t')[0]] = list(np.arange(0., int(interval), 250000.)) + [int(interval)] # [0]
        if chromCount > 22:
            break
    bedHist = BedTool('\n'.join('\n'.join(
        '\t'.join([key] + [str(int(x)) for x in [histInterval[key][i], histInterval[key][i + 1]]]) for i in
        range(len(histInterval[key]) - 1)) for key in histInterval.keys()), from_string=True).saveas('BedHist.bed')

bedTrans = bedHist.coverage(BedTool(transposonGFF)).sort().saveas('coverage.bed')

transposonDensity = []
for line in str(bedTrans).splitlines():
    if line:
        lineList = line.split('\t')
        transposonDensity.append(lineList[0:3]+[float(lineList[-1])])
transposonDensityArr = np.array(transposonDensity)
print transposonDensityArr
window = 21
for chrom in histInterval:
    print chrom
    if len(histInterval[chrom]) > window:
        arraySubset = transposonDensityArr[transposonDensityArr[:,0] == chrom]
        print arraySubset
        filtered_density = np.vectorize(np.poly1d(np.polyfit(range(len(arraySubset[:,0])),gaussian_filter1d(savgol_filter(arraySubset[:,-1],polyorder=8,window_length=window),4),15)))(range(len(arraySubset[:,0])))#scipy.ndimage.filters., hilbert np.abs
        idxs_peaks = argrelextrema(filtered_density, np.greater)[0]
        idxs_valleys = argrelextrema(filtered_density, np.less)[0]
        maxIndex = np.argmax(filtered_density, axis=0)
        fig = plt.figure()
        plt.plot(maxIndex,filtered_density[maxIndex],marker='*',color='g')
        plt.plot(arraySubset[:,-1],color='r')
        plt.plot(filtered_density,color = 'b')
        plt.savefig('transposonDensity_%s.png'%(chrom))
