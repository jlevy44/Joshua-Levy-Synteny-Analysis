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
from pyfaidx import Fasta
import sys

"""python generate"""
try:
    try:
        fastaFile, transposonGFF, numberChr, histInt, window, poly_degree, savgol_order, gaussian_bandwidth, option = sys.argv[1:]
    except:
        fastaFile, transposonGFF = sys.argv[1:3]
    # search for fai file and transposonDensity file
    """
    for file in os.listdir('.'):
        if 'cds' not in file and (file.endswith('.fa') or file.endswith('.fasta')):
            Fasta(file)
    faiFile = [file for file in os.listdir('.') if file.endswith('.fai')][0]
    transposonGFF = 'transposonDensity.gff'
    """
    try:
        numberChr = int(numberChr) + 1
        histInt = float(histInt)
        window = int(window)
        poly_degree = int(poly_degree)
        gaussian_bandwidth = int(gaussian_bandwidth)
        savgol_order = int(savgol_order)
    except:
        numberChr = 10
        histInt = 250000.
        window = 21
        poly_degree = 15
        savgol_order = 8
        gaussian_bandwidth = 4

    print """python extractCentroTelomere.py path_to_fasta repeat_GFF_file number_Chromosomes[default=10] histogramIntervalLength[default=250000.] Window_Size[default=21] Polynomial_Fit_Degree[default=15] Savgol_Filter_Order[default=8] Gaussian_Bandwidth[default=4] Option[centromere|telomere]"""


    if os.path.isfile(fastaFile + '.fai') and os.stat(fastaFile + '.fai').st_size > 0:
        faiFile = fastaFile + '.fai'
    else:
        Fasta(fastaFile)
        faiFile = fastaFile + '.fai'

    histInterval = defaultdict(list)
    chrList = []
    with open(faiFile, 'r') as f:
        chromCount = 0
        for line in f:
            chromCount += 1
            # bedread = '\n'.join('\t'.join('%s\t%d\t%d'%tuple([line.split('\t')[0]]+sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3])))))
            #interval = sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3]))
            chrList.append((line.split('\t')[0],int(line.split('\t')[1])))
    chromosomes = np.array(chrList)
    if numberChr < len(chromosomes):
        chromosomes = chromosomes[np.vectorize(int)(chromosomes[:,1]).argsort()[::-1][0:numberChr],:]
    #print chromosomes


    for i in range(len(chromosomes)):
        histInterval[chromosomes[i,0]] = list(np.arange(0., int(chromosomes[i,1]), histInt)) + [int(chromosomes[i,1])] # [0]

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
    #print transposonDensityArr
    centromereRegions = ''
    nonCentromereRegions = ''
    for chrom in histInterval:
        #print chrom
        #FIXME what if there are multiple main peaks??? eg. BhS9
        if len(histInterval[chrom]) > window:
            arraySubset = transposonDensityArr[transposonDensityArr[:,0] == chrom]
            #print arraySubset
            filtered_density = np.vectorize(np.poly1d(np.polyfit(range(len(arraySubset[:,0])),gaussian_filter1d(savgol_filter(arraySubset[:,-1],polyorder=savgol_order,window_length=window),gaussian_bandwidth),poly_degree)))(range(len(arraySubset[:,0])))#scipy.ndimage.filters., hilbert np.abs
            idxs_peaks = argrelextrema(filtered_density, np.greater)[0]
            #print len(filtered_density)
            for i,pk in enumerate(idxs_peaks):
                if pk >= len(filtered_density):
                    idxs_peaks[i] = len(filtered_density) - 1
                if pk < 0:
                    idxs_peaks[i] = 0
            #print chrom,len(filtered_density), idxs_peaks
            idxs_valleys = argrelextrema(filtered_density, np.less)[0]
            for i,pk in enumerate(idxs_valleys):
                if pk >= len(filtered_density):
                    idxs_valleys[i] = len(filtered_density) - 1
                if pk < 0:
                    idxs_valleys[i] = 0
            #print chrom,len(filtered_density),idxs_valleys,idxs_peaks
            maxIndex = np.argmax(filtered_density, axis=0)
            if any([pk > maxIndex for pk in idxs_valleys]) == 0:
                idxs_valleys = list(idxs_valleys) + [len(filtered_density) - 1]
            if any([pk < maxIndex for pk in idxs_valleys]) == 0:
                idxs_valleys = list(idxs_valleys) + [0]
            halfWidth = np.min(np.abs(idxs_valleys-maxIndex))
            centromereInterval = [maxIndex - halfWidth, maxIndex + halfWidth]
            fig = plt.figure()
            plt.plot(arraySubset[:,-1],color='r',zorder = 0)
            plt.plot(filtered_density,color = 'b', zorder = 1)
            for pk in centromereInterval:
                try:
                    plt.scatter(pk,filtered_density[pk],marker='o',color='g',s=30, zorder = pk)
                except:
                    print chrom, pk
                    plt.scatter(pk, 0., marker='o', color='y', s=30, zorder=pk)
            plt.savefig('transposonDensity_%s.png'%(chrom))
            centromereRegions += '%s\t%s\t%s\tCentromere\n'%(chrom,arraySubset[centromereInterval[0],1],arraySubset[centromereInterval[1],2])
            nonCentromereRegions += '%s\t%s\t%s\n%s\t%s\t%s\tNoCentromere\n'%(chrom,arraySubset[0,1],arraySubset[centromereInterval[0]-1,2],chrom,arraySubset[centromereInterval[0]+1,1],arraySubset[len(arraySubset)-1,2])

    BedTool(centromereRegions,from_string=True).sort().saveas('centromere.bed')
except:
    print """python extractCentroTelomere.py path_to_fasta repeat_GFF_file number_Chromosomes[default=10] histogramIntervalLength[default=250000.] Window_Size[default=21] Polynomial_Fit_Degree[default=15] Savgol_Filter_Order[default=8] Gaussian_Bandwidth[default=4] Option[centromere|telomere]"""
