import os, sys
import subprocess, shutil
import pybedtools
from pyfaidx import Fasta
import pandas as pd
import multiprocessing as mp
import numpy as np
from time import time, clock
from collections import Counter, defaultdict
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn import metrics
import subprocess
from pyfaidx import Fasta
from pybedtools import BedTool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
#import peakutils
from sklearn.neighbors import KernelDensity
from collections import defaultdict
from scipy.signal import argrelextrema
#from sklearn.feature_selection import SelectKBest
from sklearn.cluster import FeatureAgglomeration
from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import KernelPCA
import plotly.graph_objs as go
import plotly.offline as py
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import SpectralClustering



def writeKmercount(args):
    """Takes list of fasta files and runs kmercountexact.sh to generate with only one column for pos and converts them to proper bedgraph format to be sorted"""
    fastaPath, kmercountPath = args
    for fastaFile in os.listdir(fastaPath):
        if fastaFile.endswith('.fa') or fastaFile.endswith('.fasta'):
            f = fastaFile.rstrip()
            print f
            scriptName = f[:f.rfind('.')] + '.sh'
            outFileName = f[:f.rfind('.')]+'.kcount'
            lineOutputList = [fastaPath, fastaFile, kmercountPath, outFileName]
            bbtoolsI = open(scriptName, 'w')
            bbtoolsI.write('#!/bin/bash\nmodule load bbtools\nkmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s/%s out=%s/%s k=23 -Xmx100g\n' % tuple(lineOutputList))
            bbtoolsI.close()
            try:
                subprocess.call('nohup sh %s' % scriptName, shell=True)
            except:
                print 'Unable to run %s via command line..' % outFileName

def kmer2Fasta(args):
    kmercountPath = args[0]
    for kmer in os.listdir(kmercountPath):
        if kmer.endswith('.kcount'):
            with open(kmercountPath+kmer,'r') as f, open(kmercountPath+kmer+'.fa','w') as f2:
                for line in f:
                    if line and int(line.split('\t')[-1]) >= 100:
                        f2.write('>%s\n%s\n'%tuple([line.split('\t')[0]]*2))


def writeBlast(args):
    """make blast database for whole genome assembly"""
    genome, blastPath, kmercountPath, fastaPath = args
    genomeName = genome[:genome.rfind('.')]
    dbscriptName = genomeName + '.db.sh'
    blastdb = open(dbscriptName, 'w')
    database_list = [fastaPath+genome, genomeName] #FIXME add genome path
    blastdb.write(
        '#!/bin/bash\nmodule load blast+\nmakeblastdb -in %s -dbtype nucl -out %s.blast_db\n' % tuple(database_list))
    blastdb.close()
    try:
        subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName
    """ blast differential kmers to whole genome assembly"""
    for file in [file2 for file2 in os.listdir(kmercountPath) if file2.endswith('.fa') or file2.endswith('.fasta')]:
        inputFile = kmercountPath+file#'%s.higher.kmers.txt' % (file.split('.')[0])
        f = file.rstrip()
        print f
        blastscriptName = '%s.blast.sh' % (file.split('.')[0])
        outFileName = f[:f.rfind('.')]+'.BLASTtsv.txt'
        lineOutputList = [genomeName, inputFile, blastPath, outFileName]
        blast = open(blastscriptName, 'w')
        blast.write('#!/bin/bash\nmodule load blast+/2.6.0\nblastn -db ./%s.blast_db -query %s -task "blastn-short" -outfmt 6 -out %s/%s -num_threads 8 -evalue 1e-2\n' % tuple(lineOutputList))
        blast.close()
        """ send blast jobs"""
        try:
            subprocess.call('nohup sh %s' % blastscriptName, shell=True)
        except:
            print 'Unable to run %s via command line..' % blastscriptName

def findScaffolds(args):
    with open('correspondence.bed','r') as f:
        lineList = sorted([line.split('\t')[-1].strip('\n') for line in f.readlines() if line])
    return lineList

def findKmerNames(args):
    kmercountPath, genome = args
    for file in os.listdir(kmercountPath):
        if file.startswith(genome[:genome.find('.fa')]) and (file.endswith('.fa') or file.endswith('.fasta')):
            print file
            with open(kmercountPath + file,'r') as f:
                listKmer = sorted([line.strip('\n') for line in f.readlines()[1::2] if line])
            return listKmer

def blast2bed(args):
    blastFile = args[0]
    with open(blastFile,'r') as f, open('blasted.bed','w') as f2:
        for line in f:
            if line:
                l1 = line.split('\t')[1]
                f2.write('\t'.join([l1] + ['0',str(int(l1.split('_')[-1])-int(l1.split('_')[-2]))] + [line.split('\t')[0]])+'\n')

    b = pybedtools.BedTool('blasted.bed').sort().merge(c=4,o='collapse',)
    b.saveas('blasted_merged.bed')
    #return b

def kmerRelatedHistogram(args):
    save = int(args[0])
    if save == 1:
        load = 0
    else:
        load = 1
    if save:
        with open('kmerPrevalence.txt', 'r') as f:
            kmerDict = {line.split('\t')[0]: line.count(',') for line in f if line}
        pickle.dump(kmerDict, open('kmerDict.p', 'wb'))
    if load:
        kmerDict = pickle.load(open('kmerDict.p', 'rb'))
    histData = np.array(kmerDict.values())[:, np.newaxis]
    reverseLookup = defaultdict(list)
    for k in kmerDict:
        reverseLookup[kmerDict[k]].append(k)
    if save:
        pickle.dump(reverseLookup, open('reverseLook.p', 'wb'))
    if load:
        reverseLookup = pickle.load(open('reverseLook.p', 'rb'))
    xplot = np.linspace(0, np.max(histData), 300000)[:, np.newaxis]
    del kmerDict

    kde = KernelDensity(kernel='gaussian', bandwidth=500).fit(histData)
    exp_log_dens = np.exp(kde.score_samples(xplot))

    idxs_peaks = argrelextrema(exp_log_dens, np.greater)[0]  # peakutils.indexes(exp_log_dens)
    idxs_valleys = argrelextrema(exp_log_dens, np.less)[0]  # peakutils.indexes(-exp_log_dens)

    fig = plt.figure()
    plt.plot(xplot[:], exp_log_dens, '-',
             label="Envelope Kernel '{0}' Density Function".format('Gaussian'))
    plt.plot(xplot[idxs_peaks], exp_log_dens[idxs_peaks], '*', color='r', label='Peaks')
    plt.plot(xplot[idxs_valleys], exp_log_dens[idxs_valleys], 'o', color='r', label='Valleys')
    plt.hist(histData[:, 0], bins=50, label='Histogram of Counts', normed=True)
    plt.xlabel('Number of Related Kmers (Normalized)')
    plt.ylabel('Count')
    plt.legend(loc='upper left')
    plt.title('Histogram of Number of Related Kmers (Normalized)')
    plt.savefig('KmerHistogram.png')

    intervals = [(0., xplot[idxs_valleys[0]][0])] + [(xplot[idxs_valleys[i]][0], xplot[idxs_valleys[i + 1]][0]) for i in
                                                     range(len(idxs_valleys) - 1)] + [
                    (xplot[idxs_valleys[-1]][0], xplot[-1][0])]
    counts = np.array(reverseLookup.keys())
    #print reverseLookup
    for interval in enumerate(intervals):
        try:
            #print interval
            keys = counts[np.where((counts <= np.floor(interval[1][1])) & (counts >= np.ceil(interval[1][0])))]
            #print keys
            # np.vectorize(lambda x: x <= interval[1] and x >= interval[0])(counts)
            with open('Peak%d_CountInt_%d_%d.fa' % tuple([interval[0] + 1] + map(int, list(interval[1]))), 'w') as f:
                for key in keys:
                    f.write('\n'.join('>%s\n%s' % (val, val) for val in reverseLookup[key]) + '\n')
            print 'Peak%d_CountInt_%d_%d.fa' % tuple([interval[0] + 1] + map(int, list(interval[1])))
        except:
            with open('ErrFile.txt', 'a') as f:
                f.write(str(interval[0]) + '\t%s' % str(interval[1]) + '\n')

def splitFasta(args):
    fastaFile = Fasta(args[0])

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
            """try:
                print posFinal[0:4]
                print key
            except:
                print 'length scaffold very small'"""
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
    subprocess.call('bedtools getfasta -fi %s -fo %s -bed %s -name'%(args[0],args[1] + args[0].split('.fa')[0]+'_split.fa','correspondence.bed'),shell=True)
    Fasta((args[1] + args[0]).split('.fa')[0]+'_split.fa')
    print args[0].split('.fa')[0]+'_split.fa'

def generateClusteringMatrixAndKmerPrevalence(args):
    kmercountPath, genome, save = args
    save = int(save)
    kmers = findKmerNames((kmercountPath, genome))
    scaffolds = findScaffolds(1)
    kmerDict = {kmer: [kmer] for kmer in kmers}
    # print kmerDict
    #d = defaultdict(list)
    dfMatrix = pd.DataFrame(columns=kmers,index=scaffolds)
    with open('blasted_merged.bed', 'r') as f:
        for line in f:
            if line:
                # print scaffold
                listLine = line.rstrip('\n').split('\t')
                # scaffold = listLine[0]
                counts = Counter(listLine[-1].split(','))
                interval = (abs(float(listLine[2]) - float(listLine[1]))) / 5000.
                # print scaffold, counts
                for key in counts:
                    try:
                        kmerDict[key] += counts.keys()
                        kmerDict[key] = list(set(kmerDict[key]))
                    except:
                        with open('keyErrors.txt', 'a') as f2:
                            f2.write(listLine[0] + ' ' + key + '\n')

                    # print kmerDict[key]
                    counts[key] /= interval
                try:
                    dfMatrix.loc[listLine[0]] = pd.Series(counts)
                except:
                    for key in counts:
                        if key not in kmers:
                            del counts[key]
                    dfMatrix.loc[listLine[0]] = pd.Series(counts)
                #d[listLine[0]] = counts
                # dfMatrix.set_value(scaffold, key, float(counts[key])/interval)

    with open('kmerPrevalence.txt', 'w') as f:
        for key in kmerDict:
            f.write('%s\t%s\t%d\n' % (key, ','.join(kmer for kmer in kmerDict[key]), len(kmerDict[key])))

    del kmerDict

    #dfMatrix = pd.DataFrame(d).fillna(0.).T
    dfMatrix = dfMatrix.fillna(0.)
    dfMatrix = dfMatrix.reset_index()
    dfMatrix.to_feather('clusteringMatrix.feather')
    # dfMatrix.to_csv('clusteringMatrix3.csv', index=True)
    kmers = list(dfMatrix.axes[1])
    scaffolds = list(dfMatrix.axes[0])
    with open('rowNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), scaffolds[i]]) for i in range(len(scaffolds))))
    with open('colNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), kmers[i]]) for i in range(len(kmers))))
    pickle.dump(scaffolds,open('scaffolds.p','wb'),protocol=2)
    pickle.dump(kmers,open('kmers.p','wb'),protocol=2)


def peakClusteringMatrix(args):
    kmercountPath, peakFasta, blastFile, save = args
    with open(blastFile, 'r') as f, open('%s_blasted.bed'%peakFasta.split('_')[0], 'w') as f2:
        for line in f:
            if line:
                l1 = line.split('\t')[1]
                f2.write('\t'.join(
                    [l1] + ['0', str(int(l1.split('_')[-1]) - int(l1.split('_')[-2]))] + [line.split('\t')[0]]) + '\n')

    b = pybedtools.BedTool('%s_blasted.bed'%peakFasta.split('_')[0]).sort().merge(c=4, o='collapse', )
    b.saveas('%s_blasted_merged.bed'%peakFasta.split('_')[0])
    with open(peakFasta,'r') as f:
        kmers = f.readlines()[1::2]
    scaffolds = pickle.load(open('scaffolds.p','rb'))
    save = int(save)
    # print kmerDict
    #d = defaultdict(list)
    dfMatrix = pd.DataFrame(columns=kmers,index=scaffolds)
    with open('%s_blasted_merged.bed'%peakFasta.split('_')[0], 'r') as f:
        for line in f:
            if line:
                # print scaffold
                listLine = line.rstrip('\n').split('\t')
                # scaffold = listLine[0]
                counts = Counter(listLine[-1].split(','))
                interval = (abs(float(listLine[2]) - float(listLine[1]))) / 5000.
                # print scaffold, counts
                for key in counts:
                    counts[key] /= interval
                try:
                    dfMatrix.loc[listLine[0]] = pd.Series(counts)
                except:
                    for key in counts:
                        if key not in kmers:
                            del counts[key]
                    dfMatrix.loc[listLine[0]] = pd.Series(counts)
                #d[listLine[0]] = counts
                # dfMatrix.set_value(scaffold, key, float(counts[key])/interval)


    dfMatrix = dfMatrix.fillna(0.)
    dfMatrix = dfMatrix.reset_index()
    dfMatrix.to_feather('%s_clusteringMatrix.feather'%peakFasta.split('_')[0])
    # dfMatrix.to_csv('clusteringMatrix3.csv', index=True)
    kmers = list(dfMatrix.axes[1])
    scaffolds = list(dfMatrix.axes[0])
    with open('rowNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), scaffolds[i]]) for i in range(len(scaffolds))))
    with open('colNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), kmers[i]]) for i in range(len(kmers))))
    transform_plot('1')


def transform_plot(args):
    peak = args[0]
    if peak == '1':
        df = pd.read_feather('clusteringMatrix.feather')
        peak = 'main'
    else:
        df = pd.read_feather('%s_clusteringMatrix.feather'%peak)
    df = df.set_index(['index'])
    scaffolds = list(df.axes[0])
    dimensionalityReducers = {'kpca':KernelPCA(n_components=3),'factor':FactorAnalysis(n_components=3),'feature':FeatureAgglomeration(n_clusters=3)}
    data = StandardScaler().fit_transform(df)
    for model in dimensionalityReducers:
        try:
            transformed_data = dimensionalityReducers[model].fit_transform(df)
            np.save('%s_%s_transformed3D.npy'%(peak,model), transformed_data)
            plots = []
            plots.append(
                go.Scatter3d(x=transformed_data[:, 0], y=transformed_data[:, 1], z=transformed_data[:, 2], name='Data',
                             mode='markers',
                             marker=dict(color='b', size=2), text=scaffolds))
            fig = go.Figure(data=plots)
            py.plot(fig, filename=peak + '_' + model + 'Reduction.html')
        except:
            pass

def cluster(args):
    file = args[0]
    global scaffolds

    Tname = file.split('transformed3D')[0]
    transformed_data = np.load(file)

    transformed_data = StandardScaler().fit_transform(transformed_data)

    clustering_names = ['SpectralClustering']
    n_clusters = 3
    spectral = SpectralClustering(n_clusters=n_clusters,
                                          eigen_solver='arpack',
                                          affinity="nearest_neighbors", gamma=0.3)

    clustering_algorithms = [spectral]

    for name, algorithm in zip(clustering_names, clustering_algorithms):
        try:
            algorithm.fit(transformed_data)
            if hasattr(algorithm, 'labels_'):
                y_pred = algorithm.labels_.astype(np.int)
            else:
                y_pred = algorithm.predict(transformed_data)
            N = len(set(y_pred))
            c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
            # plot
            plots = []
            clusterSize = defaultdict(list)
            for key in set(y_pred):
                # print key
                cluster_scaffolds = np.array(scaffolds)[y_pred == key]
                if list(cluster_scaffolds):
                    clusterSize[key] = len(cluster_scaffolds)
                    if clusterSize[key] == max(clusterSize.values()):
                        testCluster = key
                    plots.append(
                        go.Scatter3d(x=transformed_data[y_pred == key, 0], y=transformed_data[y_pred == key, 1],
                                     z=transformed_data[y_pred == key, 2],
                                     name='Cluster %d, %d points' % (key, len(cluster_scaffolds)), mode='markers',
                                     marker=dict(color=c[key], size=2), text=cluster_scaffolds))

            if hasattr(algorithm, 'cluster_centers_'):
                centers = algorithm.cluster_centers_
                plots.append(
                    go.Scatter3d(x=centers[:, 0], y=centers[:, 1], z=centers[:, 2], mode='markers',
                                 marker=dict(color='purple', symbol='circle', size=12),
                                 opacity=0.4,
                                 name='Centroids'))

            fig = go.Figure(data=plots)
            py.plot(fig, filename=name + Tname + 'n%d' % n_clusters + 'ClusterTest.html')

        except:
            print 'Unable to cluster completely using ' + name + ' for ' + Tname


#os.chdir('../../..')

funct = sys.argv[1]
arguments = sys.argv[2:]

options = {
    'writeKmerCount': writeKmercount,
    'kmer2Fasta':kmer2Fasta,
    'writeBlast':writeBlast,
    'findScaffolds':findScaffolds,
    'findKmerNames':findKmerNames,
    'blast2bed':blast2bed,
    'kmerRelatedHistogram':kmerRelatedHistogram,
    'splitFasta':splitFasta,
    'genClusterKmer':generateClusteringMatrixAndKmerPrevalence,
    'genPeakMatrix': peakClusteringMatrix,
    'transform_plot': transform_plot,
    'cluster': cluster
}

def main():
    options[funct](tuple(arguments))

main()