import os, sys
import subprocess, shutil
import pybedtools
from pyfaidx import Fasta
import pandas as pd
import multiprocessing as mp
import numpy as np
from time import clock
from collections import Counter, defaultdict
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn import metrics
import subprocess
from pyfaidx import Fasta
from pybedtools import BedTool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import pyamg
#import peakutils
from sklearn.feature_selection import SelectKBest
from sklearn.pipeline import Pipeline
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
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
from sklearn.feature_selection import chi2
import scipy.sparse as sps
import networkx as nx
from itertools import combinations
from networkx.algorithms.components import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA



def writeKmercount(args):
    """Takes list of fasta files and runs kmercountexact.sh to generate with only one column for pos and converts them to proper bedgraph format to be sorted"""
    fastaPath, kmercountPath, kmerLength, blastMem = args
    kmerLengths = kmerLength.split(',')
    blastMemStr = "export _JAVA_OPTIONS='-Xms5G -Xmx%sG'"%(blastMem)
    for fastaFile in os.listdir(fastaPath):
        if (fastaFile.endswith('.fa') or fastaFile.endswith('.fasta')) and '_split' in fastaFile:
            kmerFiles = []
            f = fastaFile.rstrip()
            print f
            outFileNameFinal = f[:f.rfind('.')] + '.kcount'
            open(kmercountPath + '/' + outFileNameFinal, 'w').close()
            for kmerL in kmerLengths:
                scriptName = f[:f.rfind('.')] + kmerL + '.sh'
                outFileName = f[:f.rfind('.')]+kmerL+'.kcount'
                kmerFiles.append(kmercountPath+'/'+outFileName)
                lineOutputList = [fastaPath, fastaFile, kmercountPath, outFileName, kmerL]
                bbtoolsI = open(scriptName, 'w')
                bbtoolsI.write('#!/bin/bash\nmodule load bbtools\n'+blastMemStr+'\nkmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s/%s out=%s/%s k=%s -Xmx100g\n' % tuple(lineOutputList))
                bbtoolsI.close()
                try:
                    subprocess.call('sh %s' % scriptName, shell=True)#nohup
                except:
                    print 'Unable to run %s via command line..' % outFileName
                with open(kmercountPath + '/' + outFileName,'r') as f1, open(kmercountPath + '/' + outFileNameFinal, 'a') as f2:
                    f2.write(f1.read()+'\n')
            for fname in kmerFiles:
                os.remove(fname)
def kmer2Fasta(args):
    kmercountPath, kmer_low_count = args
    try:
        kmer_low_count = int(kmer_low_count)
    except:
        kmer_low_count = 100
    for kmer in os.listdir(kmercountPath):
        if kmer.endswith('.kcount'):
            with open(kmercountPath+kmer,'r') as f, open(kmercountPath+kmer+'.fa','w') as f2:
                for line in f:
                    if line and int(line.split('\t')[-1]) >= kmer_low_count:
                        f2.write('>%s\n%s\n'%tuple([line.split('\t')[0]]*2))

def writeBlast(args):
    """make blast database for whole genome assembly"""
    genome, blastPath, kmercountPath, fastaPath, BB, blastMem = args
    genomeName = genome[:genome.rfind('.')]
    blastMemStr = "export _JAVA_OPTIONS='-Xms5G -Xmx%sG'" % (blastMem)
    #dbscriptName = genomeName + '.db.sh'
    #blastdb = open(dbscriptName, 'w')
    #database_list = [fastaPath+genome, genomeName] #FIXME add genome path
    #blastdb.write(
    #    '#!/bin/bash\nmodule load blast+/2.6.0\nmakeblastdb -in %s -dbtype nucl -out %s.blast_db\n' % tuple(database_list))
    #blastdb.close()
    #try:
    #    subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    #except:
    #    print 'Unable to run %s via command line..' % dbscriptName
    """ blast differential kmers to whole genome assembly"""
    for file in [file2 for file2 in os.listdir(kmercountPath) if 'higher.kmers' in file2 and (file2.endswith('.fa') or file2.endswith('.fasta'))]:
        inputFile = kmercountPath+file#'%s.higher.kmers.txt' % (file.split('.')[0])
        f = file.rstrip()
        outFileName = f[:f.rfind('.')]+'.BLASTtsv.txt'
        lineOutputList = [genomeName, inputFile, blastPath, outFileName]
        if BB:
            subprocess.call(blastMemStr + ' && bbmap.sh vslow=t ambiguous=all noheader=t secondary=t perfectmode=t threads=8 maxsites=2000000000 outputunmapped=f ref=%s in=%s path=%s/ outm=%s'%(fastaPath+genome,inputFile,blastPath,blastPath+'/'+f[:f.rfind('.')]+'.sam'),shell=True)
        else:
            subprocess.call(blastMemStr + ' && module load blast+/2.6.0 && blastn -db ./%s.blast_db -query %s -task "blastn-short" -outfmt 6 -out %s/%s -num_threads 8 -evalue 1e-2' % tuple(lineOutputList),shell=True)

def findScaffolds(args=''):
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
    blastFile, BB ,lowMemory = args
    try:
        BB = int(BB)
        lowMemory = int(lowMemory)
    except:
        BB = 0
        lowMemory = 0
    if BB:
        with open(blastFile,'r') as f, open('blasted.bed','w') as f2:
            for line in f:
                if line:
                    l1 = line.split('\t')[2].split('::')[0]
                    f2.write('\t'.join([l1] + ['0', str(int(l1.split('_')[-1]) - int(l1.split('_')[-2]))] + [line.split('\t')[0]]) + '\n')
    else:
        with open(blastFile,'r') as f, open('blasted.bed','w') as f2:
            for line in f:
                if line:
                    l1 = line.split('\t')[1].split('::')[0] # FIXME .split('::')[0] added for blast+
                    f2.write('\t'.join([l1] + ['0',str(int(l1.split('_')[-1])-int(l1.split('_')[-2]))] + [line.split('\t')[0]])+'\n')
    if lowMemory == 0:
        b = pybedtools.BedTool('blasted.bed').sort().merge(c=4,o='collapse',)
        b.saveas('blasted_merged.bed')
    #return b

def kmerRelatedHistogram(args):
    kmercountName, save = args
    if int(save) == 1:
        load = 0
    else:
        load = 1
    if int(save):
        with open('kmerPrevalence.txt', 'r') as f:
            kmerDict = {line.split('\t')[0]: int(line.split('\t')[1]) for line in f if line}
            #kmerDict = {line.split('\t')[0]: line.count(',') for line in f if line}
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
    with open('PeaksOutNames.txt','w') as f2:
        for interval in enumerate(intervals):
            try:
                #print interval
                keys = counts[np.where((counts <= np.floor(interval[1][1])) & (counts >= np.ceil(interval[1][0])))]
                #print keys
                # np.vectorize(lambda x: x <= interval[1] and x >= interval[0])(counts)
                with open('Peak%d_CountInt_%d_%d.fa' % tuple([interval[0] + 1] + map(int, list(interval[1]))), 'w') as f:
                    for key in keys:
                        f.write('\n'.join('>%s\n%s' % (val, val) for val in reverseLookup[key]) + '\n')
                f2.write('Peak%d_CountInt_%d_%d.fa' % tuple([interval[0] + 1] + map(int, list(interval[1]))) + '\n')
            except:
                with open('ErrFile.txt', 'a') as f:
                    f.write(str(interval[0]) + '\t%s' % str(interval[1]) + '\n')

def splitFasta(args):
    fastaFile = Fasta(args[1] + args[0])
    try:
        splitLength = int(args[2])
    except:
        splitLength = 75000

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
            #if __name__ == '__main__':
            #p=mp.Pool(processes=8)
            splitLines = map(grabLine,posFinal)
            #p.close()
            #p.join()
            #print wrappedLines[0:10]
            return splitLines
        else:
            return ''
    bedText = []
    for key in fastaFile.keys():
        inputStr = fastaFile[key][:].seq
        bedText += split(inputStr,splitLength)
    with open('correspondence.bed','w') as f:
        f.write('\n'.join(bedText))
    BedTool('correspondence.bed').sort().saveas('correspondence.bed')
    subprocess.call('bedtools getfasta -fi %s -fo %s -bed %s -name'%(args[1] + args[0],args[1] + args[0].split('.fa')[0]+'_split.fa','correspondence.bed'),shell=True)
    Fasta((args[1] + args[0]).split('.fa')[0]+'_split.fa')
    print args[0].split('.fa')[0]+'_split.fa'

def generateClusteringMatrixAndKmerPrevalence(args):
    kmercountPath, save, genome, chunkSize, minChunkSize, removeNonChunk, minChunkThreshold, lowMemory = args
    try:
        removeNonChunk = int(removeNonChunk)
        chunkSize = int(chunkSize)
        minChunkSize = int(minChunkSize)
        minChunkThreshold = int(minChunkThreshold)
        lowMemory = int(lowMemory)
        save = int(save)
    except:
        chunkSize = 75000
        minChunkSize = 0
        minChunkThreshold = 0
        lowMemory = 0
        removeNonChunk = 0
        save = 0

    kmers = findKmerNames((kmercountPath, genome))
    scaffolds = findScaffolds(1)
    if removeNonChunk:
        scaffolds = list(filter(lambda scaffold: abs(int(scaffold.split('_')[-1]) - int(scaffold.split('_')[-2])) == chunkSize,scaffolds))
    elif minChunkThreshold:
        scaffolds = list(filter(lambda scaffold: abs(int(scaffold.split('_')[-1]) - int(scaffold.split('_')[-2])) >= minChunkSize,scaffolds))
    scaffoldIdx = {scaffold:i for i,scaffold in enumerate(scaffolds)}
    kmerIdx = {kmer:i for i,kmer in enumerate(kmers)} # np.uint32(
    #kmerCount = {kmer:[0.,0] for kmer in kmers}
    # print kmerDict
    #d = defaultdict(list)
    #dfMatrix = pd.DataFrame(columns=kmers,index=scaffolds)

    #open('kmerPrevalence.bed','w').close()
    #open('Progress.txt', 'w').close()
    #data = np.zeros((len(scaffolds),len(kmers)),dtype=float)
    #row_idx = []
    #column_idx = []
    #values = []
    #kmer_SparseMatrix = sps.dok_matrix((len(kmers),len(kmers)),dtype=np.int8)
    #kmer_SparseMatrix = [[i] for i in range(len(kmers))]
    #f2 = open('kmerPrevalence.bed','a')
    data = sps.dok_matrix((len(scaffolds), len(kmers)),dtype=np.float32)
    #f3 = open('Progress.txt', 'w')
    #kmerGraph = nx.Graph()
    #c = 0
    #d=0
    start = clock()
    #f3.write(str(start) + '\n')
    if lowMemory:
        with open('blasted.bed', 'r') as f:  # , open(classifyFolder + 'blasted.bed', 'w') as f2:
            for line in f:
                if line:
                    lineL = line.strip('\n').split('\t')
                    data[scaffoldIdx[lineL[0]], kmerIdx[lineL[-1]]] += 1.
    else:
        with open('blasted_merged.bed', 'r') as f:
            for line in f:
                if line:
                    listLine = line.rstrip('\n').split('\t')
                    if listLine[0] in scaffolds:
                        counts = Counter(listLine[-1].split(','))
                        interval = (abs(float(listLine[2]) - float(listLine[1]))) / 5000.
                        #scalingFactor = interval * 5000. / len(kmers)
                        #kmerText = ','.join(counts.keys())
                        #kmerCountsNum = len(counts.keys())
                        #f3.write(listLine[0] + '\n')
                        #kmerCountsNumbers = [kmerIdx[kmer] for kmer in kmerCountsNames]
                        #kmerGraph.add_edges_from(list(combinations(kmerCountsNumbers, 2)))
                        for key in counts:
                            #f2.write('%s\t0\t100\t%s\n'%(key,kmerText))
                            try:
                                #row_idx.append(scaffoldIdx[listLine[0]])
                                #column_idx.append(kmerIdx[key])
                                #values.append(counts[key] / interval)
                                data[scaffoldIdx[listLine[0]], kmerIdx[key]] = counts[key] / interval
                                #kmerCount[key][0] += kmerCountsNum #/ scalingFactor
                                #kmerCount[key][1] += 1
                            except:
                                pass
                    #[(np.uint32(kmerIdx[key]),np.uint32(kmerIdx[kmername])) for kmername in kmerCountsNames])
                        #for kmername in kmerCountsNames:
                        #    kmer_SparseMatrix[kmerIdx[key],kmerIdx[kmername]] = 1
                    #c += 1
                    #if c > 100:
                    #    d += c
                    #f3.write(str(d) + '\t' + str(clock()-start) + '\n')
                    #    c=0
                    #    print line

                        #f2.close()
                        #f3.write(str(d) + '\n' + line + '\n' + str(counts) + '\n')
                        #BedTool('kmerPrevalence.bed').sort().merge(c=4, o='distinct', delim=',').saveas(
                        #    'kmerPrevalence.bed')
                        #f2 = open('kmerPrevalence.bed','a')
                    #    c = 0

                        #data[scaffoldIdx[listLine[0]],kmerIdx[key]] = counts[key] / interval
    #f2.close()
    #f3.close()
    #with open('kmerPrevalence.txt','w') as f:
    #    for kmer, count in {key:(int(np.ceil(float(kmerCount[key][0])/kmerCount[key][1])) if kmerCount[key][1] > 0 else 0) for key in kmerCount}.items():#[(kmers[i],kmerGraph.degree(i)) for i in kmerGraph.nodes()]:#zip([kmers[i] for i in kmerGraph.nodes()],[len(kmerGraph.neighbors(kmer2)) for kmer2 in kmerGraph.nodes()]): #list(np.vectorize(int)(np.asarray(np.squeeze((kmer_SparseMatrix.sum(axis=1))))[0]))
    #        f.write('%s\t%d\n'%(kmer,count))
    #data = data[:,~np.all(data==0.,axis=0)]
    #BedTool('kmerPrevalence.bed').sort().merge(c=4,o='distinct',delim=',').saveas('kmerPrevalence.bed')

    #with open('kmerPrevalence.bed') as f, open('kmerPrevalence.txt', 'w') as f2:
    #    for line in f:
    #        if line:
    #            lineList = line.strip('\n').split('\t')
    #            f2.write('%s\t%s\n'%(lineList[0],','.join(set(lineList[-1].split(',')))))
        #for key in kmerDict:
        #    f.write('%s\t%s\t%d\n' % (key, ','.join(kmer for kmer in kmerDict[key]), len(kmerDict[key])))

    #del scaffoldIdx,kmerIdx
    #data = sps.coo_matrix((values, (row_idx, column_idx)), (len(scaffolds), len(kmers)))
    del scaffoldIdx, kmerIdx#, values, row_idx, column_idx
    #dfMatrix = pd.DataFrame(d).fillna(0.).T
    #dfMatrix = pd.DataFrame(data,index=scaffolds)
    #del data
    #dfMatrix = dfMatrix.reset_index()
    #dfMatrix.to_feather('clusteringMatrix.feather')
    # dfMatrix.to_csv('clusteringMatrix3.csv', index=True)
    #kmers = list(dfMatrix.axes[1])
    #scaffolds = list(dfMatrix.axes[0])
    sps.save_npz('clusteringMatrix.npz',data.tocsc())
    with open('rowNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), scaffolds[i]]) for i in range(len(scaffolds))))
    with open('colNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), kmers[i]]) for i in range(len(kmers))))
    pickle.dump(scaffolds,open('scaffolds.p','wb'),protocol=2)
    pickle.dump(kmers,open('kmers.p','wb'),protocol=2)

def transform_main(args): #FIXME
    try:
        main, reclusterFolder, model, n_subgenomes, metric  = args
    except:
        main, reclusterFolder, model,n_subgenomes, metric = ('1','null','kpca','2','straight')
    transform_plot((main,reclusterFolder,model,n_subgenomes,metric))

def peakClusteringMatrix(args):
    kmercountPath, peakFasta, blastFile, save = args
    with open(blastFile, 'r') as f, open('%s_blasted.bed'%peakFasta.split('_')[0], 'w') as f2:
        for line in f:
            if line:
                l1 = line.split('\t')[1]
                f2.write('\t'.join(
                    [l1] + ['0', str(int(l1.split('_')[-1]) - int(l1.split('_')[-2]))] + [line.split('\t')[0]]) + '\n')
    b = pybedtools.BedTool('%s_blasted.bed'%peakFasta.split('_')[0]).sort().merge(c=4, o='collapse' )
    b.saveas('%s_blasted_merged.bed'%peakFasta.split('_')[0])
    with open(peakFasta,'r') as f:
        kmers = f.readlines()[1::2]
    scaffolds = pickle.load(open('scaffolds.p','rb'))
    scaffoldIdx = {scaffold: i for i, scaffold in enumerate(scaffolds)}
    kmerIdx = {kmer: i for i, kmer in enumerate(kmers)}
    save = int(save)
    # print kmerDict
    #d = defaultdict(list)
    #dfMatrix = pd.DataFrame(columns=kmers,index=scaffolds)
    #data = np.zeros((len(scaffolds), len(kmers)), dtype=float)
    row_idx = []
    column_idx = []
    values = []
    with open('%s_blasted_merged.bed'%peakFasta.split('_')[0], 'r') as f:
        for line in f:
            if line:
                listLine = line.rstrip('\n').split('\t')
                counts = Counter(listLine[-1].split(','))
                interval = (abs(float(listLine[2]) - float(listLine[1]))) / 5000.
                for key in counts:
                    row_idx.append(scaffoldIdx[listLine[0]])
                    column_idx.append(kmerIdx[key])
                    values.append(counts[key] / interval)
                    #data[scaffoldIdx[listLine[0]],kmerIdx[key]] = counts[key] / interval

    #data = data[:,~np.all(data==0.,axis=0)]
    data = sps.coo_matrix((values,(row_idx,column_idx)),(len(scaffolds),len(kmers)))
    del scaffoldIdx, kmerIdx, values, row_idx, column_idx
    # dfMatrix = pd.DataFrame(d).fillna(0.).T
    #dfMatrix = pd.SparseDataFrame(data, index=scaffolds)
    sps.save_npz('%s_clusteringMatrix.npz'%peakFasta.split('_')[0],data)
    #dfMatrix = dfMatrix.fillna(0.)
    #dfMatrix = dfMatrix.reset_index()
    #dfMatrix.to_feather('%s_clusteringMatrix.feather'%peakFasta.split('_')[0])
    # dfMatrix.to_csv('clusteringMatrix3.csv', index=True)

    #kmers = list(dfMatrix.axes[1])
    #scaffolds = list(dfMatrix.axes[0])
    with open('rowNames.txt', 'w') as f:
        f.write('\n'.join('\t'.join([str(i), scaffolds[i]]) for i in range(len(scaffolds))))
    #with open('colNames.txt', 'w') as f:
    #    f.write('\n'.join('\t'.join([str(i), kmers[i]]) for i in range(len(kmers))))


def transform_plot(args):
    peak, reclusterFolder, model, n_subgenomes, metric = args
    try:
        n_subgenomes = int(n_subgenomes)
    except:
        n_subgenomes = 2
    if peak == '1':
        peak = 'main'
    if os.path.exists(peak + '_' + model + '_%d'%(n_subgenomes) + 'Reduction.html') == 0:#isfile
        if peak == 'main':
            # df = pd.read_feather('clusteringMatrix.feather')
            data = sps.load_npz('clusteringMatrix.npz')
        else:
            # df = pd.read_feather('%s_clusteringMatrix.feather'%peak)
            data = sps.load_npz('%s/%s_clusteringMatrix.npz' % (reclusterFolder, peak))
        scaffolds = pickle.load(open('scaffolds.p', 'rb'))
        # df = df.set_index(['index'])
        # scaffolds = list(df.axes[0])
        N = n_subgenomes + 1 #FIXME change the kernel
        dimensionalityReducers = {'kpca': KernelPCA(n_components=N,kernel=metric), 'factor': FactorAnalysis(n_components=N),
                                  'feature': FeatureAgglomeration(n_clusters=N)}
        data = StandardScaler(with_mean=False).fit_transform(data)
        # for model in dimensionalityReducers:
        if model != 'kpca' and peak.startswith('recluster') == 0:
            data = KernelPCA(n_components=499).fit_transform(data)
        if peak.startswith('recluster'):
            data = data.toarray()
        transformed_data = dimensionalityReducers[model].fit_transform(data)
        np.save('%s_%s_%d_transformed3D.npy'%(peak,model,n_subgenomes), transformed_data)
        if n_subgenomes > 2:
            transformed_data = KernelPCA(n_components=3,kernel=metric).fit_transform(transformed_data)
        plots = []
        plots.append(
            go.Scatter3d(x=transformed_data[:, 0], y=transformed_data[:, 1], z=transformed_data[:, 2], name='Data',
                         mode='markers',
                         marker=dict(color='b', size=2), text=scaffolds))
        fig = go.Figure(data=plots)
        py.plot(fig, filename=peak + '_' + model + '_%d'%(n_subgenomes) + 'Reduction.html')
    else:
        subprocess.call('touch %s'%(peak + '_' + model + '_%d'%(n_subgenomes) + 'Reduction.html'),shell=True)

def cluster(args):
    file, reclusterFolder, kmer500Path, clusterMethod, n_subgenomes, metric, n_neighbors = args
    print clusterMethod
    try:
        n_subgenomes = int(n_subgenomes)
    except:
        n_subgenomes = 2
    try:
        n_neighbors = int(n_neighbors)
    except:
        n_neighbors = 10
    n_clusters = n_subgenomes + 1
    clustering_algorithms = {'SpectralClustering': SpectralClustering(n_clusters=n_clusters, eigen_solver='amg', affinity= 'precomputed', random_state=42),#,gamma=1),arpack#amg,affinity="nearest_neighbors")#, n_neighbors=30, gamma=1),# nearestneighbors
                             'KMeans': MiniBatchKMeans(n_clusters=n_clusters)}
    name, algorithm = clusterMethod , clustering_algorithms[clusterMethod]
    #metric = 'cosine'
    if 'recluster' not in file:
        dataOld = sps.load_npz('clusteringMatrix.npz')
        scaffolds = pickle.load(open('scaffolds.p', 'rb'))
        kmers = pickle.load(open('kmers.p', 'rb'))
        #kmerIdx = {i : kmer for i, kmer in enumerate(kmers)}
        Tname = file.split('transformed3D')[0]
        transformed_data = np.load(file)

        transformed_data = StandardScaler().fit_transform(transformed_data)

        #clustering_names = ['SpectralClustering','KMeans']
        #n_clusters = 3

        #clustering_algorithms = [SpectralClustering(n_clusters=n_clusters,eigen_solver='amg',affinity="nearest_neighbors", gamma=1),MiniBatchKMeans(n_clusters=3)]

        #for name, algorithm in zip(clustering_names, clustering_algorithms):
            #try:
        if os.path.exists(name + Tname + 'n%d' % n_clusters + 'ClusterTest.html') == 0:
            try:
                os.mkdir('analysisOutputs/' + name + Tname + 'n%d' % n_clusters)
            except:
                pass

            if clusterMethod == 'SpectralClustering':
                #try:
                neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm = 'brute' , metric=metric)
                #except:
                #    print 'NN not working'
                neigh.fit(transformed_data)
                fit_data = neigh.kneighbors_graph(transformed_data)
                connected = sps.csgraph.connected_components(fit_data)
                if connected[0] > 1:
                    #D_num = fit_data.sum(axis=0,dtype=np.uint8)
                    #D = sps.dok_matrix((len(scaffolds), len(scaffolds)), dtype=np.uint8)
                    #for i in range(len(fit_data)):
                    #    D[i,i] = D_num[i]
                    #D = D.tocsc()
                    #L = D-fit_data
                    counts = Counter(connected[1])
                    subgraph_idx = max(counts.iteritems(), key=lambda x: x[1])[0]
                    scaffBool = connected[1] == subgraph_idx
                    #G = nx.Graph(fit_data)#nx.from_scipy_sparse_matrix(fit_data)
                    #mapping = {i: scaffolds[i] for i in range(len(scaffolds))}
                    #print G.nodes()
                    #G = nx.relabel_nodes(G, mapping, copy=False)
                    #print scaffolds[0:10], G.nodes()[0:10]
                    #print list(connected_component_subgraphs(G))[0]
                    #print sorted(nx.connected_components(G), key=len, reverse=True)
                    #gs = list(connected_component_subgraphs(G))[0]
                    #scaffolds2 = [scaffold for scaffold in scaffolds if scaffold in gs.nodes()]
                    #print scaffolds2
                    #s2 = gs.nodes()
                    #scaffBool = np.vectorize(lambda scaffold: scaffold in s2)(scaffolds)
                    #fit_data = nx.to_scipy_sparse_matrix(gs)
                    #fit_data = fit_data[scaffBool].T[scaffBool]
                    dataOld = dataOld[scaffBool]
                    print transformed_data[0:10,:],scaffolds[0:10]
                    transformed_data = transformed_data[scaffBool,:]
                    scaffolds_noconnect = list(np.array(scaffolds)[scaffBool == False])
                    #print scaffolds[0:10],np.array(scaffolds2)[0:10]#, np.array(scaffolds2)[0:10]
                    scaffolds = list(np.array(scaffolds)[scaffBool])
                    print transformed_data[0:10,:],scaffolds[0:10]
                    n_connected = connected[0]
                    while(n_connected > 1):
                        neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm='brute', metric=metric)
                        neigh.fit(transformed_data)
                        fit_data = neigh.kneighbors_graph(transformed_data)
                        connected = sps.csgraph.connected_components(fit_data)
                        counts = Counter(connected[1])
                        subgraph_idx = max(counts.iteritems(), key=lambda x: x[1])[0]
                        scaffBool = connected[1] == subgraph_idx
                        #fit_data = fit_data[scaffBool].T[scaffBool]
                        if connected[0] > 1:
                            dataOld = dataOld[scaffBool]
                            transformed_data = transformed_data[scaffBool, :]
                            scaffolds_noconnect += list(np.array(scaffolds)[scaffBool == False])
                            scaffolds = list(np.array(scaffolds)[scaffBool])
                        n_connected = connected[0]

                else:
                    scaffolds_noconnect = []

                # FIXME do something with thrown out points in future
                if n_subgenomes > 2:
                    t_data = KernelPCA(n_components=3).fit_transform(transformed_data)
                else:
                    t_data = transformed_data
                np.save('analysisOutputs/' + name + Tname + 'n%d' % n_clusters +'/graphInitialPositions.npy', t_data)
                del t_data
                sps.save_npz('analysisOutputs/' + name + Tname + 'n%d' % n_clusters +'/spectralGraph.npz', fit_data.tocsc())
                pickle.dump(scaffolds,open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/scaffolds_connect.p', 'wb'))
                pickle.dump(scaffolds_noconnect,open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/scaffolds_noconnect.p', 'wb'))
                #clusterGraph(('analysisOutputs/' + name + Tname + 'n%d' % n_clusters +'/spectralGraph.npz','analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/scaffolds_connect.p','analysisOutputs/' + name + Tname + 'n%d' % n_clusters))
                # above takes too much time!!! fix please :)
                # sps.save_npz('spectralGraph.npz',fit_data)
                # print fit_data.toarray()
            else:
                fit_data = transformed_data
            algorithm.fit(fit_data)
            if n_subgenomes > 2:
                reduction = KernelPCA(n_components=3)
                reduction.fit(transformed_data)
                reductionT = reduction.transform(transformed_data)
                scaledfit = StandardScaler()
                scaledfit.fit(reductionT)
                transformed_data2 = scaledfit.transform(reductionT)
            else:
                transformed_data2 = transformed_data
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
                print key, y_pred[0:10], cluster_scaffolds[0:10], scaffolds[0:10]
                if list(cluster_scaffolds):
                    clusterSize[key] = np.mean(np.apply_along_axis(lambda x: np.linalg.norm(x),1,transformed_data[y_pred == key,:]))#len(cluster_scaffolds)
                    if clusterSize[key] == min(clusterSize.values()):
                        testCluster = key
                    plots.append(
                        go.Scatter3d(x=transformed_data2[y_pred == key, 0], y=transformed_data2[y_pred == key, 1],
                                     z=transformed_data2[y_pred == key, 2],
                                     name='Cluster %d, %d points, %f distance' % (key, len(cluster_scaffolds),clusterSize[key]), mode='markers',
                                     marker=dict(color=c[key], size=2), text=cluster_scaffolds))

            if hasattr(algorithm, 'cluster_centers_'):
                if n_subgenomes > 2:
                    centers = scaledfit.transform(reduction.transform(algorithm.cluster_centers_)) #FIXME modify!!! need scaler fitted model from transformeddata
                else:
                    centers = algorithm.cluster_centers_
                plots.append(
                    go.Scatter3d(x=centers[:, 0], y=centers[:, 1], z=centers[:, 2], mode='markers',
                                 marker=dict(color='purple', symbol='circle', size=12),
                                 opacity=0.4,
                                 name='Centroids'))
            for key in set(y_pred)-{testCluster}:
                with open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/subgenome_%d.txt' % key, 'w') as f:
                    f.write('\n'.join(np.array(scaffolds)[y_pred == key]))

            fig = go.Figure(data=plots)

            #trainData = transformed_data[y_pred != testCluster]
            #trainData_scaffolds = np.array(scaffolds)[y_pred != testCluster]


            trainLabels = y_pred[y_pred != testCluster]
            trainData = dataOld[y_pred != testCluster]
            kbest = SelectKBest(chi2,k = 500)
            kbest.fit(trainData,trainLabels)
            bestFeatures = kbest.pvalues_.argsort()[:500]
            best_500_kmers = [kmers[i] for i in bestFeatures]
            sps.save_npz('%s/recluster%s_clusteringMatrix.npz' %(reclusterFolder, name + Tname + 'n%d' % n_clusters), dataOld[:,bestFeatures])#.tocsc()#.tocoo())
            with open('%s/kmer500Best_%s.fa'%(kmer500Path,name + Tname + 'n%d' % n_clusters),'w') as f:
                f.write('\n'.join('>%s\n%s'%(kmer,kmer) for kmer in best_500_kmers))
            subprocess.call('touch ' + 'analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '.txt',
                            shell=True)

            py.plot(fig, filename=name + Tname + 'n%d' % n_clusters + 'ClusterTest.html')
            #except:
                #print 'Unable to cluster completely using ' + name + ' for ' + Tname
    else:
        print 'Exception!!!'
        if clusterMethod != 'SpectralClustering':
            scaffolds = pickle.load(open('scaffolds.p', 'rb'))
            # kmerIdx = {i : kmer for i, kmer in enumerate(kmers)}
            Tname = file.split('transformed3D')[0]
            transformed_data = np.load(file)

            transformed_data = StandardScaler().fit_transform(transformed_data)

            #clustering_names = ['SpectralClustering', 'KMeans']
            #n_clusters = 3

            #clustering_algorithms = {'SpectralClustering':SpectralClustering(n_clusters=n_clusters, eigen_solver='amg', affinity="nearest_neighbors", gamma=1),
            #                         'KMeans':MiniBatchKMeans(n_clusters=3)}


            #for name, algorithm in zip(clustering_names, clustering_algorithms):
                #try:
            if os.path.exists(name + Tname + 'n%d' % n_clusters + '500kmerreclusterTest.html') == 0: #.isfile
                algorithm.fit(transformed_data)
                if n_subgenomes > 2:
                    reduction = KernelPCA(n_components=3)
                    reduction.fit(transformed_data)
                    reductionT = reduction.transform(transformed_data)
                    scaledfit = StandardScaler()
                    scaledfit.fit(reductionT)
                    transformed_data2 = scaledfit.transform(reductionT)
                else:
                    transformed_data2 = transformed_data
                if hasattr(algorithm, 'labels_'):
                    y_pred = algorithm.labels_.astype(np.int)
                else:
                    y_pred = algorithm.predict(transformed_data)
                N = len(set(y_pred))
                c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
                # plot
                plots = []
                clusterSize = defaultdict(list)
                try:
                    os.mkdir('analysisOutputs/' + name + Tname + 'n%d' % n_clusters)
                except:
                    pass
                for key in set(y_pred):
                    # print key
                    cluster_scaffolds = np.array(scaffolds)[y_pred == key]
                    if list(cluster_scaffolds):
                        clusterSize[key] = np.mean(np.apply_along_axis(lambda x: np.linalg.norm(x), 1,
                                                                       transformed_data[y_pred == key,
                                                                       :]))  # len(cluster_scaffolds)
                        if clusterSize[key] == min(clusterSize.values()):
                            testCluster = key
                        plots.append(
                            go.Scatter3d(x=transformed_data2[y_pred == key, 0],
                                         y=transformed_data2[y_pred == key, 1],
                                         z=transformed_data2[y_pred == key, 2],
                                         name='Cluster %d, %d points, %f distance' % (key, len(cluster_scaffolds),clusterSize[key]),
                                         mode='markers',
                                         marker=dict(color=c[key], size=2), text=cluster_scaffolds))

                if hasattr(algorithm, 'cluster_centers_'):
                    if n_subgenomes > 2:
                        centers = scaledfit.transform(reduction.transform(algorithm.cluster_centers_)) #FIXME modify
                    else:
                        centers = algorithm.cluster_centers_
                    plots.append(
                        go.Scatter3d(x=centers[:, 0], y=centers[:, 1], z=centers[:, 2], mode='markers',
                                     marker=dict(color='purple', symbol='circle', size=12),
                                     opacity=0.4,
                                     name='Centroids'))
                for key in set(y_pred)-{testCluster}:
                    with open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/recluster_subgenome_%d.txt' % key, 'w') as f:
                        f.write('\n'.join(np.array(scaffolds)[y_pred == key]))

                fig = go.Figure(data=plots)
                subprocess.call('touch ' + 'analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '.txt',shell=True)
                py.plot(fig, filename=name + Tname + 'n%d' % n_clusters + '500kmerreclusterTest.html')
                #except:
                #    print 'Unable to cluster completely using ' + name + ' for ' + Tname

def fai2bed(args):
    genome = args[0]
    Fasta(genome)
    bedFastaDict = defaultdict(list)
    with open(genome+'.fai','r') as f, open(genome + '.bed','w') as f2:
        for line in f:
            if line:
                lineList = line.split('\t')
                bedline = '%s\t%s\t%s\t%s\n'%(lineList[0],'0',lineList[1],lineList[0])
                f2.write(bedline)
                bedFastaDict[lineList[0]] = [bedline]
    return bedFastaDict

def writeKmerCountSubgenome(args):
    subgenomeFolder, kmerLength, blastMem, diff_kmer_threshold  = args
    blastMemStr = "export _JAVA_OPTIONS='-Xms5G -Xmx%sG'" % (blastMem)
    kmerLengths = kmerLength.split(',')
    try:
        os.mkdir(subgenomeFolder+'/kmercount_files/')
    except:
        pass
    kmercountPath = subgenomeFolder+'/kmercount_files/'
    for fastaFile in os.listdir(subgenomeFolder):
        if 'higher.kmers' not in fastaFile and '_split' not in fastaFile and (fastaFile.endswith('.fa') or fastaFile.endswith('.fasta')):
            kmerFiles = []
            f = fastaFile.rstrip()
            outFileNameFinal = f[:f.rfind('.')] + '.kcount'
            open(kmercountPath + '/' + outFileNameFinal, 'w').close()
            for kmerL in kmerLengths:
                print f
                outFileName = f[:f.rfind('.')] + kmerL + '.kcount'
                kmerFiles.append(kmercountPath + '/' + outFileName)
                lineOutputList = [subgenomeFolder+'/', fastaFile, kmercountPath, outFileName,kmerL]
                subprocess.call(blastMemStr + ' && module load bbtools && kmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s/%s out=%s/%s k=%s -Xmx60g' % tuple(
                        lineOutputList),shell=True)
                with open(kmercountPath + '/' + outFileName, 'r') as f1, open(kmercountPath + '/' + outFileNameFinal,'a') as f2:
                    f2.write(f1.read() + '\n')
            for fname in kmerFiles:
                os.remove(fname)

    compareKmers(diff_kmer_threshold, kmercountPath)

def kmercounttodict(kmercount2fname,kmercountPath):
    """ kmercounttodict function creates kmer : count key value pairs, takes path and file name of a kmer count file"""
    inputFile = open(kmercountPath+kmercount2fname,'r')
    print 'my input to kcount to dict is: %s' % inputFile
    dictConverted = {}
    for line in inputFile:
        if line and len(line.split()) == 2:
            lineList = line.split()
            dictConverted[lineList[0]]=(int(lineList[1].strip('\n')))
    inputFile.close()
    return dictConverted

def compareKmers(diff_kmer_threshold,kmercountPath):
    try:
        ratio_threshold = int(diff_kmer_threshold)
    except:
        ratio_threshold = 20
    dictOfGenes = {}
    #end_dinucleotide = 'GG'
    kmercountFiles = [file for file in os.listdir(kmercountPath) if file.endswith('.kcount') and '_split' not in file]
    for file in kmercountFiles:
        # creates a dictionary that associates a species to its dictionary of the kmer : count key value pairs
        # kmercounttodict function is called to create the kmer : count key value pairs
        dictOfGenes[file[:file.rfind('.')]] = kmercounttodict(file,kmercountPath) #.split('.')[0]

    # we now have two dictionaries (or more if we add more than two kmer count files to the kmercount_files path
    # now compare the two dictionaries in both directions to find kmers that are high in kmer dict 1 and low in kmer dict 2 and vice versa
    # this gets the dictionary name for the first kmer dict
    #kmerDicts = [dictOfGenes[kmercountFiles[i][:kmercountFiles[i].rfind('.')]] for i in range(len(kmercountFiles))]
    #dict1 = dictOfGenes[kmercountFiles[0][:kmercountFiles[0].rfind('.')]]#.split('.')[0]]
    #dict2 = dictOfGenes[kmercountFiles[1][:kmercountFiles[1].rfind('.')]]

    # output kmers and counts for differential kmers
    # output file names
    outFileNames = defaultdict(list)
    for file in kmercountFiles:
        outFileNames[kmercountPath + "/%s.higher.kmers.fa" % (file.split('.')[0])] = dictOfGenes[file[:file.rfind('.')]]
        #outFileNames.append(kmercountPath + "/%s.higher.kmers.fa" % (file.split('.')[0]))
    # {file.higherkmer : {kmer:count}}
    # create files for writing
    for filename in outFileNames:
        open(filename, 'w').close()
        print 'creating %s' % filename

    for outfilename, dict1 in outFileNames.iteritems():
        # check dict 1 against dict 2
        out1 = open(outfilename, 'w')
        # iterate through the keys of dict1 and identify kmers that are at least 10 fold higher in dict1
        for key, value in dict1.iteritems():
            val1 = value
            values = []
            for outfilename2 in outFileNames:
                if outfilename2 != outfilename:
                    values.append(outFileNames[outfilename2].get(key,3))
            # require at least 30 fold higher kmers in dict1
            if any([(val1 / val2) > ratio_threshold for val2 in values]): # all()
                out1.write('>%s.%d.%s\n%s\n' % (key, val1, '.'.join(map(str,values)), key))
        out1.close()
    """
    # do same for other direction of query, # check dict 2 against dict 1
    out2 = open(outFileNames[1], 'w')
    for key, value in dict2.iteritems():
        key1 = str(key)
        val1 = value
        val2 = dict1.get(key, 3)
        if (val1 / val2) > ratio_threshold:
            out2.write('>%s.%d.%d\n%s\n' % (key, val1, val2, key))
    out2.close()"""

def blast2bed3(subgenomeFolder,blastPath, bedPath, sortPath, genome,BB):
    """Takes a list of genotype files with only one column for pos and converts them to proper bedgraph format to be sorted"""
    print 'blast files contains'
    blastFiles = os.listdir(blastPath)
    print('\n'.join('{}: {}'.format(*k) for k in enumerate(blastFiles)))
    if BB:
        endname = '.sam'
    else:
        endname = '.BLAST.tsv.txt'
    for blastFile in blastFiles:
        if blastFile.endswith(endname):
            f = blastFile.rstrip()
            outFileName = f[:f.rfind('.')]+'.bed3'
            input_list = [blastPath, f]
            inpath = '%s/%s' % tuple(input_list)
            inputFile = open(inpath, 'r')
            outpath = os.path.join(bedPath, outFileName)
            bo = open(outpath, 'w')
            if BB:
                for line in inputFile:
                    lineInList = line.split()
                    lineOutputList = [lineInList[2],int(lineInList[3]),int(lineInList[3])+1]
                    bo.write('%s\t%d\t%d\n' % tuple(lineOutputList))
                    """blast:ATATGTTGTAATATTTGAGCACT.322.13	Nt01_118425000_118500000	100.000	23	23	24631	24609	2.35e-04	46.1"""
                    """sam:ATATGTTGTAATATTTGAGCACT.322.13  16      Nt01_118425000_118500000        24609   3       23=     *       0       0       AGTGCTCAAATATTACAACATAT *       XT:A:R  NM:i:0  AM:i:3"""
            else:
                for line in inputFile:
                    lineInList = line.split()
                    lineOutputList = [lineInList[1], int(lineInList[8]), int(lineInList[8])+1]
                    bo.write('%s\t%d\t%d\n' % tuple(lineOutputList))
            inputFile.close()
            bo.close()
            sortedName = f[:f.rfind('.')] + '.sorted.bed3'
            si = os.path.join(bedPath, outFileName)
            so = os.path.join(sortPath, sortedName)
            coveragename = subgenomeFolder + '/' + f[:f.rfind('.')] + '.sorted.cov'
            if not os.path.exists(sortPath):
                os.makedirs(sortPath)
            b = BedTool(si)
            # if chromosomes, then use the 1Mb window
            if not os.path.exists(genome.replace('.fai','')+'.bed'):
                fai2bed((genome,))
            shutil.copy(genome.replace('.fai','')+'.bed',subgenomeFolder)
            windows = '%s.bed' % genome
            a = BedTool(windows)
            b.sort().saveas(so)
            a.coverage(b).saveas(coveragename)
            bedgname = f[:f.rfind('.')] + '.sorted.cov.bedgraph'
            open(subgenomeFolder + '/' + bedgname, 'w').close()
            bedgo = open(subgenomeFolder + '/' + bedgname, 'w')
            covFile = open(coveragename, 'r')
            for line in covFile:
                lineInList = line.split()
                lineOutputList = [lineInList[0], int(lineInList[1]), int(lineInList[2]), int(lineInList[5]) ]#int(lineInList[3]) ]
                bedgo.write('%s\t%d\t%d\t%d\n' % tuple(lineOutputList))
            covFile.close()
            bedgo.close()

def bed2unionBed(genome, subgenomeFolder, bedPath):
    # for blastFile in blastFiles:
    #     f = blastFile.rstrip()
    #     outFileName = f[:f.rfind('.')]+'.bed3'
    bedGraphFiles = [file for file in os.listdir(subgenomeFolder) if file.endswith('.sorted.cov.bedgraph')]
    #bedFile = blastFiles[0]
    #f = bedFile.rstrip()
    #print f
    inputName = 'subgenomes'#f[:f.rfind('.')]
    outputFileName = subgenomeFolder + '/' +inputName + '.union.bedgraph'
    genomeName = genome[:genome.rfind('.')]
    subprocess.call('cut -f 1-2 %s.fai > %s.genome'%(genome,subgenomeFolder + '/genome'),shell=True)
    genomeFile = subgenomeFolder + '/genome.genome'
    #a_blast = blastFiles[0]
    #b_blast = blastFiles[1]
    #a_name = a_blast[:a_blast.rfind('.')] + '.sorted.cov.bedgraph'
    #b_name = b_blast[:b_blast.rfind('.')] + '.sorted.cov.bedgraph'
    #a = pybedtools.BedTool(a_name)
    #a_sort = a.sort()
    #b = pybedtools.BedTool(b_name)
    #b.sort()
    #b_sort = b.sort()
    bedGraphBedSortFn = [BedTool(subgenomeFolder+'/' + file).sort().fn for file in bedGraphFiles]
    x = BedTool()
    result = x.union_bedgraphs(i=bedGraphBedSortFn, g=genomeFile, empty=True)
    result.saveas(outputFileName)

def make_plots(genome, bedFiles):
    f = bedFiles[0]
    inputName = f[:f.rfind('.')]
    outputFileName = inputName + '.union.bedgraph'
    lineOutputList = [outputFileName, genome]
    dbscriptName = outputFileName + '.plot.sh'
    c = open(dbscriptName, 'w')
    c.write('module load bedtools/2.25.0 && python ./GenomeKmerComp/multiplot_gRNA_chrIter_2colBedGraphv2.py %s %s\n' % tuple(lineOutputList))
    c.close()
    try:
        subprocess.call('sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName

def kmerratio2scaffasta(subgenomePath, originalSubgenomePath, fastaPath, genomeName, originalGenome, BB, bootstrap, iteration, kmerLength,runFinal,kmer500Path, metric, originalStr, blastMem, kmer_low_count, diff_kmer_threshold, unionbed_threshold):
    try:
        original = int(originalStr)
    except:
        original = 0
    try:
        absolute_threshold, ratio_threshold = tuple(map(int,unionbed_threshold.split(',')))
    except:
        absolute_threshold, ratio_threshold = 10, 2
    if original:
        genomeName = originalGenome
    genome = fastaPath + genomeName
    blastMemStr = "export _JAVA_OPTIONS='-Xms5G -Xmx%sG'" % (blastMem)
    a = subgenomePath + '/subgenomes.union.bedgraph'#[file for file in os.listdir(subgenomePath) if 'union.bedgraph' in file][0]
    ubedg = open(a, 'r')
    genomeFastaObj = Fasta(genome)
    extractPath = subgenomePath+'/extractedSubgenomes/'
    try:
        os.mkdir(extractPath)
    except:
        pass
    ### define output filenames
    genomeprefix = genome[genome.rfind('/')+1:genome.rfind('.')]
    outputSubgenomes = [extractPath + genomeprefix + '.subgenome' + chr(i+65) + '.fasta' for i in range(len(ubedg.readline().split('\t')[3:]))]
    ubedg.seek(0)
    scaffoldsOut = [[] for subgenome in outputSubgenomes]#{subgenome:[] for subgenome in outputSubgenomes}
    ambiguousScaffolds = []

    # parse the unionbed file to subset
    for line in ubedg:
        if line:#'Chr' in line or 'chr' in line or 'B' in line or 'Sc' in line:
            # print line
            lineList = line.split('\t')
            scaff = str((lineList[0]).rstrip())
            ambiguous = 1
            x = [float((lineList[i]).rstrip()) for i in range(3,len(lineList))]
            for i in range(len(x)):
                x_others = x[:i] + x[i+1:] #FIXME maybe add a break below, or find max of list to speed up algorithm, not a big deal
                if all([(x_i == 0 and x[i] > absolute_threshold) or (x_i > 0 and (x[i]/x_i) > ratio_threshold) for x_i in x_others]):
                    scaffoldsOut[i].append(scaff)
                    ambiguous = 0
            if ambiguous:
                ambiguousScaffolds.append(scaff)
    no_kill = all([len(subgenome) > 0 for subgenome in scaffoldsOut])
    ubedg.close()
    for subgenome, scaffolds in zip(outputSubgenomes,scaffoldsOut):
        with open(subgenome,'w') as f:
            for scaff in scaffolds:
                f.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
        subprocess.call(blastMemStr + ' && reformat.sh in=%s out=%s fastawrap=60'%(subgenome,subgenome.replace('.fasta','_wrapped.fasta')),shell=True)
    with open(extractPath+'ambiguousScaffolds.fasta','w') as f:
        for scaff in ambiguousScaffolds:
            f.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
    subprocess.call(blastMemStr + ' && reformat.sh in=%s out=%s fastawrap=60' % (extractPath+'ambiguousScaffolds.fasta', extractPath+'ambiguousScaffolds_wrapped.fasta'), shell=True)

    if iteration == 0:
        subgenomePath = subgenomePath + '/bootstrap_1'
    if iteration < bootstrap:
        iteration += 1
        subgenomePath = subgenomePath[:subgenomePath.rfind('/')] + '/bootstrap_%d' % (iteration)
        try:
            os.mkdir(subgenomePath)
        except:
            pass
        for i, scaffolds in enumerate(scaffoldsOut):
            with open(subgenomePath + '/subgenome_%d.txt'%i,'w') as f:
                f.write('\n'.join(scaffolds))
    elif iteration == bootstrap:
        iteration += 1
        try:
            os.mkdir(originalSubgenomePath + '/classify')
        except:
            pass
        for i, scaffolds in enumerate(scaffoldsOut):
            with open(originalSubgenomePath + '/classify/subgenome_%d.txt' % i, 'w') as f:
                f.write('\n'.join(scaffolds))
    if no_kill == 0:
        return
    if runFinal == 0:
        subgenomeExtraction((subgenomePath, originalSubgenomePath, fastaPath, genomeName, originalGenome, BB, bootstrap, iteration, kmerLength, runFinal,kmer500Path,metric, originalStr, blastMem, kmer_low_count, diff_kmer_threshold, unionbed_threshold))
    else:
        return
    """
            x3 = float((line.split()[3]).rstrip())
            #print 'my x3 is %s\n' % x3
            x4 = float((line.split()[4]).rstrip())
            #print 'my x4 is %s\n' % x4
            if x4 > 0 and (x3/x4) > 2:
                o1h.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
            elif x3 > 0 and (x4/x3) > 2:
                o2h.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
            elif x4 == 0 and x3 > 10:
                o1h.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
            elif x3 == 0 and x4 > 10:
                o2h.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
    o1h.close()
    o2h.close()"""


def subgenomeExtraction(args):
    subgenome_folder, originalSubgenomePath, fastaPath, genomeName, originalGenome, BB, bootstrap, iteration, kmerLength, runFinal,kmer500Path, metric, originalStr, blastMem, kmer_low_count, diff_kmer_threshold, unionbed_threshold = args
    bedDict = fai2bed((fastaPath+genomeName,))
    try:
        bootstrap = int(bootstrap)
        iteration = int(iteration)
        runFinal = int(runFinal)
    except:
        bootstrap = 0
        iteration = 0
        runFinal = 0
    try:
        BB = int(BB)
    except:
        BB = 0
    blastPath = subgenome_folder + '/blast_files/'
    bedPath = subgenome_folder + '/bed_files/'
    sortPath = subgenome_folder + '/sortedbed_files/'
    try:
        os.mkdir(blastPath)
    except:
        pass
    try:
        os.mkdir(bedPath)
    except:
        pass
    try:
        os.mkdir(sortPath)
    except:
        pass
    if bootstrap >= iteration or runFinal == 1:
        for file in os.listdir(subgenome_folder):
            if file and file.endswith('.txt'):
                with open(subgenome_folder+'/'+file,'r') as f,open(subgenome_folder+'/'+file.replace('.txt','.bed'),'w') as f2:
                    for line in f:
                        if line:
                            f2.write(bedDict[line.strip('\n')][0])
                subprocess.call('bedtools getfasta -fi %s -fo %s -bed %s -name'%(fastaPath+genomeName,subgenome_folder + '/%s_'%('model')+file.replace('.txt','.fa'),subgenome_folder+'/'+file.replace('.txt','.bed')),shell=True)
        writeKmerCountSubgenome((subgenome_folder,kmerLength,blastMem,diff_kmer_threshold))
        writeBlast((originalGenome,blastPath,subgenome_folder+'/kmercount_files/',fastaPath,BB,blastMem))
        blast2bed3(subgenome_folder, blastPath, bedPath, sortPath, fastaPath+originalGenome, BB)
        bed2unionBed(fastaPath+originalGenome, subgenome_folder, bedPath)
        kmerratio2scaffasta(subgenome_folder, originalSubgenomePath, fastaPath, genomeName, originalGenome, BB, bootstrap, iteration, kmerLength,runFinal,kmer500Path,metric,originalStr, blastMem, kmer_low_count, diff_kmer_threshold, unionbed_threshold)
    if bootstrap < iteration and runFinal == 0:
        subgenome_folder = originalSubgenomePath
        if 0: #FIXME classifier has been removed because it has been crashing... Will not be implemented in v1.0, needs work
            classifyFolder = subgenome_folder + '/classify/'
            subgenome_folder, runFinal = classify(classifyFolder,fastaPath,genomeName,kmerLength, originalSubgenomePath.split('/')[-1],kmer500Path,metric,originalStr,originalGenome, blastMem, kmer_low_count)
        subgenome_folder, runFinal = 'null', 0
        iteration += 1
        if runFinal:
            subgenomeExtraction((subgenome_folder, originalSubgenomePath, fastaPath, genomeName, originalGenome, BB, bootstrap, iteration, kmerLength,runFinal,kmer500Path,metric,originalStr, blastMem, kmer_low_count, diff_kmer_threshold, unionbed_threshold))

def classify(classifyFolder, fastaPath, genomeName, kmerLength,model,kmer500Path,metric,originalStr,originalGenome,blastMem,kmer_low_count):
    subgenome_files = [np.vectorize(lambda line: line.strip('\n'))(open(classifyFolder + file, 'r').readlines()) for file in os.listdir(classifyFolder) if file.endswith('.txt') and os.stat(classifyFolder+file).st_size]
    blastMemStr = "export _JAVA_OPTIONS='-Xms5G -Xmx%sG'" %(blastMem)
    try:
        kmer_low_count = int(kmer_low_count)
    except:
        kmer_low_count = 100
    if len(subgenome_files) > 1:
        # FIXME maybe filter out missing sequences...
        scaffoldLabel_dict = defaultdict(lambda: 0)
        for i in range(len(subgenome_files)):
            for scaffold in subgenome_files[i]:
                scaffoldLabel_dict[scaffold] = i
        n_subgenomes = len(subgenome_files)
        total_subgenome_scaffolds = np.concatenate(tuple(subgenome_files),axis=0)
        subgenomeFolder=classifyFolder + 'finalClassifiedOutputs'
        try:
            os.mkdir(subgenomeFolder)
        except:
            pass
        try:
            original = int(originalStr)
        except:
            original = 0
        if original:
            with open(fastaPath+originalGenome+'.fai', 'r') as f:
                scaffolds = np.array([line.split()[0] for line in f if line])
            genomeName = originalGenome
        else:
            with open(fastaPath+genomeName+'.fai', 'r') as f:
                scaffolds = np.array([line.split()[0] for line in f if line])
        print fastaPath + genomeName
        print classifyFolder + 'ambiguous.kcount.fa'
        print classifyFolder + 'ambiguous.sam'
        print blastMemStr
        print blastMemStr + ' && bbmap.sh vslow=t ambiguous=all noheader=t secondary=t perfectmode=t threads=8 maxsites=2000000000 outputunmapped=f ref=%s path=%s/ in=%s outm=%s' % (
                fastaPath + genomeName, classifyFolder, classifyFolder + 'ambiguous.kcount.fa',
                classifyFolder + 'ambiguous.sam')
            #scaffolds = np.array(pickle.load(open('scaffolds.p', 'rb')))
        runFinal = 0
        scaffolds_unchecked = np.setdiff1d(scaffolds, total_subgenome_scaffolds)
        print 'scaffolds_unchecked', scaffolds_unchecked
        if list(scaffolds_unchecked):
            kmerLengths = kmerLength.split(',')
            genomeFastaObj = Fasta(fastaPath + genomeName)
            with open(classifyFolder + 'ambiguous.fa', 'w') as f:
                for scaff in total_subgenome_scaffolds:
                    f.write('>%s\n%s\n' % (scaff, str(genomeFastaObj[scaff][:])))
            subprocess.call(blastMemStr + ' && reformat.sh in=%s out=%s fastawrap=60' % (classifyFolder + 'ambiguous.fa', classifyFolder + 'ambiguous_wrapped.fa'), shell=True)
            kmerFiles = []
            open(classifyFolder + 'ambiguous.kcount','w').close()
            for kmerL in kmerLengths:
                kmerFiles.append(classifyFolder + 'ambiguous%s.kcount'%kmerL)
                subprocess.call(blastMemStr + ' && kmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s out=%s k=%s -Xmx60g' % (classifyFolder + 'ambiguous_wrapped.fa', classifyFolder + 'ambiguous%s.kcount'%kmerL,kmerL), shell=True)
                with open(classifyFolder + 'ambiguous%s.kcount'%kmerL,'r') as f1, open(classifyFolder + 'ambiguous.kcount','a') as f2:
                    f2.write(f1.read()+'\n')
            for fname in kmerFiles:
                os.remove(fname)
            with open(classifyFolder+'ambiguous.kcount','r') as f, open(classifyFolder+'ambiguous.kcount.fa','w') as f2:#kmer2Fasta((classifyFolder)) #ambiguous.kcount.fa
                for line in f:
                    if line and int(line.split('\t')[-1]) >= kmer_low_count:
                        f2.write('>%s\n%s\n' % tuple([line.split('\t')[0]] * 2))
            #with open(classifyFolder+'test.txt','w') as f:
            #    f.write('\n'.join([]))
            subprocess.call(blastMemStr + ' && bbmap.sh vslow=t ambiguous=all noheader=t secondary=t perfectmode=t threads=8 maxsites=2000000000 outputunmapped=f ref=%s nodisk in=%s outm=%s' % (fastaPath + genomeName, classifyFolder + 'ambiguous.kcount.fa', classifyFolder + 'ambiguous.sam'), shell=True)
            #print 'ALLCLEAR'
            kmerIdx = {line[1].split('\t')[0]: line[0] for line in
                       enumerate(open(classifyFolder + 'ambiguous.kcount', 'r').readlines())}
            scaffoldIdx = {scaffold[1]: scaffold[0] for scaffold in enumerate(scaffolds)}
            data = sps.dok_matrix((len(scaffolds), len(kmerIdx.values())), dtype=np.float32)
            scaffoldL = np.array([map(int,scaffold.split('_')[-2:]) for scaffold in scaffolds])
            scaffoldLengths = abs(scaffoldL[:,1]-scaffoldL[:,0])*5000.
            with open(classifyFolder + 'ambiguous.sam', 'r') as f:#, open(classifyFolder + 'blasted.bed', 'w') as f2:
                for line in f:
                    if line:
                        data[scaffoldIdx[line.split('\t')[2].split('::')[0]],kmerIdx[line.split('\t')[0]]] += 1.
                        #f2.write('\t'.join([l1] + ['0', str(int(l1.split('_')[-1]) - int(l1.split('_')[-2]))] + [
                        #    line.split('\t')[0]]) + '\n')
            # FIXME pca only for now, maybe use supervised LDA in future!!!!!!!!!
            data = data.tocsc()
            # divide every row by scaffold length
            for i in range(len(scaffoldLengths)):
                data[i,:]/=scaffoldLengths[i]
            #intermediateTransform = KernelPCA(n_components=int(len(kmerIdx.keys())/10)).fit_transform(StandardScaler(with_mean=False).fit_transform(data))
            #lda = LDA(n_components=n_subgenomes + 1)
            scaffBool = np.vectorize(lambda scaffold: scaffold in total_subgenome_scaffolds)(scaffolds)
            trainLabels = np.vectorize(lambda scaffold: scaffoldLabel_dict[scaffold])(total_subgenome_scaffolds)
            #lda.fit(intermediateTransform[scaffBool],trainLabels)
            # FIXME change kernel type
            transformed_data = StandardScaler().fit_transform(KernelPCA(n_components=n_subgenomes+1,kernel=metric).fit_transform(StandardScaler(with_mean=False).fit_transform(data)))#lda.transform(intermediateTransform)nt(len(kmerIdx.keys())/10)
            trainData = transformed_data[scaffBool]
            testData = transformed_data[np.vectorize(lambda scaffold: scaffold in scaffolds_unchecked)(scaffolds)]
            knn = KNeighborsClassifier()
            knn.fit(trainData,trainLabels)
            testLabels = knn.predict(testData)
            kbest = SelectKBest(chi2, k=500)
            kbest.fit(data[scaffBool], trainLabels)
            bestFeatures = kbest.pvalues_.argsort()[:500]
            kmers = kmerIdx.keys()
            best_500_kmers = [kmers[i] for i in bestFeatures]
            with open('%s/kmer500Best_%s.fa' % (kmer500Path, model+'_preClassify'), 'w') as f:
                f.write('\n'.join('>%s\n%s' % (kmer, kmer) for kmer in best_500_kmers))
            for i in range(len(subgenome_files)):
                subgenome_files[i] = np.concatenate((subgenome_files[i],scaffolds_unchecked[testLabels == i]))
                with open(subgenomeFolder + '/subgenome_%d.txt'%i,'w') as f:
                    f.write('\n'.join(subgenome_files[i]))
            runFinal = 1
        else:
            #print 'heyHEY'
            return 'null', 0
    else:
        #print 'HOHO'
        return 'null', 0
    #print 'HEHE'
    return subgenomeFolder, runFinal

def generateKmerGraph(args):
    kmerPath, kmerName, n_subgenomes, BB = args
    try:
        n_subgenomes = int(n_subgenomes)
    except:
        n_subgenomes = 2
    try:
        BB = int(BB)
    except:
        BB = 0
    if BB:
        k = 2
    else:
        k = 1
    blastPath = kmerPath + '/' + kmerName + '/'
    with open(blastPath+ kmerName + '.blast.txt', 'r') as f, open(blastPath + 'blasted.bed', 'w') as f2:
        for line in f:
            if line:
                l1 = line.split('\t')[k].split('::')[0]
                f2.write('\t'.join([l1] + ['0', str(int(l1.split('_')[-1]) - int(l1.split('_')[-2]))] + [
                    line.split('\t')[0]]) + '\n')

    b = pybedtools.BedTool(blastPath + 'blasted.bed').sort().merge(c=4, o='collapse', )
    b.saveas(blastPath + 'blasted_merged.bed')
    kmerGraph = nx.Graph()
    totalKmerCount = {}
    totalKmerCount = defaultdict(lambda:0, totalKmerCount)
    with open(blastPath + 'blasted_merged.bed', 'r') as f:
        for line in f:
            if line:
                #kmerCountsNumbers = [kmerIdx[kmer] for kmer in set(line.rstrip('\n').split('\t').split(','))]
                kmerz = line.rstrip('\n').split('\t')[-1].split(',')
                kmerCount = Counter(kmerz)
                for key in kmerCount:
                    totalKmerCount[key] += kmerCount[key]
                kmerGraph.add_edges_from(list(combinations(kmerCount.keys(),2)))#kmerCountsNumbers, 2)))

    histData = np.array([kmerGraph.degree(kmer) for kmer in kmerGraph.nodes()])[:, np.newaxis]
    xplot = np.linspace(0, np.max(histData), 300000)[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=np.max(histData)/18.).fit(histData)
    exp_log_dens = np.exp(kde.score_samples(xplot))

    idxs_peaks = argrelextrema(exp_log_dens, np.greater)[0]  # peakutils.indexes(exp_log_dens)
    idxs_valleys = argrelextrema(exp_log_dens, np.less)[0]  # peakutils.indexes(-exp_log_dens)

    plt.figure()
    plt.plot(xplot[:], exp_log_dens, '-',
             label="Envelope Kernel '{0}' Density Function".format('Gaussian'))
    plt.plot(xplot[idxs_peaks], exp_log_dens[idxs_peaks], '*', color='r', label='Peaks')
    plt.plot(xplot[idxs_valleys], exp_log_dens[idxs_valleys], 'o', color='r', label='Valleys')
    plt.hist(histData[:, 0], bins=50, label='Histogram of Counts', normed=True)
    plt.xlabel('Number of Related Kmers (Normalized)')
    plt.ylabel('Count')
    plt.legend(loc='upper left')
    plt.title('Histogram of Number of Related Kmers (Normalized)')
    plt.savefig(blastPath + 'KmerHistogram.png')

    plt.figure()
    plt.hist(totalKmerCount.values(),bins = 50)
    plt.xlabel('Total Kmer Count')
    plt.ylabel('Count')
    plt.legend(loc='upper left')
    plt.title('Histogram of Total Kmer Count Throughout Genome')
    plt.savefig(blastPath + 'TotalKmerCountHistogram.png')


    #2d graph
    plt.figure()
    plt.axis('off')
    nx.draw_networkx(kmerGraph, pos=nx.spring_layout(kmerGraph), edge_color='b', with_labels=0)#nodecolor='r'
    plt.savefig(blastPath + 'kmerGraph_%s.png'%(kmerName), bbox_inches="tight")

    #3d graph
    plt.figure()
    pos = nx.spring_layout(kmerGraph,dim=3)#,k=0.7,iterations=50,scale=5,dim=3)
    Xed = []
    Yed = []
    Zed = []
    for edge in kmerGraph.edges():
        Xed += [pos[edge[0]][0], pos[edge[1]][0], None]
        Yed += [pos[edge[0]][1], pos[edge[1]][1], None]
        Zed += [pos[edge[0]][2], pos[edge[1]][2], None]
    plots = []
    if len(idxs_valleys) > 0:
        intervals = [(0., xplot[idxs_valleys[0]][0])] + [(xplot[idxs_valleys[i]][0], xplot[idxs_valleys[i + 1]][0]) for i in
                                                         range(len(idxs_valleys) - 1)] + [
                        (xplot[idxs_valleys[-1]][0], xplot[-1][0])]
    else:
        intervals = [(0.,np.max(histData))]
    N = len(intervals)
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N+1)]

    kmerCountData = np.array([[kmer, int(kmerGraph.degree(kmer))] for kmer in kmerGraph.nodes()])
    print intervals
    allCounts = np.vectorize(int)(kmerCountData[:,1])
    #idxs = np.arange(len(kmerGraph.nodes()))
    #print kmerCountData
    for i,interval in enumerate(intervals):
        try:
            #print np.where((kmerCountData[:, 1] <= np.floor(interval[1])) & (kmerCountData[:, 1] >= np.ceil(interval[0])))
            idx_slice = np.where(( allCounts <= np.floor(interval[1])) & (allCounts >= np.ceil(interval[0])))
            nodesData = kmerCountData[:,0][idx_slice]
            counts = allCounts[idx_slice]
            nodesData = nodesData[np.argsort(counts)]
            counts = counts[np.argsort(counts)]
            print nodesData
            nodes = list(nodesData)
            G = kmerGraph.copy()
            G.remove_nodes_from(nodes)
            plt.figure()
            plt.axis('off')
            nx.draw_networkx(G, pos=nx.spring_layout(G), edge_color='b',with_labels=False)  # nodecolor='r'
            plt.savefig(blastPath + 'kmerGraph_Interval_%d_%d.png' % (int(np.ceil(interval[0])),int(np.floor(interval[1]))), bbox_inches="tight")
            with open(blastPath + 'kmers_Interval_%d_%d.txt' % (int(np.ceil(interval[0])),int(np.floor(interval[1]))),'w') as f:
                f.write('\n'.join('%s\t%d'%(nodesData[j],counts[j]) for j in range(len(nodesData))))# nodesData[np.argsort(nodesData[:,1])]))
            Xv = [pos[k][0] for k in nodes]
            Yv = [pos[k][1] for k in nodes]
            Zv = [pos[k][2] for k in nodes]
            nodesText = ['%s, %d connections, %d total kmer count'%(kmer, int(kmerGraph.degree(kmer)), totalKmerCount[kmer]) for kmer in nodes]
            plots.append(go.Scatter3d(x=Xv,
                                  y=Yv,
                                  z=Zv,
                                  mode='markers',
                                  name='Peak Interval %d %d'%(int(np.ceil(interval[0])),int(np.floor(interval[1]))),
                                  marker=go.Marker(symbol='dot',
                                                   size=5,
                                                   color=c[i],
                                                   line=go.Line(color='rgb(50,50,50)', width=0.5)
                                                   ),
                                  text=nodesText,
                                  hoverinfo='text'
                                  ))
        except OSError as e:
            print e.errno
            print e.filename
            print e.strerror

    plots.append(go.Scatter3d(x=Xed,
                              y=Yed,
                              z=Zed,
                              mode='lines',
                              line=go.Line(color='rgb(210,210,210)', width=1),
                              hoverinfo='none'
                              ))

    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''
                )

    layout = go.Layout(
        title="Graph of 500 Best Kmers",
        width=1000,
        height=1000,
        showlegend=True,
        scene=go.Scene(
            xaxis=go.XAxis(axis),
            yaxis=go.YAxis(axis),
            zaxis=go.ZAxis(axis),
        ),
        margin=go.Margin(
            t=100
        ),
        hovermode='closest',
        annotations=go.Annotations([
            go.Annotation(
                showarrow=False,
                text="",
                xref='paper',
                yref='paper',
                x=0,
                y=0.1,
                xanchor='left',
                yanchor='bottom',
                font=go.Font(
                    size=14
                )
            )
        ]), )

    data1 = go.Data(plots)
    fig1 = go.Figure(data=data1, layout=layout)
    py.plot(fig1, filename=blastPath + 'KmerPeakGraph3D.html')

    affinity_matrix = nx.to_numpy_matrix(kmerGraph)
    s = SpectralClustering(n_clusters=n_subgenomes, affinity='precomputed')
    labels = s.fit_predict(affinity_matrix)
    plots = []

    plots.append(go.Scatter3d(x=Xed,
                              y=Yed,
                              z=Zed,
                              mode='lines',
                              line=go.Line(color='rgb(210,210,210)', width=1),
                              hoverinfo='none'
                              ))

    N = len(set(labels))
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N + 1)]
    for i,label in enumerate(list(set(labels))):
        nodes = list(kmerCountData[:,0][labels == label])
        counts = list(allCounts[labels == label])
        Xv = [pos[k][0] for k in nodes]
        Yv = [pos[k][1] for k in nodes]
        Zv = [pos[k][2] for k in nodes]
        nodesText = [
            '%s, %d connections, %d total kmer count' % (kmer, int(kmerGraph.degree(kmer)), totalKmerCount[kmer]) for
            kmer in nodes]
        plots.append(go.Scatter3d(x=Xv,
                                  y=Yv,
                                  z=Zv,
                                  mode='markers',
                                  name='Cluster %d' %i,
                                  marker=go.Marker(symbol='dot',
                                                   size=5,
                                                   color=c[i],
                                                   line=go.Line(color='rgb(50,50,50)', width=0.5)
                                                   ),
                                  text=nodesText,
                                  hoverinfo='text'
                                  ))
        with open(blastPath + 'kmers_cluster_%d.txt' %(i),
                  'w') as f:
            f.write('\n'.join('%s\t%d' % (nodes[j], counts[j]) for j in range(len(nodes))))
    data1 = go.Data(plots)
    fig1 = go.Figure(data=data1, layout=layout)
    py.plot(fig1, filename=blastPath + 'KmerClusteredGraph3D.html')

def clusterGraph(args): #FIXME under development
    helpStr = """python subgenomeClusteringInterface.py clusterGraph graph_Matrix_path[.npz] connectedScaffoldsPath[.p] outputDir initialPosition=random|spectral|path_to_.npy iterations features_bed=Path_to_Bed\npython subgenomeClusteringInterface.py clusterGraph help/-h"""
    try:
        if args[0] == 'help' or args[0] == '-h':
            print helpStr
            quit()
    except:
        print helpStr
        quit()
    try:
        sparse_matrix_file, scaffoldsFile, outDir, initialPos, iteration, bedFeatures = args
    except:
        print helpStr
        quit()
    try:
        bedFeaturesFile = bedFeatures.split('=')[1]
        if bedFeaturesFile.endswith('.bed') == 1:
            featureMap = 1
        else:
            featureMap = 0
    except:
        featureMap = 0
    G = nx.from_scipy_sparse_matrix(sps.load_npz(sparse_matrix_file))
    scaffolds = pickle.load(open(scaffoldsFile, 'rb'))
    N = 2
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N + 1)]
    mapping = {i:scaffolds[i] for i in range(len(scaffolds))}
    G=nx.relabel_nodes(G,mapping, copy=False)
    nodes = G.nodes()
    if featureMap:
        scaffoldsDict = {scaffold : '\t'.join(['_'.join(scaffold.split('_')[0:-2])]+scaffold.split('_')[-2:]) for scaffold in scaffolds}
        outputFeatures = defaultdict(list)
        scaffoldIdx = {scaffolds[i] : i for i in range(len(scaffolds))}
        scaffoldsBed = BedTool('\n'.join(scaffoldsDict.values()),from_string=True)
        featureBed = BedTool(bedFeaturesFile)
        featuresDict = {scaffold: '' for scaffold in scaffolds}
        finalBed = scaffoldsBed.intersect(featureBed,wa=True,wb=True).sort().merge(d=-1,c=7,o='distinct')
        #print finalBed.head()
        finalBed.saveas(outDir+'/finalBed.bed')
        omittedRegions = scaffoldsBed.intersect(featureBed,v=True,wa=True)
        omittedRegions.saveas(outDir+'/ommitted.bed')
        #print omittedRegions.head()
        for line in str(finalBed).splitlines()+[line2+'\tunlabelled' for line2 in str(omittedRegions).splitlines()]:
            lineList = line.strip('\n').split('\t')
            feature = lineList[-1]
            scaffold = '_'.join(lineList[0:-1])
            if ',' in feature:
                featuresDict[scaffold] = '|'.join(feature.split(','))#[ft for ft in feature.split(',')])#.split('_')[0]
                feature = 'ambiguous'
            else:
                featuresDict[scaffold] = feature#.split('_')[0]
            #idx = scaffoldIdx[scaffold]
            outputFeatures[scaffold] = feature#tuple(feature)#.split('_'))
        mainFeatures = set(outputFeatures.values())
        N = len(mainFeatures)
        c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N + 1)]
        #print c
        featuresColors = {feature : c[i] for i, feature in enumerate(mainFeatures)}
        #print featuresColors
        outputFeaturesArray = np.array([outputFeatures[scaffold] for scaffold in nodes])
        names = np.vectorize(lambda name: 'Scaffolds: ' + name)(outputFeaturesArray)#[:,0]
        #print names
        colors = np.vectorize(lambda feature: featuresColors[feature])(outputFeaturesArray)#[:,1]
        nodesText = np.array(['%s, %d connections, feature= %s' % (scaffold, int(G.degree(scaffold)),featuresDict[scaffold]) for scaffold in nodes])  # , Related: %s'%(kmer, int(G.degree(kmer)), ' '.join(G[kmer].keys())) for kmer in nodes]
        print np.column_stack((scaffolds,outputFeaturesArray,colors,nodesText,names))
        #print [names[names==feature] for feature in mainFeatures]
        #print { feature : np.vectorize(lambda feature: featuresColors[feature])(list(set(names[names==feature]))) for feature in mainFeatures}
    else:
        names = 'Scaffolds'
        colors = c[0]
        nodesText = ['%s, %d connections, ' % (scaffold, int(G.degree(scaffold))) for scaffold in nodes]  # , Related: %s'%(kmer, int(G.degree(kmer)), ' '.join(G[kmer].keys())) for kmer in nodes]
    try:
        initialPos_name = initialPos.split('=')[1]
    except:
        print 'Invalid initial position'
        print helpStr
        quit()
    if initialPos_name == 'spectral':
        pos_i = nx.spectral_layout(G,dim=3)
    elif initialPos_name == 'random':
        pos_i = nx.random_layout(G, dim=3)
    elif initialPos_name.endswith('.npy'):
        pos_i = defaultdict(list)
        t_data = np.load(initialPos_name)
        for i in range(len(scaffolds)):
            pos_i[scaffolds[i]] = tuple(t_data[i,:])
    else:
        print 'Invalid initial position name: ' + initialPos_name
        print helpStr
        quit()
    try:
        iterations = map(int,iteration.split(','))
    except:
        print 'Proper iterations input is comma delimitted'
        print helpStr
        quit()
    masterData = [] #FIXME get slider to work
    for idx,i in enumerate(iterations):
        if i != 0: #FIXME change computation route, takes too much time to calculate all of those iterations, maybe do 1 iteration at a time
            pos = nx.spring_layout(G,dim=3,iterations=i,pos=pos_i)
        else:
            pos = pos_i
        plots = []
        Xv = np.array([pos[k][0] for k in nodes])
        Yv = np.array([pos[k][1] for k in nodes])
        Zv = np.array([pos[k][2] for k in nodes])
        Xed = []
        Yed = []
        Zed = []
        for edge in G.edges():
            Xed += [pos[edge[0]][0], pos[edge[1]][0], None]
            Yed += [pos[edge[0]][1], pos[edge[1]][1], None]
            Zed += [pos[edge[0]][2], pos[edge[1]][2], None]
        #print names
        if featureMap:
            for name in mainFeatures:
                plots.append(go.Scatter3d(x=Xv[outputFeaturesArray == name],
                                      y=Yv[outputFeaturesArray == name],
                                      z=Zv[outputFeaturesArray == name],
                                      mode='markers',
                                      name= name,
                                      marker=go.Marker(symbol='dot',
                                                       size=5,
                                                       color=featuresColors[name],
                                                       line=go.Line(color='rgb(50,50,50)', width=0.5)
                                                       ),
                                      text=nodesText[outputFeaturesArray == name],
                                      hoverinfo='text'
                                      ))
        else:
            plots.append(go.Scatter3d(x=Xv,
                                      y=Yv,
                                      z=Zv,
                                      mode='markers',
                                      name=names,
                                      marker=go.Marker(symbol='dot',
                                                       size=5,
                                                       color=colors,
                                                       line=go.Line(color='rgb(50,50,50)', width=0.5)
                                                       ),
                                      text=nodesText,
                                      hoverinfo='text'
                                      ))
        plots.append(go.Scatter3d(x=Xed,
                                  y=Yed,
                                  z=Zed,
                                  mode='lines',
                                  line=go.Line(color='rgb(210,210,210)', width=1),
                                  hoverinfo='none'
                                  ))

        if idx == 0:
            sliders_dict = {
                'active': 0,
                'yanchor': 'top',
                'xanchor': 'left',
                'currentvalue': {
                    'font': {'size': 20},
                    'prefix': 'Frame:',
                    'visible': True,
                    'xanchor': 'right'
                },
                'transition': {'duration': 300, 'easing': 'cubic-in-out'},
                'pad': {'b': 10, 't': 50},
                'len': 0.9,
                'x': 0.1,
                'y': 0,
                'steps': []
            }
        """
        layout = go.Layout(
            title="Graph of Scaffolds",
            width=1000,
            height=1000,
            showlegend=True,
            scene=go.Scene(
                xaxis=go.XAxis(axis),
                yaxis=go.YAxis(axis),
                zaxis=go.ZAxis(axis),
            ),
            margin=go.Margin(
                t=100
            ),
            hovermode='closest',
            annotations=go.Annotations([
                go.Annotation(
                    showarrow=False,
                    text="",
                    xref='paper',
                    yref='paper',
                    x=0,
                    y=0.1,
                    xanchor='left',
                    yanchor='bottom',
                    font=go.Font(
                        size=14
                    )
                )
            ]), )"""
        slider_step = {'args': [
            [str(i)],
            {'frame': {'duration': 300, 'redraw': False},
             'mode': 'immediate',
             'transition': {'duration': 300}}
            ],
            'label': str(i),
            'method': 'animate'}
        sliders_dict['steps'].append(slider_step)
        masterData.append({'data' : go.Data(plots),'name' : str(i)})#, 'layout': layout})
    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''
                )
    masterLayout = dict(
        title="Graph of Scaffolds",
        updatemenus=[{'direction': 'left',
                      'pad': {'r': 10, 't': 87},
                      'showactive': False,
                      'type': 'buttons',
                      'x': 0.1,
                      'xanchor': 'right',
                      'y': 0,
                      'yanchor': 'top', 'buttons': [
                {
                    'args': [None, {'frame': {'duration': 500, 'redraw': False},
                                    'fromcurrent': True,
                                    'transition': {'duration': 300, 'easing': 'quadratic-in-out'}}],
                    'label': 'Play',
                    'method': 'animate'
                },
                {
                    'args': [[None], {'frame': {'duration': 0, 'redraw': False}, 'mode': 'immediate',
                                      'transition': {'duration': 0}}],
                    'label': 'Pause',
                    'method': 'animate'
                }
            ]}],
        sliders=[sliders_dict],
        width=1000,
        height=1000,
        showlegend=True,
        scene=go.Scene(
            xaxis=go.XAxis(axis),
            yaxis=go.YAxis(axis),
            zaxis=go.ZAxis(axis),
        ),
        margin=go.Margin(
            t=100
        ),
        hovermode='closest',
        annotations=go.Annotations([
            go.Annotation(
                showarrow=False,
                text="",
                xref='paper',
                yref='paper',
                x=0,
                y=0.1,
                xanchor='left',
                yanchor='bottom',
                font=go.Font(
                    size=14
                )
            )
        ]), )
    """
    masterLayout['sliders']= {
            'args': [
                'transition', {
                    'duration': 400,
                    'easing': 'cubic-in-out'
                }
            ],
            'initialValue': str(iterations[0]),
            'plotlycommand': 'animate',
            'values': map(str,iterations),
            'visible': True
        }
    masterLayout['sliders'] = [sliders_dict]"""
    fig1 = go.Figure(data=masterData[0]['data'], layout=masterLayout, frames=masterData)
    py.plot(fig1, filename=outDir + '/OutputGraph_frames_%s.html'%iteration)



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
    'cluster': cluster,
    'transform_main': transform_main,
    'subgenomeExtraction': subgenomeExtraction,
    'generateKmerGraph': generateKmerGraph,
    'clusterGraph': clusterGraph
}

def main():
    options[funct](tuple(arguments))

main()