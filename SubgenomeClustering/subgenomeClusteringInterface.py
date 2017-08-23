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
                subprocess.call('sh %s' % scriptName, shell=True)#nohup
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
        subprocess.call('module load blast+/2.6.0 && blastn -db ./%s.blast_db -query %s -task "blastn-short" -outfmt 6 -out %s/%s -num_threads 8 -evalue 1e-2' % tuple(lineOutputList),shell=True)

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
    kmercountPath, save, genome = args
    save = int(save)
    kmers = findKmerNames((kmercountPath, genome))
    scaffolds = findScaffolds(1)
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
    with open('blasted_merged.bed', 'r') as f:
        for line in f:
            if line:
                listLine = line.rstrip('\n').split('\t')
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
        main, reclusterFolder, model  = args
    except:
        main, reclusterFolder, model = ('1','null','kpca')
    transform_plot((main,reclusterFolder,model))

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
    peak, reclusterFolder, model = args
    if peak == '1':
        peak = 'main'
    if os.path.exists(peak + '_' + model + 'Reduction.html') == 0:#isfile
        if peak == 'main':
            # df = pd.read_feather('clusteringMatrix.feather')
            data = sps.load_npz('clusteringMatrix.npz')
        else:
            # df = pd.read_feather('%s_clusteringMatrix.feather'%peak)
            data = sps.load_npz('%s/%s_clusteringMatrix.npz' % (reclusterFolder, peak))
        scaffolds = pickle.load(open('scaffolds.p', 'rb'))
        # df = df.set_index(['index'])
        # scaffolds = list(df.axes[0])
        dimensionalityReducers = {'kpca': KernelPCA(n_components=3), 'factor': FactorAnalysis(n_components=3),
                                  'feature': FeatureAgglomeration(n_clusters=3)}
        data = StandardScaler(with_mean=False).fit_transform(data)
        # for model in dimensionalityReducers:
        if model != 'kpca' and peak.startswith('recluster') == 0:
            data = KernelPCA(n_components=499).fit_transform(data)
        if peak.startswith('recluster'):
            data = data.toarray()
        transformed_data = dimensionalityReducers[model].fit_transform(data)
        np.save('%s_%s_transformed3D.npy'%(peak,model), transformed_data)
        plots = []
        plots.append(
            go.Scatter3d(x=transformed_data[:, 0], y=transformed_data[:, 1], z=transformed_data[:, 2], name='Data',
                         mode='markers',
                         marker=dict(color='b', size=2), text=scaffolds))
        fig = go.Figure(data=plots)
        py.plot(fig, filename=peak + '_' + model + 'Reduction.html')
    else:
        subprocess.call('touch %s'%(peak + '_' + model + 'Reduction.html'),shell=True)

def cluster(args):
    file, reclusterFolder, kmer500Path, clusterMethod = args
    n_clusters = 3
    clustering_algorithms = {'SpectralClustering': SpectralClustering(n_clusters=n_clusters, eigen_solver='lobpcg', affinity= 'nearest_neighbors', n_neighbors=30,random_state=42),#,gamma=1),arpack#amg,affinity="nearest_neighbors")#, n_neighbors=30, gamma=1),# nearestneighbors
                             'KMeans': MiniBatchKMeans(n_clusters=3)}
    name, algorithm = clusterMethod , clustering_algorithms[clusterMethod]
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
            try:
                os.mkdir('analysisOutputs/' + name + Tname + 'n%d' % n_clusters)
            except:
                pass
            for key in set(y_pred):
                # print key
                cluster_scaffolds = np.array(scaffolds)[y_pred == key]
                if list(cluster_scaffolds):
                    clusterSize[key] = np.mean(np.apply_along_axis(lambda x: np.linalg.norm(x),1,transformed_data[y_pred == key,:]))#len(cluster_scaffolds)
                    if clusterSize[key] == min(clusterSize.values()):
                        testCluster = key
                    plots.append(
                        go.Scatter3d(x=transformed_data[y_pred == key, 0], y=transformed_data[y_pred == key, 1],
                                     z=transformed_data[y_pred == key, 2],
                                     name='Cluster %d, %d points, %f distance' % (key, len(cluster_scaffolds),clusterSize[key]), mode='markers',
                                     marker=dict(color=c[key], size=2), text=cluster_scaffolds))

            if hasattr(algorithm, 'cluster_centers_'):
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
                        go.Scatter3d(x=transformed_data[y_pred == key, 0],
                                     y=transformed_data[y_pred == key, 1],
                                     z=transformed_data[y_pred == key, 2],
                                     name='Cluster %d, %d points, %f distance' % (key, len(cluster_scaffolds),clusterSize[key]),
                                     mode='markers',
                                     marker=dict(color=c[key], size=2), text=cluster_scaffolds))

            if hasattr(algorithm, 'cluster_centers_'):
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
    subgenomeFolder, null  = args
    try:
        os.mkdir(subgenomeFolder+'/kmercount_files/')
    except:
        pass
    kmercountPath = subgenomeFolder+'/kmercount_files/'
    for fastaFile in os.listdir(subgenomeFolder):
        if 'higher.kmers' not in fastaFile and (fastaFile.endswith('.fa') or fastaFile.endswith('.fasta')):
            f = fastaFile.rstrip()
            print f
            outFileName = f[:f.rfind('.')] + '.kcount'
            lineOutputList = [subgenomeFolder+'/', fastaFile, kmercountPath, outFileName]
            subprocess.call('module load bbtools && kmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s/%s out=%s/%s k=23 -Xmx60g' % tuple(
                    lineOutputList),shell=True)

    compareKmers(subgenomeFolder, kmercountPath)

def kmercounttodict(kmercount2fname,kmercountPath):
    """ kmercounttodict function creates kmer : count key value pairs, takes path and file name of a kmer count file"""
    inputFile = open(kmercountPath+kmercount2fname,'r')
    print 'my input to kcount to dict is: %s' % inputFile
    dictConverted = {}
    for line in inputFile:
        if line:
            lineList = line.split()
            dictConverted[lineList[0]]=(int(lineList[1].strip('\n')))
    inputFile.close()
    return dictConverted

def compareKmers(subgenomeFolder,kmercountPath):
    dictOfGenes = {}
    ratio_threshold = 20
    end_dinucleotide = 'GG'
    kmercountFiles = os.listdir(kmercountPath)
    for file in kmercountFiles:
        # creates a dictionary that associates a species to its dictionary of the kmer : count key value pairs
        # kmercounttodict function is called to create the kmer : count key value pairs
        dictOfGenes[file[:file.rfind('.')]] = kmercounttodict(file,kmercountPath) #.split('.')[0]

    # we now have two dictionaries (or more if we add more than two kmer count files to the kmercount_files path
    # now compare the two dictionaries in both directions to find kmers that are high in kmer dict 1 and low in kmer dict 2 and vice versa
    # this gets the dictionary name for the first kmer dict
    dict1 = dictOfGenes[kmercountFiles[0][:kmercountFiles[0].rfind('.')]]#.split('.')[0]]
    dict2 = dictOfGenes[kmercountFiles[1][:kmercountFiles[1].rfind('.')]]

    # output kmers and counts for differential kmers
    # output file names
    outFileNames = []
    for file in kmercountFiles:
        outFileNames.append(kmercountPath + "/%s.higher.kmers.fa" % (file.split('.')[0]))

    # create files for writing
    for filename in outFileNames:
        open(filename, 'w').close()
        print 'creating %s' % filename


    # check dict 1 against dict 2
    out1 = open(outFileNames[0], 'w')
    # iterate through the keys of dict1 and identify kmers that are at least 10 fold higher in dict1
    for key, value in dict1.iteritems():
        key1 = str(key)
        val1 = value
        val2 = dict2.get(key, 3)
        # require at least 30 fold higher kmers in dict1
        if (val1 / val2) > ratio_threshold:
            out1.write('>%s.%d.%d\n%s\n' % (key, val1, val2, key))
    out1.close()

    # do same for other direction of query, # check dict 2 against dict 1
    out2 = open(outFileNames[1], 'w')
    for key, value in dict2.iteritems():
        key1 = str(key)
        val1 = value
        val2 = dict1.get(key, 3)
        if (val1 / val2) > ratio_threshold:
            out2.write('>%s.%d.%d\n%s\n' % (key, val1, val2, key))
    out2.close()

def blast2bed3(subgenomeFolder,blastPath, bedPath, sortPath, genome):
    """Takes a list of genotype files with only one column for pos and converts them to proper bedgraph format to be sorted"""
    print 'blast files contains'
    blastFiles = os.listdir(blastPath)
    print('\n'.join('{}: {}'.format(*k) for k in enumerate(blastFiles)))
    for blastFile in blastFiles:
        f = blastFile.rstrip()
        outFileName = f[:f.rfind('.')]+'.bed3'
        input_list = [blastPath, f]
        inpath = '%s/%s' % tuple(input_list)
        inputFile = open(inpath, 'r')
        outpath = os.path.join(bedPath, outFileName)
        bo = open(outpath, 'w')
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
        if not os.path.exists(genome.replace('.fai','.bed')):
            fai2bed((genome))
        shutil.copy(genome.replace('.fai','.bed'),subgenomeFolder)
        windows = '%s.bed' % genome
        a = BedTool(subgenomeFolder + '/' + windows)
        b.sort().saveas(so)
        a.coverage(b).saveas(coveragename)
        bedgname = f[:f.rfind('.')] + '.sorted.cov.bedgraph'
        open(subgenomeFolder + '/' + bedgname, 'w').close()
        bedgo = open(subgenomeFolder + '/' + bedgname, 'w')
        covFile = open(coveragename, 'r')
        for line in covFile:
            lineInList = line.split()
            lineOutputList = [lineInList[0], int(lineInList[1]), int(lineInList[2]), int(lineInList[3]) ]
            bedgo.write('%s\t%d\t%d\t%d\n' % tuple(lineOutputList))
        covFile.close()
        bedgo.close()

def bed2unionBed(genome, subgenomeFolder, bedPath):
    # for blastFile in blastFiles:
    #     f = blastFile.rstrip()
    #     outFileName = f[:f.rfind('.')]+'.bed3'
    blastFiles = os.listdir(bedPath)
    bedFile = blastFiles[0]
    f = bedFile.rstrip()
    print f
    inputName = f[:f.rfind('.')]
    outputFileName = subgenomeFolder + '/' +inputName + '.union.bedgraph'
    genomeName = genome[:genome.rfind('.')]
    subprocess.call('cut -f 1-2 %s.fai > %s.genome'%(genome,subgenomeFolder + '/genome'),shell=True)
    genomeFile = subgenomeFolder + '/genome.genome'
    a_blast = blastFiles[0]
    b_blast = blastFiles[1]
    a_name = a_blast[:a_blast.rfind('.')] + '.sorted.cov.bedgraph'
    b_name = b_blast[:b_blast.rfind('.')] + '.sorted.cov.bedgraph'
    a = pybedtools.BedTool(a_name)
    a_sort = a.sort()
    b = pybedtools.BedTool(b_name)
    b.sort()
    b_sort = b.sort()
    x = pybedtools.BedTool()
    result = x.union_bedgraphs(i=[a_sort.fn, b_sort.fn], g=genomeFile, empty=True)
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

def kmerratio2scaffasta(subgenomePath, genome):
    a = subgenomePath + '/' + [file for file in os.listdir(subgenomePath) if 'union.bedgraph' in file][0]
    ubedg = open(a, 'r')
    genomeFastaObj = Fasta(genome)

    ### define output filenames
    genomeprefix = genome[:genome.rfind('.')]
    o1 = subgenomePath + '/' + genomeprefix + '.subgenomeA.fasta'
    o2 = subgenomePath + '/' + genomeprefix + '.subgenomeB.fasta'
    o1h = open(o1, 'w')
    o2h = open(o2, 'w')

    # parse the unionbed file to subset
    for line in ubedg:
        if line:#'Chr' in line or 'chr' in line or 'B' in line or 'Sc' in line:
            # print line
            scaff = str((line.split()[0]).rstrip())
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
    o2h.close()
    ubedg.close()

def subgenomeExtraction(args):
    subgenome_folder, fastaPath, genomeName, originalGenome = args
    bedDict = fai2bed((fastaPath+genomeName,))

    for file in os.listdir(subgenome_folder):
        if file and file.endswith('.txt'):
            with open(subgenome_folder+'/'+file,'r') as f,open(subgenome_folder+'/'+file.replace('.txt','.bed'),'w') as f2:
                for line in f:
                    if line:
                        f2.write(bedDict[line.strip('\n')][0])
            subprocess.call('bedtools getfasta -fi %s -fo %s -bed %s -name'%(fastaPath+genomeName,subgenome_folder + '/%s_'%('model')+file.replace('.txt','.fa'),subgenome_folder+'/'+file.replace('.txt','.bed')),shell=True)
    blastPath = subgenome_folder+'/blast_files/'
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
    writeKmerCountSubgenome((subgenome_folder,0))
    writeBlast((originalGenome,blastPath,subgenome_folder+'/kmercount_files/',fastaPath))
    blast2bed3(subgenome_folder, blastPath, bedPath, sortPath, fastaPath+originalGenome)
    bed2unionBed(fastaPath+originalGenome, subgenome_folder, bedPath)
    kmerratio2scaffasta(subgenome_folder, fastaPath+originalGenome)
    #make_plots(genome, bedFiles)

def generateKmerGraph(args):
    kmerPath, kmerName = args
    blastPath = kmerPath + '/' + kmerName + '/'
    with open(blastPath+ kmerName + '.blast.txt', 'r') as f, open(blastPath + 'blasted.bed', 'w') as f2:
        for line in f:
            if line:
                l1 = line.split('\t')[1]
                f2.write('\t'.join([l1] + ['0', str(int(l1.split('_')[-1]) - int(l1.split('_')[-2]))] + [
                    line.split('\t')[0]]) + '\n')

    b = pybedtools.BedTool(blastPath + 'blasted.bed').sort().merge(c=4, o='collapse', )
    b.saveas(blastPath + 'blasted_merged.bed')
    kmerGraph = nx.Graph()
    with open(blastPath + 'blasted_merged.bed', 'r') as f:
        for line in f:
            if line:
                #kmerCountsNumbers = [kmerIdx[kmer] for kmer in set(line.rstrip('\n').split('\t').split(','))]
                kmerGraph.add_edges_from(list(combinations(set(line.rstrip('\n').split('\t')[-1].split(',')),2)))#kmerCountsNumbers, 2)))

    histData = np.array([kmerGraph.degree(kmer) for kmer in kmerGraph.nodes()])[:, np.newaxis]
    xplot = np.linspace(0, np.max(histData), 300000)[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=25).fit(histData)
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
    plt.axis('off')
    nx.draw(kmerGraph, pos=nx.spectral_layout(kmerGraph), nodecolor='r', edge_color='b')
    plt.savefig(blastPath + 'kmerGraph_%s.png'%(kmerName), bbox_inches="tight")




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
    'generateKmerGraph': generateKmerGraph
}

def main():
    options[funct](tuple(arguments))

main()