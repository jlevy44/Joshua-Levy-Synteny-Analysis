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
    dbscriptName = genomeName + '.db.sh'
    blastdb = open(dbscriptName, 'w')
    database_list = [fastaPath+genome, genomeName] #FIXME add genome path
    blastdb.write(
        '#!/bin/bash\nmodule load blast+/2.6.0\nmakeblastdb -in %s -dbtype nucl -out %s.blast_db\n' % tuple(database_list))
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
    scaffolds = pickle.load(open('scaffolds.p', 'rb'))
    if peak == '1':
        #df = pd.read_feather('clusteringMatrix.feather')
        data = sps.load_npz('clusteringMatrix.npz')
        peak = 'main'
    else:
        #df = pd.read_feather('%s_clusteringMatrix.feather'%peak)
        data = sps.load_npz('%s/%s_clusteringMatrix.npz'%(reclusterFolder,peak))
    #df = df.set_index(['index'])
    #scaffolds = list(df.axes[0])
    dimensionalityReducers = {'kpca':KernelPCA(n_components=3),'factor':FactorAnalysis(n_components=3),'feature':FeatureAgglomeration(n_clusters=3)}
    data = StandardScaler(with_mean=False).fit_transform(data)
    #for model in dimensionalityReducers:

    if os.path.isfile(peak + '_' + model + 'Reduction.html') == 0:
        if model != 'kpca':
            data = KernelPCA(n_components=49).fit_transform(data)
        transformed_data = dimensionalityReducers[model].fit_transform(data)
        np.save('%s_%s_transformed3D.npy'%(peak,model), transformed_data)
        plots = []
        plots.append(
            go.Scatter3d(x=transformed_data[:, 0], y=transformed_data[:, 1], z=transformed_data[:, 2], name='Data',
                         mode='markers',
                         marker=dict(color='b', size=2), text=scaffolds))
        fig = go.Figure(data=plots)
        py.plot(fig, filename=peak + '_' + model + 'Reduction.html')


def cluster(args):
    file, reclusterFolder, kmer50Path = args
    if 'recluster' not in file:
        dataOld = sps.load_npz('clusteringMatrix.npz')
        scaffolds = pickle.load(open('scaffolds.p', 'rb'))
        kmers = pickle.load(open('kmers.p', 'rb'))
        #kmerIdx = {i : kmer for i, kmer in enumerate(kmers)}
        Tname = file.split('transformed3D')[0]
        transformed_data = np.load(file)

        transformed_data = StandardScaler().fit_transform(transformed_data)

        clustering_names = ['SpectralClustering','KMeans']
        n_clusters = 3

        clustering_algorithms = [SpectralClustering(n_clusters=n_clusters,eigen_solver='arpack',affinity="nearest_neighbors", gamma=1),MiniBatchKMeans(n_clusters=3)]

        for name, algorithm in zip(clustering_names, clustering_algorithms):
            try:
                if os.path.isfile(name + Tname + 'n%d' % n_clusters + 'ClusterTest.html') == 0:
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
                    os.mkdir('analysisOutputs/' + name + Tname + 'n%d' % n_clusters)
                    for key in set(y_pred):
                        # print key
                        cluster_scaffolds = np.array(scaffolds)[y_pred == key]
                        if list(cluster_scaffolds):
                            clusterSize[key] = np.mean(np.apply_along_axis(lambda x: np.linalg.norm(x),1,transformed_data[y_pred == key,:]))#len(cluster_scaffolds)
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
                    for key in set(y_pred)-{testCluster}:
                        with open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/subgenome_%d.txt' % key) as f:
                            f.write('\n'.join(np.array(scaffolds)[y_pred == key]))

                    fig = go.Figure(data=plots)
                    py.plot(fig, filename=name + Tname + 'n%d' % n_clusters + 'ClusterTest.html')

                    #trainData = transformed_data[y_pred != testCluster]
                    #trainData_scaffolds = np.array(scaffolds)[y_pred != testCluster]


                    trainLabels = y_pred[y_pred != testCluster]
                    trainData = dataOld[y_pred != testCluster]
                    kbest = SelectKBest(chi2,k = 50)
                    kbest.fit(trainData,trainLabels)
                    bestFeatures = kbest.pvalues_.argsort()[:50]
                    best_50_kmers = [kmers[i] for i in bestFeatures]
                    sps.save_npz('%s/recluster%s_clusteringMatrix.npz' %(reclusterFolder, name + Tname + 'n%d' % n_clusters), dataOld[:,bestFeatures])#.tocsc()#.tocoo())
                    with open('%s/kmer50Best_%s.fa'%(kmer50Path,name + Tname + 'n%d' % n_clusters),'w') as f:
                        f.write('\n'.join('>%s\n%s'%(kmer,kmer) for kmer in best_50_kmers))

            except:
                print 'Unable to cluster completely using ' + name + ' for ' + Tname
    else:
        scaffolds = pickle.load(open('scaffolds.p', 'rb'))
        # kmerIdx = {i : kmer for i, kmer in enumerate(kmers)}
        Tname = file.split('transformed3D')[0]
        transformed_data = np.load(file)

        transformed_data = StandardScaler().fit_transform(transformed_data)

        clustering_names = ['SpectralClustering', 'KMeans']
        n_clusters = 3

        clustering_algorithms = [
            SpectralClustering(n_clusters=n_clusters, eigen_solver='arpack', affinity="nearest_neighbors", gamma=1),
            KMeans(n_clusters=3)]

        for name, algorithm in zip(clustering_names, clustering_algorithms):
            try:
                if os.path.isfile(name + Tname + 'n%d' % n_clusters + '50kmerreclusterTest.html') == 0:
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
                    os.mkdir('analysisOutputs/' + name + Tname + 'n%d' % n_clusters)
                    for key in set(y_pred):
                        # print key
                        cluster_scaffolds = np.array(scaffolds)[y_pred == key]
                        with open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/recluster_subgenome_%d.txt'%key) as f:
                            f.write('\n'.join(cluster_scaffolds))
                        if list(cluster_scaffolds):
                            clusterSize[key] = np.mean(np.apply_along_axis(lambda x: np.linalg.norm(x), 1,
                                                                           transformed_data[y_pred == key,
                                                                           :]))  # len(cluster_scaffolds)
                            if clusterSize[key] == max(clusterSize.values()):
                                testCluster = key
                            plots.append(
                                go.Scatter3d(x=transformed_data[y_pred == key, 0],
                                             y=transformed_data[y_pred == key, 1],
                                             z=transformed_data[y_pred == key, 2],
                                             name='Cluster %d, %d points' % (key, len(cluster_scaffolds)),
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
                        with open('analysisOutputs/' + name + Tname + 'n%d' % n_clusters + '/subgenome_%d.txt' % key) as f:
                            f.write('\n'.join(np.array(scaffolds)[y_pred == key]))

                    fig = go.Figure(data=plots)
                    py.plot(fig, filename=name + Tname + 'n%d' % n_clusters + '50kmerreclusterTest.html')
            except:
                print 'Unable to cluster completely using ' + name + ' for ' + Tname

def fai2bed(args):
    genome = args[0]
    Fasta(genome)
    bedFastaDict = defaultdict(list)
    with open(genome+'.fai','r') as f, open(genome + '.bed','w') as f2:
        for line in f:
            if line:
                lineList = line.split('\t')
                bedline = '%s\t%s\t%s\n'%(lineList[0],'0',lineList[1])
                f2.write(bedline)
                bedFastaDict[lineList[0]] = [bedline]
    return bedFastaDict

def writeKmerCountSubgenome((args)):
    subgenomeFolder,  = args
    try:
        os.mkdir(subgenomeFolder+'kmercount_files/')
    except:
        pass
    kmercountPath = subgenomeFolder+'kmercount_files/'
    for fastaFile in os.listdir(subgenomeFolder):
        if fastaFile.endswith('.fa') or fastaFile.endswith('.fasta'):
            f = fastaFile.rstrip()
            print f
            scriptName = f[:f.rfind('.')] + '.sh'
            outFileName = f[:f.rfind('.')] + '.kcount'
            lineOutputList = [subgenomeFolder, fastaFile, kmercountPath, outFileName]
            subprocess.call('module load bbtools && kmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s/%s out=%s/%s k=23 -Xmx60g' % tuple(
                    lineOutputList),shell=True)



def subgenomeExtraction(args):
    subgenome_folder, fastaPath, genomeName, model, originalKmerCount = args
    bedDict = fai2bed((fastaPath+genomeName,))

    for file in os.listdir(subgenome_folder):
        if file and file.endswith('.txt'):
            with open(subgenome_folder+file,'r') as f,open(subgenome_folder+file.replace('.txt','.bed'),'w') as f2:
                for line in f:
                    if line:
                        f2.write(bedDict[line.strip('\n')][0])
            subprocess.call('bedtools getfasta -fi %s -fo %s -bed %s -name'%(fastaPath+genomeName,subgenome_folder + '%s_'%(model)+file.replace('.txt','.fa'),subgenome_folder+file.replace('.txt','.bed')),shell=True)


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
    'transform_main': transform_main
}

def main():
    options[funct](tuple(arguments))

main()