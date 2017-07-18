import os
import subprocess, shutil
import pybedtools
from pyfaidx import Fasta
import pandas as pd
import multiprocessing as mp
import numpy as np
from time import time, clock
from collections import Counter
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn import metrics

# modify to run kmercount on a set of fasta files:
def writeKmercount(fastaPath, fastaFiles, kmercountPath):
    """Takes list of fasta files and runs kmercountexact.sh to generate with only one column for pos and converts them to proper bedgraph format to be sorted"""
    for fastaFile in fastaFiles:
        f = fastaFile.rstrip()
        print f
        scriptName = f[:f.rfind('.')] + '.sh'
        outFileName = f[:f.rfind('.')]+'.kcount'
        lineOutputList = [fastaPath, fastaFile, kmercountPath, outFileName]
        bbtoolsI = open(scriptName, 'w')
        bbtoolsI.write('#!/bin/bash\nmodule load bbtools\nkmercountexact.sh overwrite=true fastadump=f mincount=3 in=%s/%s out=%s/%s k=23 -Xmx20g\n' % tuple(lineOutputList))
        bbtoolsI.close()
        try:
            subprocess.call('nohup sh %s' % scriptName, shell=True)
        except:
            print 'Unable to run %s via command line..' % outFileName


def kmercounttodict(kmercount2fname,kmercountPath):
    """ kmercounttodict function creates kmer : count key value pairs, takes path and file name of a kmer count file"""
    inputFile = open(kmercountPath+kmercount2fname,'r')
    dictConverted = {}
    for line in inputFile:
        if line:
            lineList = line.split()
            dictConverted[lineList[0]]=(int(lineList[1].strip('\n')))
    inputFile.close()
    return dictConverted

def kmer2Fasta(kmercountPath):
    for kmer in os.listdir(kmercountPath):
        if kmer.endswith('.kcount'):
            with open(kmercountPath+kmer,'r') as f, open(kmercountPath+'kmerFasta.fa','w') as f2:
                for line in f:
                    if line and int(line.split('\t')[-1]) >= 100:
                        f2.write('>%s\n%s\n'%tuple([line.split('\t')[0]]*2))

def compareKmers(kmercountPath, kmercountFiles):
    dictOfGenes = {}
    ratio_threshold = 20
    end_dinucleotide = 'GG'
    for file in kmercountFiles:
        # creates a dictionary that associates a species to its dictionary of the kmer : count key value pairs
        # kmercounttodict function is called to create the kmer : count key value pairs
        dictOfGenes[file.split('.')[0]] = kmercounttodict(file,kmercountPath)

    # we now have two dictionaries (or more if we add more than two kmer count files to the kmercount_files path
    # now compare the two dictionaries in both directions to find kmers that are high in kmer dict 1 and low in kmer dict 2 and vice versa
    # this gets the dictionary name for the first kmer dict
    dict1 = dictOfGenes[kmercountFiles[0].split('.')[0]]
    dict2 = dictOfGenes[kmercountFiles[1].split('.')[0]]

    # output kmers and counts for differential kmers
    # output file names
    outFileNames = []
    for file in kmercountFiles:
        outFileNames.append("%s.higher.kmers.txt" % (file.split('.')[0]))

    # create files for writing
    for filename in outFileNames:
        open(filename, 'w').close()

    # check dict 1 against dict 2
    out1 = open(outFileNames[0], 'w')
    # iterate through the keys of dict1 and identify kmers that are at least 10 fold higher in dict1
    for key, value in dict1.iteritems():
        key1 = str(key)
        val1 = value
        val2 = dict2.get(key, 3)
        # require at least 30 fold higher kmers in dict1
        if (val1 / val2) > ratio_threshold:
            # require PAM site on kmer end
            if key1[-2:] == end_dinucleotide:
                out1.write('>%s.%d.%d\n%s\n' % (key, val1, val2, key))
    out1.close()

    # do same for other direction of query, # check dict 2 against dict 1
    out2 = open(outFileNames[1], 'w')
    for key, value in dict2.iteritems():
        key1 = str(key)
        val1 = value
        val2 = dict1.get(key, 3)
        if (val1 / val2) > ratio_threshold:
            # require PAM site on kmer end
            if key1[-2:] == end_dinucleotide:
                out2.write('>%s.%d.%d\n%s\n' % (key, val1, val2, key))
    out2.close()


# modify to run blast on fasta files output from compareKmers
"""Takes differential kmercount fasta files and BLASTs them to genome"""
def writeBlast(genome, blastPath, kmercountPath, kmercountFiles):
    """make blast database for whole genome assembly"""
    genomeName = genome[:genome.rfind('.')]
    dbscriptName = genomeName + '.db.sh'
    blastdb = open(dbscriptName, 'w')
    database_list = [genome, genomeName]
    blastdb.write(
        '#!/bin/bash\nmodule load blast+\nmakeblastdb -in %s -dbtype nucl -out %s.blast_db\n' % tuple(database_list))
    blastdb.close()
    try:
        subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName
    """ blast differential kmers to whole genome assembly"""
    for file in [file2 for file2 in os.listdir(kmercountPath) if file2.endswith('.fa')]:
        inputFile = kmercountPath+file#'%s.higher.kmers.txt' % (file.split('.')[0])
        f = file.rstrip()
        print f
        blastscriptName = '%s.blast.sh' % (file.split('.')[0])
        outFileName = f[:f.rfind('.')]+'.BLASTtsv.txt'
        lineOutputList = [genomeName, inputFile, blastPath, outFileName]
        blast = open(blastscriptName, 'w')
        blast.write('#!/bin/bash\nmodule load blast+\nblastn -db ./%s.blast_db -query %s -task "blastn-short" -outfmt 6 -out %s/%s -num_threads 8 -evalue 1e-2\n' % tuple(lineOutputList))
        blast.close()
        """ send blast jobs"""
        try:
            subprocess.call('nohup sh %s' % blastscriptName, shell=True)
        except:
            print 'Unable to run %s via command line..' % blastscriptName


def blast2bed3(blastFiles, bedPath, sortPath, genome):
    """Takes a list of genotype files with only one column for pos and converts them to proper bedgraph format to be sorted"""
    print 'blast files contains'
    print('\n'.join('{}: {}'.format(*k) for k in enumerate(blastFiles)))
    for blastFile in blastFiles:
        f = blastFile.rstrip()
        outFileName = f[:f.rfind('.')]+'.bed3'
        inputFile = open(f, 'r')
        # open(outFileName, 'w').close()
        # outputFile = open(outFileName, 'w')
        outpath = os.path.join(bedPath, outFileName)
        if not os.path.exists(bedPath):
            os.makedirs(bedPath)
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
        coveragename = f[:f.rfind('.')] + '.sorted.cov.bedgraph'
        if not os.path.exists(sortPath):
            os.makedirs(sortPath)
        b = pybedtools.BedTool(si)
        # genomeName = genome[:genome.rfind('.')]
        windows = '%s.win.1Mb.bed' % genome
        a = pybedtools.BedTool(windows)
        b.sort().saveas(so)
        b.coverage(a).saveas(coveragename)


# def bedSort(bedFiles):
#     for bedFile in bedFiles:
#         f = bedFile.rstrip()
#         print f
#         outFileName = f[:f.rfind('.')] + '.sorted.bed3'
#         b = pybedtools.BedTool(bedFile)
#         b.sort().saveas(outFileName)

# >> > a = pybedtools.example_bedtool('a.bed')
# >> > b = pybedtools.example_bedtool('b.bed')
# >> > c = b.coverage(a)


def genomewindows(genome):
    """make genome file and fixed windows for whole genome assembly"""
    genomeName = genome[:genome.rfind('.')]
    dbscriptName = genomeName + '.genome.sh'
    g = open(dbscriptName, 'w')
    lineOutputList = [genome, genome]
    refGenomeFastaObj = Fasta(genome)
    refGenomeFastaObj.close()
    g.write('#!/bin/bash\ncut -f 1-2 %s.fai > %s.genome\n' % tuple(lineOutputList))
    g.write('bedtools makewindows -g %s.genome -w 1000000 > %s.win.1Mb.bed\n' % tuple(lineOutputList))
    g.close()
    try:
        subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName

    # """calculate blast coverage over the whole genome assembly"""
    # for bedFile in bedFiles:
    #     f = bedFile.rstrip()
    #     print f
    #     inputName = f[:f.rfind('.')]
    #     inputFileName = inputName + '.sorted.bed3'
    #     dbscriptName = inputName + '.cov.sh'
    #     c = open(dbscriptName, 'w')
    #     lineOutputList = [inputFileName, genomeName, genomeName]
    #     c.write('#!/bin/bash\nmodule load bedtools\ncoverageBed -a %s -b GENE_LENGTH_%s.win.1Mb.bed | cut -f 1-4 | sortBed -i stdin > %s.cov\n' % tuple(lineOutputList))
    #     c.close()
    #     try:
    #         subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    #     except:
    #         print 'Unable to run %s via command line..' % dbscriptName


def bed2unionBed(genome, bedFiles):
    bedFile = bedFiles[0]
    f = bedFile.rstrip()
    print f
    inputName = f[:f.rfind('.')]
    outputFileName = inputName + '.union.bedgraph'
    genomeName = genome[:genome.rfind('.')]
    genomeFile = 'GENE_LENGTH_%s' % genome
    sortedbedFiles = []
    #FIXME
    a = sortedbedFiles[0]
    b = sortedbedFiles[1]
    x = pybedtools.BedTool()
    result = x.union_bedgraphs(i=[a.fn, b.fn], g=genomeFile)
    result.saveas(outputFileName)

# scrap code:
#  make union bedgraph from the two bed files
# fall back
# module load bedtools
# unionBedGraphs -i Nitab_v45_genome_Chr_S.cov Nitab_v45_genome_Chr_T.cov -empty -g GENE_LENGTH_Nitab_v4.5_genome_Chr_Edwards2017.fasta > Seans.tobacco.gRNA.bedgraph

# finalBedIntersect = BedTool(bedList[0])
# for bed in bedList[1:]:
#         finalBedIntersect = finalBedIntersect.intersect(BedTool(bed)).sort()
#
# a = pybedtools.BedTool('BhS99ends.bed')
# a.intersect('BhS99ends.bed').saveas('snps-in-exons.bed')


# make plots
# python ../GenomeKmerComp/multiplot_gRNA_chrIter_2colBedGraphv2.py Seans.tobacco.gRNA.bedgraph
def make_plots(genome, bedFiles):
    for bedFile in bedFiles:
        f = bedFile.rstrip()
        print f
        inputName = f[:f.rfind('.')]
    outputFileName = inputName + '.union.bedgraph'
    dbscriptName = outputFileName + '.plot.sh'
    c = open(dbscriptName, 'w')
    c.write('#!/bin/bash\nmodule load bedtools\npython ../GenomeKmerComp/multiplot_gRNA_chrIter_2colBedGraphv2.py\n' % outputFileName)
    c.close()
    try:
        subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName

def findScaffolds():
    with open('correspondence.bed','r') as f:
        lineList = sorted([line.split('\t')[-1].strip('\n') for line in f.readlines() if line])
    return lineList

def findKmerNames(kmercountPath):
    for file in os.listdir(kmercountPath):
        if file.endswith('.fa'):
            print file
            with open(kmercountPath + file,'r') as f:
                listKmer = sorted([line.strip('\n') for line in f.readlines()[1::2] if line])
            return listKmer
def blast2bed(blastFile):
    with open(blastFile,'r') as f, open('blasted.bed','w') as f2:
        for line in f:
            if line:
                l1 = line.split('\t')[1]
                f2.write('\t'.join([l1] + ['0',str(int(l1.split('_')[-1])-int(l1.split('_')[-2]))] + [line.split('\t')[0]])+'\n')

    b = pybedtools.BedTool('blasted.bed').sort().merge(c=4,o='collapse',)
    return b

global dfMatrix
global count
global target
def populateKMatrix(line):
    global dfMatrix
    global count
    global target
    #print line
    if line:
        dfMatrix[line.split('\t')[1]][line.split('\t')[0]] += 1
        count += 1

        if float(count) / 14487787. > target:
            print target
            target += 1

def real_main():
    # parse config file
    configFile = open('kmerCountCompare.txt','r') #generate in integrated analysis
    def parseConfigFindPath(stringFind,configFile):
        """findPath will find path of associated specified string"""
        for line in configFile:
            if stringFind in line: # if find string specified, return pathname/info
                configFile.seek(0)
                return line.split()[-1].strip('\n')

    # paths parsed from config file:
    path = parseConfigFindPath('systemPath',configFile)
    fastaPath = parseConfigFindPath('fastaPath', configFile)
    fastaFiles = filter(None, str(subprocess.Popen(['ls', '%s' % fastaPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
    genome = parseConfigFindPath('genome', configFile)
    print fastaFiles
    # make fixed windows across the genome in bed format
    #genomewindows(genome)

    # first analysis step is to split whole genome into subgenome fasta files
    # use something like python ../GenomeKmerComp/fetch.tobacco_subgenome_sequences.py
    # we should automate this soon

    # run kmercountexact.sh on each subgenome to yield kmercount_files
    # write the shell scripts for generating kmercounts
    # pass fastaFiles to bbtools for running kmerCounts
    kmercountPath = parseConfigFindPath('kmercountPath',configFile)
    #writeKmercount(fastaPath, fastaFiles, kmercountPath)

    kmercountFiles = filter(None,str(subprocess.Popen(['ls', '%s' % kmercountPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
    # kmercount_files are input to compareKmers function
    # output are fasta files of 23mers ending with PAM
    #compareKmers(kmercountPath, kmercountFiles)
    #kmer2Fasta(kmercountPath)

    # blast fasta files against the whole genome
    # could filter BLAST output, if needed
    blastPath = parseConfigFindPath('blastPath',configFile)
    #writeBlast(genome, blastPath, kmercountPath, kmercountFiles)

    # check for BLAST output
    blastFiles = filter(None,str(subprocess.Popen(['ls', '%s' % blastPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
    # then parse blast output to bed3 format:
    #bedPath = parseConfigFindPath('bedPath',configFile)
    #sortPath = parseConfigFindPath('sortPath',configFile)
    #blast2bed3(blastFiles, bedPath, sortPath, genome)

    kmers = findKmerNames(kmercountPath)

    scaffolds = findScaffolds()

    print len(scaffolds), len(kmers)
    #d = {scaffold: {kmer: 0 for kmer in kmers} for scaffold in scaffolds}
    #print d.keys()[0:10]
    global count, target, dfMatrix
    start = clock()
    dfMatrix = pd.DataFrame(np.full((len(scaffolds),len(kmers)),0,dtype=float))
    dfMatrix.columns = kmers
    dfMatrix.rows = scaffolds
    print clock()-start
    #dfMatrix.to_csv('emptyMatrix.csv')
    count = 0
    target = 1
    #with open(blastPath + 'kmerFasta.BLASTtsv.txt','r') as f:
    #    lines = f.readlines()
    """

    """
    start = clock()
    labels = dfMatrix.rows

    #dfMatrix.to_csv('clusteringMatrix.csv')
    #b = blast2bed(blastPath + 'kmerFasta.BLASTtsv.txt')
    #print b.head()
    #b.saveas('blasted_merged.bed')
    #print time.clock() - start
    def addToKMatrix(line):
        global dfMatrix
        scaffold = line.split('\t')[0]
        #print scaffold
        print scaffold
        counts = Counter(line.strip('\n').split('\t')[-1].split(','))
        for key in counts.keys():
            dfMatrix.set_value(scaffold, key, counts[key])
        print scaffold
    print dfMatrix.rows
    """
    if __name__ == '__main__':
        p = mp.Pool(8)
        with open('blasted_merged.bed', 'r') as f:
            for line in f:
                if line:
                    p.apply_async(addToKMatrix,(line,))
        #p.map(populateKMatrix, lines)
        p.close()
        p.join()"""
    with open('blasted_merged.bed', 'r') as f:
        for line in f:
            if line:
                scaffold = line.split('\t')[0]
                # print scaffold
                counts = Counter(line.strip('\n').split('\t')[-1].split(','))
                interval = abs(float(line.split('\t')[2]) - float(line.split('\t')[1]))
                print scaffold, counts
                for key in counts.keys():
                    dfMatrix.set_value(scaffold, key, float(counts[key])/interval)
    """
    with open('blasted_merged.bed','r') as f:
        for line in f:
            if line:
                scaffold = line.split('\t')[0]
                counts = Counter(line.strip('\n').split('\t')[-1].split(','))
                #print counts
                for key in counts.keys():
                    dfMatrix.set_value(scaffold,key,counts[key])
                    """
    dfMatrix.to_csv('clusteringMatrix.csv')
    sample_size = len(scaffolds)

    def bench_k_means(estimator, name, data):
        t0 = time.time()
        estimator.fit(data)
        print('% 9s   %.2fs    %i   %.3f   %.3f   %.3f   %.3f   %.3f    %.3f'
              % (name, (time.time() - t0), estimator.inertia_,
                 metrics.homogeneity_score(labels, estimator.labels_),
                 metrics.completeness_score(labels, estimator.labels_),
                 metrics.v_measure_score(labels, estimator.labels_),
                 metrics.adjusted_rand_score(labels, estimator.labels_),
                 metrics.adjusted_mutual_info_score(labels, estimator.labels_),
                 metrics.silhouette_score(data, estimator.labels_,
                                          metric='euclidean',
                                          sample_size=sample_size)))

    pca = PCA(n_components=min(len(scaffolds),len(kmers)))
    pca = pca.fit_transform(dfMatrix.as_matrix())
    bench_k_means(KMeans(init=pca.components_, n_clusters=3, n_init=1),
                  name="PCA-based",
                  data=dfMatrix.as_matrix())
    # turn above into function and run clustering model on this...


    #bedFiles = filter(None, str(subprocess.Popen(['ls', '%s' % bedPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))

    # don't think we need to find sorted files
    #sortedbedFiles = filter(None, str(subprocess.Popen(['ls', '%s' % bedPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))


    # make union between bedfiles
    #bed2unionBed(genome, bedFiles)

    # make plots
    #make_plots(genome, bedFiles)


if ( __name__ == '__main__' ):
    real_main()

