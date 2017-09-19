import os
import subprocess, shutil
import pybedtools
from pyfaidx import Fasta

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
    print 'my input to kcount to dict is: %s' % inputFile
    dictConverted = {}
    for line in inputFile:
        if line:
            lineList = line.split()
            dictConverted[lineList[0]]=(int(lineList[1].strip('\n')))
    inputFile.close()
    return dictConverted

def compareKmers(kmercountPath, kmercountFiles, gRNA_search):
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
        print 'creating %s' % filename

    """this version of kmerCompare checks for GG at end of Kmer assuming a search for gRNA"""
    if int(gRNA_search):
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
    else:
        """this version of kmerCompare does NOT check for PAM (GG) at end of Kmer"""
        print 'kmerCompare without check for PAM site'
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


# modify to run blast on fasta files output from compareKmers
"""Takes differential kmercount fasta files and BLASTs them to genome"""
def writeBlast(genome, blastPath, kmercountFiles):
    """make blast database for whole genome assembly"""
    genomeName = genome[:genome.rfind('.')]
    dbscriptName = genomeName + '.db.sh'
    blastdb = open(dbscriptName, 'w')
    database_list = [genome, genomeName]
    blastdb.write('#!/bin/bash\nmodule load blast+\nmakeblastdb -in %s -dbtype nucl -out %s.blast_db\n' % tuple(database_list))
    blastdb.close()
    try:
        subprocess.call('nohup sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName
    """ blast differential kmers to whole genome assembly"""
    blast_script_suffix = '.blast.sh'
    for file in kmercountFiles:
        inputFile = '%s.higher.kmers.txt' % (file.split('.')[0])
        f = file.rstrip()
        print f
        blastscriptPrefix = file.split('.')[0]
        blastscriptName = '%s%s' % (blastscriptPrefix, blast_script_suffix)
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

# subprocess.call('qsub -P %s -N %s -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e -l h_rt=26:00:00'
#                  '  %s%s\n'%(projectName,blastscriptPrefix,nerscUsername,'runAllmaps2.sh'),
#                             shell=True)
# qsub -P plant-assembly.p -N BhD1 -cwd -b yes -now no -j yes -m abes -M sgordon@lbl.gov -w e -l exclusive.c -l h_rt=76:00:00 ./runAllmaps.sh


def blast2bed3(blastPath, blastFiles, bedPath, sortPath, genome, chromosomes):
    """Takes a list of genotype files with only one column for pos and converts them to proper bedgraph format to be sorted"""
    print 'blast files contains'
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
        coveragename = f[:f.rfind('.')] + '.sorted.cov'
        if not os.path.exists(sortPath):
            os.makedirs(sortPath)
        b = pybedtools.BedTool(si)
        # if chromosomes, then use the 1Mb window
        if int(chromosomes):
            windows = '%s.win.1Mb.bed' % genome
        else:
            windows = '%s.bed' % genome
        a = pybedtools.BedTool(windows)
        b.sort().saveas(so)
        a.coverage(b).saveas(coveragename)
        bedgname = f[:f.rfind('.')] + '.sorted.cov.bedgraph'
        open(bedgname, 'w').close()
        bedgo = open(bedgname, 'w')
        covFile = open(coveragename, 'r')
        for line in covFile:
            lineInList = line.split()
            lineOutputList = [lineInList[0], int(lineInList[1]), int(lineInList[2]), int(lineInList[3]) ]
            bedgo.write('%s\t%d\t%d\t%d\n' % tuple(lineOutputList))
        covFile.close()
        bedgo.close()

def genomewindows(genome, chromosomes):
    """make genome file and fixed windows for whole genome assembly"""
    if int(chromosomes):
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
    else:
        """ modified version of genomewindows that makes intervals that cover the entire scaffold rather than fixed intervals"""
        genomeName = genome[:genome.rfind('.')]
        faifile = genome + '.fai'
        scaffintf = genome + '.bed'
        open(scaffintf, 'w').close()
        oh = open(scaffintf, 'w')
        inh = open(faifile, 'r')
        for line in inh:
            lineInList = line.split()
            lineOutputList = [lineInList[0], int(0), int(lineInList[1])]
            oh.write('%s\t%d\t%d\n' % tuple(lineOutputList))
        inh.close()
        oh.close()

def bed2unionBed(genome, blastFiles):
    # for blastFile in blastFiles:
    #     f = blastFile.rstrip()
    #     outFileName = f[:f.rfind('.')]+'.bed3'
    bedFile = blastFiles[0]
    f = bedFile.rstrip()
    print f
    inputName = f[:f.rfind('.')]
    outputFileName = inputName + '.union.bedgraph'
    genomeName = genome[:genome.rfind('.')]
    genomeFile = '%s.genome' % genome
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
    f = bedFiles[0]
    inputName = f[:f.rfind('.')]
    outputFileName = inputName + '.union.bedgraph'
    lineOutputList = [outputFileName, genome]
    dbscriptName = outputFileName + '.plot.sh'
    c = open(dbscriptName, 'w')
    c.write('#!/bin/bash\nmodule load bedtools\npython ./GenomeKmerComp/multiplot_gRNA_chrIter_2colBedGraphv2.py %s %s\n' % tuple(lineOutputList))
    c.close()
    try:
        subprocess.call('sh %s' % dbscriptName, shell=True)
    except:
        print 'Unable to run %s via command line..' % dbscriptName



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
    kmercountPath = parseConfigFindPath('kmercountPath',configFile)
    blastPath = parseConfigFindPath('blastPath',configFile)
    bedPath = parseConfigFindPath('bedPath',configFile)
    sortPath = parseConfigFindPath('sortPath',configFile)

    fastaFiles = filter(None, str(subprocess.Popen(['ls', '%s' % fastaPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
    genome = parseConfigFindPath('genome', configFile)
    gRNA_search = parseConfigFindPath('gRNA_search', configFile)
    chromosomes= parseConfigFindPath('chromosomes', configFile)

    # make subdirs
    if not os.path.exists(bedPath):
        os.makedirs(bedPath)
    if not os.path.exists(kmercountPath):
        os.makedirs(kmercountPath)
    if not os.path.exists(blastPath):
        os.makedirs(blastPath)
    if not os.path.exists(bedPath):
        os.makedirs(bedPath)
    if not os.path.exists(sortPath):
        os.makedirs(sortPath)

    # make fixed windows across the genome in bed format
    genomewindows(genome, chromosomes)


    # first analysis step is to split whole genome into subgenome fasta files
    # use something like python ../GenomeKmerComp/fetch.tobacco_subgenome_sequences.py
    # we should automate this soon


    # run kmercountexact.sh on each subgenome to yield kmercount_files
    # write the shell scripts for generating kmercounts
    # pass fastaFiles to bbtools for running kmerCounts
    # writeKmercount(fastaPath, fastaFiles, kmercountPath)

    # get kmer count files
    kmercountFiles = filter(None,str(subprocess.Popen(['ls', '%s' % kmercountPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
    # kmercount_files are input to compareKmers function
    # output are fasta files of 23mers ending with PAM (or not, depending on gRNA_search value in config. set to 1 to do gRNA search.
    # compareKmers(kmercountPath, kmercountFiles, gRNA_search)

    # blast fasta files against the whole genome
    # could filter BLAST output, if needed
    writeBlast(genome, blastPath, kmercountFiles)

    # get list of BLAST files
    blastFiles = filter(None,str(subprocess.Popen(['ls', '%s' % blastPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))

    # then parse blast output to bed3 format:
    blast2bed3(blastPath, blastFiles, bedPath, sortPath, genome, chromosomes)

    # get list of bedfiles
    bedFiles = filter(None, str(subprocess.Popen(['ls', '%s' % bedPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))

    # make union between bedfiles
    bed2unionBed(genome, bedFiles)

    # make plots
    make_plots(genome, bedFiles)


if ( __name__ == '__main__' ):
    real_main()
