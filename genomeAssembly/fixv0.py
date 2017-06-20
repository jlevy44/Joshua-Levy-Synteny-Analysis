import subprocess, os, sys
from collections import defaultdict
import numpy as np
from multiprocessing import Pool, Queue
import subprocess,shutil
from jcvi.formats import gff
from pyfaidx import Fasta

references = [folder for folder in os.listdir('referenceGenomes') if '.' not in folder and folder]

buildSamples = np.vectorize(lambda x: 'Bdist_%s_v0'%(x))(sys.argv[1:])

root = os.getcwd()+'/'
CDSOld = root+'referenceGenomes/'+'314'+'/'+[file for file in os.listdir(root+'referenceGenomes/'+'314') if file.endswith('.fa') and '.cds' in file.lower()][0]
runCommand = lambda x: subprocess.call(x,shell=True)

def formatSamplev0(sample):
    global root
    commands = ['python %sformatBed.py s %s v0'%(root,sample),'python %sformatCDS.py s %s v0'%(root,sample)]
    for command in commands:
        runCommand(command)
        os.chdir(root)

def buildSamplesv0(sample): #sample = Bdist_xxx_v0.fa
    global root
    global  CDSOld
    print sample
    os.chdir('v0/'+sample)
    print 'c ' + os.getcwd()
    fastaNew = sample+'.fa'
    geneNaming = sample.replace('_','')
    writeCommands = ['samtools faidx %s' %fastaNew,'python -m jcvi.formats.gff bed --type=CDS --key=Name %s -o %s' % (
                     geneNaming + '.gff3', sample + '.CDS.bed'),'bedtools getfasta -name -fi ./%s -bed %s.CDS.bed -fo %s.cds'%(fastaNew,sample,sample) ]#'gffread -E %s -o- > %s' % (geneNaming + '.gff3', geneNaming + '.cufflinks.gff'),
                     #'python -m jcvi.formats.gff load %s %s --feature=CDS --id_attribute=Name -o %s' % (geneNaming + '.cufflinks.gff', fastaNew,sample + '.cds')]
    writeCommands2 = ['samtools faidx %s' % fastaNew,
                     'gmap_build --dir=. -d %s %s' % (geneNaming, fastaNew),
                     'gmap --dir=. -d %s -B 5 -A --format=gff3_gene -n 1 -t 8 %s > %s 2> %s' % (
                         geneNaming, CDSOld, geneNaming + '.gff3', geneNaming + '.log'),
                     'python %srenameGenes.py %s %s %s' % (root, geneNaming + '.gff3', 'Bradi', geneNaming),
                     'python -m jcvi.formats.gff bed --type=mRNA --key=Name %s -o %s' % (
                     geneNaming + '.gff3', sample + '.bed'),
                     'python -m jcvi.formats.gff load %s %s --feature=CDS --id_attribute=Name -o %s' % (
                     geneNaming + '.gff3', fastaNew, sample + '.cds')]

    for command in writeCommands:
        runCommand(command)
        print command
    os.chdir(root)



if __name__ == '__main__':
    with open('output.txt', 'a') as f:
        f.write('Outv1')
    p = Pool(processes=6)
    p.map(func=buildSamplesv0, iterable=buildSamples)
    p.map(func=formatSamplev0, iterable=buildSamples)
    #for sample in buildSamples:
    #    os.chdir(root)
    #    buildSamplesv0(sample)
    #    os.chdir(root)
    #    formatSamplev0(sample)