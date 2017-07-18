from pybedtools import *
from pyfaidx import Fasta
import numpy as np
import subprocess, os, shutil, sys
import glob
from shutil import copyfile
import stat

# tobacco
# analysisPath = '/global/projectb/scratch/sgordon/plum_DNASeq/LSR2.mode2.121.meraculous_run/chr_separated/'
analysisPath = '/global/projectb/scratch/sgordon/crispr/tobacco_genomes_to_analyze/'
refGenome = 'Nitab_v4.5_genome_Chr_Edwards2017.fasta'

nerscUsername = 'sgordon'
# projectName = 'plant-analysis.p'
# projectName = 'plant-annotation.p'
projectName = 'plant-assembly.p'

cpus = 3


# # try to run samtools to generate final.scaffolds.fa.unfiltered.fa.fai file:
# samtoolsI = open('runSamIndex.sh', 'w')
# samtoolsI.write('#!/bin/bash\nmodule load samtools\nsamtools faidx GCF_000390325.2_Ntom_v01_genomic.fna\n')
# samtoolsI.write('samtools faidx GCF_000393655.1_Nsyl_genomic.fna\n')
# samtoolsI.close()
#
# subprocess.call('nohup sh runSamIndex.sh', shell=True)
#
refGenomeFastaObj = Fasta(refGenome)

# output hybridum chromosomes to individual fasta files:
faiFile = open(refGenome + '.fai','r')
s = ""
for line in faiFile:
    if 'Chr' in line or 'B' in line or 'chr' in line or 'Nt' in line:
        chrom = line.split()[0]
        chrFastaOutFile = open(chrom + '.fasta', 'w')
        chrFastaOutFile.write('>%s\n%s\n' % (chrom, str(refGenomeFastaObj[chrom][:])))
        chrFastaOutFile.close()
        s += "{0}.fasta,".format(chrom)
faiFile.close()
s = s.rstrip(",")
print s

# write a shell script to run seal
sealI = open('runSeal.sh', 'w')
# sealI.write('#!/bin/bash\nmodule load bbtools\nseal.sh ref='+s+', in=GCF_000393655.1_Nsyl_genomic.fna pattern=out_%.seal.fasta outu=unmatched.seal.fasta ambig=all refnames overwrite=t refstats=match.stats threads=16\n')
sealI.write('#!/bin/bash\nmodule load bbtools\nseal.sh ref='+s+', in=GCF_000390325.2_Ntom_v01_genomic.fna pattern=out_%.seal.fasta outu=unmatched.seal.fasta ambig=all refnames overwrite=t refstats=match.stats threads=16\n')
sealI.close()
st = os.stat('runSeal.sh')
os.chmod('runSeal.sh', st.st_mode | stat.S_IEXEC)
# if you want to run it immediately
try:
    subprocess.call('nohup sh runSeal.sh', shell=True)
except:
    print 'Unable to run seal via command line..'

# # once the scaffolds are partitioned to chromosomes then we need to create a map for allmaps for each chromosome and output to a new dir
# sealFiles = []
# for file in os.listdir("./"):
#     if file.endswith(".seal.fasta"):
#         sealFiles.append(file)
#
# s = ""
# for file in sealFiles:
#     # make a dir and move the files there
#     directory = file + '_map'
#     if not os.path.exists(directory):
#         os.makedirs(directory)
#     try:
#         copypath = './' + directory + '/' + file
#         print copypath
#
#         copyfile(file, copypath)
#         os.chdir(directory)
#         #make index for seq in seal chrom file
#         chromFasta = Fasta(file)
#         # write a shell script to run seal
#         mapI = open('runbbMap.sh', 'w')
#         mapI.write(
#             '#!/bin/bash\nmodule load bbtools\nrm -R ./ref\nSHRED=' + file + '\n'
#             'shred.sh in=$SHRED out=$SHRED.shreds.overlapping.fa length=500 minlength=500 overlap=25 threads=16\n'
#             'shred.sh in=$SHRED out=$SHRED.shreds.sparse.fa length=500 minlength=500 overlap=-500 threads=16\n'
#              #ref should be probably just be only the chromosome from reference by which scaffolds were binned
#             'REF=' + analysisPath + refGenome + '\n'
#             'bbmap.sh ref=$REF threads=16\n'
#             'bbmap.sh idtag=t outputunmapped=f in=$SHRED.shreds.overlapping.fa outm=$SHRED.mapped.overlapping.sam minid=0.97 maxindel=3 ef=0.01 threads=16\n'
#             'bbmap.sh idtag=t outputunmapped=f in=$SHRED.shreds.sparse.fa outm=$SHRED.mapped.sparse.sam minid=0.97 maxindel=3 ef=0.01 threads=16\n')
#         mapI.close()
#
#         st = os.stat('runbbMap.sh')
#         os.chmod('runbbMap.sh', st.st_mode | stat.S_IEXEC)
#         subprocess.call('nohup sh runbbMap.sh', shell=True)
#         os.chdir("..")
#     except:
#         print 'Unable to run mapping via command line for' + file + '..'
#
# os.chdir(analysisPath)
#
# for file in sealFiles:
#     # make a dir and move the files there
#     directory = file + '_map'
#     os.chdir(directory)
#
#     myFile = file + '.fai'
#     faiFile = open(myFile, 'r')
#
#     ### use peach
#     sparseFile = file + '.mapped.sparse.sam'
#     denseFile = file + '.mapped.overlapping.sam'
#     samfileSparse = open(sparseFile, 'r')
#     samfileDense = open(denseFile, 'r')
#
#     ### create output files for the subset of the samfile into D or S relevant parts:
#     open('new_v1.sam', 'w').close()
#     samFileNew = open('new_v1.sam', 'w')
#
#     ### modify to take in two different sam files, a dense SAM marker set for small scaffolds and a less dense SAM for scaffolds > 25kb
#
#     ### find the sam entries for the D subgenome:
#     chrDictSparse = {}
#     chrDictDense = {}
#
#     for line in faiFile:
#         # print line
#         scaffoldLength = line.split()[1]
#         if scaffoldLength > 3000:
#             chrDictSparse[line.split()[0]] = {}
#         else:
#             chrDictDense[line.split()[0]] = {}
#     faiFile.close()
#
#     for line in samfileSparse:
#         if '@' not in line:
#             key = line.split()[0].split('_')[0]
#             if chrDictSparse.has_key(key):
#                 samFileNew.write(line)
#     samfileSparse.close()
#
#     for line in samfileDense:
#         if '@' not in line:
#             key = line.split()[0].split('_')[0]
#             if chrDictDense.has_key(key):
#                 samFileNew.write(line)
#     samfileDense.close()
#
#     samFileNew.close()
#
#     ########## end of initial analysis
#
#
#     open('allmaps.queryGene99.bed', 'w').close()
#     open('allmaps.targetGene99.bed', 'w').close()
#     open('allmaps.queryGene98.bed', 'w').close()
#     open('allmaps.targetGene98.bed', 'w').close()
#
#     open('allmaps.syntelogs99.txt', 'w').close()
#     open('allmaps.syntelogs98.txt', 'w').close()
#
#     fileBhD99 = open('BhD99.bed', 'w')
#     fileBhD98 = open('BhD98.bed', 'w')
#
#     allmapsBhDqueryGene99 = open('allmaps.queryGene99.bed', 'w')
#     allmapsBhDtargetGene99 = open('allmaps.targetGene99.bed', 'w')
#
#     allmapsBhDqueryGene98 = open('allmaps.queryGene98.bed', 'w')
#     allmapsBhDtargetGene98 = open('allmaps.targetGene98.bed', 'w')
#
#     syntelogsBhD99 = open('allmaps.syntelogs99.txt', 'w')
#     syntelogsBhD98 = open('allmaps.syntelogs98.txt', 'w')
#
#     D = 0
#
#     samfile = 'new_v1.sam'
#     inputFile = open(samfile, 'r')
#     for line in inputFile:
#         if '_' in line and '-' in line:
#             global D
#             D += 1
#             Schr = line.split()[0].split('_')[0]
#             Sstart = int(line.split()[0].split('_')[1].split('-')[0])
#             Send = int(line.split()[0].split('_')[1].split('-')[1])
#             Chrom = line.split()[2]
#             Cstart = int(line.split()[3])
#             Cend = Cstart + 500
#             queryGeneD = str('T0000')
#             targetGeneD = str('Q_0000')
#             if 'YI' in line:
#                 if float(line.split()[-1].split(':')[-1]) > 99:
#                     identity = float(line.split()[-1].split(':')[-1])
#                     fileBhD99.write('%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n' % (
#                     Schr, 'JGI', 'shreds', Sstart, Send, Chrom, Cstart, Cend, identity))
#                     syntelogsBhD99.write('%s%d\t%s%d\t%d\n' % (targetGeneD, D, queryGeneD, D, 1000))
#                     allmapsBhDqueryGene99.write(
#                         '%s\t%d\t%d\t%s%d\t%d\t%s\n' % (Schr, Sstart, Send, queryGeneD, D, 1000, '+'))
#                     allmapsBhDtargetGene99.write(
#                         '%s\t%d\t%d\t%s%d\t%d\t%s\n' % (Chrom, Cstart, Cend, targetGeneD, D, 1000, '+'))
#                 else:
#                     fileBhD98.write('%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n' % (
#                     Schr, 'JGI', 'shreds', Sstart, Send, Chrom, Cstart, Cend))
#                     syntelogsBhD98.write('%s%d\t%s%d\t%d\n' % (targetGeneD, D, queryGeneD, D, 1000))
#                     allmapsBhDqueryGene98.write(
#                         '%s\t%d\t%d\t%s%d\t%d\t%s\n' % (Schr, Sstart, Send, queryGeneD, D, 1000, '+'))
#                     allmapsBhDtargetGene98.write(
#                         '%s\t%d\t%d\t%s%d\t%d\t%s\n' % (Chrom, Cstart, Cend, targetGeneD, D, 1000, '+'))
#     inputFile.close()
#
#     os.chdir("..")
#
# os.chdir(analysisPath)
#
# for file in sealFiles:
#
#     # make syntelog and associated files for allmaps
#     directory = file + '_map'
#     os.chdir(directory)
    # allmap = open('setupAllmaps.sh', 'w')
    # allmap.write('#!/bin/bash\n')
    # allmap.write('python -m jcvi.assembly.syntenypath bed allmaps.syntelogs99.txt --switch --scale=10000 --qbed=allmaps.targetGene99.bed --sbed=allmaps.queryGene99.bed -o map99.bed\n')
    # allmap.write('python -m jcvi.assembly.syntenypath bed allmaps.syntelogs98.txt --switch --scale=10000 --qbed=allmaps.targetGene98.bed --sbed=allmaps.queryGene98.bed -o map98.bed\n')
    # allmap.write('cat map98.bed | sed \'s/5:/25:/g\' | sed \'s/4:/24:/g\' | sed \'s/3:/23:/g\' | sed \'s/2:/22:/g\' | sed \'s/1:/21:/g\' | sed \'s/6:/26:/g\' | sed \'s/7:/27:/g\' | sed \'s/8:/28:/g\' | sed \'s/9:/29:/g\' | sed \'s/10:/31:/g\' > map98tr.bed\n')
    # allmap.write('python -m jcvi.assembly.allmaps mergebed map99.bed map98tr.bed -o mapMer.bed\n')
    # allmap.close()
    # st = os.stat('setupAllmaps.sh')
    # os.chmod('setupAllmaps.sh', st.st_mode | stat.S_IEXEC)
    # try:
    #     # test
    #     print 'running allmaps setup for ' + file
    #     subprocess.call('nohup sh setupAllmaps.sh', shell=True)
    # except:
    #     print 'Unable to run setup allmaps files for' + file + '..'
    # run allmaps
    # runAllmap = open('runAllmaps.sh', 'w')
    # runAllmap.write('#!/bin/bash\n')
    # runAllmap.write('module unload python/2.7.4\n')
    # runAllmap.write('module load python/2.7-anaconda_4.2.0\n')
    # runAllmap.write('source activate allmaps_env_genepool\n')
    # runAllmap.write('python -m jcvi.assembly.allmaps path mapMer.bed ' + file + ' --skipconcorde --cpus=%d --weightsfile=sweight.txt\n' % (cpus))
    # runAllmap.close()
    # st = os.stat('runAllmaps.sh')
    # os.chmod('runAllmaps.sh', st.st_mode | stat.S_IEXEC)
    # # #make proper weights file
    # # weightFile = open('sweight.txt', 'w')
    # # weightFile.write('map99 8\n')
    # # weightFile.write('map98tr 1\n')
    # # weightFile.close()
    # try:
    #     # submitting to high.q queue, change cpus variable to 30 cores at top of file
    #     print 'submitting qsub for ' + file
    #     subprocess.call('qsub -q high_excl.q -P %s -N %s_allmaps -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
    #                     ' -l exclusive.c -l h_rt=56:00:00 %s\n' % (projectName, file, nerscUsername, './runAllmaps.sh'), shell=True)
    # except:
    #     print 'Unable to run allmaps via qsub command for' + file + '..'

    #     # submitting to normal queue, change cpus variable to 14 cores at top of file
    #     print 'submitting qsub for ' + file
    #     subprocess.call('qsub -P %s -N %s_allmaps -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
    #                     ' -l exclusive.c -l h_rt=76:00:00 %s\n' % (projectName, file, nerscUsername, './runAllmaps.sh'), shell=True)
    # except:
    #     print 'Unable to run allmaps via qsub command for' + file + '..'
    #
    #     # submitting to normal queue, change cpus variable to 14 cores at top of file
    #     print 'submitting qsub for ' + file
    #     subprocess.call('qsub -P %s -N %s_allmaps -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
    #                     ' -pe pe_slots 5 -l h_rt=146:00:00 %s\n' % (projectName, file, nerscUsername, './runAllmaps.sh'), shell=True)
    # except:
    #     print 'Unable to run allmaps via qsub command for' + file + '..'

    # # run on the gpint
    # st = os.stat('runAllmaps.sh')
    # os.chmod('runAllmaps.sh', st.st_mode | stat.S_IEXEC)
    # subprocess.call('nohup sh runAllmaps.sh', shell=True)
    #
    #
    # os.chdir("..")



