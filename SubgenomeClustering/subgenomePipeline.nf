#!/usr/bin/env nextflow

String findValue(String inputStr) {
configFile = file('kmerCountCompare.txt')
allLines = configFile.readLines()
for( line in allLines ) {
    if(line.startsWith(inputStr)){
        println (line - inputStr - '= ')
        return (line - inputStr - '= ');
    }
}

}

// parse config

pythonPath = findValue('pythonPath ' );
systemPath = findValue('systemPath ' );
blastPath = findValue('blastPath ' );
kmercountPath = findValue('kmercountPath ' );
fastaPath = findValue('fastaPath ' );
bedPath = findValue('bedPath ' );
sortPath = findValue('sortPath ' );
sortedbedPath = findValue('sortedbedPath ' );
genome = findValue('genome ' );
preprocess = findValue('preprocess ').asType(Integer);
save = findValue('save ');
old = findValue('old ');
peakPath = findValue(' peakPath')
splitFast = findValue('splitFasta ').asType(Integer);
writeKmer = findValue('writeKmer ').asType(Integer);
fromFasta = findValue('fromFasta ').asType(Integer);


genomeSplitName = genome - '_split' - '.fa' + '_split.fa'

genomeFullPath = fastaPath + genomeSplitName

kmercountName = genomeSplitName - '.fa' + '.kcount' + '.fa'

blastName = kmercountName - '.fa' + '.BLASTtsv.txt'

workingDir = new File('').getAbsolutePath()

genomeChan = Channel.from(genome)

genomeChan2 = Channel.create()

process splitFasta {

clusterOptions = { splitFast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan

output:
    stdout genomeChan2

script:
if(splitFast == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py splitFasta ${genomeName} ${fastaPath}
    """
else
    """
    echo ${genomeSplitName}
    """

}

genomeChan3 = genomeChan2.last()

genomeChan4 = Channel.create()

process writeKmerCount {

clusterOptions = { writeKmer == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan3

output:
    val genomeName into genomeChan4


script:
if(writeKmer == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py writeKmerCount ${fastaPath} ${kmercountPath}
    """
else
    """
    touch done
    """

}

genomeChan5 = Channel.create()

process kmer2Fasta {

clusterOptions = { fromFasta == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan4

output:
    val genomeName into genomeChan5


script:
if(fromFasta == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py kmer2Fasta ${kmercountPath}
    """
else
    """
    touch done
    """

}


genomeChan6 = Channel.create()

process BlastOff {

clusterOptions = { writeBlast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan4

output:
    val genomeName into genomeChan5


script:
if(writeBlast == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    makeblastdb -in ${genomeFullPath} -dbtype nucl -out ${genomeName}.blast_db
    blastn -db ./${genomeName}.blast_db -query ${kmercountName} -task "blastn-short" -outfmt 6 -out ${blastPath}/${blastName} -num_threads 8 -evalue 1e-2
    """
else
    """
    touch done
    """

}


genomeChan7 = Channel.create()

process blast2bed {

clusterOptions = { b2b == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan6

output:
    val genomeName into genomeChan7


script:
if(b2b == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py blast2bed ${blastPath}${blastName}
    """
else
    """
    touch done
    """

}

genomeChan8 = Channel.create()

process genClusterMatrix_kmerPrevalence {

clusterOptions = { genMat == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan7

output:
    val genomeName into genomeChan8


script:
if(genMat == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py genClusterKmer ${kmercountPath} ${genomeName} ${save}
    """
else
    """
    touch done
    """

}

peaks = Channel.create()//.watchPath(peakPath+'*.fa')

process kmerHist {

clusterOptions = { kHist == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan7

output:
    stdout peaks


script:
if(kHist == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py kmerRelatedHistogram ${save}
    """
else
    """
    echo ${kmercountName}
    """

}

peaksWait = peaks.last()

peaksFinal = Channel.watchPath(peakPath+'*.fa','create')
                    .flatMap{peak -> [peak, peak - '.fa' + '.BLASTtsv.txt']}

peaks2 = Channel.create()

process blastPeaks {

clusterOptions = { bPeak == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    set(peakName,peakBlastName) from peaksFinal

output:
    set(peakName,peakBlastName) into peaks2


script:
if(bPeak == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    blastn -db ./${genomeName}.blast_db -query ${peakName} -task "blastn-short" -outfmt 6 -out ${blastPath}/${peakBlastName} -num_threads 8 -evalue 1e-2
    """
else
    """
    touch done
    """

}

peaks3 = Channel.create()

process genPeakCluster {

clusterOptions = { bPeak2 == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    set(peakName,peakBlastName) from peaks2

output:
    //val peakBlastName into peaks3


script:
if(bPeak2 == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py peakClusteringMatrix ${kmercountPath} ${peakName} ${peakBlastName} ${save}
    """
else
    """
    touch done
    """

}

// Add clustering algorithms HERE!!!

