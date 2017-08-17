//peaks = Channel.create()//genomeChan8.last()
        //            .subscribe {println it}
        //            .fromPath()
                    //.watchPath(peakPath+'*.fa')


process kmerHist {

clusterOptions = { kHist == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan8

output:
    file '${workingDir}PeaksOutNames.txt' into peaks


script:
if(kHist == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py kmerRelatedHistogram ${kmercountName} ${save}
    cd -
    rsync ${workingDir}/PeaksOutNames.txt .
    touch PeaksOutNames.txt
    """
else
    """
    #!/bin/bash
    rsync ${workingDir}/PeaksOutNames.txt .
    touch PeaksOutNames.txt
    """

}

//peaksWait = peaks.last()

peaksFinal = peaks.splitText()
                  .map {it -> file(it)}
                  .splitFasta(by: 30000, file: true)
                  .set { blast_peaks }
                //.flatMap{it -> tuple(it, it - '.fa' + '.BLASTtsv.txt')}
                //.flatMap{file -> tuple(file.name, file.name - '.fa' + '.BLASTtsv.txt')}

                    //.last()
                    //.subscribe {println it}
                    //.fromPath(peakPath+'*.fa')
                    //.watchPath(peakPath+'*.fa','create')
                    //.flatMap{file -> tuple(file, file.name, file.name - '.fa' + '.BLASTtsv.txt')}
                    //.splitText(by: 5000, file: True)
                    //.set { chunksChannel }

peaks2 = Channel.create()

process blastPeaks {

clusterOptions = { bPeak == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    each 'query.fa' from blast_peaks.toList()
    //set peakName, peakBlastName from peaksFinal

output:
    file blast_output
    //set peakName, peakBlastName into peaks2

//-out ${blastPath}/${peakBlastName} >>
script:
if(bPeak == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    blastn -db ./${genomeName}.blast_db -query query.fa -task "blastn-short" -outfmt 6 -num_threads 8 -evalue 1e-2 -out > blast_output
    """//${peakName} ${blastPath}/${peakBlastName}
else
    """
    touch done
    """

}

blast_output.collectFile() { item -> }

// can try .collectFile().subscribe { merged_file -> merged_file.copyTo(out_dir) }



peaks3 = Channel.create()

process genPeakCluster {

clusterOptions = { bPeak2 == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    set peakName, peakBlastName from peaks2

output:
    val peakName into peaks3


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


