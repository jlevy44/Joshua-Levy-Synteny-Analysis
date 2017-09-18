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

String findPeak(String peakfname){

    //String[] parts = peakfname.split("_");
    //String peakfinalname = parts[0];
    int finalPosition = peakfname.lastIndexOf('_');
    String peakfinalname = peakfname.substring(0,finalPosition);
    return peakfinalname;

}

// parse config

pythonPath = findValue('pythonPath ' );
systemPath = findValue('systemPath ' );
blastPath = findValue('blastPath ' );
kmercountPath = findValue('kmercountPath ' );
blastMemory = findValue('blastMemory ' );
reclusterPath = findValue('reclusterPath ');
best500kmerPath = findValue('kmer500BestPath ');
fastaPath = findValue('fastaPath ' );
bedPath = findValue('bedPath ' );
sortPath = findValue('sortPath ' );
sortedbedPath = findValue('sortedbedPath ' );
genome = findValue('genome ' );
peakPath = findValue(' peakPath');
BBstr = findValue('BB ');
n_subgenomes = findValue('n_subgenomes ');
splitLength = findValue('splitFastaLineLength ');
bootstrap = findValue('bootstrap ');
kmerLength = findValue('kmerLength ');
transformMetric = findValue('transformMetric ');
n_neighbors = findValue('n_neighbors ');
metric = findValue('metric ');
kmer_low_count = findValue('kmer_low_count ')
diff_kmer_threshold = findValue('diff_kmer_threshold ')
unionbed_threshold = findValue('unionbed_threshold ')
minChunkSize = findValue('minChunkSize ')
removeNonChunk = findValue('removeNonChunk ')
minChunkThreshold = findValue('minChunkThreshold ')
lowMemory = findValue('lowMemory ')



reduction_techniques = findValue('reduction_techniques ').split(',')
clusterModels = findValue('clusterMethods ').split(',')

//clusterModels = ['KMeans','SpectralClustering']

//preprocess = findValue('preprocess ').asType(Integer);


save = findValue('save ');
old = findValue('old ');

slurm = findValue('slurm ').asType(Integer);
BB = findValue('BB ').asType(Integer);
splitFast = findValue('splitFasta ').asType(Integer);
writeKmer = findValue('writeKmer ').asType(Integer);
fromFasta = findValue('kmer2Fasta ').asType(Integer);
original = findValue('original ').asType(Integer);
writeBlast = findValue('writeBlast ').asType(Integer);
runBlastParallel = findValue('runBlastParallel ').asType(Integer);
b2b = findValue('blast2bed ').asType(Integer);
genMat = findValue('generateClusteringMatrix ').asType(Integer);
kHist = findValue('kmerHistogram ').asType(Integer);
bPeak = findValue('blastPeaks ').asType(Integer);
bPeak2 = findValue('peakClusterMatrix ').asType(Integer);
trans = findValue('transformData ').asType(Integer);
trans2 = findValue('transformDataChildren ').asType(Integer);
clust = findValue('ClusterAll ').asType(Integer);
extract = findValue('extract ').asType(Integer);
kmerBlast = findValue('kmerBlast ').asType(Integer);
kmerGraph = findValue('kmerGraph ').asType(Integer);


genomeSplitName = genome - '_split' - '.fa' + '_split.fa';
blastDBName = genomeSplitName - '.fa';
genomeFullPath = fastaPath + genomeSplitName;
n_clusters = (n_subgenomes.asType(Integer) + 1).asType(String);
originalStr = original.asType(String);
blastMemStr = "export _JAVA_OPTIONS='-Xms5G -Xmx" + blastMemory + "G'"
//genome2 = { original == 1 ? genome - '\n' : genomeSplitName };
//println genome2
//originalGenome = { original == 1 ? fastaPath + genome : genomeFullPath };
kmercountName = genomeSplitName - '.fa' + '.kcount' + '.fa';
blastName = kmercountName - '.fa' + '.BLASTtsv.txt';
workingDir = new File('').getAbsolutePath();
//blastDBName2 = {original == 1 ? genome - '.fa' : blastDBName };
if ( original == 1){
    genome2 = genome - '\n';
    blastDBName2 = genome - '\n' -'.fa';
    originalGenome = fastaPath + genome;
}
else {
    genome2 = genomeSplitName;
    blastDBName2 = blastDBName;
    originalGenome = genomeFullPath;
}
println genome2

check = Channel.from(genomeSplitName,genomeFullPath,kmercountName,blastName,workingDir)
                .subscribe{println it}

//gCh = Channel.from(genome)
//            .subscribe {println it + 'aaaa'}

genomeChan = Channel.from(genome - '\n')

genomeChan2 = Channel.create()

process splitFastaProcess {
executor = 'local'
cpus = { splitFast == 1 ? 2 : 1 }
//'-N 2 -p regular -D . '}

input:
    val genomeName from genomeChan

output:
    stdout genomeChan2

script:
if(splitFast == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py splitFasta ${genomeName} ${fastaPath} ${splitLength}
    """
else
    """
    echo ${genomeSplitName}
    """

}


genomeChan3 = genomeChan2.last()

gCh = genomeChan3
                .subscribe{println it -'\n' + 'aaniofs'}

genomeChan4 = Channel.create()

process writeKmerCount {

//clusterOptions = { slurm == 0 ? { writeKmer == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 9 -p regular -D . '}
cpus = { writeKmer == 1 ? 9 : 1 }

input:
    val genomeName from genomeChan3

output:
    val genomeName into genomeChan4


script:
if(writeKmer == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py writeKmerCount ${fastaPath} ${kmercountPath} ${kmerLength} ${blastMemory}
    """
else
    """
    touch done
    """

}

genomeChan5 = Channel.create()
genomeChan55 = Channel.create()

process kmer2Fasta {
executor = 'local'
//clusterOptions = { slurm == 0 ? { fromFasta == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
//clusterOptions = { slurm == 1 ? '-D .' : '-P plant-analysis.p -cwd' }
cpus = { fromFasta == 1 ? 2 : 1 }
memory = { fromFasta == 1 ? '100 GB' : '60 GB' }


input:
    val genomeName from genomeChan4

output:
    val genomeName into genomeChan5
    val genomeName into genomeChan55


script:
if(fromFasta == 1 && BB == 0)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    python subgenomeClusteringInterface.py kmer2Fasta ${kmercountPath} ${kmer_low_count}
    ${blastMemStr} && makeblastdb -in ${genomeFullPath} -dbtype nucl -out ${blastDBName}.blast_db
    """
else if(fromFasta == 1 && BB == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py kmer2Fasta ${kmercountPath} ${kmer_low_count}
    ${blastMemStr} && bbmap.sh ref=${genomeFullPath}
    """
else
    """
    touch done
    """

}

process createOrigDB {
executor = 'local'
//clusterOptions = { slurm == 0 ? { original == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { original == 1 ? 2 : 1 }


input:
    val genomeName from genomeChan55

script:
if(original == 1 && BB == 0)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    ${blastMemStr} && makeblastdb -in ${originalGenome} -dbtype nucl -out ${blastDBName2}.blast_db
    """
else if(original == 1 && BB == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    ${blastMemStr} && bbmap.sh ref=${originalGenome}
    """
else
    """
    touch done
    """

}


//Channel.watchPath(kmercountPath+kmercountName, 'create,modify')
if(runBlastParallel == 1){
    genomeChan5.map{it -> file(kmercountPath+kmercountName)}
               .splitFasta(by: 50000,file: true)
               .into {kmerFasta; kFast2}
    }
else{
    genomeChan5.map{it -> file(kmercountPath+kmercountName)}
               .into {kmerFasta; kFast2}
}

                    //.set { kmerFasta }

kFast2.subscribe {println it}

blast_result = Channel.create()

process BlastOff {

//clusterOptions = {slurm == 0 ? { writeBlast == 1 ? '-P plant-analysis.p -cwd -l h_rt=24:00:00 -pe pe_slots 16 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : { writeBlast == 1 ? '-N 16 -p regular -D . ' : '-N 1 -p regular -D . ' }}
cpus = { writeBlast == 1 ? 16 : 1 }


input:
    file 'query.fa' from kmerFasta
    //each genomeName from genomeChan55.toList()

output:
    file blast_result
    //val genomeName into genomeChan6
    //${kmercountPath}${kmercountName}

script:
if(writeBlast == 1 && BB == 0)
    """
    #!/bin/bash
    module load blast+/2.6.0
    #cd ${workingDir}
    ${blastMemStr} && blastn -db ${workingDir}/${blastDBName}.blast_db -query query.fa -task "blastn-short" -outfmt 6 -num_threads 15 -evalue 1e-2 > blast_result
    """
else if(writeBlast == 1 && BB == 1)
    """
    #!/bin/bash
    ${blastMemStr} && bbmap.sh vslow=t ambiguous=all noheader=t secondary=t perfectmode=t threads=15 maxsites=2000000000 outputunmapped=f ref=${workingDir}/${genomeFullPath} in=query.fa outm=result.sam
    mv result.sam blast_result
    """
else
    """
    touch blast_result
    """

}

if (writeBlast == 1){
blast_result.collectFile(name: blastPath + blastName)
            .map {file -> genomeSplitName}
            .set {genomeChan6}
            }
else {
blast_result.collectFile()
            .map {file -> genomeSplitName}
            .set {genomeChan6}
}

genomeChan7 = Channel.create()

process blast2bed {
executor = 'local'
//clusterOptions = {slurm == 0 ? { b2b == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { b2b == 1 ? 2 : 1 }


input:
    val genomeName from genomeChan6

output:
    val genomeName into genomeChan7


script:
if(b2b == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py blast2bed ${blastPath}${blastName} ${BBstr} ${lowMemory}
    """
else
    """
    touch done
    """

}

genomeChan8 = Channel.create()
//genomeChan85 = Channel.create()

process genClusterMatrix_kmerPrevalence {

//clusterOptions = {slurm == 0 ? { genMat == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { genMat == 1 ? 2 : 1 }


input:
    val genomeName from genomeChan7

output:
    val genomeName into genomeChan8
    //val genomeName into genomeChan85


script:
if(genMat == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py genClusterKmer ${kmercountPath} ${save} ${genomeName} ${splitLength} ${minChunkSize} ${removeNonChunk} ${minChunkThreshold} ${lowMemory}
    """
else
    """
    touch done
    """

}

//reduction_techniques = ['factor','kpca']
//,'feature']

/*
if(trans == 0){
genomeChan8 = genomeChan8.take(1)
}
*/

process transform_main {

executor = 'local'

//clusterOptions = { slurm == 0 ? { trans == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { trans == 1 ? 2 : 1 }


input:
    val genomeName from genomeChan8
    each technique from reduction_techniques

script:
if(trans == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py transform_main 1 ${reclusterPath} ${technique} ${n_subgenomes} ${transformMetric}
    """
else
    """
    #!/bin/bash
    cd ${workingDir}
    touch main_${technique}_${n_subgenomes}_transformed3D.npy
    """

}


//peaks4 = peaks3.map { it -> findPeak(it) }
bestKmerMatrices = Channel.watchPath(reclusterPath+'*.npz','create,modify')
                          .map { file -> findPeak(file.name) }
                          .unique()
if(trans2 == 1) {
process transform {

//clusterOptions = { slurm == 0 ? { trans2 == 1 ? '-P plant-analysis.p -cwd -l high.c -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}

cpus = { trans2 == 1 ? 2 : 1 }

input:
    val kmerMat from bestKmerMatrices
    each technique from reduction_techniques

//output:
//    val peakName into peaks3


script:
if(trans2 == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py transform_plot ${kmerMat} ${reclusterPath} ${technique} ${n_subgenomes} ${transformMetric}
    """
else
    """
    #!/bin/bash
    cd ${workingDir}
    touch ${kmerMat}_${technique}_${n_subgenomes}_transformed3D.npy
    """

}
}

transformedData = Channel.watchPath('*transformed3D.npy','create,modify')
                         .unique()
                         .map {file -> file.name - 'transformed3D.npy'}
                         //.filter((file.name).startsWith('main'))


kmerBest500Files = Channel.create()
subgenomeFoldersRaw = Channel.create()
process cluster {

//clusterOptions = { slurm == 0 ? { clust == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { clust == 1 ? 2 : 1 }


input:
    val transformedData
    each model from clusterModels

output:
    file 'test.txt' into kmerBest500Files//stdout into kmerBest500Files//val 'kmer500Best_${model}${transformedData}n3.fa' into kmerBest500Files
    file 'test2.txt' into subgenomeFoldersRaw
//${best500kmerPath}/

script:
if(clust == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py cluster ${transformedData}transformed3D.npy ${reclusterPath} ${best500kmerPath} ${model} ${n_subgenomes} ${metric} ${n_neighbors}
    echo ${model}
    cd -
    echo kmer500Best_${model}${transformedData}n${n_clusters}.fa > test.txt
    echo ${model}${transformedData}n${n_clusters} > test2.txt
    """
else
    """
    echo kmer500Best_${model}${transformedData}n${n_clusters}.fa > test.txt
    echo ${model}${transformedData}n${n_clusters} > test2.txt
    """

}
subgenomeFoldersRaw.splitText()
                .filter {it.toString().size() > 1}
                .set {subgenomeFolders}
subgenomeFolders.map { it -> it - '\n' }
                .set { subgenomeFoldersFinal }
                          //Channel.watchPath('analysisOutputs/*.txt','create,modify')
                          //.map {file -> file.name - '.txt'}

//kmerBest500Files_preClassify = Channel.create()

process subgenomeExtraction {


//clusterOptions = { slurm == 0 ? { extract == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 9 -p regular -D . '}
cpus = { extract == 1 ? 9 : 1 }
//memory = {{ extract == 1 ? '30 GB' : '10 MB'}//65.GB * task.attempt : '10 MB' }
//errorStrategy = 'retry' //{ task.exitStatus == 1 ? 'retry' : 'terminate' }
//maxRetries = 0//1//2


input:
    val subgenomeFolder from subgenomeFoldersFinal

output:
    file 'test.txt' into kmerBest500Files_preClassify


script:
if(extract == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py subgenomeExtraction ./analysisOutputs/${subgenomeFolder} ./analysisOutputs/${subgenomeFolder} ${fastaPath} ${genomeSplitName} ${genome2} ${BBstr} ${bootstrap} 0 ${kmerLength} 0 ${best500kmerPath} ${transformMetric} ${originalStr} ${blastMemory} ${kmer_low_count} ${diff_kmer_threshold} ${unionbed_threshold}
    cd -
    echo kmer500Best_${subgenomeFolder}_preClassify.fa > test.txt
    touch test.txt
    """
else
    """
    echo kmer500Best_${subgenomeFolder}_preClassify.fa > test.txt
    touch test.txt
    """
}

//kmerBest500Files = Channel.watchPath(best500kmerPath + '*.fa','create,modify')
//                          .unique()
//                         .flatMap { file -> tuple(file.name, file.name - '.fa') }



kmerBest500Files_preClassify.splitText()
                .filter {it.toString().size() > 1}
                .set {best_kmer1}


kmerBest500Files.splitText()
                .filter {it.toString().size() > 1}
                .set {best_kmer}

          //.filter {( it =~ /recluster/ ) == 0 &&  it.isEmpty() == 0}
/*
best_kmer.map {it -> file(best500kmerPath + it.toString() - '\n')}
          .set {best_kmer2}

best_kmer2.filter { file -> file.exists() }
          .set {best_kmer3}

best_kmer3.map {file -> tuple(file.name - best500kmerPath , file.name - best500kmerPath - '.fa')}
          .into {best_kmerFinal; printkmerFinal}

*/
//best_kmer.map { it -> tuple(it - '\n', it - '\n' - '.fa') } //.filter({ it })//( it =~ /recluster/ ) == 0 &&  it.isEmpty() == 0})//it.contains('recluster') == 0 && it.isEmpty() == 0})
//         .set {best_kmer2}

best_kmer2 = best_kmer.mix(best_kmer1)
           .set {best_kmer3}

best_kmer3.map {it -> file(best500kmerPath + it.toString() - '\n')}
          .set {best_kmer4}

best_kmer4.filter { file -> file.exists() }
          .set {best_kmer5}

best_kmer5.map {file -> tuple(file.name - best500kmerPath , file.name - best500kmerPath - '.fa')}
          .into {best_kmerFinal; printkmerFinal}

                //.filter {exist(file(best500kmerPath + it.toString()))}
                //.into {best_kmer2; printFlat}
                //

//best_kmerSemi = best_kmer2.concat(best_kmer6)
//                .into { best_kmerFinal; printkmerFinal }


a = Channel.from(1..100)
            .randomSample(10)
b = Channel.from(300..400)
            .randomSample(13)
c = a.concat(b)
     .println()

printkmerFinal.println()

//printkmer.subscribe {println it}

//printFlat.subscribe {println it}

kmer_blasted = Channel.create()

process kmerBlastOff {

//clusterOptions = { slurm == 0 ? { kmerBlast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 16 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 16 -p regular -D . '}
cpus = { kmerBlast == 1 ? 16 : 1 }


input:
    set query, queryFolder from best_kmerFinal // kmerBest500Files

output:
    val queryFolder into kmer_blasted

script:
if(kmerBlast == 1 && BB == 0)
    """
    #!/bin/bash
    module load blast+/2.6.0
    cd ${workingDir}
    mkdir ${best500kmerPath}/${queryFolder}
    ${blastMemStr} && blastn -db ${workingDir}/${blastDBName}.blast_db -query ${best500kmerPath}/${query} -task "blastn-short" -outfmt 6 -num_threads 15 -evalue 1e-2 > ${best500kmerPath}/${queryFolder}/${queryFolder}.blast.txt
    """
else if(kmerBlast == 1 && BB == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    mkdir ${best500kmerPath}/${queryFolder}
    echo ${best500kmerPath}/${queryFolder}/${queryFolder}.blast.sam
    ${blastMemStr} && bbmap.sh vslow=t ambiguous=all noheader=t secondary=t perfectmode=t threads=15 maxsites=2000000000 outputunmapped=f ref=${genomeFullPath} in=${best500kmerPath}/${query} outm=${best500kmerPath}/${queryFolder}/${queryFolder}.blast.sam
    mv ${best500kmerPath}/${queryFolder}/${queryFolder}.blast.sam ${best500kmerPath}/${queryFolder}/${queryFolder}.blast.txt"""
else
    """
    touch blast_result
    """

}

process kmerGraphs {

//clusterOptions = { slurm == 0 ? { kmerGraph == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { kmerGraph == 1 ? 2 : 1 }

input:
    val queryFolder from kmer_blasted

script:
if(kmerGraph == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py generateKmerGraph ${best500kmerPath} ${queryFolder} ${n_subgenomes} ${BBstr}
    """
else
    """
    touch blast_result
    """

}