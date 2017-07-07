#!/usr/bin/env nextflow

String findValue(String inputStr) {
configFile = file('scaffoldConfig.txt')
allLines = configFile.readLines()
for( line in allLines ) {
    if(line.startsWith(inputStr)){
        println (line - inputStr - '= ')
        return (line - inputStr - '= ');
    }
}

}
//String workingDir = System.getProperty("user.dir");
writeSh = findValue('writeSh ').asType(Integer);
println 'writeSh ' + writeSh.asType(String)
buildRef = findValue('buildRef ').asType(Integer);
version = findValue('version ' );
CDS = findValue('CDS ' );
CDSFasta = findValue('CDSFasta ' );
geneNameOld = findValue('geneNameOld ');
buildSamp = findValue('buildSample ').asType(Integer);
constructSample = findValue('constructSample ').asType(Integer);

chanBuildSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true).toList()
//chanBuildSamples.subscribe { println "value: $it" }
workingDir = new File('').getAbsolutePath()
chanBuildRef = Channel.fromPath('referenceGenomes' + '/*',type: 'dir', relative: true).toList()
//= Channel
                    //
                    /* .subscribe(onNext: { println it}, onComplete: { println 'Done.' }) */


if(writeSh) {
linkageChannel = Channel.create()
process writeShFiles {

output:
    file 'done' into linkageChannel

"""
#!/bin/bash
touch done
python ${workingDir}/writeShFiles.py
"""
}

}
else {
println "Not writing analysis files..."
linkageChannel = Channel.from('done')
}

if(buildRef) {

linkageChannel1 = Channel.create()

process buildReference {

clusterOptions = '-P plant-analysis.p -cwd -l high.c -pe pe_slots 6'

input:
file 'done' from linkageChannel
each ref from chanBuildRef

output:
file 'done' into linkageChannel1


script:
"""
#!/bin/bash
touch done
cd ${workingDir}/referenceGenomes/${ref}
pwd
sh buildRef.sh
"""
}
}
else {
linkageChannel1 = Channel.from('done')
}



if(buildSamp) {
linkageChannel2 = Channel.create()

process buildSample{
clusterOptions = '-P plant-analysis.p -cwd -l high.c -pe pe_slots 6'

input:
each sample from chanBuildSamples
file 'done' from linkageChannel1

output:
file 'done' into linkageChannel2

//when:
//file('${version}/${sample}/${sample}.cds').length() == 0 || file('${version}/${sample}/${sample}.bed').length() == 0

script:
"""
#!/bin/bash
touch done
cd ${workingDir}/${version}/${sample}
pwd
sh build.sh
"""
}
}
else {
linkageChannel2 = Channel.from('done')
//Channel.fromPath(workingDir + '/referenceGenomes' + '/*',type: 'dir', relative: true)
//Channel.from('done')
}



if(constructSample) {
linkageChannel3 = Channel.create()
linkageChannel4 = Channel.create()
linkageChannel5 = Channel.create()
process nucmerfy {

cpus 6

input:
file 'done' from linkageChannel2
each sample from chanBuildSamples
//val chanBuildSamples
//file 'done'


output:
file 'done' into linkageChannel3
//file 'done'

when:
File('${version}/${sample}/${CDS}nuc.tiling').length() == 0

script:
"""
#!/bin/bash
touch done
cd ${workingDir}/${version}/${sample}
sh nucCommand.sh
"""

}

List references = chanBuildRef

process com_1_2 {

input:
file 'done' from linkageChannel3
each sample from chanBuildSamples

output:
file 'done' into linkageChannel4


script:
"""
#!/bin/python
import os, subprocess, sys
open('done','w').close()
os.chdir('../../..')
from pipelineFunctions import *
os.chdir(${version}/${sample})
if all([os.path.isfile('%s.%s.lifted.anchors' %(${sample}, ref)) and os.stat('%s.%s.lifted.anchors' %(${sample}, ref)).st_size > 0 for ref in ${references}]) == 0:
    subprocess.call('sh constructv1_1.sh',shell=True)
    sampleCount = 0
    for ref in weights.keys():
        sampleCount = replaceGeneNames(${sample}, ref, sampleCount)
try:
    tiling2bed('%snuc.tiling'%${CDS}, ${CDS}, ${sample}, ${sample}+'_%ssyn'%${CDS}+'.bed')
except:
    print sys.exc_info()[0]
replaceGeneNames(${sample},${CDS},0,1)
subprocess.call('sh constructv1_2.sh',shell=True)
"""
}

process allmaps {
clusterOptions = '-P plant-analysis.p -cwd -l h_rt=12:00:00 -pe pe_slots 32 -e OutputFile.txt'
queue 'long'

input:
file 'done' from linkageChannel4
each sample from chanBuildSamples

output:
file 'done' into linkageChannel5



script:
if(File('${version}/${sample}/multipleMapping.bed').length() == 0)
"""
#!/bin/bash
touch done
cd ../../../${version}/${sample}
sh qsub_build.sh
"""
else
    error "Unable to build new genome: ${sample}"
}
}