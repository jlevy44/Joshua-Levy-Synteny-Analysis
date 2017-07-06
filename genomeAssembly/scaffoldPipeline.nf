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

chanBuildSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
chanBuildSamples.subscribe { println "value: $it" }
workingDir = Channel.fromPath('../genomeDistAssemblies',type: 'dir')
chanBuildRef = Channel
                    .fromPath('referenceGenomes' + '/*',type: 'dir', relative: true)
                    /* .subscribe(onNext: { println it}, onComplete: { println 'Done.' }) */

if(writeSh) {
process writeShFiles {

input:
val workingDir

"""
#!/bin/bash
python ${workingDir}/writeShFiles.py
"""
}

}
else {
println "Not writing analysis files..."
}

if(buildRef) {
clusterOptions = '-P plant-analysis.p -cwd -l high.c -pe pe_slots 6'

process buildReference {

input:
val chanBuildRef

script:
"""
#!/bin/bash
cd ../../../referenceGenomes/${chanBuildRef}
pwd
sh buildReference.sh
"""
}
}

if(buildSamp) {
process buildSample{

input:
val chanBuildSamples

when:
File('${version}/${chanBuildSamples}/${chanBuildSampels}.cds').length() == 0 || File('${version}/${chanBuildSamples}/${chanBuildSampels}.bed').length() == 0

script:
"""
#!/bin/bash
cd ../../../${version}/${chanBuildSamples}
pwd
sh buildSample.sh
"""
}
}

if(constructSample) {
process nucmerfy {

cpus 6

input:
val chanBuildSamples

when:
File('${version}/${chanBuildSamples}/${CDS}nuc.tiling').length() == 0

script:
"""
#!/bin/bash
cd ../../../${version}/${chanBuildSamples}
sh nucCommand.sh
"""

}

List references = chanBuildRef.toList()

process com_1_2 {

input:
val chanBuildSamples

script:
"""
#!/bin/python
import os, subprocess, sys
os.chdir('../../..')
from pipelineFunctions import *
os.chdir(${version}/${chanBuildSamples})
if all([os.path.isfile('%s.%s.lifted.anchors' %(${chanBuildSamples}, ref)) and os.stat('%s.%s.lifted.anchors' %(${chanBuildSamples}, ref)).st_size > 0 for ref in ${references}]) == 0:
    subprocess.call('sh constructv1_1.sh',shell=True)
    sampleCount = 0
    for ref in weights.keys():
        sampleCount = replaceGeneNames(${chanBuildSamples}, ref, sampleCount)
try:
    tiling2bed('%snuc.tiling'%${CDS}, ${CDS}, ${chanBuildSamples}, ${chanBuildSamples}+'_%ssyn'%${CDS}+'.bed')
except:
    print sys.exc_info()[0]
replaceGeneNames(${chanBuildSamples},${CDS},0,1)
subprocess.call('sh constructv1_2.sh',shell=True)
"""
}

process allmaps {
clusterOptions = '-P plant-analysis.p -cwd -l h_rt=12:00:00 -pe pe_slots 32 -e OutputFile.txt'
queue 'long'

input:
val chanBuildSamples


script:
if(File('${version}/${chanBuildSamples}/multipleMapping.bed').length() == 0)
"""
#!/bin/bash
cd ../../../${version}/${chanBuildSamples}
sh qsub_build.sh
"""
else
    error "Unable to build new genome: ${chanBuildSamples}"
}
}