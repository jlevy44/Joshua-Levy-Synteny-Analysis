#!/bin/bash

module load bedtools/2.25.0
module load graphviz/2.38.0
nextflow run -process.echo true subgenomePipeline.nf -with-dag flowchart.png -with-timeline timeline.html
#nohup python subgenomeClustering.py &