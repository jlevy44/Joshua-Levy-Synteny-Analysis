#!/bin/bash
module load bedtools/2.25.0
nextflow run -process.echo true scaffoldPipeline.nf
