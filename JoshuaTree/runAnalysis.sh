#!/bin/bash
module load circos/0.64
module load bedtools/2.25.0
python integratedAnalysisPipeline.py
rsync FastaOut/*.fasta CactusRun/FastaFiles/
