#!/bin/bash
module load PhyML
module load bedtools/2.25.0
nohup python CNSAnalysisv1.py
module load circos/0.64
nohup python circosGeneration.py
