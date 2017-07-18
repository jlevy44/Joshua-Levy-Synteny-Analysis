__author__ = 'sgordon007'

"""
plot the genotype data with multiple subplots, one per lib.
"""

# BTNCA.11436.8.207067.ACTCGCT-TATCCTC.fastq_calling.sorted.bedgraph BTNCX.11436.8.207067.GGAGCTA-CGTCTAA.fastq_calling.sorted.bedgraph BTNGT.11436.8.207067.CGGAGCC-GTAAGGA.fastq_calling.sorted.bedgraph

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def geno_multi_read(a):

    # read in the first lib
    print 'working on file:' + a
    df1 = pd.read_csv(a, sep='\t', header=None, lineterminator='\n')
    # df1.columns = ['chr', 'start', 'stop', 'Dgenome_gRNA1', 'Sgenome_gRNA1']
    df1.columns = ['chr', 'start', 'stop', 'Sgenome_gRNA', 'Tgenome_gRNA']
    test_data = df1.head()
    print test_data

    # would be nice to automatically subset the chromosomes
    # set a list of known chromosomes, iterate through them and use an if statement to query if that is the current chromosome, if so, a
    # print(df.loc[df['B'].isin(['one', 'three'])])
    # chr_list = ['BhD1', 'BhD2', 'BhD3', 'BhD4', 'BhD5', 'BhS1', 'BhS2', 'BhS3', 'BhS4', 'BhS5', 'BhS6', 'BhS7', 'BhS8', 'BhS9', 'BhS10']
    # chr_list = ['Chr01K', 'Chr02K', 'Chr03K', 'Chr04K', 'Chr05K', 'Chr06K', 'Chr07K', 'Chr08K', 'Chr09K', 'Chr01N', 'Chr02N', 'Chr03N', 'Chr04N','Chr05N', 'Chr06N', 'Chr07N','Chr08N', 'Chr09N']
    # chr_list = ['Nt01', 'Nt07', 'Nt08', 'Nt03', 'Nt10', 'Nt11', 'Nt05', 'Nt06', 'Nt20', 'Nt16','Nt09', 'Nt02', 'Nt04', 'Nt14', 'Nt15', 'Nt19', 'Nt23', 'Nt24', 'Nt12', 'Nt13']
    chr_list = ['Nt01', 'Nt02', 'Nt03', 'Nt04', 'Nt05', 'Nt06', 'Nt07', 'Nt06', 'Nt20', 'Nt16', 'Nt09', 'Nt02', 'Nt04','Nt14', 'Nt15', 'Nt19', 'Nt23', 'Nt24', 'Nt12', 'Nt13']

    to_subset = ['Sgenome_gRNA', 'Tgenome_gRNA']
    for c in chr_list:
        chrl = [c]
        chr = c
        dfc = df1.loc[df1['chr'].isin(chrl)]
        # subset columns to plot
        dfc = dfc[to_subset]
        # Create plot for each chr
        fig = plt.figure(figsize=(7, 11))
        ax = fig.add_subplot(111)
        dfc.plot(ax=ax); ax.set_title(chr); plt.xlabel('physical distance (x10 kbp)'); plt.ylabel('gRNA cut site freq')
        f_name = 'chromosome' + chr + '.pdf'
        plt.savefig(f_name)

if __name__ == "__main__":
    geno_multi_read(a=str(sys.argv[1]))

