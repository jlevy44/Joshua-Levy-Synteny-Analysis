import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import shutil
if 0:
    plt.style.use('ggplot')
    chroms = {'chr'+str(i):[] for i in range(1,6)}
    gff = [file for file in os.listdir('.') if file.endswith('.gff3')]
    unout = [file for file in os.listdir('.') if file.endswith('.unout')]
    files = []
    for file in os.listdir('.'):
        if file.endswith('.txt') and (file.endswith('.correspondence.txt') == 0 ) and file.startswith('karyotype.') and (any([file.split('.')[1] in file2 for file2 in gff]) and any([file.split('.')[1] in file2 for file2 in unout])):
            files.append(file.split('.')[1])
            with open(file,'r') as f:
                lines = f.readlines()[0:5]
            for line in lines:
                if line:
                    chroms[line.split()[-1]].append(int(line.split()[-2]))
    for chrom in chroms:
        fig = plt.figure()
        sns.distplot(chroms[chrom])
        fig.savefig(chrom + '_sizes.png')

    df = pd.DataFrame(chroms,index=files)
    df.to_csv('chromosome_info.csv')
if 1:
    gff = [file for file in os.listdir('.') if file.endswith('.gff3')]
    unout = [file for file in os.listdir('.') if file.endswith('.unout')]
    df = pd.read_csv('chromosome_info.csv')
    fig = plt.figure()
    genome_size = np.array(df.sum(axis=1))
    sns.distplot(genome_size)
    fig.savefig('total' + '_sizes.png')
    fig = plt.figure()
    x,y = sns.kdeplot(genome_size, shade=True).get_lines()[0].get_data()
    cdf = scipy.integrate.cumtrapz(y, x, initial=0)
    nearest_05 = np.abs(cdf - 0.1).argmin()
    x_median = x[nearest_05]
    y_median = y[nearest_05]
    plt.vlines(x_median, 0, y_median)
    plt.savefig('bottom10.png')
    keep_bool = genome_size > x[nearest_05]
    keep_genomes = np.array(df.axes[0])[keep_bool]
    print keep_genomes
    print df['Unnamed: 0'][keep_bool == False]#np.array(df.axes[0])[keep_bool == False]
    try:
        os.mkdir('keep_files')
    except:
        pass
    df2 = df[keep_bool]
    df2.to_csv('keep_genomes.csv')
    for genome in keep_genomes:
        if str(genome) != '314':
            for file in gff:
                if str(genome) in file:
                    shutil.copy(file,'keep_files')
            for file in unout:
                if str(genome) in file:
                    shutil.copy(file,'keep_files')

