import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from collections import defaultdict


"""
ABRD_CNSElements_Intergenic.bed
TotalSequenceAmount=2514493bps
NumberOfElements=67480
SpeciesDistribution=BdTRi:65378,BdtwD:66021,BdTRgD:52560,ABRS:64893,BhybzeroS:63930,Bdistachyon:66882,BdTRc:67127,BdEighteen:67181,Sbicolor:67060,PobD:56250,Bstacei:64704,BhybzeroD:64583,ABRD:67479,BdTRgS:53262,BdtwS:64087,BhybsixD:63067,BhybsixS:62371,PobS:57329,Osativa:54151

ABRD_CNSElements_Intronic.bed
TotalSequenceAmount=6407201bps
NumberOfElements=193218
SpeciesDistribution=BdTRi:188935,BdtwD:190145,BdTRgD:148764,ABRS:189466,BhybzeroS:186190,Bdistachyon:192986,BdTRc:192907,BdEighteen:192902,Sbicolor:193079,PobD:156834,Bstacei:190058,BhybzeroD:186192,ABRD:193217,BdTRgS:157643,BdtwS:186254,BhybsixD:184088,BhybsixS:182712,PobS:169842,Osativa:165512

ABRD_ConservedElements.bed
"""

def plotStacked(directory):
    rawDataT = {'Species': [], 'Conserved CDS': [], 'CNS Intronic': [], 'CNS Intergenic': []}
    rawData = defaultdict(list)
    with open(directory+'CNSStatistics.txt','r') as f:
        for segment in str(f.read()).split('\n\n'):
            species = segment[:segment.find('_')]
            if species not in rawData.keys():
                rawData[species] = {'Conserved CDS': 0, 'CNS Intronic': 0, 'CNS Intergenic': 0}
            if 'CNSElements_Intergenic.bed' in segment:
                rawData[species]['CNS Intergenic'] = int(segment[segment.find('=')+1:segment.find('bps')])
            if 'CNSElements_Intronic.bed' in segment:
                rawData[species]['CNS Intronic'] = int(segment[segment.find('=')+1:segment.find('bps')])
            if 'Conserved_CDS.bed' in segment:
                rawData[species]['Conserved CDS'] = int(segment[segment.find('=')+1:segment.find('bps')])
    for species in sorted(rawData.keys()):
        if species and '\n' not in species:
            rawDataT['Species'].append(species)
            rawDataT['Conserved CDS'].append(rawData[species]['Conserved CDS'])
            rawDataT['CNS Intronic'].append(rawData[species]['CNS Intronic'])
            rawDataT['CNS Intergenic'].append(rawData[species]['CNS Intergenic'])
    df = pd.DataFrame(rawDataT,columns = ['Species','Conserved CDS','CNS Intronic','CNS Intergenic'])
    # Create the general blog and the "subplots" i.e. the bars
    f, ax1 = plt.subplots(1, figsize=(10, 5))

    # Set the bar width
    bar_width = 0.75

    bar_l = [i+1 for i in range(len(df['Species']))]

    tick_pos = [i + (bar_width / 2) for i in bar_l]

    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the pre_score data
            df['Conserved CDS'],
            # set the width
            width=bar_width,
            # with the label pre score
            label='Conserved CDS',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#F4561D')

    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the mid_score data
            df['CNS Intergenic'],
            # set the width
            width=bar_width,
            # with pre_score on the bottom
            bottom=df['Conserved CDS'],
            # with the label mid score
            label='CNS Intergenic',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#F1911E')

    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the post_score data
            df['CNS Intronic'],
            # set the width
            width=bar_width,
            # with pre_score and mid_score on the bottom
            bottom=[i + j for i, j in zip(df['Conserved CDS'], df['CNS Intergenic'])],
            # with the label post score
            label='CNS Intronic',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#F1BD1A')

    # set the x ticks with names
    plt.xticks(tick_pos, df['Species'],rotation= 'vertical')

    # Set the label and legends
    ax1.set_ylabel("Total Amountof Conserved Sequence (Bps)")
    ax1.set_xlabel("Species")
    plt.legend(loc='upper left')
    plt.title('Breakdown of Conserved Elements Total Sequence Length for Each Species',y=1.08)

    # Set a buffer around the edge
    plt.xlim([min(tick_pos) - bar_width, max(tick_pos) + bar_width])

    plt.savefig(directory+'BreakdownOfCE.pdf',bbox_inches='tight')
    print df


if __name__ == '__main__':
    p = Pool(5)
    p.map(plotStacked, ['./20Species/', './27Species/'])