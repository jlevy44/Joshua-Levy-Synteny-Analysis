import numpy as np, os
from pybedtools import BedTool
import matplotlib.pyplot as plt

listSubsetSpeciesFiles = np.intersect1d([file for file in os.listdir('20Species') if file.endswith('AllCNSElements.bed')],[file for file in os.listdir('27Species') if file.endswith('AllCNSElements.bed')])
findCoverage = np.vectorize(lambda line: float(line.split('\t')[-1]) if len(line.split('\t'))==7 else np.nan)
with open('CoverageFile.txt','w') as f:
    f.write('Percent of that one data set is a subset of the other:\n\n')
    for species in listSubsetSpeciesFiles:
        a = BedTool('./20Species/' + species).sort().merge().intersect(BedTool('./27Species/' + species).sort().merge(),wao= True).sort().merge(c=7 ,o= 'sum')
        b = BedTool('./27Species/' + species).sort().merge().intersect(BedTool('./20Species/' + species).sort().merge(),wao = True).sort().merge(c=7, o= 'sum')
        c1 = 0  # total sequence
        c2 = 0  # coverage sequence
        for line in str(a).split('\n'):  # .saveas('CoverageSmalltoLarge_'+species)
            if line:
                c1 += int(line.split('\t')[2])-int(line.split('\t')[1])
                c2 += int(line.split('\t')[-1])
        print c1,c2
        f.write(species + '\n[Small/large]: ' + str(float(c2) / float(c1) * 100.) + '\n')
        c1 = 0
        c2 = 0
        # np.vectorize()(str(BedTool('./27Species/'+species).coverage(BedTool('./20Species/'+species))))#.saveas('CoverageLargetoSmall_'+species)
        for line in str(b).split('\n'):  # .saveas('CoverageSmalltoLarge_'+species)
            if line:
                c1 += int(line.split('\t')[2]) - int(line.split('\t')[1])
                c2 += int(line.split('\t')[-1])
        f.write('[Large/small]: ' + str(float(c2) / float(c1) * 100.) + '\n\n')
        print c1, c2
        a.saveas('CoverageSmalltoLarge_' + species)
        b.saveas('CoverageLargetoSmall_' + species)
        c=BedTool('./20Species/' + species).sort().merge().coverage(BedTool('./27Species/' + species).sort().merge())
        d=BedTool('./27Species/' + species).sort().merge().coverage(
            BedTool('./20Species/' + species).sort().merge())
        d1,d2 = (findCoverage(str(c).split('\n')),findCoverage(str(d).split('\n')))
        plt.figure()
        plt.hist(d1[~np.isnan(d1)],40)
        plt.title(species+' Coverage Results for small to large')
        plt.xlabel('Coverage')
        plt.ylabel('Count')
        plt.savefig(species+'_CoverageSmallLargefig.png')
        plt.figure()
        plt.hist(d2[~np.isnan(d2)], 40)
        plt.title(species+' Coverage Results for large to small')
        plt.xlabel('Coverage')
        plt.ylabel('Count')
        plt.savefig(species + '_CoverageLargeSmallfig.png')
    """
    for species in listSubsetSpeciesFiles:
        a = BedTool('./20Species/' + species).sort().merge().coverage(BedTool('./27Species/' + species).sort().merge())
        b = BedTool('./27Species/'+species).sort().merge().coverage(BedTool('./20Species/'+species).sort().merge())
        c1 = 0 # total sequence
        c2 = 0 # coverage sequence
        for line in str(a).split('\n'):#.saveas('CoverageSmalltoLarge_'+species)
            if line:
                c1 += int(line.split('\t')[-2])
                c2 += int(line.split('\t')[-3])
        f.write(species+'\n[Small/large]: '+str(float(c2)/float(c1)*100.)+'\n')
        c1 = 0
        c2 = 0
        #np.vectorize()(str(BedTool('./27Species/'+species).coverage(BedTool('./20Species/'+species))))#.saveas('CoverageLargetoSmall_'+species)
        for line in str(b).split('\n'):#.saveas('CoverageSmalltoLarge_'+species)
            if line:
                c1 += int(line.split('\t')[-2])
                c2 += int(line.split('\t')[-3])
        f.write('[Large/small]: ' + str(float(c2) / float(c1)*100.) + '\n\n')
        a.saveas('CoverageSmalltoLarge_'+species)
        b.saveas('CoverageLargetoSmall_'+species)
        """

