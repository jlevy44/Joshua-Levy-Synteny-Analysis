import numpy as np
import matplotlib.pylab as plt
from pybedtools import BedTool
import os
with open('MASeg.bed','w') as f:
    for file in os.listdir('.'):
        if file.endswith('ConservedElements.bed'):
            with open(file,'r') as f2:
                f.writelines(np.vectorize(lambda line: line.split('\t')[-1].split(';')[-1].strip('\n') + '\t0\t1\t' +line.split('\t')[-1].split(';')[-2]+'\n')(f2.readlines()))
    f.close()
plotmat = []
for line in str(BedTool('MASeg.bed').sort().merge(c=4,o='distinct',delim='|')).split('\n'):
    if line:
        try:
            plotmat.append([int(speciesCount.split(':')[1]) for speciesCount in line.split('\t')[-1].split('|')[0].split(',')])
        except:
            print line
matrix=np.matrix(plotmat)
print matrix
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.imshow(matrix, interpolation='none', cmap=plt.cm.afmhot,extent=(0.5,np.shape(matrix)[0]+0.5,0.5,np.shape(matrix)[1]+0.5))
plt.colorbar()
plt.show()