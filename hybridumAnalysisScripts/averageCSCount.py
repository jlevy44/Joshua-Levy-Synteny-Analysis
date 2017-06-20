import os, numpy as np
import matplotlib.pyplot as plt#, plotly.plotly as py
mafFiles = [files for files in os.listdir('.') if files.endswith('.maf')]
means = []
for file in mafFiles:
    with open(file,'r') as f:
        means += [np.mean([int(line.split('\t')[3]) for line in segment.split('\n') if line.startswith('s')][:-1]) for segment in str(f.read()).split('\n\n')[1:]]
        #print means
means = [mean for mean in means if str(mean) != 'nan']
print len(means),np.sqrt(np.mean(means)) , np.mean(means), np.std(means)
plt.hist(means,bins=np.arange(0,100,4))
plt.axvline(x=16,color='r',linewidth = 5)
plt.title('')
plt.show()