import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import peakutils
from sklearn.neighbors import KernelDensity
from collections import defaultdict
from scipy.signal import argrelextrema
save=0
load=1
if save:
    with open('kmerPrevalence.txt','r') as f:
        kmerDict = {line.split('\t')[0]:line.count(',') for line in f if line}
    pickle.dump(kmerDict,open('kmerDict.p','wb'))
if load:
    kmerDict = pickle.load(open('kmerDict.p','rb'))
histData = np.array(kmerDict.values())[:,np.newaxis]
reverseLookup = defaultdict(list)
for k in kmerDict:
    reverseLookup[kmerDict[k]].append(k)
if save:
    pickle.dump(reverseLookup, open('reverseLook.p', 'wb'))
if load:
    reverseLookup = pickle.load(open('reverseLook.p', 'rb'))
xplot = np.linspace(0,np.max(histData),300000)[:,np.newaxis]
del kmerDict

kde = KernelDensity(kernel = 'gaussian', bandwidth=500).fit(histData)
exp_log_dens = np.exp(kde.score_samples(xplot))

idxs_peaks = argrelextrema(exp_log_dens,np.greater)[0] #peakutils.indexes(exp_log_dens)
idxs_valleys = argrelextrema(exp_log_dens,np.less)[0] #peakutils.indexes(-exp_log_dens)


fig = plt.figure()
plt.plot(xplot[:], exp_log_dens, '-',
            label="Envelope Kernel '{0}' Density Function".format('Gaussian'))
plt.plot(xplot[idxs_peaks],exp_log_dens[idxs_peaks],'*',color='r',label='Peaks')
plt.plot(xplot[idxs_valleys],exp_log_dens[idxs_valleys],'o',color='r',label='Valleys')
plt.hist(histData[:,0],bins = 50,label='Histogram of Counts',normed=True)
plt.xlabel('Number of Related Kmers (Normalized)')
plt.ylabel('Count')
plt.legend(loc='upper left')
plt.title('Histogram of Number of Related Kmers (Normalized)')
plt.savefig('KmerHistogram.png')


intervals = [(0,xplot[0][0])] + [(xplot[idxs_valleys[i]][0],xplot[idxs_valleys[i+1]][0]) for i in range(len(idxs_valleys)-1)] + [(xplot[idxs_valleys[-1]][0],xplot[-1][0])]
counts = np.array(reverseLookup.keys())
idx = 1
print reverseLookup
for interval in enumerate(intervals):
    try:
        print interval
        keys = counts[np.where((counts <= interval[1]) & (counts >= interval[0]))]
        print keys
        #np.vectorize(lambda x: x <= interval[1] and x >= interval[0])(counts)
        with open('Peak%d_CountInt_%d_%d.fa'%tuple([idx+1]+map(int,interval)),'w') as f:
            for key in keys:
                f.write('\n'.join('>%s\n%s'%(val,val) for val in reverseLookup[key])+'\n')
        idx += 1
    except:
        with open('ErrFile.txt','a') as f:
            f.write(str(idx)+'\t%s'%str(interval)+'\n')
