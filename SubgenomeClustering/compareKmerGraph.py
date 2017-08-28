import pandas as pd, numpy as np
import os
from collections import defaultdict

subDict = defaultdict(list)
peakDict = defaultdict(list)
finalDict = defaultdict(list)


for file in os.listdir('.'):
    with open(file,'r') as f:
        if file.endswith('.txt'):
            peakDict[file] = set([line[:line.find('\t')] for line in f.readlines()])
        if file.endswith('.fa'):
            subDict[file] = set([line.strip('\n') for line in f.readlines()[1::2]])

for key in peakDict:
    finalDict[key] = {key2: len(peakDict[key].intersection(subDict[key2])) for key2 in subDict}

#print peakDict
#print subDict

df = pd.DataFrame(finalDict)

print df

df.to_csv('SubgenomeKmers.csv')