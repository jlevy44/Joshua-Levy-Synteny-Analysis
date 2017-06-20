from collections import Counter

fileText = open('analyzeBd.txt','r').read().split()
BdBPs = Counter()
for chunk in fileText:
    if '|314' in chunk:
        BdSegs = [segment for segment in chunk.split('|') if segment.startswith('314-')]
        for segment in BdSegs:
            BdBPs[segment.split('-')[1]] += int(segment.split('-')[-1])-int(segment.split('-')[-2])

print BdBPs