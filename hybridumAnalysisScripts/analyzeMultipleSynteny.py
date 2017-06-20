from collections import defaultdict
dict1 = defaultdict(list)
dict2 = defaultdict(list)
with open('finalSyntenyMultipleSpecies2.txt','r') as f:
    for line in f:
        linelist = line.split('\t')[0:-1]
        dict1[linelist[0]].append((int(linelist[1]),int(linelist[2])))
print dict1
for key in dict1.keys():
    if len(dict1[key]) > 1:
        for i in range(len(dict1[key])-1):
            dict2[key].append(dict1[key][i+1][0]-dict1[key][i][1])

print dict2