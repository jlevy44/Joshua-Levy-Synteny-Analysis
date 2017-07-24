import matplotlib.pyplot as plt
histData = []
with open('kmerPrevalence.txt','r') as f:
    for line in f:
        if line:
            histData.append(line.count(','))


fig = plt.figure()
plt.hist(histData,bins = 50)
plt.xlabel('Number of Related Kmers')
plt.ylabel('Count')
plt.title('Histogram of Number of Related Kmers')
plt.savefig('KmerHistogram.png')


