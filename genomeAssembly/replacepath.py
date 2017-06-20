import sys

inputList = sys.argv

for e in inputList:
	if e.endswith('.bed'):
		bedFile = e

with open(bedFile,'r') as f:
	text = f.read().replace('.path1','').replace('.path2','')

with open(bedFile,'w') as f:
	f.write(text)
