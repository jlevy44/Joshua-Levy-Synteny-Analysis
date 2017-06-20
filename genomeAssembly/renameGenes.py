import sys

args = sys.argv[1:]

gffFile = args[0]

CDSgeneNaming = args[1]

newGeneName = args[2] 

with open(gffFile,'r') as f:
        text = f.read().replace(CDSgeneNaming,newGeneName)

with open(gffFile,'w') as f:
        f.write(text)
