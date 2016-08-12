import os, sys

path1 = '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/DNA_Segment_Extraction/UnOutFiles'
path2= '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/DNA_Segment_Extraction/TestFiles'

for file in os.listdir(path1):
    fileInputOpen = open(path1+'/'+file,'r')
    fileOutputName=file[:file.rfind('.')]+'Test'+'.unout'
    i=0
    open(path2+'/'+fileOutputName,'w').close()
    fileOutputOpen = open(path2+'/'+fileOutputName,'w')
    for line in fileInputOpen:
        fileOutputOpen.write(line)
        if '[' in line:
            i+=1
        if i>100:
            break
    fileInputOpen.close()
    fileOutputOpen.close()


