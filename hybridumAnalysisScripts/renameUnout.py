import subprocess
listPathUnOut = str(subprocess.Popen(['ls', '%s' % '/global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/UnOutFiles'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()).split('\n')
proteumDict = {'590':'003','601':'005','598':'007','596':'009','593':'011'}

for file in listPathUnOut:
    if file and file.endswith('.unout'):
        file = file.strip('\n')
        renameFile = file
        renameFile = renameFile.replace('564','001').replace('.pre.BhD','')
        for key in proteumDict.keys():
            if key in renameFile:
                renameFile = renameFile.replace(key,proteumDict[key])
        subprocess.call(['mv %s %s'%(file, renameFile)],shell=True)

