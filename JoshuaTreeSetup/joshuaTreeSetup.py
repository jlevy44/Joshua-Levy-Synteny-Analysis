import Tkinter as tk
from Tkinter import *
import os, shutil
import numpy as np
import threading
import subprocess


maindir = os.getcwd()+'/'

def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of relevant filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list of files or list information
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path or value of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or specific value trying to find
            configFile.seek(0)
            return line.split()[-1].strip('\n')

class integratedAnalysis():
    def __init__(self,IntconfigFile,master):
        self.Intconfig = IntconfigFile
        frame = Frame(master)
        frame.pack(side = BOTTOM)
        self.readConfigButton = Button(frame, text='Read Integrated Config', command=self.readConfig)
        self.readConfigButton.pack()


    def readConfig(self):
        # open master configuration file
        open(self.Intconfig, 'r').close()
        masterConfigFile = open(self.Intconfig, 'r')

        # grab the following information from the configuration file
        weightsList = parseConfigFindList('Weights info', masterConfigFile)
        findInfoList = np.vectorize(lambda x: parseConfigFindPath(x,masterConfigFile))(['performSynteny', 'performCircos', 'performALLMAPS', 'querySpecies', 'NameAnalysis',
                        'writeFastaOut', 'Loci_Threshold',
                        'pathPython', 'pathSystem', 'pathALLMAPS', 'BdPath', 'pathUnOut', 'pathGFF', 'pathSort',
                        'BPsMergeDist', 'genomePath',
                        'karyotypesFilesPath', 'circosConfigFilesPath', 'LinkPath', 'circosOutPath', 'BPsThreshold',
                        'multipleSeqAlignFastasPath', 'fastaOutputName', 'allMAPImageOutputPath', 'online',
                        'projectName',
                        'nerscUsername', 'nohup', 'cactusRun', 'cactusFolder'])  # find the following information
        # assign values
        self.IntConfigInfo = findInfoList


    def correctFolders(self):
        maindir = os.getcwd()+'/'
        (self.pathALLMAPS,self.BdPath,self.pathUnOut,self.pathGFF,self.pathSort,self.genomePath,self.karyotypesFilesPath,
        self.circosConfigFilesPath,self.LinkPath,self.circosOutPath,self.multipleSeqAlignFastasPath,self.allMAPImageOutputPath,self.cactusFolder) = np.vectorize(lambda x: maindir + x + '/')(['allMAPSAnalysis',
                        'BdFiles','UnOutFiles','GFFfiles','SortFiles','Genomes','CircosInputs/Data/Karyotypes','CircosInputs/ConfigFiles',
                                                    'CircosInputs/Data/Links','CircosOutputs','FastaOut','allMAPSAnalysis/ALLMAPS_Chromosome_Images','CactusRun'])

    def setupFiles(self):
        self.correctFolders()
        for file in os.listdir('.'):
            if file.endswith('.gff') or file.endswith('.gff3'):
                shutil.move(file,self.pathGFF)
            if file.endswith('.fa') or file.endswith('.fai'):
                shutil.move(file,self.genomePath)
            if file.endswith('.unout'):
                shutil.move(file,self.pathUnOut)

class CNSAnalysis():
    def __init__(self,CNSconfigFile,master):
        self.CNSconfig = CNSconfigFile
        frame = Frame(master)
        frame.pack(side = BOTTOM)
        self.readConfigButton = Button(frame,text = 'Read CNS Config', command = self.readConfig)
        self.readConfigButton.pack()

    def readConfig(self):
        configFile = open(self.CNSconfig, 'r')
        self.CNSconfigInfo = np.append(np.vectorize(lambda x: parseConfigFindPath(x,configFile))(['root_folder', 'pathPython',
                                 'checkValidity','conservedFastaPath','pickleSkip','pickleName','fasta2phylip','PhyML',
                                  'bootstrap','treeFile','treeOut','ratioCopy','outputTreeImages']),
                                       np.vectorize(lambda x: parseConfigFindList(x,configFile))(['masterListSpecies','intragenus','intergenus','subgenome']))
        configFile.close()

#app = App()
def printHi():
    print 'hi'


def grabFolders():
    if 'CNSAnalysis' not in os.listdir('.'):
        #commands = ['mkdir joshuaTree', 'git init joshuaTree', 'cd joshuaTree',  'git remote add origin https://github.com/jlevy44/Joshua-Levy-Synteny-Analysis.git',
                    #'git config core.sparseCheckout true', 'echo \"Joshua-Levy-Synteny-Analysis/IntegratedAnalysisFolderStructure/*\" >> .git/info/sparse-checkout',
                    #'git pull origin master']
        commands = ['svn checkout https://github.com/jlevy44/Joshua-Levy-Synteny-Analysis.git/trunk/JoshuaTree','mv JoshuaTree/* .']# --no-auth-cache --non-interactive --trust-server-cert']
        for command in commands:
            subprocess.call(command,shell=True)
        #subprocess.call('git clone --depth 1 https://github.com/jlevy44/Joshua-Levy-Synteny-Analysis.git joshuaTree',shell = True)
        #subprocess.call('git filter-branch --prune-empty --subdirectory-filter joshuaTree/JoshuaTree HEAD',shell = True)

root = Tk()

masterSetup = integratedAnalysis('masterConfig.txt',root)
CNSAnalysisSetup = CNSAnalysis('configCNSAnalysis.txt',root)
frame = Frame(root)
frame.pack(side = TOP)
label_1 = Label(frame,text = 'masterConfig Name:')
label_2 = Label(frame, text='CNSConfig Name:')

entry_1 = Entry(frame)
entry_2 = Entry(frame)

label_1.grid(row=0,sticky=E)
label_2.grid(row=1,sticky=E)

entry_1.grid(row=0,column=1)
entry_2.grid(row=1,column=1)

options = Label(frame,text = 'Options')
options.grid(columnspan = 2)

b1 = Button(frame, text = 'Grab Folders', command = grabFolders)
b1.grid(row = 3)


root.mainloop()

print 'hi'
"""Goal: Incorporate integrated analysis, CNS analysis, """