import numpy as np
import multiprocessing as mp
import sys

global inputStr

def grabLine(positions):
    global inputStr
    if positions[-1] != 'end':
        return inputStr[positions[0]:positions[1]]
    else:
        #print inputStr[positions[0]:]
        return inputStr[positions[0]:]
def wrap(inputStr='',width=60):
    if inputStr:
        positions = np.arange(0,len(inputStr),width)
        posFinal = []
        for i in range(len(positions[:-1])):
            posFinal += [(positions[i],positions[i+1])]
        posFinal += [(positions[-1],'end')]
        print posFinal[0:10]
        if __name__ == '__main__':
            p=mp.Pool(processes=8)
            wrappedLines = p.map(grabLine,posFinal)
            p.close()
            p.join()
            #print wrappedLines[0:10]
        return '\n'.join(wrappedLines)
    else:
        return ''

e = sys.argv

try:
    inputStr=sys.argv[1]
    try:
        width = sys.argv[2]
    except:
        width = 60
except:
    inputStr = ''
    width = 60

wrap(inputStr,width)
