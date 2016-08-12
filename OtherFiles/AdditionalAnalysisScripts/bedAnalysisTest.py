from pybedtools import *

a = BedTool('PAC4GC.524-PAC2_0.308_5.bed').sort()
b = BedTool('PAC4GC.524-PAC2_0.312_5.bed').sort()
c = BedTool('PAC4GC.524-PAC2_0.383_5.bed').sort()

# figure this out

#a.sort().merge(delim='|')
#abc2=a.cat(b,postmerge=True,force_truncate=False)
abc=a.cat(b,postmerge=False).cat(c,postmerge=False).sort().merge(o='distinct',c=4,delim='|',d=100000)
"""
abc3=a.multi_intersect(i=[b.fn,c.fn])#,distinct=True)#,wa=True,wb=True)
ab = a.intersect(b, wa = True,wb=True).sort()
abc = ab.intersect(c , wa = True, wb = True).sort()
ac = a.intersect(c, wa = True, wb = True).sort()
bc = b.intersect(c, wa = True, wb = True).sort()

"""

print abc
#524_Chr08K	29916188	29975311	308_Chr04_7714154_7689415|312_scaffold_4_35628858_35653273|383_Chr04K_49408420_49424873
for line in str(abc).split('\n'):
    if line:
        lineList = line.split('\t')
        output1 = lineList[0:3]
        print '%s %s %s %s' %(output1[0].split('_')[0],output1[0].split('_')[1],output1[1],output1[2])
        if '|' in lineList[3]:
            output2=lineList[3].split('|')
        else:
            output2 = [lineList[3]]
        for output in output2:
            outputList = output.split('_')
            outputSpecies = outputList[0]
            outputChromosome = output[output.find('_')+1:output.find(outputList[-2])-1]
            outputTuple = (outputSpecies, outputChromosome, outputList[-2],outputList[-1])
            print '%s %s %s %s'%outputTuple


#print (abc).sort().merge(c=)
#print (a+b+c).sort()

#print a.sort() == (abc).sort().merge()

