import sys,os
from collections import defaultdict, Counter
from fai2karyotypeMod import fai2karyotype
from pybedtools import BedTool
import subprocess
import numpy as np
import shutil

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
    configFile.seek(0)
    return listOfItems

confDict = {'ticks':"""show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = dims(ideogram,radius_outer)
orientation      = out
label_multiplier = 1e-6
color            = black
size             = 20p
thickness        = 3p
label_offset     = 5p
format           = %d

<tick>
spacing        = 1u
show_label     = no
size           = 10p
</tick>

<tick>
spacing        = 5u
show_label     = yes
label_size     = 20p
size           = 15p
</tick>

<tick>
spacing        = 10u
show_label     = yes
label_size     = 24p
</tick>

</ticks>""",
            'label':"""show_label       = yes
label_font       = default
label_radius     = dims(image,radius)-30p
label_size       = 24
label_parallel   = yes
label_case       = lower
label_format     = eval(sprintf("chr%s",var(label)))
""",
            'position':"""radius           = 0.775r
thickness        = 30p
fill             = yes
fill_color       = black
stroke_thickness = 2
stroke_color     = black""",
            'r0r1':"""# set track radius values based on track counter
r1  = eval(sprintf("%fr",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
r0  = eval(sprintf("%fr",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))""",
            'ideogram':"""<ideogram>

<spacing>
default = 0.01r
break   = 0.5r
</spacing>

<<include position.conf>>
<<include label.conf>>
<<include bands.conf>>

radius*       = 0.95r

</ideogram>""",
            'bands':"""show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0"""}

for conf in confDict.keys():
    with open(conf+'.conf','w') as f:
        f.write(confDict[conf])
        f.close()

with open('configCNSAnalysis.txt','r') as f:
    inputSpecies = parseConfigFindList('masterListSpecies', f)
protId = defaultdict(list)
for line in inputSpecies:
    if line:
        print line
        protId[line.split('_')[0]] = line.split('_')[1]

#protId = {'Bdistachyon':'314','Bstacei':'316','Osativa':'323','Phallii':'308','Pvirgatum':'383','Sbicolor':'313','Sitalica':'312'}
#inputList = sys.argv
inputList = protId.keys()#['Bdistachyon','Bstacei','Osativa','Phallii','Pvirgatum','Sbicolor','Sitalica']
listFiles = os.listdir('.')
speciesDict = defaultdict(list)
"""files = [[file for file in listFiles if species in file and 'Conserved_CDS' in file][0],
                                [file for file in listFiles if species in file and 'CNSElements_Intergenic' in file][0],
                                [file for file in listFiles if species in file and 'CNSElements_Intronic' in file][0]]"""

for species in inputList:
    if species:
        speciesDict[species] = [[file for file in listFiles if species in file and 'ConservedElements' in file][0],
                                [file for file in listFiles if protId[species] in file and 'Genes' in file and file.endswith('.bed3')][0],
                                 fai2karyotype([file for file in listFiles if protId[species] in file and file.endswith('.fai')][0],species,300000),
                                'heatmap.%s.txt'%species,open('heatmap.%s.txt'%species,'w')]
generateRadii = np.linspace(0.40,0.80,len(speciesDict.keys())+1)
totalNumberSpecies = float(len(speciesDict.keys()))
for species in speciesDict.keys():
    histInterval = defaultdict(list)
    with open([file for file in listFiles if protId[species] in file and file.endswith('.fai')][0],'r') as f:
        for line in f:
            #bedread = '\n'.join('\t'.join('%s\t%d\t%d'%tuple([line.split('\t')[0]]+sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3])))))
            interval = sorted(np.vectorize(lambda x: int(x))(line.split('\t')[1:3]))
            histInterval[line.split('\t')[0]] = list(np.arange(0.,interval[-1],250000.)) + [interval[-1]]
        bedHist = BedTool('\n'.join('\n'.join('\t'.join([key] + [str(int(x)) for x in [histInterval[key][i],histInterval[key][i+1]]]) for i in range(len(histInterval[key])-1)) for key in histInterval.keys()),from_string=True)
    print speciesDict[species][1]
    with open(speciesDict[species][1],'r') as f:
        bedGenes = BedTool(f.read(),from_string=True).sort().merge()
    bedHistGeneFinal = bedHist.intersect(bedGenes,wao=True).sort().merge(c=7,o='sum',d=-1)
    with open('%s_geneDensity.txt'%(species),'w') as f:
        for line in str(bedHistGeneFinal).split('\n'):
            if line:
                lineList = line.split('\t')
                f.write('\t'.join(lineList[0:3])+'\t%f'%(float(lineList[-1])/(float(lineList[2])-float(lineList[1])))+'\n')
    for species2 in speciesDict.keys():
        speciesDict[species2][-1].close()
        speciesDict[species2][-1] = open(speciesDict[species2][-2],'w')
    with open('histogramCount%s.txt'%species,'w') as f:
        for i in range(1):#3
            file = speciesDict[species][i]
            #file in speciesDict[species][0:3]:
            with open(file,'r') as bedfile:
                for line in bedfile:
                    if line:
                        try:
                            countSeq = Counter()
                            for countOfSpecies in line.split('\t')[-1].split(';')[1].split(','):
                                countSeq[countOfSpecies.split(':')[0]] = int(countOfSpecies.split(':')[1])
                            f.write(line[:line.rfind('\t')+1]+str((float(line.split('\t')[2])-float(line.split('\t')[1]))*float(len(set(countSeq.elements())))/totalNumberSpecies)+'\n')
                            for species2 in countSeq.keys():
                                if countSeq[species2] > 0:
                                    output = str(i+2)
                                    speciesDict[species2][-1].write(line[:line.rfind('\t')] + '\n')
                                else:
                                    output = '1'

                        except:
                            pass
        f.close()
    with open('histogramCount%s.txt'%species,'r') as f:
        bedHist2 = BedTool(f.read(), from_string=True).sort().merge(c=4,o='mean').saveas('histogramCount%s.txt'%species)
        bedSpeciesHist = bedHist.intersect(bedHist2, wao=True).sort().merge(c=[7,8], o=['sum','sum'], d=-1)
    with open('histogramCount%s.txt'%species,'w') as f:
        for line in str(bedSpeciesHist).split('\n'):
            if line:
                #if not float(line.split('\t')[3]):
                #    print line
                if line.split('\t')[4] != '0':
                    f.write('\t'.join(line.split('\t')[0:3]+[str(float(line.split('\t')[3])/float(line.split('\t')[4]))])+'\n')
                else:
                    f.write('\t'.join(
                        line.split('\t')[0:3] + [str(float(line.split('\t')[3]))]) + '\n')
    for species2 in speciesDict.keys():
        speciesDict[species2][-1].close()
        with open(speciesDict[species2][-2],'r') as f:
            reads = f.read()
        with open(speciesDict[species2][-2],'w') as f2:
            for line in str(bedHist.intersect(BedTool(reads,from_string=True).sort().merge(),wao=True).sort().merge(c=7,o='sum',d=-1)).split('\n'):
                if line:
                    try:
                        f2.write('\t'.join(line.split('\t')[0:3]+[str(float(line.split('\t')[-1])/250000.)])+'\n')
                    except:
                        pass

    # now configure circos files
    print speciesDict[species][2], [(os.getcwd()+'/',species2) for species2 in speciesDict.keys()], os.getcwd()+'/'+'histogramCount%s.txt' %species,os.getcwd()+'/'+'%s_geneDensity.txt'%species
    circosconf = """
show_histogram = yes
show_heatmap   = yes
use_rules      = yes

<<include colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>
<<include bands.conf>>
<<include position.conf>>
<<include label.conf>>



<image>
<<include etc/image.conf>>
</image>

karyotype         = %s
chromosomes_units = 1000000
chromosomes_display_default = yes

# to see how reversing ideograms work - comment out the chromosomes
# line below
#chromosomes       = hs2

# and uncomment the two definitions below
# - first split hs2 into three ideograms
# - now reverse the ideogram with tag "b"
#chromosomes       = hs2[a]:0-60;hs2[b]:70-140;hs2[c]:150-)
#chromosomes_reverse = b

#chromosomes        = hs2[a]:0-30;hs2[b]:50-80;hs2[c]:100-130;hs2[d]:150-180;hs2[e]:190-200;hs2[f]:210-)
#chromosomes_radius = a:0.95r;b:0.9r;c:0.85r;d:0.8r;e:0.75r;f:0.7r

<plots>

show = no


%s

<plot>
show         = conf(show_histogram)
type         = histogram
file         = %s
thickness    = 2
#color        = black
fill_under   = yes
fill_color   = blue
r0           = 0.80r
r1           = 0.90r
orientation = out
max_gap      = 5u
z = 10
</plot>

<plot>
show         = conf(show_histogram)
type         = histogram
file         = %s
orientation  = out
thickness    = 1
#color        = black
fill_under   = yes
fill_color   = green
r0           = 0.90r
r1           = 1r
max_gap      = 5u
z = 10
</plot>
</plots>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim"""%(os.getcwd()+'/'+speciesDict[species][2],'\n'.join("""<plot>
show             = conf(show_heatmap)
type             = heatmap
layers      = 1
margin      = 0.02u
#orientation = out
color = white, spectral-7-div, black
color_mapping = 1
thickness   = 1
padding     = 1
min = 0
max = 0.45
#color            = black
#fill_color = yellow
stroke_thickness = 5
#scale_log_base   = 0.25
#stroke_color     = black
file             = %s
r0             = %fr
r1             = %fr
#<rules>
#use = conf(use_rules)
#<rule>
#condition     = var(value) == 1
#color         = white
#</rule>
#<rule>
#condition     = var(value) > 1
#color         = black
#</rule>
#</rules>
</plot>"""%(os.getcwd()+'/'+'heatmap.'+speciesDict.keys()[i]+'.txt',generateRadii[i],generateRadii[i+1]) for i in range(len(speciesDict.keys()))),os.getcwd()+'/'+'histogramCount%s.txt'%(species),os.getcwd()+'/'+species+'_geneDensity.txt')
    with open('circos.conf','w') as f:
        f.write(circosconf)
        f.close()
    print os.getcwd()+'/'+'circos.conf'
    print ['circos','-conf',os.getcwd()+'/'+'circos.conf','-outputfile',
                     species,'-outputdir',os.getcwd()]
    subprocess.call(['circos','-conf',os.getcwd()+'/'+'circos.conf','-outputfile',species,'-outputdir',os.getcwd()])


"""<plot>
<<include r0r1.conf>>
file             = data/6/variation.heatmap.txt
stroke_thickness = 0
min              = 2000
max              = 250000
</plot>

<plot>
<<include r0r1.conf>>
scale_log_base   = 0.5
</plot>

<plot>
<<include r0r1.conf>>
scale_log_base   = 1   # this is the default value
</plot>

<plot>
<<include r0r1.conf>>
scale_log_base   = 2
</plot>

<plot>
<<include r0r1.conf>>
scale_log_base   = 3
</plot>

<plot>
<<include r0r1.conf>>
scale_log_base   = 5
</plot>

<plot>
<<include r0r1.conf>>
color            = conf(plots,color_alt)
file             = data/6/heatmap.step.txt
pattern          = hline,vline
color_mapping    = 0  # default
min              = 0
max              = 10
stroke_thickness = 0
</plot>

<plot>
<<include r0r1.conf>>
color            = conf(plots,color_alt)
file             = data/6/heatmap.step.txt
pattern          = hline,solid,vline
color_mapping    = 1
min              = 0
max              = 10
stroke_thickness = 0
</plot>

<plot>
<<include r0r1.conf>>
color            = conf(plots,color_alt)
file             = data/6/heatmap.step.txt
pattern          = hline,solid,vline
color_mapping    = 2
min              = 0
max              = 10
stroke_thickness = 0
</plot>

<plot>
<<include r0r1.conf>>
color            = conf(plots,color_alt)
file             = data/6/heatmap.step.txt
pattern          = hline,checker,vline
color_mapping    = 2
min              = 2
max              = 8
stroke_thickness = 0
</plot>"""

"""<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

karyotype   = %s

chromosomes_units  = 1000000
#chromosomes        = hs1;hs2
#chromosomes_breaks = -hs1:120-140
chromosomes_display_default = yes

track_width = 0.05
track_pad   = 0.02
track_start = 0.95

<plots>

type    = heatmap

<rules>
<rule>
condition     = var(value) = 0
color         = white
</rule>
<rule>
condition     = var(value) = 1
color         = black
</rule>
#<rule>
#condition     = var(value) = 2
#color         = green
#</rule>
#<rule>
#condition     = var(value) = 3
#color         = blue
#</rule>
</rules>

# default file for all tracks
#file             = data/6/snp.number.1mb.txt

# a 9 color diverging spectral palette specified using a color list name
color  = spectral-9-div

# referenced via conf(plots,color_alt)
color_alt = black,spectral-8-div,grey

# or the reverse list
#color = spectral-9-div-rev

# or you can even combine lists
# color = ylorrd-9-seq-rev,ylgnbu-9-seq

stroke_thickness = 1
stroke_color     = black
min              = 1000
max              = 5000

%s
<\plots>
<plots>
<plot>

# The type sets the format of the track.

type = histogram
file = %s

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the
# ideogram.

r1   = 0.75r
r0   = 0.80r

# Histograms can have both a fill and outline. The default outline is 1px thick black.

fill_color = vdgrey

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.

#thickness = 0p

# Do not join histogram bins that do not abut.

extend_bin = no

# Like for links, rules are used to dynamically alter formatting of
# each data point (i.e. histogram bin). Here, I include the <rule>
# block from a file, which contains the following
#
# <rule>
# condition = on(hs1)
# show      = no
# </rule>
#
# to avoid displaying any data on hs1. The rule is included from a
# file because it is reused again in the track below.

<rules>
</rules>
</plot>
<plot>

# The type sets the format of the track.

type = histogram
file = %s

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the
# ideogram.

r1   = 0.65r
r0   = 0.75r

# Histograms can have both a fill and outline. The default outline is 1px thick black.

fill_color = vdgrey

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.

#thickness = 0p

# Do not join histogram bins that do not abut.

extend_bin = no

# Like for links, rules are used to dynamically alter formatting of
# each data point (i.e. histogram bin). Here, I include the <rule>
# block from a file, which contains the following
#
# <rule>
# condition = on(hs1)
# show      = no
# </rule>
#
# to avoid displaying any data on hs1. The rule is included from a
# file because it is reused again in the track below.

<rules>
</rules>
</plot>

</plots>
<<include etc/housekeeping.conf>>
data_out_of_range* = trim"""#%(os.getcwd()+'/'+speciesDict[species][2],'\n'.join("""<plot>
#<<include r0r1.conf>>
#file             = %sheatmap.%s.txt
#</plot>"""%(os.getcwd()+'/',species2) for species2 in speciesDict.keys()),os.getcwd()+'/'+'histogramCount%s.txt' %species,os.getcwd()+'/'+'%s_geneDensity.txt'%species)
