
"""Sbicolor_313_Sb
Osativa_323_Os
Bdistachyon_314_BdD
Bstacei_316_Bs
ABRD_001_BhD
BdtwD_003_BhtD
BdTRgD_005_BhTgD
PobD_007_BhPD
BhybsixD_009_BhSD
BhybzeroD_011_BhZD
ABRS_002_BhS
BdtwS_004_BhtS
BdTRgS_006_BhTgS
PobS_008_BhPS
BhybsixS_010_BhSS
BhybzeroS_012_BhZS
BdTRi_348_BdTi
BdEighteen_362_BdE
BdTRc_354_BdTc"""
if 0:
    with open('prot_dict','r') as f:
        species = [line.split('     ')[1].strip('\n') for line in f if line]

    masterSpeciesList = ['_'.join([species_Name,species_Name.replace('Bdist',''),'|'+species_Name.replace('ist','')+'~']) for species_Name in species]
    shortName = set(['|'+species_Name.replace('ist','')+'~' for species_Name in species])
    intra = set([name for name in shortName if '460' in name or '314' in name])
    extra = shortName - intra
from ete2 import Tree

t = Tree('grassnewick_091117temp7.nh')
species = [node.name for node in t.traverse('postorder') if node.name]

with open('prot_dict_hyb','w') as f:
    f.write('\n'.join(['%s     %s'%(s.replace('Bhyb','').replace('Bdist',''),s) for s in species]))

masterSpeciesList = ['_'.join([species_Name,species_Name.replace('Bdist','').replace('Bhyb',''),'|'+species_Name.replace('ist','').replace('yb','')+'~']) for species_Name in species]
shortName = set(['|'+species_Name.replace('ist','').replace('yb','')+'~' for species_Name in species])
intra = set([name for name in shortName if '001' in name or '003' in name or '011' in name or '009' in name])
extra = shortName - intra









configText = """root_folder =
checkValidity = 0
pathPython = /global/u2/j/jlevy/python_modules
pathSystem = /global/u2/j/jlevy/virtualBin
conservedFastaPath = /global/projectb/scratch/jlevy/Spring2017/Synteny/BhybridumFinal/Dist110HybD/CNSAnalysis/ConservedFasta/
pickleName = MAFWork.p
ratioCopy = 0
pickleSkip = 0
fasta2phylip = 1
PhyML = 1
bootstrap = 1
treeOut = 1
treeFile = 0
short_fasta_output_name = 1
masterListSpecies =
%s
Stop
intragenus =
%s
Stop
intergenus =
%s
Stop
subgenome =
xxxxx
Stop"""%tuple(map(lambda x: '\n'.join(x),[masterSpeciesList,intra,extra]))

with open('configCNSAnalysis_hyb.txt','w') as f:
    f.write(configText)