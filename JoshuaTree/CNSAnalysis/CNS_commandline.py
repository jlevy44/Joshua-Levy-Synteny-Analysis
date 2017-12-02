import begin, ete3
from ete3 import PhyloTree,Tree,TreeStyle,NodeStyle
import subprocess, os
import re
from Bio import Phylo
from Bio.Phylo import Consensus as CS

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

@begin.subcommand
def runCNSAnalysis():
    subprocess.call('python runCNSAnalysis.py')

@begin.subcommand
def output_Tree(aln_output_txt,aln_file, out_fname):
    if aln_file.endswith('.phylip'):
        print 'Input must be fasta file for now'
        quit()
    elif aln_file.endswith('.fasta') or aln_file.endswith('.fa'):
        subprocess.call("awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' %s > aln.fasta"%aln_file,shell=True)
        t = PhyloTree(aln_output_txt,alignment='aln.fasta',alg_format='fasta')
    else:
        t = Tree(aln_output_txt)
    ts = TreeStyle()
    ns = NodeStyle()
    ns['size']=0
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = True
    for n in t.traverse():
        n.set_style(ns)
    #t.show(tree_style=ts)
    t.render(out_fname,tree_style = ts)

def change_of_coordinates(in_file,out_file): # FIXME USE PYTABIX FIXES NEEDED, maffilter change coordinates??
    with open(in_file,'r') as f, open(out_file,'w') as f2:
        for line in f:
            if line.startswith('#') == 0:
                break
            else:
                f2.write(line)
        #    offset = f.tell()
        #f.seek(offset) #FIXME TEST
        for line in f:
            lineList = line.split()
            lineList[1] = str(int(lineList[0].split('_')[-2]) + int(lineList[1]))
            lineList[0] = lineList[0][lineList[0].find('.')+1:[m.start(0) for m in re.finditer('_',lineList[0])][-2]]
            f2.write('\t'.join(lineList)+'\n')

def check_vcf_empty(vcf_in):
    with open(vcf_in,'r') as f:
        for line in f:
            if line.startswith('#') == 0:
                break
            offset = f.tell()
        f.seek(offset)
        if f.readline():
            return False
        else:
            return True

def filter_and_format():
    return 'in dev'

def sort_vcf(vcf_in,vcf_out):
    subprocess.call("""bgzip -c {0} > vcfs/out.vcf.gz
                (zcat vcfs/out.vcf.gz | head -300 | grep ^#;
                zcat vcfs/out.vcf.gz | grep -v ^# | sort -k1,1d -k2,2n;) \
                | bgzip -c > {1}.gz
                bcftools index {1}.gz""".format(vcf_in,vcf_out), shell=True)

@begin.subcommand
def maf2vcf(cns_config, reference_species, change_coordinates):
    change_coordinates = int(change_coordinates)
    mafFiles = [file for file in os.listdir('.') if file.endswith('.maf') and file.startswith('out_maf') == 0]
    try:
        os.mkdir('vcfs')
    except:
        pass
    with open(cns_config,'r') as f:
        master_species = parseConfigFindList("masterListSpecies",f)
    all_species = set([species.split('_')[0] for species in master_species])
    all_species_but_one = all_species - {reference_species}
    finalOutVCFFiles = []
    for i,maf in enumerate(mafFiles):
        subprocess.call("sed '/Anc/d' %s > out_maf.maf"%maf,shell=True)
        with open('maf_filter_config.bpp','w') as f:
            f.write("""
    input.file=./out_maf.maf
    input.format=Maf
    output.log=out.log
    maf.filter=\
        Subset(\
                strict=yes,\
                keep=no,\
                species=(%s),\
                remove_duplicates=yes),\
        Subset(\
                strict=yes,\
                keep=yes,\
                species=(%s),\
                remove_duplicates=yes),\
        VcfOutput(\
                file=vcfs/Out%d_all.vcf,\
                genotypes=(%s),\
                all=no,\
                reference=%s)
            """%(','.join(all_species),reference_species,i,','.join(all_species_but_one),reference_species))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
        if check_vcf_empty('vcfs/Out%d_all.vcf'%i) == 0:
        # see if there is anything in file
        # if yes, change coordinates and sort
            finalOutVCFFiles.append('vcfs/Out%d_all.vcf.gz'%i)
            change_file = 'vcfs/Out%d_all.vcf'%i
            if change_coordinates:
                change_of_coordinates('vcfs/Out%d_all.vcf'%i,'vcfs/Out%d_new_coord_unsorted.vcf'%i)
                change_file = 'vcfs/Out%d_new_coord_unsorted.vcf'%i
            sort_vcf(change_file,'vcfs/Out%d_all.vcf'%i)
            #tabix -p vcf %s.gz
        # change the coordinate system and sort the file, also remove empty vcf files, fix merge/concat
    subprocess.call('bcftools concat -O v -o vcfs/final_all.vcf %s'%(' '.join(finalOutVCFFiles)),shell=True)
    subprocess.call("sed 's/##source=Bio++/##INFO=<ID=AC>/g' vcfs/final_all.vcf > vcfs/final_all_edit.vcf && mv vcfs/final_all_edit.vcf vcfs/final_all.vcf && rm vcfs/final_all_edit.vcf",shell=True)
    sort_vcf('vcfs/final_all.vcf','vcfs/final_all_sorted.vcf')
    subprocess.call('gzip -d vcfs/final_all_sorted.vcf.gz',shell=True)
    #subprocess.call('bcftools sort -o vcfs/final_all_sorted.vcf vcfs/final_all.vcf',shell=True)
    #FIXME add sort function

@begin.subcommand
def estimate_phylogeny(cns_config, consensus_algorithm, major_cutoff, min_block_length):
    min_block_length = int(min_block_length)
    try:
        os.mkdir('./maf_trees')
    except:
        pass
    mafFiles = [file for file in os.listdir('.') if file.endswith('.maf') and file.startswith('out_maf') == 0]
    try:
        os.mkdir('vcfs')
    except:
        pass
    with open(cns_config,'r') as f:
        master_species = parseConfigFindList("masterListSpecies",f)
    all_species = set([species.split('_')[0] for species in master_species])
    fileOutTrees = []
    for i,maf in enumerate(mafFiles):
        subprocess.call("sed '/Anc/d' %s > out_maf.maf"%maf,shell=True)
        with open('maf_filter_config.bpp','w') as f:
            f.write("""
    input.file=./out_maf.maf
    input.format=Maf
    output.log=out.log
    maf.filter=\
        MinBlockLength(min_length=%d), \
        Subset(\
                strict=yes,\
                keep=no,\
                species=(%s),\
                remove_duplicates=yes),\
        DistanceEstimation(\
                method=ml,\
                model=GTR,\
                extended_names=no),\
        DistanceBasedPhylogeny(\
                method=bionj,\
                dist_mat=MLDistance),\
        OutputTrees(\
                tree=BioNJ,\
                file=./maf_trees/trees_%d.nwk,\
                compression=none,\
                strip_names=yes)
            """%(min_block_length,','.join(all_species),i))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
        fileOutTrees.append('./maf_trees/trees_%d.nwk'%i)
    with open('./maf_trees/final_trees.nwk','w') as f:
        for tree in fileOutTrees:
            with open(tree,'r') as f2:
                for line in f2:
                    if line:
                        f.write(line)
    #subprocess.call('for f in ./maf_trees/trees_*.nwk; do cat "$f"; echo "\\newline"; done > out && mv out ./maf_trees/final_trees.nwk',shell=True)
    trees = list(Phylo.parse('./maf_trees/final_trees.nwk', 'newick'))
    if consensus_algorithm == 'strict':
        tree = CS.strict_consensus(trees)
    elif consensus_algorithm == 'majority':
        tree = CS.majority_consensus(trees,float(major_cutoff))
    else:
        tree = CS.adam_consensus(trees)
    Phylo.write(tree,'./maf_trees/output_tree_consensus.nh','newick')


@begin.subcommand
def find_Tree():
    return 'in dev'

@begin.subcommand
def intersect_vcf(vcf_in, bed_regions, vcf_out):
    subprocess.call('bedtools intersect -wa -a %s -b %s > %s'%(vcf_in, bed_regions, vcf_out),shell=True)



"""
sed '/Anc/d' FastaOut1.maf > out_maf.maf

input.file=./out_maf.maf
input.format=Maf
output.log=out.log
maf.filter=\
        Subset(\
                strict=yes,\
                keep=no,\
                species=(Sbicolor,Osativa,Bdistachyon,Bstacei,ABRD,BdtwD,BhybsixD,BhybzeroD,ABRS,BdtwS,BhybsixS,BhybzeroS,Phallii,PvirgatumK,PvirgatumN,ZmaysAGPv3),\
                remove_duplicates=yes),\
        Subset(\
                strict=yes,\
                keep=yes,\
                species=(Bdistachyon),\
                remove_duplicates=yes),\
        VcfOutput(\
                file=vcfs/all.vcf,\

                genotypes=(Sbicolor,Osativa,Bstacei,ABRD,BdtwD,BhybsixD,BhybzeroD,ABRS,BdtwS,BhybsixS,BhybzeroS,Phallii,PvirgatumK,PvirgatumN,ZmaysAGPv3),\
                all=no,\
                reference=Bdistachyon)


                CAN ADD TO ABOVE with .gz extension included
                compression=gzip,\
nohup ./maffilter param=filterTest.bpp &
"""

@begin.start
def main():
    pass