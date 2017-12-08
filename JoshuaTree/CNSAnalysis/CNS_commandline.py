import begin, ete3
from ete3 import PhyloTree,Tree,TreeStyle,NodeStyle,EvolTree
import subprocess, os
import re
from Bio import Phylo
from Bio.Phylo import Consensus as CS
import dill as pickle
import numpy as np
import sys
import pandas as pd
from collections import Counter

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
def maf2vcf(cns_config, reference_species, change_coordinates, out_all_species, overlaps):
    change_coordinates = int(change_coordinates)
    out_all_species = int(out_all_species)
    overlaps = int(overlaps)
    mafFiles = [file for file in os.listdir('.') if file.endswith('.maf') and file.startswith('FastaOut')]
    try:
        os.mkdir('vcfs')
    except:
        pass
    if cns_config.endswith('.txt'):
        with open(cns_config,'r') as f:
            master_species = parseConfigFindList("masterListSpecies",f)
    else:
        master_species = cns_config.split(',')
    all_species = set([species.split('_')[0] for species in master_species])
    all_species_but_one = all_species - {reference_species}
    if out_all_species:
        species = all_species
    else:
        species = all_species_but_one
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
            """%(','.join(all_species),reference_species,i,','.join(species),reference_species))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
        if check_vcf_empty('vcfs/Out%d_all.vcf'%i) == 0:
        # see if there is anything in file
        # if yes, change coordinates and sort
            finalOutVCFFiles.append('vcfs/Out%d_all.vcf.gz'%i)
            subprocess.call("sed 's/##source=Bio++/##INFO=<ID=AC>/g' %s > vcfs/temp.vcf && mv vcfs/temp.vcf %s && rm vcfs/temp.vcf"%('vcfs/Out%d_all.vcf'%i,'vcfs/Out%d_all.vcf'%i),shell=True)
            change_file = 'vcfs/Out%d_all.vcf'%i
            if change_coordinates:
                change_of_coordinates('vcfs/Out%d_all.vcf'%i,'vcfs/Out%d_new_coord_unsorted.vcf'%i)
                change_file = 'vcfs/Out%d_new_coord_unsorted.vcf'%i
            sort_vcf(change_file,'vcfs/Out%d_all.vcf'%i)
            #tabix -p vcf %s.gz
        # change the coordinate system and sort the file, also remove empty vcf files, fix merge/concat
    # FIXME problem with overlaps, --allow overlaps works with some analyses but not others, may want to just throw in contig list!!! below
    subprocess.call('bcftools concat%s -O v -o vcfs/final_all.vcf %s'%(' --allow-overlaps' if overlaps else '',' '.join(finalOutVCFFiles)),shell=True)
    #subprocess.call("sed 's/##source=Bio++/##INFO=<ID=AC>/g' vcfs/final_all.vcf > vcfs/final_all_edit.vcf && mv vcfs/final_all_edit.vcf vcfs/final_all.vcf && rm vcfs/final_all_edit.vcf",shell=True)
    sort_vcf('vcfs/final_all.vcf','vcfs/final_all_sorted.vcf')
    subprocess.call('rm vcfs/final_all_sorted.vcf && gzip -d vcfs/final_all_sorted.vcf.gz',shell=True)
    #subprocess.call('bcftools sort -o vcfs/final_all_sorted.vcf vcfs/final_all.vcf',shell=True)
    #FIXME add sort function

@begin.subcommand
def estimate_phylogeny(cns_config, consensus_algorithm, major_cutoff, min_block_length, concat_size, consensus_tree):
    min_block_length = int(min_block_length)
    consensus_tree = int(consensus_tree)
    concat_size = int(concat_size)
    try:
        os.mkdir('./maf_trees')
    except:
        pass
    mafFiles = [file for file in os.listdir('.') if file.endswith('.maf') and (file.startswith('merged') + file.startswith('out_maf') == 0)]
    try:
        os.mkdir('vcfs')
    except:
        pass
    with open(cns_config,'r') as f:
        master_species = parseConfigFindList("masterListSpecies",f)
    all_species = set([species.split('_')[0] for species in master_species])
    if consensus_tree:
        fileOutTrees = []
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
            MaskFilter(species=(%s)),\
            MinBlockLength(min_length=%d), \
            Merge(species=(%s)),\
            Concatenate(minimum_size=%d),\
            RemoveEmptySequences(),\
            DistanceEstimation(\
                    method=ml,\
                    model=GTR,\
                    gap_option=no_gap,\
                    parameter_estimation=initial,\
                    gaps_as_unresolved=yes,\
                    unresolved_as_gap=yes,\
                    extended_names=no),\
            DistanceBasedPhylogeny(\
                    method=bionj,\
                    dist_mat=MLDistance),\
            OutputTrees(\
                    tree=BioNJ,\
                    file=./maf_trees/trees_%d.nwk,\
                    compression=none,\
                    strip_names=yes)
                """%(','.join(all_species),','.join(all_species),min_block_length,','.join(all_species),concat_size,i))
            subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
            with open('./maf_trees/trees_%d.nwk'%i,'r') as f:
                if f.readline():
                    fileOutTrees.append('./maf_trees/trees_%d.nwk'%i)
        with open('./maf_trees/final_trees.nwk','w') as f: #FIXME maf filter throws exception while calculating distance matrix at empty sites
            for tree in fileOutTrees:
                with open(tree,'r') as f2:
                    for line in f2:
                        if line and line.count('(') == line.count(')'):
                            f.write(line)
        #subprocess.call('for f in ./maf_trees/trees_*.nwk; do cat "$f"; echo "\\newline"; done > out && mv out ./maf_trees/final_trees.nwk',shell=True)
        trees_gen = Phylo.parse('./maf_trees/final_trees.nwk', 'newick')
        trees = []
        with open('./maf_trees/final_trees.nwk','r') as f:
            for i in range(len(f.readlines())):
                try:
                    trees.append(trees_gen.next())
                except:
                    pass
        if consensus_algorithm == 'strict':
            tree = CS.strict_consensus(trees)
        elif consensus_algorithm == 'majority':
            tree = CS.majority_consensus(trees,float(major_cutoff))
        else:
            tree = CS.adam_consensus(trees)
        Phylo.write(tree,'./maf_trees/output_tree_consensus.nh','newick')
    else:
        subprocess.call('rm merged.maf',shell=True)
        for file in mafFiles:
            subprocess.call("sed -e '/Anc/d;/#/d' %s >> merged.maf"%file,shell=True)
        with open('merged.maf','r') as f: #FIXME
            txt = f.read().replace('\n\n\n','\n\n')
        with open('merged.maf','w') as f:
            f.write(txt)
        del txt
        #subprocess.call("sed '/#/d' merged.maf > out_maf.maf",shell=True)
        with open('maf_filter_config.bpp','w') as f:
                f.write("""
        input.file=./merged.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\
            Subset(\
                    strict=yes,\
                    keep=no,\
                    species=(%s),\
                    remove_duplicates=yes),\
            MaskFilter(species=(%s)),\
            MinBlockLength(min_length=%d), \
            Merge(species=(%s)),\
            Concatenate(minimum_size=1000000000000),\
            DistanceEstimation(\
                    method=ml,\
                    model=GTR,\
                    gap_option=no_gap,\
                    parameter_estimation=initial,\
                    gaps_as_unresolved=yes,\
                    unresolved_as_gaps=yes,\
                    extended_names=no),\
            DistanceBasedPhylogeny(\
                    method=bionj,\
                    dist_mat=MLDistance),\
            OutputTrees(\
                    tree=BioNJ,\
                    file=./maf_trees/output_tree_consensus.nh,\
                    compression=none,\
                    strip_names=yes)
                """%(','.join(all_species),','.join(all_species),min_block_length,','.join(all_species)))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)

def maf_change_coordinates(segment,ref_species):
    aln_lines = segment.splitlines()
    for i,line in enumerate(aln_lines):
        if line.startswith('s'):
            lineList = line.split()
            orientation = lineList[4]
            lineList2 = lineList[1].split('.')
            lineList3 = lineList2[-1].split('_')[-2:]
            lineList2[2] = lineList2[2].replace('_'+'_'.join(lineList3),'')
            if orientation == '-':
                lineList[2] = str(int(lineList3[-1])-int(lineList[2]))#-int(lineList[3]))
            else:
                lineList[2] = str(int(lineList3[-2]) + int(lineList[2]))
            lineList[1] = '.'.join(lineList2[::2])
            aln_lines[i] = '\t'.join(lineList)
            if lineList2[0] == ref_species:
                chrom = lineList2[2]
                position = int(lineList[2])
    return chrom,position,'\n'.join(filter(None,aln_lines))+'\n\n',

@begin.subcommand
def selective_pressure_statistics(cns_config,reference_species, min_block_length, dist_max, window_size, root_species): # FIXME add ingroup outgroup options
    min_block_length = int(min_block_length)
    window_size = int(window_size)
    dist_max = int(dist_max)
    try:
        os.mkdir('./maf_trees')
    except:
        pass
    if cns_config.endswith('.txt'):
        with open(cns_config,'r') as f:
            master_species = parseConfigFindList("masterListSpecies",f)
    else:
        master_species = cns_config.split(',')
    master_species = set([species.split('_')[0] for species in master_species])
    print master_species
    mafFiles = [file for file in os.listdir('.') if file.startswith('FastaOut') and file.endswith('.maf')]
    if 0: #FIXME tests
        subprocess.call('rm merged.maf',shell=True)
        for file in mafFiles:
            subprocess.call("sed -e '/Anc/d;/#/d' %s >> merged.maf"%file,shell=True)
        with open('merged.maf','r') as f: #FIXME
            txt = f.read().replace('\n\n\n','\n\n')
        with open('merged.maf','w') as f:
            f.write(txt)
        del txt
        with open('maf_filter_config.bpp','w') as f:
                    f.write("""
            input.file=./merged.maf
            input.format=Maf
            output.log=out.log
            maf.filter=\
                Subset(\
                        strict=yes,\
                        keep=no,\
                        species=(%s),\
                        remove_duplicates=yes),\
                MaskFilter(species=(%s)),\
                MinBlockLength(min_length=%d),\
                Output(\
                        file=merged.filtered.maf,\
                        compression=none,\
                        mask=yes)
                    """%(','.join(master_species),','.join(master_species),min_block_length))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
        # FIXME sort them by coordinates then output, and keep order of species
        maf_sort_structure = []
        with open('merged.filtered.maf','r') as f:
            for segment in f.read().split('\n\n'): # FIXME can turn this to generator in future, especially with large maf size
                if segment:
                    maf_sort_structure.append(maf_change_coordinates(segment,reference_species))
        maf_sort_structure = pd.DataFrame(maf_sort_structure).sort_values([0,1])
        with open('merged.filtered.new_coords.maf','w') as f2:
            for seg in maf_sort_structure.itertuples():
                f2.write(seg[3])
        del maf_sort_structure #FIXME may be bad when maf file is hundreds of gb large
        with open('maf_filter_config.bpp','w') as f:
            f.write("""
            input.file=./merged.filtered.new_coords.maf
            input.format=Maf
            output.log=out.log
            maf.filter=\
                Merge(\
                        species=(%s),\
                        dist_max=%d),\
                WindowSplit(\
                        preferred_size=%d,\
                        align=adjust,\
                        keep_small_blocks=yes),\
                RemoveEmptySequences(),\
                Output(\
                        file=merged.filtered.new_coords.syntenic.maf,\
                        compression=none,\
                        mask=no)
                    """%(reference_species,dist_max,window_size))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
    maf_configs = []
    maf_configs.append("""
        input.file=./merged.filtered.new_coords.syntenic.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\
            OutputCoordinates(\
                    file=coordinates.txt,\
                    compression=none,\
                    species=(%s),\
                    output_src_size=yes)
            """%reference_species)
    maf_configs.append("""
        input.file=./merged.filtered.new_coords.syntenic.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\
            DistanceEstimation(\
                    method=ml,\
                    model=GTR,\
                    gap_option=no_gap,\
                    parameter_estimation=initial,\
                    gaps_as_unresolved=yes,\
                    unresolved_as_gaps=yes,\
                    extended_names=no),\
            DistanceBasedPhylogeny(\
                    method=bionj,\
                    dist_mat=MLDistance),\
            NewOutgroup(\
                    tree_input=BioNJ,\
                    tree_output=BioNJ,\
                    outgroup=%s),\
            SequenceStatistics(\
                statistics=(\
                    BlockCounts,\
                    AlnScore,\
                    BlockLength,\
                    CountClusters(\
                        tree=bionj,\
                        threshold=0.001),\
                    SiteFrequencySpectrum(\
                        bounds=(-0.5,0.5,1.5,2.5,3.5,4.5),\
                        ingroup=(%s),\
                        outgroup=%s),\
                    SiteStatistics(\
                        species=%s),\
                        ),\
                    DiversityStatistics(ingroup=(%s)),\
                ref_species=%s,\
                file=data.statistics.csv)"""%(root_species,','.join(master_species),root_species,reference_species,','.join(master_species),reference_species))
    maf_configs.append("""
        input.file=./merged.filtered.new_coords.syntenic.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\
            DistanceEstimation(\
                    method=ml,\
                    model=GTR,\
                    gap_option=no_gap,\
                    parameter_estimation=initial,\
                    gaps_as_unresolved=yes,\
                    unresolved_as_gaps=yes,\
                    extended_names=no),\
        DistanceBasedPhylogeny(\
                    method=bionj,\
                    dist_mat=MLDistance),\
        NewOutgroup(\
                    tree_input=BioNJ,\
                    tree_output=BioNJ,\
                    outgroup=%s),\
        OutputTrees(\
                    tree=BioNJ,\
                    file=./maf_trees/output_trees.nwk,\
                    compression=none,\
                    strip_names=yes)"""%(root_species))
    # FIXME master species on 4th entry
    for maf_config in maf_configs:
        with open('maf_filter_config.bpp','w') as f:
            f.write(maf_config)
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
    # output will be trees and statistics, use pandas and grab tree lengths to pandas, encode array, then add polymorphism density from vcf possibly, and then pca the results and cluster the regions
    # neural network model???


@begin.subcommand
def evaluate_selective_pressure(maf_structure_pickle, neutral_tree_nwk):
    # FIXME codeml only works in protein coding regions https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1470900/ need to add for CNS!!! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3725466/
    # FIXME branch specific rate changes??
    # FIXME also degree of conservation when looking at selective pressure density
    # FIXME add new metric? Conserved pressure metric, conserved ratio (# conserved/total conserved weighted) * selective pressure density?? for different elements??
    codeMLPath = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin'
    try:
        os.mkdir('./selectivePressureTests')
    except:
        pass
    maf = pickle.load(open('MAFWork.p','rb'))
    onesCondition = maf[3].keys()[np.argmax([sum(val.values()) if type(val) != type([]) else 0 for val in maf[3].values()])]
    for maseg in maf[1][onesCondition]['CS']:
        model_tree = EvolTree(open(neutral_tree_nwk,'rb').read(),binpath=codeMLPath) #FIXME is this the right tree to be using???
        aln_lines = maf[1][onesCondition]['CS'][maseg].splitlines()
        coord_info = {}
        for i in range(len(aln_lines))[::2]:
            coord_info[aln_lines[i].find('>')+1:aln_lines[i].find('.')] = aln_lines[i][aln_lines[i].find('.'):] #FIXME output these coordinates to bed file
            aln_lines[i] = aln_lines[i][:aln_lines[i].find('.')]
            alignment = '\n'.join(aln_lines)
        model_tree.link_to_alignment(alignment)
        model_tree.workdir = './selectivePressureTests'
        model_tree.run_model('M1')
        model_tree.run_model('M2')
        pval = model_tree.get_most_likely('M2','M1')
        if pval < 0.05:
            for s in range(len(model2.sites['BEB']['aa'])):
                print 'positively selected site %s at position: %s, with probability: %s' % (model2.sites['BEB']['aa'][s], s+1, model2.sites['BEB']['p2'][s])
        # FIXME first fix codeml, then print out selective pressure density after removing --, append that to masegment, bed file with coordinates, maSeg, and selective pressure density, or could just output one-base coordinate and calculate density through bedtools
        else:
            print 'fail to reject null hypothesis'
        """
        best_model = None
        best_lnl = float('-inf')
        for starting_omega in [0.2, 0.7, 1.2]: # FIXME optimize and choose parameters
            model_tree.run_model('b_free.'+str(starting_omega))
            current_model = model_treetree.get_evol_model('b_free.'+str(starting_omega))
            if current_model.lnL > best_lnl:
                best_lnl = current_model.lnL
                best_model = current_model
        """
        #print best_model
        tree.render(maseg+'.png')
        break # FIXME remove after testing
        "FIXME UNDER DEVELOPMENT"

@begin.subcommand
def reroot_Tree(tree_in,root_species,tree_out):
    t = PhyloTree(open(tree_in,'r').read())
    t.set_outgroup(root_species)
    t.write(tree_out)

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

        MaskFilter(\
                species=(%s),\
                window.size=10,\
                window.step=1,\
                max.masked=3),\
        MinBlockLength(min_length=%d), \
            XFullGap(species=(%s)),\
            MinBlockLength(min_length=25),\

                CAN ADD TO ABOVE with .gz extension included
                compression=gzip,\
nohup ./maffilter param=filterTest.bpp &
"""

@begin.start
def main():
    pass