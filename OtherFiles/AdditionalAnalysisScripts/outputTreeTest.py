from ete3 import Tree,TreeStyle,NodeStyle
import os
for tree in [file for file in os.listdir('.') if file and file.endswith('_phyml_tree.txt')]:
    #print tree
    try:
        with open( tree, 'r') as f:  # FIXME boot_trees verus phyml_tree
            t = Tree(open(tree,'r').read())
            ts = TreeStyle()
            ns = NodeStyle()
            ns['size']=0
            ts.show_leaf_name = True
            ts.show_branch_length = False
            ts.show_branch_support = True
            for n in t.traverse():
                n.set_style(ns)
            #t.show(tree_style=ts)
            t.render( os.getcwd()+'/'+tree.replace('_phyml_tree.txt', '.png'),tree_style = ts)
    except:
        pass