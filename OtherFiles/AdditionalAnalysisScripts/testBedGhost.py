from pybedtools import *

a = BedTool('a.bed').sort()
b = BedTool('b.bed').sort()

print b.intersect(a)

print BedTool('333-pseudomolecule_1\t46605380\t46781712',from_string=True)

print BedTool('333-psuedomolecule_1\t46739946\t46742497',from_string=True)

print BedTool('333-pseudomolecule_1\t46605380\t46781712',from_string=True).intersect(BedTool('333-pseudomolecule_1\t46739946\t46742497',from_string=True)).count()


a = BedTool('a\t5\t12',from_string=True)
b = BedTool('a\t6\t10',from_string=True)

print a.intersect(b)