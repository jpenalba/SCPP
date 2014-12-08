#!/usr/bin/home python

# This script uses the output of generate references to make a composite
# reference file to recover as much loci as possible
# The input references for this script should all be from closely related
# individuals.

import argparse
import os
import sys

comp = argparse.ArgumentParser(description = 'Generates a composite reference from the output from the first iteration of the generateReferences.py script. The idea behind this script is to generate a reference file that will contain the maximum amount of loci from the different output references')

comp.add_argument('-r', dest='ref', help='path to directory of references')

if len(sys.argv) == 1:
    comp.print_help()
    sys.exit(1)
args = comp.parse_args()


### SCRIPT STARTS HERE ###
allfiles = os.listdir(args.ref)
refs = []
for each in allfiles:
    if each.endswith('.fa'): refs.append(each)
    else: pass
complib = {}
lenlib = {}
for ref in refs:
    reffile = open(args.ref+ref, 'r')
    for lines in reffile:
        info = lines.strip().split()
        if info[0].startswith('>'): locus = info[0]
        else: 
            seq = info[0]
            if locus not in complib: 
                if locus not in lenlib:
                    lenlib[locus] = len(seq)
                    complib[locus] = seq
                elif locus in lenlib:
                    if len(seq) > lenlib[locus]: 
                        complib[locus] = seq
                        lenlib[locus] = len(seq)
                    else: pass
            else: pass
    reffile.close()
compout = open(args.ref+'composite.fa', 'w')
for loc in complib:
    compout.write('%s\n%s\n' % (loc, complib[loc]))
compout.close()
    
