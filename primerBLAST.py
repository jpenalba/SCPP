#!/usr/bin/home python

# Primer BLAST

# This script is written to recover assemblies from SCPP bait where only the primer sequences are available.
#

########################
# SETTING UP ARGUMENTS #
########################

import argparse
import os
import sys

pb = argparse.ArgumentParser(description='Runs BLAST for primer sequences to look for associated assemblies. This is for SCPP targets where only the primer sequences are available')
pb.add_argument('-p', dest='prim', type=str, help='Fasta file of primer sequences')
pb.add_argument('-a', dest='assembly', type=str, help='Path to final assemblies')
pb.add_argument('-o', dest='out', type=str, help='Output directory')

if len(sys.argv) == 1:
    pb.print_help()
    sys.exit(1)
args = pb.parse_args()

try: os.mkdir(args.out)
except OSError: pass

########################
# SETTING UP FUNCTIONS #
########################

def fa_parser(fasta):
    fa_file = open(fasta, 'r')
    fa_dic = {}
    for lines in fa_file:
        info = lines.strip()
        if info.startswith('>'): 
            name = info[1:]
            fa_dic[name] = ''
        else: fa_dic[name] = fa_dic[name] + info
    return fa_dic

######################
# MAKING INPUT FILES #
######################
'''
print '==> SETTING UP BLAST INPUT...\n'
primer_dict = fa_parser(args.prim)
primer_out = open(os.path.join(args.out, "primers.in"), 'w')

loci = []
for primer in primer_dict:
    locus = primer.split('-')[0]
    if locus not in loci:
        loci.append(locus)
        primer_out.write('>%s\n%s%s%s\n' % (locus, primer_dict[locus+'-F'], 'N'*20, primer_dict[locus+'-R']))
primer_out.close()


cat_out = open(os.path.join(args.out, "blastIN.fa"), 'w')
all_files = os.listdir(args.assembly)
for files in all_files:
    if files.endswith('.fa.final'):
        path = args.assembly+files
        lib_name = files.split('.')[0]
        for contig in fa_parser(path):
            cat_out.write('>%s_%s\n' % (lib_name, contig))
            cat_out.write(fa_parser(path)[contig]+'\n')
cat_out.close()
'''
#################
# RUNNING BLAST #
#################

print '==> RUNNING PRIMER BLAST...\n'
os.system('makeblastdb -in %s -dbtype nucl -out %s' % (os.path.join(args.out, 'blastIN.fa'), os.path.join(args.out, 'blastIN')))
os.system('blastn -task blastn-short -db %s -query %s -outfmt 6 -out %s' % (os.path.join(args.out, 'blastIN'), os.path.join(args.out, 'primers.in'), os.path.join(args.out, 'primerBLAST.out')))

