#!/usr/bin/home python

###############################################
#                                             #
#   Clustering Reciprocal BLAST Hit (cRBH)    # 
#           Reference Generator               #
#                                             #
#  Author: Joshua Penalba                     #
#  Date Written: 15 Aug 2014                  #
#  Date Last Modified: 29 Sep 2014            #
#                                             #
###############################################

###############################
# ESTABLISING ARGUMENT PARSER #
###############################
import argparse
import os
import sys

rbh_args = argparse.ArgumentParser(description='Uses a clustering Reciprocal BLAST hit algorithm to recover assembled contigs that corresponds to each target for a SCPP project')

rbh_args.add_argument('-r', dest='ref', type=str, help='fasta file of target sequences', required=True)
rbh_args.add_argument('-a', dest='assembly', type=str, help='directory of final assemblies (ends in .fa.final)')
rbh_args.add_argument('-o', dest='out', type=str, help='output directory')
rbh_args.add_argument('-e', dest='evalue', type=str, help='e-value for BLAST search [1e-5]', default='1e-5')

if len(sys.argv) == 1:
    rbh_args.print_help()
    sys.exit(1)
args = rbh_args.parse_args()

if not args.out.endswith('/'): args.out = args.out + '/'
if not args.assembly.endswith('/'): args.assembly = args.assembly + '/'


try: os.mkdir(args.out)
except OSError: pass

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

def rev_comp(seq):
    seq = seq[::-1]
    seq = seq.replace('A','X').replace('T','A').replace('X','T')
    seq = seq.replace('C','X').replace('G','C').replace('X','G')
    return seq

#######################
# CONCATENATE CONTIGS #
#######################
# Renames each of the contigs with the library name in the beginning.
# Concatenates each of those sets of contigs; including the reference.
print '==> SETTING UP BLAST FILE ...\n'
cat_out = open(args.out+'blastIN.fa', 'w')

for locus in fa_parser(args.ref):
    cat_out.write('>REF_%s\n' % locus)
    cat_out.write(fa_parser(args.ref)[locus]+'\n')

all_files = os.listdir(args.assembly)
for files in all_files:
    if files.endswith('.fa.final'):
        path = args.assembly+files
        lib_name = files.split('.')[0]
        for contig in fa_parser(path):
            cat_out.write('>%s_%s\n' % (lib_name, contig))
            cat_out.write(fa_parser(path)[contig]+'\n')

cat_out.close()

###############
# RECIP BLAST #
###############
# Calls BLAST+ on the input file made earlier.

print '==> RUNNING BLAST ...'
os.system('makeblastdb -in %s -dbtype nucl -out %s' % (args.out+'blastIN.fa', args.out+'blastIN'))
os.system('blastn -db %s -query %s -outfmt 6 -out %sblast.out' % (args.out+'blastIN', args.out+'blastIN.fa', args.out))


##########
# FILTER #
##########
# Creates a smaller file with necessary data.
# Filters overlap of < 100bp.
# Filters self matches.
print ''
print '==> FILTERING RECIPROCAL BLASTS ...\n'
blast_out = open(args.out+'blast.out', 'r')
filt_out = open(args.out+'blast_filt.out', 'w')
for lines in blast_out:
    info = lines.strip().split()
    (sp1, sp2) = (info[0].split('_')[0], info[1].split('_')[0])
    if sp1 != sp2 and int(info[3]) > 100:
        if sp1 == 'REF': sp1 = info[0]
        if sp2 == 'REF': sp2 = info[1]
        rev = 'N'
        if int(info[8]) > int(info[9]): rev = 'Y'
        filt_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sp1, info[0], sp2, info[1], info[3], info[11], rev))

blast_out.close()
filt_out.close()

###########
# CLUSTER #    
###########

# Keep only contigs with reciprocal hits
hits = {}
filt_out = open(args.out+'blast_filt.out','r')
pairs = set()
keep_pairs = set()
for lines in filt_out:
    info = lines.strip().split()
    pairs.add((info[1],info[3]))
for pair in pairs:
    recip = (pair[1],pair[0])
    if recip in pairs: 
        keep_pairs.add(pair)    
filt_out.close()

#SETS UP TWO DICTIONARIES FOR RECIPROCAL HITS
#############################################
ref_dict_main = {} # NOTE: locus -> lib -> contigs -> contig info
lib_dict_main = {} # NOTE: lib1 -> lib2 -> contig_lib1 -> contig_lib2 -> contig info
filt_out = open(args.out+'blast_filt.out','r')
for lines in filt_out:
    info = lines.strip().split()
    if (info[1],info[3]) not in keep_pairs: continue
    else: 
        (sp1, sp2, contig1, contig2) = (info[0], info[2], info[1], info[3])
        if sp2.startswith('REF'): continue
        if sp1.startswith('REF'):
            if sp1 not in ref_dict_main: ref_dict_main[sp1] = {}
            if sp2 not in ref_dict_main[sp1]: 
                ref_dict_main[sp1][sp2] = {contig2:info}
            elif sp2 in ref_dict_main[sp1]:
                if contig2 in ref_dict_main[sp1][sp2]:
                    if int(info[4]) > int(ref_dict_main[sp1][sp2][contig2][4]): ref_dict_main[sp1][sp2][contig2] = info
                elif contig2 not in ref_dict_main[sp1][sp2]:
                    ref_dict_main[sp1][sp2][contig2] = info
        else: 
            if sp1 not in lib_dict_main: 
                lib_dict_main[sp1] = {}
            if sp2 not in lib_dict_main[sp1]:
                lib_dict_main[sp1][sp2] = {}
            elif sp2 in lib_dict_main[sp1]:
                if contig1 not in lib_dict_main[sp1][sp2]:
                    lib_dict_main[sp1][sp2][contig1] = {contig2:info}
                elif contig1 in lib_dict_main[sp1][sp2]:
                    lib_dict_main[sp1][sp2][contig1][contig2] = info

def filter_contigs(ref_dict, lib_dict):
    for locus in ref_dict:
        for lib in ref_dict[locus]:
            if type(ref_dict[locus][lib]) is dict:
                if len(ref_dict[locus][lib].keys()) > 1:
                    (hi_hits, lib_key, ref_bitscore, hit_dict) = (0, ref_dict[locus][lib], 0, {})
                    #lib_key = ref_dict[locus][lib]
                    #ref_bitscore = 0
                    #hit_dict = {}
                    for contig in lib_key:
                        if float(lib_key[contig][5]) > ref_bitscore:
                            ref_bitscore = float(lib_key[contig][5])
                            best_ref_hit = contig
                        hits = 0
                        for other_lib in ref_dict[locus]: #If the other library is in the reference dictionary for that locus
                            if other_lib != lib: # If the other lib is not the same as the lib we are looking at 
                                for other_contig in ref_dict[locus][other_lib]: # For each contig for a different library for that locus
                                    lib_bitscore = 0
                                    if other_contig in lib_dict[other_lib][lib]:
                                        for each_contig in lib_dict[other_lib][lib][other_contig]:
                                            if float(lib_dict[other_lib][lib][other_contig][each_contig][5]) > lib_bitscore:
                                                best_lib_hit = each_contig
                                                lib_bitscore = float(lib_dict[other_lib][lib][other_contig][each_contig][5])
                                    else: pass
                                    if best_lib_hit == contig: hits += 1
                        hit_dict[contig] = hits
                    hit_dict[best_ref_hit] += 2
                    top_hit = 0
                    for each_hit in hit_dict:
                        if hit_dict[each_hit] > top_hit: top_hit = hit_dict[each_hit]
                    for each_hit in hit_dict:
                        if hit_dict[each_hit] < top_hit: del ref_dict[locus][lib][each_hit]
            else: 
                if len(ref_dict[locus][lib]) > 1:
                    hit_dict = {}
                    for contig in ref_dict[locus][lib]:
                        hits = 0
                        for other_lib in ref_dict[locus]:
                            if other_lib != lib:
                                for other_contig in ref_dict[locus][other_lib]:
                                    lib_bitscore = 0
                                    if other_contig in lib_dict[other_lib][lib]:
                                        for each_contig in lib_dict[other_lib][lib][other_contig]:
                                            if float(lib_dict[other_lib][lib][other_contig][each_contig][5]) > lib_bitscore:
                                                best_lib_hit = each_contig
                                                lib_bitscore = float(lib_dict[other_lib][lib][other_contig][each_contig][5])
                                    else: pass
                                    if best_lib_hit == contig: hits += 1
                        hit_dict[contig] = hits
                    top_hit = 0
                    for each_hit in hit_dict:
                        if hit_dict[each_hit] > top_hit: top_hit = hit_dict[each_hit]
                    for each_hit in hit_dict:
                        if hit_dict[each_hit] < top_hit: del ref_dict[locus][lib][ref_dict[locus][lib].index(each_hit)]
    return ref_dict

def cluster_contigs(ref_dict, lib_dict):
    libraries = lib_dict.keys()
    for locus in ref_dict:
        no_contig = []
        yes_contig = []
        for lib in libraries:
            if lib in ref_dict[locus]: yes_contig.append(lib)
            else: no_contig.append(lib)
        for Ncontig_lib in no_contig:
            hit_dict = {}
            for Ycontig_lib in yes_contig: #For each library with a representative contig
                for ref_contig in ref_dict[locus][Ycontig_lib]: #For each representative contig in the library
                    if ref_contig in lib_dict[Ycontig_lib][Ncontig_lib]: #If that contig exists in the lib_dict heirarchy
                        alt_contig_dict = lib_dict[Ycontig_lib][Ncontig_lib][ref_contig]
                        alt_bitscore = 0
                        for alt_contig in alt_contig_dict: #For all the alternative contigs, look for the highest bitscore
                            if alt_contig_dict[alt_contig][5] > alt_bitscore:
                                alt_bitscore = alt_contig_dict[alt_contig][5]
                                best_contig = alt_contig
                        if best_contig not in hit_dict: 
                            hit_dict[best_contig] = 1
                        else: hit_dict[best_contig] += 1
            if len(hit_dict.keys()) > 0:
                highest = 0
                for each_hit in hit_dict:
                    if hit_dict[each_hit] > highest: highest = hit_dict[each_hit]
                for each_hit in hit_dict:
                    if hit_dict[each_hit] == highest: 
                        if Ncontig_lib not in ref_dict[locus]: ref_dict[locus][Ncontig_lib] = [each_hit]
                    elif Ncontig_lib in ref_dict[locus]: ref_dict[locus][Ncontig_lib].append(each_hit)
    return ref_dict

def count_libs(ref_dict):
    total_libs = 0
    for locus in ref_dict:
        total_libs += len(ref_dict[locus].keys())
    return total_libs

def count_contigs(ref_dict):
    contig_count = 0
    for locus in ref_dict:
        for lib in ref_dict[locus]:
            if type(ref_dict[locus][lib]) is dict: 
                contig_count += len(ref_dict[locus][lib].keys())
            elif type(ref_dict[locus][lib]) is list:
                contig_count += len(ref_dict[locus][lib])
    return contig_count

# Call clustering
##################
print ''
print '==> STARTING CONTIG RECOVERY ...\n'
print 'Initial numbers - Libraries: %i, Contigs: %i\n' % (count_libs(ref_dict_main), count_contigs(ref_dict_main))
ref_dict_main = filter_contigs(ref_dict_main, lib_dict_main)
print 'Post filtering numbers - Libraries: %i, Contigs: %i\n' % (count_libs(ref_dict_main), count_contigs(ref_dict_main))

lib_limit = 0
cluster_round = 0
while count_libs(ref_dict_main) > lib_limit:
    cluster_round += 1
    lib_limit = count_libs(ref_dict_main)
    print '... Cluster round %i ...' % cluster_round
    ref_dict_main = cluster_contigs(ref_dict_main, lib_dict_main)
    print 'Post cluster recovery - Libraries: %i, Contigs: %i\n' % (count_libs(ref_dict_main), count_contigs(ref_dict_main))

ref_dict_main = filter_contigs(ref_dict_main, lib_dict_main)
print 'Post final filtering - Libraries: %i, Contigs: %i\n' % (count_libs(ref_dict_main), count_contigs(ref_dict_main))

# Check reverse compliment
##########################
print '==> CREATING REV-COMP DICT...'
print 'make sure it doesn\'t get stuck here\n'
rev_dict = {}
other_dict = {}
other_list = []

for locus in ref_dict_main:
    for lib in ref_dict_main[locus]:
        if type(ref_dict_main[locus][lib]) is dict: 
            for contig in ref_dict_main[locus][lib]: 
                rev_dict[contig] = ref_dict_main[locus][lib][contig][6]
        elif type(ref_dict_main[locus][lib]) is list:
            for contig in ref_dict_main[locus][lib]: 
                other_list.append(contig)
print len(rev_dict.keys())
print len(other_list)
prev_len = len(other_list)+1
while len(other_list) != prev_len:
    working_list = other_list
    prev_len = len(other_list)
    for contig in working_list:
        lib = contig.split('_')[0]
        for other_contig in rev_dict:
            other_lib = other_contig.split('_')[0]
            try:
                if lib_dict_main[lib][other_lib][contig][other_contig]:
                    print 'YES'
                    contig_info = lib_dict_main[lib][other_lib][contig][other_contig]
                    if rev_dict[other_contig] == 'Y':
                        if contig_info[6] == 'Y': other_dict[contig] = 'N'
                        elif contig_info[6] == 'N': other_dict[contig] = 'Y'
                    elif rev_dict[other_contig] == 'N':
                        other_dict[contig] = contig_info[6]
                    del other_list[other_list.index(contig)]
                    break
            except KeyError: pass
    rev_dict.update(other_dict.items())

######################
# PRINT OUTPUT FILES #
######################

print '==> PRINTING OUTPUT FILES ...\n'
# Reformat dictionary
#####################
ref_by_samp = {}
for locus in ref_dict_main:
    for lib in ref_dict_main[locus]:
        if type(ref_dict_main[locus][lib]) is dict: contig_list = ref_dict_main[locus][lib].keys()
        elif type(ref_dict_main[locus][lib]) is list: contig_list = ref_dict_main[locus][lib]
        if lib not in ref_by_samp: 
            ref_by_samp[lib] = {locus:contig_list}
        elif lib in ref_by_samp:
            ref_by_samp[lib][locus] = contig_list

# Make summary file
###################
sum_file = open(args.out+'cRBH_summary.txt', 'w')
all_loci = fa_parser(args.ref).keys()
all_loci.sort()
sum_file.write('\t%s\n' % ('\t').join(all_loci))
lib_keys = ref_by_samp.keys()
lib_keys.sort()
for lib in lib_keys:
    contig_list = []
    for locus in all_loci:
        try: contig_list.append(','.join(ref_by_samp[lib]['REF_'+locus]))
        except KeyError: contig_list.append('NA')
    sum_file.write('%s\t%s\n' % (lib, ('\t').join(contig_list)))
sum_file.close()

# Print final references
########################
for lib in lib_keys:
    all_contigs = fa_parser(args.assembly+lib+'.fa.final')
    out_file = open(args.out+lib+'.fa', 'w')
    for locus in all_loci:
        if 'REF_'+locus not in ref_by_samp[lib]: continue
        if len(ref_by_samp[lib]['REF_'+locus]) > 1:
               for each_contig in ref_by_samp[lib]['REF_'+locus]:
                   out_file.write('>%s\n' % locus+each_contig.split('_')[1])
                   if rev_dict[each_contig] == 'Y':
                       out_file.write('%s\n' % rev_comp(all_contigs[each_contig.split('_')[1]]))
                   elif rev_dict[each_contig] == 'N':
                       out_file.write('%s\n' % all_contigs[each_contig.split('_')[1]])
        elif len(ref_by_samp[lib]['REF_'+locus]) == 1:
            out_file.write('>%s\n' % locus)
            try:
                if rev_dict[ref_by_samp[lib]['REF_'+locus][0]] == 'Y':
                    out_file.write('%s\n' % rev_comp(all_contigs[ref_by_samp[lib]['REF_'+locus][0].split('_')[1]]))
                elif rev_dict[ref_by_samp[lib]['REF_'+locus][0]] == 'N':
                    out_file.write('%s\n' % all_contigs[ref_by_samp[lib]['REF_'+locus][0].split('_')[1]])
            except KeyError: 
                out_file.write('%s\n' % all_contigs[ref_by_samp[lib]['REF_'+locus][0].split('_')[1]])
                print 'Check '+ref_by_samp[lib]['REF_'+locus][0]
    out_file.close()

#os.system('rm %sblastIN* %sblast.out %sblast_filt.out' % (args.out, args.out, args.out))
print '==> COMPLETE!\n'                                  

