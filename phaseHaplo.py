#!/usr/bin/home python

# This is a script to use the output of makeConsensus.py to make haplotypes
# This uses both read-backed phasing from GATK and statistical phasing using PHASE
# This assumes a vcf format 4.1
# 31 Oct 2014
# Joshua Penalba
# Dependencies: MAFFT, PHASE


# SETTING UP ARGUMENTS
######################

import argparse
import os
import sys

ph = argparse.ArgumentParser(description='Creates haplotype sequences from the consensus of the SCPP data')
ph.add_argument('-i', dest='in_path', type=str, help='path to consensus and vcf files (files end in .fa and .vcf)')
ph.add_argument('-m', dest='meta', type=str, help='path to meta data')
ph.add_argument('-o', dest='out', type=str, help='path to output')
ph.add_argument('-DP', dest='cov', type=int, help='minimal coverage filter [10]', default=10)
ph.add_argument('-GQ', dest='qual', type=int, help='minimal quality filter[21]', default=21)
ph.add_argument('-M', dest='maf', type=float, help='minimal MAF filter [0.2]', default=0.2)
ph.add_argument('-H', dest='freq', type=float, help='minimal frequency to call homozygous sites [0.8]', default=0.8)
if len(sys.argv) == 1:
    ph.print_help()
    sys.exit(1)
args = ph.parse_args()
if args.out.endswith('/'): pass
else: args.out = args.out+'/'
try: os.mkdir(args.out)
except OSError: pass

# SETTING UP FUNCTIONS
######################

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


# ORGANIZE DATA INTO SYSTEMS
############################

meta_data = {}
meta_file = open(args.meta, 'r')
meta_file.readline()
for lines in meta_file:
    info = lines.strip().split()
    if info[3] not in meta_data:
        meta_data[info[3]] = {info[0]:info[1]}
    elif info[3] in meta_data:
        meta_data[info[3]][info[0]] = info[1]
meta_file.close()

for system in meta_data:
    out_path = args.out+system+'/'
    try: os.mkdir(out_path)
    except OSError: pass
    for lib in meta_data[system]:
        os.system('cp %s%s.vcf %s%s_consensus.fa %s' % (args.in_path, lib, args.in_path, lib, out_path))
    phased = {}
    unphased = {}
    consensus_seq = {}
    system_dict = {}

    # Read-backed phasing using the vcf
    all_files = os.listdir(out_path)
    for each_file in all_files:
        if each_file.endswith('.vcf'):
            vcf_dict = {}
            indel_dict = {}
            vcf_file = open(out_path+each_file, 'r')
            for lines in vcf_file:
                info = lines.strip().split()
                (locus, pos, ref_base, alt_base) = (info[0], info[1], info[3], info[4]) 
                (data, labels, phase_data) = (info[7], info[8].split(':'), info[9].split(':'))
                cov = int(phase_data[labels.index('DP')])
                qual = int(phase_data[labels.index('GQ')])
                if locus not in vcf_dict: # --->  Reinitiate variables with new locus
                    vcf_dict[locus] = list()
                    curr_locus = locus
                    haplo_group = 0
                    curr_group = ''
                try: phasing = phase_data[labels.index('HP')].split(',') # ---> Extract phasing information for heterozygotes
                except ValueError: phasing = 'F'
                if cov >= args.cov and qual >= args.qual: #---> Only keep those that exceed the coverage and quality threshold
                    # ---> If this is a heterozygote
                    if data.startswith('ABHet'):
                        if phasing[0].split('-')[0] != curr_group: # ---> Check for a different haplotype block
                            curr_group = phasing[0].split('-')[0]
                            haplo_group += 1
                        freq = float(data.split(';')[0].split('=')[1]) 
                        if freq >= args.freq: # ---> Keep homo reference base if the freq exceeds setting
                            vcf_dict[locus].append([pos, ref_base, ref_base])
                        elif freq < args.maf:
                            vcf_dict[locus].append([pos, alt_base, alt_base]) #---> Keep homo alt base if the freq is below setting
                        elif freq >= args.maf and freq < args.freq: # ---> Keep hetero if falls between range
                            if phasing[0].split('-')[1] == '1': # ---> Check for phasing directionality
                                vcf_dict[locus].append([pos, ref_base, alt_base, haplo_group])
                            else: 
                                vcf_dict[locus].append([pos, alt_base, ref_base, haplo_group])
                    elif data.startswith('ABHom'):
                        if phase_data[0].startswith('1/1'):
                            vcf_dict[locus].append([pos, ref_base, ref_base])
                        elif phase_data[0].startswith('0/0'):
                            vcf_dict[locus].append([pos, alt_base, alt_base])
                    elif data.startswith('AC='):
                        if alt_base.find(',') != -1:
                            vcf_dict[locus].append([pos, alt_base.split(',')[0], alt_base.split(',')[1], "?"])
                elif cov < args.cov: # ---> If below coverage threshold, turn into N
                    vcf_dict[locus].append([pos, 'N', 'N'])

            # WRITE SUMMARY FILE
            summ_file = open(out_path+each_file+'.summary', 'w')
            summ_file.write('Summary Format:\n#########\nLOCUS_NAME\nBase1\tBase2\tBase3\n')
            summ_file.write('Allele1A\tAllele2A\tAllele3A\nAllele1B\tAllele2B\tAllele3B\n')
            summ_file.write('HaploBlock1\tHaploBlock2\tHaploBlock3\n##########\n\n')
            for locus in vcf_dict:
                summ_file.write(locus+'\n')
                positions = []
                alleleA = []
                alleleB = []
                haplo_block = []
                for variant in range(len(vcf_dict[locus])):
                    var_data = vcf_dict[locus][variant]
                    positions.append(str(var_data[0]))
                    alleleA.append(str(var_data[1]))
                    alleleB.append(str(var_data[2]))
                    if len(var_data) == 4:
                        if var_data[3] > 1:
                            library = each_file.split('.')[0]
                            if locus not in unphased: 
                                unphased[locus] = {library:[var_data[3], 0]} # ---> 0 refers to single SNPs and the var_data[3] is the number of blocks
                            elif locus in unphased: 
                                if library in unphased[locus]:
                                    unphased[locus][library] = [var_data[3],1] # ---> 1 refers to haplotype blocks
                                elif library not in unphased[locus]:
                                    unphased[locus][library] = [var_data[3],0]
                        haplo_block.append(str(var_data[3]))
                    elif len(var_data) == 3: haplo_block.append('-')
                summ_file.write('\t'.join(positions)+'\n')
                summ_file.write('\t'.join(alleleA)+'\n')
                summ_file.write('\t'.join(alleleB)+'\n')
                summ_file.write('\t'.join(haplo_block)+'\n\n')
            summ_file.close()
            system_dict[each_file.split('.')[0]] = vcf_dict
            
        elif each_file.endswith('_consensus.fa'):
            consensus_seq[each_file.split('_')[0]] = fa_parser(out_path+each_file)


    # Statistical phasing using PHASE
    # ---> Call MAFFT
    unphased_singles = {}
    for locus in unphased:
        locus_out = open(out_path+locus+'.fa', 'w')
        for library in consensus_seq:
            try: locus_out.write('>%s\n%s\n' % (library, consensus_seq[library][locus]))
            except KeyError: pass
        locus_out.close()
        
        locus_dict = fa_parser(out_path+locus+'.fa')
        locus_dict_keys = locus_dict.keys()
        if len(locus_dict_keys) == 1:
            single_lib = locus_dict_keys[0]
            if single_lib not in unphased_singles:
                unphased_singles[single_lib] = {locus:system_dict[single_lib][locus]}
            elif single_lib in unphased_singles:
                unphased_singles[single_lib][locus] = system_dict[single_lib][locus]
            continue

        os.system('mafft %s%s.fa > %s%s.mafft.fa' % (out_path, locus, out_path, locus))
        mafft_dict = fa_parser(out_path+locus+'.mafft.fa')
        front_shift = {} # ---> LIB -> length
        indel_shift = {} # ---> LIB -> location, length
        samples = mafft_dict.keys()
        samples.sort()
        for each in mafft_dict:
            front_length = 0
            for nuc in mafft_dict[each]:
                if nuc == '-':
                    front_length += 1
                else: break
            front_shift[each] = front_length
            main_seq = mafft_dict[each].strip('-')
            # ---> This looks for the position and the length of indels to figure out other shifts.
            # The position is relative to the original sequence not the MAFFT alignment
            if main_seq.find('-') != -1:
                indel_shift[each] = {}
                split_seq = main_seq.split('-')
                record = 'ON'
                curr_position,prev_len = -1,0
                for fragments in split_seq:
                    print split_seq
                    curr_len = len(fragments)
                    if curr_len == 0:
                        if prev_len > 0: curr_position += 2
                        elif prev_len == 0: curr_position += 1
                        if record == 'ON':
                            record = 'OFF'
                            indel_position = curr_position
                            indel_shift[each][indel_position] = 2
                        elif record == 'OFF':
                            indel_shift[each][indel_position] += 1
                    elif curr_len > 0:
                        if prev_len > 0: 
                            indel_shift[each][curr_position] = 1
                        record = 'ON'
                        curr_position += len(fragments)
                    prev_len = len(fragments)
        print indel_shift
        # ---> This looks for fixed differences
        locus_length = len(mafft_dict[each])
        fixed_diff = {}
        for nuc in range(locus_length):
            cat_bases = ''
            for each in samples:
                cat_bases += mafft_dict[each][nuc]
            counter = 0
            for base in cat_bases:
                if base == '-' or base == 'n':
                    continue
                elif counter == 0:
                    focal_base = base
                    counter += 1
                elif counter == 1:
                    if base == focal_base: 
                        continue
                    elif base != focal_base:
                        fixed_diff[nuc+1] = cat_bases.upper()
        # ---> Add the vcf called SNPs
        for each in samples:
            try: 
                for variant in enumerate(system_dict[each][locus]):
                    to_add = 0
                    try: to_add += front_shift[each]
                    except KeyError: pass
                    try: 
                        indel_list = indel_shift[each].keys()
                        indel_list.sort()
                        curr_shift = int(variant[1][0]) + to_add
                        for indel in indel_list:
                            print variant, to_add
                            if curr_shift >= indel+front_shift[each]:
                                print locus, each, curr_shift, indel, variant
                                to_add += indel_shift[each][indel]
                                curr_shift += indel_shift[each][indel]
                            else: pass
                    except KeyError: pass
                    #print 'FINAL\t', each,variant, to_add
                    system_dict[each][locus][variant[0]][0] = int(variant[1][0]) + to_add
            except KeyError: pass

        # ---> Make directory to turn into PHASE input
        phase_info_dict = {} # ---> pos > lib > bases, haplo_group
        for position in fixed_diff: # ---> Adds fixed positions
            phase_info_dict[position] = {}
            for base in enumerate(fixed_diff[position]):
                phase_info_dict[position][samples[base[0]]] = [base[1],base[1],1]

        for each in samples: # ---> Adds vcf calls
            try:
                for variant in system_dict[each][locus]:
                    if len(variant) == 4:
                        position = int(variant[0])
                        if position not in phase_info_dict:
                            phase_info_dict[position] = {each:variant[1:]}
                        elif position in phase_info_dict:
                            phase_info_dict[position][each] = variant[1:]
            except KeyError: pass

        # ---> Adds missing calls 
        for position in phase_info_dict:
            included_libs = phase_info_dict[position].keys()
            if len(included_libs) < len(samples):
                to_do = samples[:]
                for each_lib in included_libs:
                    del to_do[to_do.index(each_lib)]
                for each_lib in to_do:
                    phase_info_dict[position][each_lib] = [mafft_dict[each_lib][position+1].upper(), mafft_dict[each_lib][position+1].upper(), 1]
        
        # ---> Makes the PHASE input files 
        all_pos = phase_info_dict.keys()
        all_pos.sort()
        phase_input = open(out_path+locus+'_phase.inp', 'w')
        phase_known = open(out_path+locus+'_phase.known', 'w')
        phase_input.write(str(len(samples))+'\n')
        phase_input.write(str(len(all_pos))+'\n')
        position_input = [str(i) for i in all_pos]
        phase_input.write('P '+' '.join(position_input)+'\n')
        phase_input.write('S'*len(all_pos)+'\n')
        input_dict = {}
        for each_pos in all_pos:
            for lib in phase_info_dict[each_pos]:
                if lib not in input_dict:
                    input_dict[lib] = {each_pos:phase_info_dict[each_pos][lib]}
                elif lib in input_dict:
                    input_dict[lib][each_pos] = phase_info_dict[each_pos][lib]
        for each_sam in samples:
            phase_input.write(each_sam+'\n')
            top_line, bottom_line, phase_line = [], [], [] 
            for each_pos in all_pos:
                top_line.append(input_dict[each_sam][each_pos][0])
                bottom_line.append(input_dict[each_sam][each_pos][1])
                try:
                    if input_dict[each_sam][each_pos][2] == 1:
                        phase_line.append('0')
                    elif input_dict[each_sam][each_pos][2] != 1:
                        phase_line.append('*')
                except IndexError: phase_line.append('0')
            phase_input.write(' '.join(top_line).replace('-','?').replace('N','?')+'\n')
            phase_input.write(' '.join(bottom_line).replace('-','?').replace('N','?')+'\n')
            phase_known.write(''.join(phase_line)+'\n')
        phase_input.close()
        phase_known.close()
        
        # ---> Run PHASE
        os.system('PHASE -k%s %s %s' % (out_path+locus+'_phase.known', out_path+locus+'_phase.inp', out_path+locus+'_phase.out'))
        os.system('rm %s%s.fa' % (out_path, locus))

        # ---> Parse PHASE output (if there is one)
        phase_out = open(out_path+locus+'_phase.out', 'r')
        first_line = phase_out.readline().strip().split()
        if len(first_line) != 0:
            curr_line = ''
            while curr_line != 'BEGIN BESTPAIRS2':
                curr_line = phase_out.readline().strip()
                split_line = curr_line.split(':')
                if len(split_line) > 0 and split_line[0] == 'Positions of loci': phased_positions = split_line[1].split()
            phased_snps = {}
            for lines in phase_out:
                info = lines.strip().split()
                if info[0] == 'END': break
                elif info[0].startswith('0'): 
                    phased_lib = info[1]
                    line = 'A'
                else:
                    if line == 'A': 
                        phase_lineA, line = info, 'B'
                    elif line == 'B':
                        phase_lineB = info
                        phased_snps[phased_lib] = [phase_lineA, phase_lineB]

            # ---> Adds the phased SNPs into the main dictionary
            for phased_lib in phased_snps:
                for pos in enumerate(phased_positions):
                    if phased_snps[phased_lib][0][pos[0]] == '=' or phased_snps[phased_lib][0][pos[0]] == '?': 
                        pass
                    else: 
                        for unphased_snps in enumerate(system_dict[phased_lib][locus]):
                            if unphased_snps[1][0] == int(pos[1]):
                                system_dict[phased_lib][locus][unphased_snps[0]] = [int(pos[1]), phased_snps[phased_lib][0][pos[0]], phased_snps[phased_lib][1][pos[0]], 1]

        # ---> Add the mafft alignments to the main library
        # This has to be done because the position has been changed in relation to the mafft position
        for mafft_lib in mafft_dict:
            consensus_seq[mafft_lib][locus] = mafft_dict[mafft_lib]

    # ---> Print the final haplotypes
    haplo_dict = {}
    for library in consensus_seq:
        haplo_dict[library] = {}
        for locus in consensus_seq[library]:
            hap_one, hap_two = consensus_seq[library][locus], consensus_seq[library][locus]
            haplo_dict[library][locus] = [hap_one, hap_two]
    for library in haplo_dict:
        out_hap = open(out_path+library+'_phased.fa', 'w')
        for locus in haplo_dict[library]:
            hap_seq = haplo_dict[library][locus]
            hap_seq[0], hap_seq[1] = list(hap_seq[0]), list(hap_seq[1])
            try: 
                for variant in system_dict[library][locus]:
                    position = int(variant[0]) - 1
                    if len(variant[1]) == 3:
                        variant[1], variant[2] = variant[1][1], variant[2][1]
                    hap_seq[0][position], hap_seq[1][position] = variant[1], variant[2]
            except KeyError: pass
            hap_seq[0], hap_seq[1] = ''.join(hap_seq[0]).replace('-',''), ''.join(hap_seq[1]).replace('-','')
            out_hap.write('>%s_hap0\n%s\n>%s_hap1\n%s\n' % (locus, haplo_dict[library][locus][0].upper(), locus, haplo_dict[library][locus][1].upper()))
        out_hap.close()
    try: os.mkdir(out_path+'tmp_files')
    except OSError: pass
    os.system('mv %s* %stmp_files' % (out_path, out_path))
    os.system('mv %stmp_files/*_phased.fa %s' % (out_path, out_path))

