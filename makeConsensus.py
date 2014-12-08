#!/usr/bin/home python

#This is a script to make consensus out of alignments
# 26 Feb 2014
# Joshua Penalba


#SETTING UP ARGUMENTS
import argparse
import os
import sys

con = argparse.ArgumentParser(description='Creates consensus sequences of alignments to SCPP data')
con.add_argument('-G', dest='gatk', type=str, help='path to GATK directory')
con.add_argument('-p', dest='picc', type=str, help='path to picard directory')
con.add_argument('-b', dest='bam', type=str, help='path to bam files')
con.add_argument('-f', dest='ref', type=str, help='path to references')
con.add_argument('-c', dest='summ', type=str, help='path to coverage summary')
con.add_argument('-o', dest='out', type=str, help='path to output directory')
con.add_argument('-D', dest='cov', type=int, help='minimal coverage filter [10]', default=10)
con.add_argument('-Q', dest='qual', type=int, help='minimal quality filter [21]', default=21)
con.add_argument('-M', dest='maf', type=float, help='minimal MAF filter [0.2]', default=0.2)
con.add_argument('-H', dest='freq', type=float, help='minimal frequency to call homozygous sites [0.8]', default=0.8)
if len(sys.argv) == 1:
    con.print_help()
    sys.exit(1)
args = con.parse_args()
if args.out.endswith('/'): pass
else: args.out = args.out+'/'
try: os.mkdir(args.out)
except OSError: pass
GATK = args.gatk+'GenomeAnalysisTK.jar'
AddOrReplace = args.picc+'AddOrReplaceReadGroups.jar'
CreateDict = args.picc+'CreateSequenceDictionary.jar'

#INITIALIZE A COMMAND
def SNP(x,y):
    z = x+y
    if z == 'AC' or z == 'CA': a = 'M'
    elif z == 'AG' or z == 'GA': a = 'R'
    elif z == 'AT' or z == 'TA': a = 'W'
    elif z == 'CG' or z == 'GC': a = 'S'
    elif z == 'CT' or z == 'TC': a = 'Y'
    elif z == 'GT' or z == 'TG': a = 'K'
    return a

#######################
# PARSE COVERAGES

cov_dict = {}
cov_summ = open(args.summ, 'r')
loci = cov_summ.readline().strip().split()
print loci
for lines in cov_summ:
    info = lines.strip().split()
    covs = info[1:]
    print covs
    cov_dict[info[0]] = {}
    for cov in enumerate(covs): 
        if cov[1] != 'NA': cov_dict[info[0]][loci[cov[0]]] = float(cov[1])
cov_summ.close()

######################
#CALLING PICARD AND GATK
allfiles = os.listdir(args.bam)
for each in allfiles:
    if each.endswith('.sorted.bam'):
        lib = each.split('.')[0]
        org_ref = args.ref+lib+'.fa'
        ref_path = args.ref+lib+'_new.fa'
        try: org_file = open(org_ref, 'r')
        except IOError: continue
        new_ref = open(ref_path, 'w')
        for lines in org_file:
            info = lines.strip().split()
            if len(info) > 0:
                if info[0].startswith('>'): new_ref.write(lines)
                else:
                    line = info[0].replace('M','A')
                    line = line.replace('R','A')
                    line = line.replace('W','A')
                    line = line.replace('S','G')
                    line = line.replace('Y','T')
                    line = line.replace('K','T')
                    new_ref.write(line+'\n')
            else: pass
        org_file.close()
        new_ref.close()
        bam_path = args.bam+each
        out_path = args.out+lib
        os.system('samtools faidx %s' % ref_path)
        index_ref = args.ref+lib+'_new.dict'
        os.system('java -Dsamjdk.try_use_intel_deflater=false -jar %s INPUT=%s OUTPUT=%s.rg.bam RGID=%s RGLB=exon_capture RGPL=illumina RGPU=lane3 RGSM=%s' %(AddOrReplace, bam_path, out_path, lib, lib))
        os.system('java -Dsamjdk.try_use_intel_deflater=false -jar %s R=%s O=%s' % (CreateDict, ref_path, index_ref))
        os.system('samtools sort %s.rg.bam %s.rg.sort' % (out_path, out_path))
        os.system('samtools index %s.rg.sort.bam' % out_path)
        os.system('java -Xmx8g -jar %s -T HaplotypeCaller -R %s -I %s.rg.sort.bam -o %s.rg.vcf' % (GATK, ref_path, out_path, out_path))
        os.system('java -Xmx8g -jar %s -T ReadBackedPhasing -R %s -I %s.rg.sort.bam --variant %s.rg.vcf --min_base_quality_score %s -o %s.raw.vcf' % (GATK, ref_path, out_path, out_path, args.qual, out_path))
        os.system('java -Xmx8g -jar %s -T VariantAnnotator -A DepthPerAlleleBySample -A AlleleBalance -A FisherStrand -A HaplotypeScore -A HardyWeinberg -R %s -I %s.rg.sort.bam --variant %s.raw.vcf -o %s.annotated.vcf' % (GATK, ref_path, out_path, out_path, out_path))
        os.system("grep -v \'#\'  %s.annotated.vcf > %s.raw.vcf_org" % (out_path, out_path))
        os.system("grep -v \'ABHet=1.00\' %s.raw.vcf_org > %s.raw2.vcf" % (out_path, out_path))
        os.system("rm  %s.raw.vcf* %s.annotated.vcf* %s.rg.bam* %s.rg.sort* %s.rg.vcf* %s%s.dict %s.fai" % (out_path, out_path, out_path, out_path, out_path, args.ref, lib, ref_path))
        

#MAKE DIR OF REFERENCE LOCI
        ref_file = open(ref_path, 'r')
        ref_loci = {}
        for lines in ref_file:
            info = lines.strip().split()
            if len(info) == 0: continue
            else:
                if info[0].startswith('>'): 
                    locus = info[0][1:]
                    ref_loci[locus] = []
                else: ref_loci[locus].append(info[0])
        for each in ref_loci:
            if len(ref_loci[each]) == 1: seq = ref_loci[each][0]
            elif len(ref_loci[each]) > 1: seq = ''.join(ref_loci[each])
            ref_loci[each] = list(seq)
        ref_file.close()

#CALL SNPS
        vcf_file = open(out_path+'.raw2.vcf', 'r')
        indels = {}
        for lines in vcf_file:
            info = lines.strip().split()
            locus = info[0]
            pos = int(info[1])
            ref_base = info[3]
            alt_base = info[4]
            data = info[7]
            if data.startswith('ABH'): 
                depth = int(data.split(';DP=')[1].split(';')[0])
                state = data[0:5]
                freq = float(data.split(';')[0].split('=')[1])
                call = info[9][0:3]
                if depth >= args.cov: 
                    if data.startswith('ABHet'): 
                        if freq >= args.freq: 
                            new_base = ref_base
                        elif freq <= args.maf: 
                            new_base = alt_base
                        else: new_base = SNP(ref_base,alt_base)
                    elif data.startswith('ABHom'): 
                        if freq >= args.freq: 
                            if call == '1/1': new_base = alt_base
                            elif call == '0/0': new_base = ref_base
                        elif freq < args.freq: new_base = 'N'
                elif depth < args.cov: new_base = 'N'
                ref_loci[locus][pos-1] = new_base
            else: 
                if locus not in indels: indels[locus] = {pos:[ref_base, alt_base]}
                else: indels[locus][pos] = [ref_base, alt_base]
        '''
#CALL INDELS
        for locus in indels:
            indel_pos = indels[locus].keys()
            indel_pos.sort()
            shift = 0
            for each in indel_pos:
                refbase = indels[locus][each][0]
                altbase = indels[locus][each][1]
                if refbase < altbase:
                    ins_size = len(altbase)-1
                    ins = 'N' * ins_size
                    ref_loci[locus].insert(each+shift,ins)
                    shift += 1
                elif refbase > altbase:
                    del_size = len(refbase)-1
                    for i in range(del_size):
                        ref_loci[locus][each+i+shift] = 'N'
        '''
#####################
#MAKE OUTPUT
        out_loci = ref_loci.keys()
        out_loci.sort()
        out_file = open(out_path+'_consensus.fa', 'w')
        for locus in out_loci:
            try:
                if cov_dict[lib][locus] >= args.cov:
                    out_file.write('>%s\n' % locus)
                    out_file.write(''.join(ref_loci[locus])+'\n')
            except KeyError: pass
        out_file.close()
        os.system('mv %s.raw2.vcf %s.vcf' % (out_path, out_path))
    else: pass
