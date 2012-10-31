#! /usr/bin/env python

"""snp_mapper.py: Mapping snp to coding regions"""

__author__ = "Chen Cao"
__copyright__ = "Copyright 2011, University of Maryland College Park"

import os,sys
abspath = os.path.dirname(__file__)
sys.path.append(abspath)
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from snp_mapper_lib import ExonInfo
from operator import attrgetter

write = sys.stdout.write
exon_seq_handle = open(abspath + "/Human_hg19_refseq_seq_exon") #load the file for the mRNA sequence from the UCSC genome browser
exon_info_handle = open(abspath + "/Human_hg19_refseq_info") #load the exon info handle the coordinate file from UCSC genome browser
prefix = "hg19_refGene_"
gene_info_list = {} #a dictionary of ExonInfo class
exon_info_lines = exon_info_handle.readlines()[1:] #neglect the header line
exon_seq_records = list(SeqIO.parse(exon_seq_handle,"fasta"))

for i,(info_line,seq_record) in enumerate(zip(exon_info_lines,exon_seq_records)):
    info_list = info_line.split('\t')
    if seq_record.id == prefix + info_list[1]:
        if not info_list[2] in gene_info_list: gene_info_list[info_list[2]]=[]
        gene_info_list[info_list[2]].append(ExonInfo(info_list))
        gene_info_list[info_list[2]][-1].seq = seq_record.seq
    else: print "Gene id not match!!" + info_list[1]
exon_info_handle.close()
exon_seq_handle.close()

for chr in gene_info_list:
    gene_info_list[chr].sort(key=attrgetter('cds_start'))

#function to map all mutation
def map_all_mutation_location(chr, pos, ori, mut,fd,gene_info_list):
    chr.lower()
    for gene in gene_info_list[chr]:
        if (pos >= (gene.tx_start + 1 ) and pos < (gene.cds_start) and gene.strand=='+') or (pos >= (gene.cds_end + 1 ) and pos < (gene.tx_end) and gene.strand=='-'):
            fd.write(chr + "\t" + str(pos) + "\t" + ori + "\t" + mut + "\t" + gene.refseq_id + "\t" + "5-UTR\n" )
            write("\t5'UTR of " + gene.refseq_id)
        if (pos >= (gene.tx_start + 1 ) and pos < (gene.cds_start) and gene.strand=='-') or (pos >= (gene.cds_end + 1 ) and pos < (gene.tx_end) and gene.strand=='+'):
            fd.write(chr + "\t" + str(pos) + "\t" + ori + "\t" + mut + "\t" + gene.refseq_id + "\t" + "3-UTR\n" )
            write("\t3'UTR of " + gene.refseq_id)
        if pos >= (gene.cds_start + 1) and pos < (gene.cds_end + 1):
            exon = 0 #whether the SNP is inside exon
            for i, (start, end) in enumerate(zip(gene.exon_starts,gene.exon_ends)):
               if pos >= (start + 1) and pos <(end + 1):
                   exon = 1
                   fd.write(chr + "\t" + str(pos) + "\t" + ori + "\t" + mut + "\t")
                   offset = pos - (start + 1)
                   hit_exon_start = 1
                   for j in range(i):
                       hit_exon_start += (gene.exon_ends[j] - gene.exon_starts[j]) 
                   mut_pos_exon = hit_exon_start + offset
                   write(" The " + str(i+1) + "th exon in " + gene.refseq_id)
                   fd.write(gene.refseq_id + "\t")
                   if not hasattr(gene,'cds_start_pos_in_exons'):continue
                   mut_pos_cds = mut_pos_exon - gene.cds_start_pos_in_exons + 1
                   write(" at the " + str(mut_pos_exon) + " nt in the exon sequence")
                   if gene.strand == '+': coding_seq = gene.seq[(gene.cds_start_pos_in_exons - 1):(gene.cds_end_pos_in_exons - 1)]
                   else: coding_seq = gene.seq.reverse_complement()[(gene.cds_start_pos_in_exons - 1):(gene.cds_end_pos_in_exons - 1)]
                   if len(coding_seq) == 0:continue
                   try:
                       if gene.strand == '+': ori_protein = coding_seq.translate(cds = True)
                       else: ori_protein = coding_seq.reverse_complement().translate(cds = True)
                   except:
                       write(" Exception! type:" + str(sys.exc_info()[0]) + " value:" + str(sys.exc_info()[1]) + " traceback:" + str(sys.exc_info()[2]))
                       fd.write("Exception\n")
                       break
                   coding_seq_mut = coding_seq.tomutable()
                   write(" at the " + str(mut_pos_cds) +
                          " nt in the coding sequence")
                   if coding_seq_mut[mut_pos_cds - 1] == ori : coding_seq_mut[mut_pos_cds - 1] = mut
                   else : 
                       write(" The original nt doesn't match with the sequence")
                       fd.write("Exception\n")
                       break 
                   if gene.strand == '+':
                       mut_protein = coding_seq_mut.toseq().translate()
                       ori_protein = coding_seq.translate()
                   else: 
                       mut_protein = coding_seq_mut.toseq().reverse_complement().translate()
                       ori_protein = coding_seq.reverse_complement().translate()
                   ns = 0 #flag for nonsynonymous mutation
                   for k, (ori_aa, mut_aa) in enumerate(zip(ori_protein,mut_protein)):
                       if not ori_aa == mut_aa : 
                           write(" residue mutation" + ori_aa + str(k+1) + mut_aa)
                           ns = 1
                           fd.write(ori_aa + str(k+1) + mut_aa + "\n")
                   if ns == 0 : 
                           write(" This is a synonymous mutation")
                           fd.write("synonymous\n")
            if exon == 0:
                fd.write(chr + "\t" + str(pos) + "\t" + ori + "\t" + mut + "\t" + gene.refseq_id + "\t" + "intron\n" )
                write("\tintron of " + gene.refseq_id)       
        else:continue

#function to map coding region mutation only
def map_mutation_location(chr, pos, ori, mut,fd,gene_info_list):
    chr.lower()
    for gene in gene_info_list[chr]:
        if pos >= (gene.cds_start + 1) and pos < (gene.cds_end + 1):
            for i, (start, end) in enumerate(zip(gene.exon_starts,gene.exon_ends)):
               if pos >= (start + 1) and pos <(end + 1):
                   fd.write(chr + "\t" + str(pos) + "\t" + ori + "\t" + mut + "\t")
                   offset = pos - (start + 1)
                   hit_exon_start = 1
                   for j in range(i):
                       hit_exon_start += (gene.exon_ends[j] - gene.exon_starts[j])
                   mut_pos_exon = hit_exon_start + offset
                   write(" The " + str(i+1) + "th exon in " + gene.refseq_id)
                   fd.write(gene.refseq_id + "\t")
                   if not hasattr(gene,'cds_start_pos_in_exons'):continue
                   mut_pos_cds = mut_pos_exon - gene.cds_start_pos_in_exons + 1
                   write(" at the " + str(mut_pos_exon) + " nt in the exon sequence")
                   if gene.strand == '+': coding_seq = gene.seq[(gene.cds_start_pos_in_exons - 1):(gene.cds_end_pos_in_exons - 1)]
                   else: coding_seq = gene.seq.reverse_complement()[(gene.cds_start_pos_in_exons - 1):(gene.cds_end_pos_in_exons - 1)]
                   if len(coding_seq) == 0:continue
                   try:
                       if gene.strand == '+': ori_protein = coding_seq.translate(cds = True)
                       else: ori_protein = coding_seq.reverse_complement().translate(cds = True)
                   except:
                       write(" Exception! type:" + str(sys.exc_info()[0]) + " value:" + str(sys.exc_info()[1]) + " traceback:" + str(sys.exc_info()[2]))
                       break
                   coding_seq_mut = coding_seq.tomutable()
                   write(" at the " + str(mut_pos_cds) +
                          " nt in the coding sequence")
                   if coding_seq_mut[mut_pos_cds - 1] == ori : coding_seq_mut[mut_pos_cds - 1] = mut
                   else :
                       write(" The original nt doesn't match with the sequence")
                       break
                   if gene.strand == '+':
                       mut_protein = coding_seq_mut.toseq().translate()
                       ori_protein = coding_seq.translate()
                   else:
                       mut_protein = coding_seq_mut.toseq().reverse_complement().translate()
                       ori_protein = coding_seq.reverse_complement().translate()
                   ns = 0 #flag for nonsynonymous mutation
                   for k, (ori_aa, mut_aa) in enumerate(zip(ori_protein,mut_protein)):
                       if not ori_aa == mut_aa :
                           write(" residue mutation" + ori_aa + str(k+1) + mut_aa)
                           ns = 1
                           fd.write(ori_aa + str(k+1) + mut_aa + "\n")
                   if ns == 0 :
                           write(" This is a sysnonymous mutation")
                           fd.write("sysnonymous\n")
        else:continue


var_handle = open(sys.argv[1])
output_handle = open(sys.argv[2],'w')
for line in var_handle:
    if line[0] == '#':continue
    if line[:4] == 'chrM': continue
    items = line.split('\t')
    write (line.rstrip() + "\t")
    map_all_mutation_location(items[0], int(items[1]), items[3],items[4],output_handle,gene_info_list)
    write("\n")       
output_handle.close()
