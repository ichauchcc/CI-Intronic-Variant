#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 12:24:42 2020

@author: yuchen
"""

##from subprocess import call
# if want to read from vcf could use 'pysam'
import pysam
from maxentpy import maxent
from itertools import islice, tee

# get the +/- 500 (len = 1000) sequence from the variant point

def reverse_seq(x):
  return x[::-1]

def sliding_window(iterable, size): 
    iterators = tee(iterable, size) 
    iterators = [islice(iterator, i, None) for i, iterator in enumerate(iterators)]  
    yield from zip(*iterators)

#list(sliding_window(range(5), 3))

f=open('test.fasta')
seq={}
for line in f:
        if line.startswith('>'):
                name=line.replace('>','').split()[0]
                seq[name]=''
        else:
                seq[name]+=line.replace('\n','').strip()
f.close()

strand = '-'
# 11	47367305	G	A - (new donor)
ref1_chr = 11
ref1_pos = 47367305
ref1_ref = 'G'
ref1_var = 'A'

ref1_flag = 'donor'

# 11	47364865	C	T - (new acceptor)
ref2_chr = 11
ref2_pos = 47364865
ref2_ref = 'C'
ref2_var = 'T'
ref2_flag = 'acceptor'

for sequence in seq.values():
    if strand == '-':
        rev_seq = reverse_seq(str(sequence))
        #print(rev_seq[int((len(rev_seq)-1)/2)])
        if ref1_flag == 'donor':
            tar_seq = rev_seq[int((len(rev_seq)-1)/2):]
            print(tar_seq,len(tar_seq))

#seq = str(pysam.faidx('/home/groups/cardio/References/Bwa_b37/hs37d5.fa','{0}:{1}-{2}'.format(ref1_chr, ref1_pos, ref1_pos)))

##for command in ("samtools faidx /home/groups/cardio/References/Bwa_b37/hs37d5.fa {chr}:{start}-{end}"\.format(chr=ref1_chr,start=)):
##    # shell = True is so you can handle redirects like output a file
##    # e.g. samtools indaxstarts filename.bam > filename.txt
##    call(command, shell = True)

#maxent.score5('cagGTAAGT')  # 3 bases in exon and 6 bases in intron

#maxent.score3('ttccaaacgaacttttgtAGgga')  # 20 bases in the intron and 3 base in the exon
