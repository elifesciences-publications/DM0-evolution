#!/usr/bin/env python

'''
grep-maeA-citT-polymorphism.py by Rohan Maddamsetti.

This script greps through all annotated.gd to report polymorphic mutations
related to maeA (sfcA in these genomes), citT, rnk
(whose promoter regulates citT in Cit+), dctA, and makes a table
for the manuscript.

Usage: python grep-maeA-citT-polymorphism.py > ../results/maeA-citT-polymorphism.csv

'''

from os.path import exists, join, basename, isdir
from os import makedirs, listdir, chdir, getcwd

def main():
    projdir = "/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/"
    genomedir = join(projdir,"genomes/polymorphism/")

    evidence_codes = {'RA','MC', 'JC', 'UN'}
    ## rnk included for its promoter.
    relevant_genes = {'rnk','citT','sfcA','dctA'}

    for g in [x for x in listdir(genomedir) if x.startswith('ZDB')]:
        gdfile = join(genomedir,g,"output/evidence/annotated.gd")
        gd_h = open(gdfile)
        for l in gd_h:
            ## skip evidence lines.
            if any([y in l for y in evidence_codes]):
                continue
            if 'frequency=1\t' in l: ## skip non-polymorphic mutations
                continue
            if any([z in l for z in relevant_genes]):
                data = l.split()
                print(','.join(['genome='+g]+data))
            #    mut_type = data[0]
            #    position = data[4]
            #    new_aa = data[6]
            #    aa_pos = data[7]
            #    old_aa = data[8]
            #    new_codon = data[9]
            #    codon_num = data[10]
            #    old_codon = data[12]
            #    freq = data[13]
            #    gene_name = data[14]
            #    if 'sfcA' in gene_name:
            #        gene_name = 'gene_name=maeA'
            #    relevant_data = [g,gene_name,freq,position,old_aa,aa_pos,new_aa,old_codon,codon_num,new_codon]
            #    print(','.join(relevant_data))

main()
