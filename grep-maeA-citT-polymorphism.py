#!/usr/bin/env python

'''
grep-maeA-citT-polymorphism.py by Rohan Maddamsetti.

This script greps through all annotated.gd to report polymorphic mutations
related to maeA (sfcA in these genomes), citT, rnk
(whose promoter regulates citT in Cit+), and dctA.

Usage: python grep-maeA-citT-polymorphism.py > ../results/maeA-citT-polymorphism.csv

'''

from os.path import exists, join, basename, isdir, dirname, realpath
from os import makedirs, listdir, chdir, getcwd

def main():
    
    srcdir = dirname(realpath(__file__))
    assert srcdir.endswith('src')
    projdir = dirname(srcdir)
    assert projdir.endswith('DM0-evolution')

    genomedir = join(projdir,"genomes","polymorphism")

    evidence_codes = {'RA','MC', 'JC', 'UN'}
    ## rnk included for its promoter.
    relevant_genes = {'rnk','citT','sfcA','dctA'}

    for g in [x for x in listdir(genomedir) if x.startswith('ZDB') or x.startswith('CZB')]:
        gdfile = join(genomedir,g,"output","evidence","annotated.gd")
        gd_h = open(gdfile)
        for l in gd_h:
            ## skip evidence lines.
            if any([y in l for y in evidence_codes]):
                continue

            if any([z in l for z in relevant_genes]):
                ## change sfcA to synonym maeA.
                l = l.replace('sfcA','maeA')
                data = l.split('\t')        
                print(','.join(['genome='+g]+data))

main()
