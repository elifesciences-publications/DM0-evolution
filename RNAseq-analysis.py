#!/usr/bin/env python

''' 

RNAseq-analysis.py by Rohan Maddamsetti. 

This script does two tasks.

1) make FASTA file of all ORFs in REL606.

2) concatenates RNAseq fastq files by genome and biological replicate.

NOTE: this will run on my office computer-- not on my laptop.

'''

from os.path import join, expanduser, exists
from os import makedirs, listdir
from Bio import SeqIO
import pandas as pd
import subprocess

def generateORF_fasta(ref_file,outf):

    with open(outf,'w') as outfh:
        rec = SeqIO.read(open(ref_file),"genbank")
        for feat in rec.features:
            ## only consider CDS for now.
            if feat.type != 'CDS':
                continue
            
            locus_tag = 'NA'
            gene = 'NA'
            Note = 'NA'
            location = str(feat.location)
            ORF_seq = feat.location.extract(rec).seq
            try:
                locus_tag = feat.qualifiers['locus_tag'][0]
            except KeyError:
                print('warning: locus_tag missing')
            try:
                gene = feat.qualifiers['gene'][0]
            except KeyError:
                print('warning: gene name missing')
            try:
                note = feat.qualifiers['Note'][0]
            except KeyError:
                print('warning: note missing')
                
            outfh.write(">locus_tag=%s|gene=%s|Note=%s|Location=%s\n%s\n" % (
                locus_tag,
                gene,
                note,
                location,
                ORF_seq
            ))

def generate_kallisto_quant_cmds(RNAseqdir,samples_f,outdir):
    cmd_strings = []    
    sample_df = pd.read_csv(samples_f)
    grouped = sample_df.groupby(['Genome','BiologicalReplicate'])
    for name, group in grouped:
        search_strs = [str(x) for x in group['Admera Health ID']]
        input_files = []
        for s in search_strs:
            in_f = [join(RNAseqdir,x) for x in listdir(RNAseqdir) if s in x]
            input_files = input_files + in_f
        ## use gzcat to speed up kallisto quant.
        input_files = ['<(gzcat ' + x + ')' for x in input_files]
        input_fstring = ' '.join(input_files)
        output_file = join(outdir,'_'.join([str(x) for x in name]))
        cmd = 'kallisto quant -i ../results/RNAseq-analysis/LCA_ORFs.idx -o ' + output_file + ' -b 100 ' + input_fstring
        cmd_strings.append(cmd)
    return cmd_strings

            
def main():
    homedir = expanduser("~")
    projdir = join(homedir,"BoxSync/DM0-evolution")
    RNAseq_analysis_dir = join(projdir,"results/RNAseq-analysis")

    if not exists(RNAseq_analysis_dir):
        makedirs(RNAseq_analysis_dir)
    
    LCAref = join(projdir,'genomes/curated-diffs/LCA.gbk')
    ORF_file = join(RNAseq_analysis_dir,"LCA_ORFs.fasta")

    ## this is surprisingly slow. Only run if the file doesn't exist.
    if not exists(ORF_file):
        generateORF_fasta(LCAref,ORF_file)

    raw_RNAseq_dir = join(projdir,"Kenyon-RNAseq")
    RNAseq_samples_f = join(projdir,"data/rohan-formatted/Kenyon-RNAseq-samples.csv")

    quant_cmds = generate_kallisto_quant_cmds(raw_RNAseq_dir,RNAseq_samples_f, RNAseq_analysis_dir)

    for qcmd in quant_cmds:
        outdir = qcmd.split()[5]
        ## only run kallisto quant if the outdir doesn't exist.
        if not exists(outdir):
            print(qcmd)
            subprocess.run(qcmd,shell=True,executable='/bin/bash')

main()
