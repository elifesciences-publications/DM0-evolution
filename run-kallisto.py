#!/usr/bin/env python

''' 
run-kallisto.py by Rohan Maddamsetti. 

This script does two tasks.
1) make FASTA file of all ORFs in REL606.
2) concatenates RNAseq fastq files by genome and biological replicate.
NOTE: this will run on my office computer-- not on my laptop.

'''

from os.path import join, exists, dirname, realpath
from os import makedirs, listdir
from Bio import SeqIO
import pandas as pd
import subprocess
## need shlex so that we don't split on strings in double quotes.
import shlex 

def quote_string(x):
    ''' add double quotes to a string, in case it contains spaces or
        other weirdness. This is required for system calls using the
        shell, but not when working within python itself. '''
    return('\"' + x + '\"')

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
                
            outfh.write(">locus_tag=%s|gene=%s\n%s\n" % (
                locus_tag,
                gene,
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
        input_files = ['<(gzcat ' + quote_string(x) + ')' for x in input_files]
        input_fstring = ' '.join(input_files)
        output_file = join(outdir,'_'.join([str(x) for x in name]))
        idx_file = join('..','results','RNAseq-analysis','LCA_ORFs.idx')
        cmd = 'kallisto quant -i ' + quote_string(idx_file) + ' -o ' + quote_string(output_file) + ' -b 100 ' + input_fstring
        cmd_strings.append(cmd)
    return cmd_strings

            
def main():

    srcdir = dirname(realpath(__file__))
    assert srcdir.endswith('src')
    projdir = dirname(srcdir)
    assert projdir.endswith('DM0-evolution')
    
    RNAseq_analysis_dir = join(projdir,"results","RNAseq-analysis")

    if not exists(RNAseq_analysis_dir):
        makedirs(RNAseq_analysis_dir)
    
    LCAref = join(projdir,'genomes","curated-diffs","LCA.gbk')
    ORF_file = join(RNAseq_analysis_dir,"LCA_ORFs.fasta")

    ## this is surprisingly slow. Only run if the file doesn't exist.
    if not exists(ORF_file):
        generateORF_fasta(LCAref,ORF_file)

    raw_RNAseq_dir = join(projdir,"Kenyon-RNAseq")
    RNAseq_samples_f = join(projdir,"data","rohan-formatted","Kenyon-RNAseq-samples.csv")

    ## run kallisto index.
    fasta_input = join(projdir, 'results','RNAseq-analysis','LCA_ORFs.fasta')
    idx_output = join(projdir, 'results','RNAseq-analysis','LCA_ORFs.idx')

    index_cmd = ' '.join(['kallisto index --make-unique -i', quote_string(idx_output), quote_string(fasta_input)])
    ## only run if the index doesn't exist.
    if not exists(idx_output):
        subprocess.run(index_cmd,shell=True,executable='bash')
    
    quant_cmds = generate_kallisto_quant_cmds(raw_RNAseq_dir,RNAseq_samples_f, RNAseq_analysis_dir)

    for qcmd in quant_cmds:
        outdir = shlex.split(qcmd)[5]
        ## only run kallisto quant if the outdir doesn't exist.
        if not exists(outdir):
            print(qcmd)
            subprocess.run(qcmd,shell=True,executable='bash')

main()
