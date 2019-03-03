#!/usr/bin/env python

''' IS150-analysis.py by Rohan Maddamsetti.

    This script generate alignments of the regions
    flanking IS150 insertions in the following experiments:

    1) Zack's DM0 and DM25 evolution experiments.
    2) LTEE
    3) Mutation Accumulation evolution experiment.

    I did not find any obvious evidence of IS150 insertion motifs.
'''

import pandas as pd
from os.path import expanduser, join
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from random import randint

def row_to_record(tup, refseq,upstream=10,downstream=10):
    ''' Input: tuple (header string, position)
    and a genome sequence.
    Output: a Seq.Record.'''
    header_string, pos = tup
    region = refseq[pos-upstream:pos+downstream]
    return SeqRecord(id=header_string,seq=region)

def df_to_records(df, refseq):
    ''' Input: Pandas Dataframe with a 'Position' column,
    and a reference sequence (Seq object).
    Output: a list of SeqRecords.'''
    positions = list(df['Position'])
    column_names = list(df)
    ''' generate a list of header strings'''
    headers = df.apply(
        lambda y:'|'.join([c+'='+str(x) for x,c in zip(y,column_names)]),
        axis=1
    )

    records = [row_to_record(t,refseq) for t in zip(headers,positions)]
    return records

def calcGC_of_all_records(rec_list):
    concatenated = Seq("", generic_dna)
    for rec in rec_list:
        concatenated += rec.seq
    return GC(concatenated)

def GC_content_statistical_test():
    ''' this code was cut and pasted from main-- needs additional
    edits to run.
    '''

    ''' calculate GC content for the insertion sites, and for the REL606 genome'''
    insertion_site_GC = calcGC_of_all_records(all_recs)
    print(insertion_site_GC)
    print(calcGC_of_all_records(DM0_IS150_recs))
    print(calcGC_of_all_records(LTEE_MAE_recs))
    print(GC(REL606seq))

    ''' difference of about 4%-- see if this is significant.
    10000 times: take 568 substrings of length 20 from the REL606 genome, and concatenate.
    then measure the GC content. How many times larger than 54%?

    I get a p-value = 0.045 -- not worth reporting.
    '''
    reps = 10000
    samples = len(all_recs)
    genomic_pos_range = len(REL606seq) - 20 ## don't take illegal substrings past the end.
    GC_samples = []
    for i in range(reps):
        concatenated = Seq("", generic_dna)
        for j in range(samples):
            randpos = randint(0,genomic_pos_range)
            my_substring = REL606seq[randpos:randpos+1]
            concatenated += my_substring
        GC_samples.append(GC(concatenated))
    passed = [x for x in GC_samples if x >= insertion_site_GC]
    print(len(passed)/len(GC_samples))


def main():

    homedir = expanduser("~")
    projdir = join(homedir,"BoxSync/DM0-evolution")
    resultsdir = join(projdir,"results/genome-analysis")

    REL606ref = join(projdir,"genomes/REL606.7.gbk")
    REL606seq = SeqIO.read(REL606ref,"genbank").seq

    LCAref = join(projdir,"genomes/curated-diffs/LCA.gbk")
    LCAseq = SeqIO.read(LCAref,"genbank").seq

    DM0_f = join(resultsdir,"IS_insertions.csv")
    DM0_df = pd.read_csv(DM0_f, skipinitialspace=True)
    DM0_df['Position'] = DM0_df['genome_start']

    LTEE_MAE_f = join(resultsdir,"LTEE_MAE_IS150_insertions.csv")
    LTEE_MAE_df = pd.read_csv(LTEE_MAE_f, skipinitialspace=True)
    ## split into MAE_df and LTEE_df to handle double-counting in LTEE_df.
    MAE_df = LTEE_MAE_df[LTEE_MAE_df.Environment == 'MAE']

    LTEE_df = LTEE_MAE_df[LTEE_MAE_df.Environment == 'LTEE']
    LTEE_df = LTEE_df.groupby(['Environment','IS_element','Position','Population']).agg('count').reset_index()
    ## This contains independent IS150 insertions across populations
    LTEE_df = LTEE_df.drop(columns=['Clone'])

    LTEE_IS150_recs = df_to_records(LTEE_df, REL606seq)
    MAE_IS150_recs = df_to_records(MAE_df,REL606seq)
    DM0_IS150_recs = df_to_records(DM0_df,LCAseq)

    LTEE_MAE_recs = LTEE_IS150_recs + MAE_IS150_recs
    all_recs = LTEE_MAE_recs + DM0_IS150_recs
    training_unaligned = MultipleSeqAlignment(LTEE_MAE_recs)
    validate_unaligned = MultipleSeqAlignment(DM0_IS150_recs)
    full_unaligned = MultipleSeqAlignment(all_recs)

    ''' write unaligned sequences to disk for sequence alignment'''

    training_unaln = join(projdir,"results/LTEE_MAE_IS150_unaligned.fasta")
    with open(training_unaln, "w") as handle:
        count = SeqIO.write(training_unaligned, handle, "fasta")

    validate_unaln = join(projdir,"results/DM0_IS150_unaligned.fasta")
    with open(validate_unaln, "w") as handle:
        count = SeqIO.write(validate_unaligned, handle, "fasta")

    full_unaln = join(projdir,"results/all_IS150_unaligned.fasta")
    with open(full_unaln, "w") as handle:
        count = SeqIO.write(full_unaligned, handle, "fasta")

    ''' run MAFFT on these unaligned FASTA files. '''
    cline_training = MafftCommandline("mafft", input=training_unaln)
    cline_validate = MafftCommandline("mafft", input=validate_unaln)
    cline_full = MafftCommandline("mafft", input=full_unaln)

    stdout, stderr = cline_training()
    training_aln = join(projdir,"results/LTEE_MAE_IS150_aligned.fasta")
    with open(training_aln, "w") as handle:
        handle.write(stdout)

    stdout, stderr = cline_validate()
    validate_aln = join(projdir,"results/DM0_IS150_aligned.fasta")
    with open(validate_aln, "w") as handle:
        handle.write(stdout)

    stdout, stderr = cline_full()
    full_aln = join(projdir,"results/all_IS150_aligned.fasta")
    with open(full_aln, "w") as handle:
        handle.write(stdout)

main()
