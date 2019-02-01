'''
deal-with-gdiffs.py by Rohan Maddamsetti.

I first create a directory structure to analyze breseq output for Zach's
DM0 and DM25 Cit+ evolution experiments. Then, I examine the evolved genomes
for possible contamination or cross-contamination. Then, I take the evolved genomes,
with ancestral mutations subtracted out, and write out a table for make-figures.R
to plot the counts and types of mutations in each population.
'''

from os.path import join, basename, exists, isfile, expanduser
from os import listdir, getcwd, makedirs, walk
import pathlib
import pandas as pd
from subprocess import run
from itertools import product,combinations
from pprint import pprint

## add genomediff parser to path.
import sys
sys.path.append("external/genomediff-python/")
import genomediff

class CallGDtools:

    def __init__(self,calltype,outfile,infiles=[],mut_type=None,conditions=None,overwrite=False):
        if calltype not in ['INTERSECT','SUBTRACT', 'REMOVE','UNION']:
            print('ERROR: NOT IMPLEMENTED.')
            quit()
        self.calltype = calltype
        self.outfile = outfile
        self.infiles = infiles
        self.overwrite = overwrite
        self.args = ['gdtools', self.calltype, '-o', self.outfile] + self.infiles
        if calltype == 'REMOVE':
            self.conditions = conditions
            if mut_type is not None:
                if mut_type not in ['SNP','DEL','MOB','SUB', 'INS', 'AMP']:
                    print('ERROR: unknown mut_type:',mut_type)
                    quit()
                self.mut_type = mut_type
                self.args = ['gdtools', self.calltype, '-o', self.outfile, '-m', mut_type] + [x for pair in zip(['-c']*len(conditions), conditions) for x in pair] + self.infiles
            else:
                self.args = ['gdtools', self.calltype, '-o', self.outfile] + [x for pair in zip(['-c']*len(conditions),conditions) for x in pair] + self.infiles
        self.cmd = ' '.join(self.args)

    def call(self):
        if not self.overwrite and isfile(self.outfile):
            print(self.outfile,"already exists-- not overwritten.")
            return None
        else:
            return run(self.cmd,shell=True,executable='/bin/bash')

    def __str__(self):
        return self.cmd

def contamination_check(pathdict, pop_clone_labels,resultsdir):
    '''
    Assert that each clone contains all mutations that occur in its immediate
    parent strain.

    ignore mutations in the known troublesome loci: 23S rRNA (rrlA), rhsA, rhsB,
    and insE-3 and insE-5. These will be filtered out anyway.
    '''
    
    ## make sure that all basenames of the paths are in pop_clone_labels.    
    assert set(pop_clone_labels.Name) == set([x for x in pathdict])

    for k in pathdict.keys():
        child_gd = join(pathdict[k],"output/evidence/annotated.gd")
        parent = pop_clone_labels[pop_clone_labels['Name']==k].ParentClone.item()
        parent_gd = join(pathdict[parent],"output/evidence/annotated.gd")
        assert isfile(child_gd)
        assert isfile(parent_gd)

        child_diff = genomediff.GenomeDiff.read(open(child_gd, 'r', encoding='utf-8'))
        parent_diff = genomediff.GenomeDiff.read(open(parent_gd, 'r', encoding='utf-8'))

        ## remove the bad loci from the diffs, following gdtools REMOVE semantics.
        parent_diff.remove("gene_name==rrlA")
        parent_diff.remove("gene_name==rhsA")
        parent_diff.remove("gene_name==rhsB")
        parent_diff.remove("genes_overlapping==insE-3")
        parent_diff.remove("genes_overlapping==insE-5")

        child_diff.remove("gene_name==rrlA")
        child_diff.remove("gene_name==rhsA")
        child_diff.remove("gene_name==rhsB")
        child_diff.remove("genes_overlapping==insE-3")
        child_diff.remove("genes_overlapping==insE-5")

        for mut in parent_diff.mutations:
            try:
                assert mut in child_diff.mutations
            except AssertionError:
                ## look at the missing mutation. In a missing region of coverage?
                print(parent,"\n",mut,"\n",k,"\n")
                
def subtract_ancestors_and_filter_mutations(pathdict,pop_clone_labels, resultsdir):
    ''' First copy annotated genomediffs into a nice directory structure, filtering out
        mutations in the ancestral strains.'''
    evol_genomes = [x for x in pathdict.keys() if x.startswith('ZDBp')]
    evol_pop_clone_labels = pop_clone_labels[pop_clone_labels['Name'].isin(evol_genomes)]
    for evol in evol_genomes:
        anc = evol_pop_clone_labels[evol_pop_clone_labels['Name']==evol].ParentClone.item()
        evol_gd = join(pathdict[evol],"output/evidence/annotated.gd")
        anc_gd = join(pathdict[anc],"output/evidence/annotated.gd")
        assert isfile(evol_gd)
        assert isfile(anc_gd)
        infiles = [evol_gd, anc_gd]
        outname = evol+"_minus_"+anc+".gd"
        my_founder = evol_pop_clone_labels[evol_pop_clone_labels['Name']==evol].Founder.item()
        my_ara_status = evol_pop_clone_labels[evol_pop_clone_labels['Name']==evol].AraStatus.item()
        my_environment = evol_pop_clone_labels[evol_pop_clone_labels['Name']==evol].Environment.item()
        outdir = join(resultsdir,my_founder,my_ara_status,my_environment)
        tempf = join(outdir,"temp_"+outname)
        outf = join(outdir,outname)
        CallGDtools('SUBTRACT',tempf,infiles).call()

        ''' Now remove all mutation calls in:
        rrlA (23S ribosomal RNA)
        rhsA and rhsB (repeat-rich proteins)
        Q48Q mutation in insE-3 or insE-5 (uncertain genomic variation in insE transposon).
        '''
        all_conditions = ['gene_name==rrlA','gene_name==rhsA','gene_name==rhsB','gene_name==insE-3','gene_name==insE-5']
        while len(all_conditions): ## run gdtools REMOVE for each condition to remove all matches.
            clist = [all_conditions.pop(0)]
            ## add a digit to the temp outfile name to keep them distinct.
            clen = len(all_conditions)
            if clen == 0: ## in final iteration: use outf.
                tempoutf = outf
            else:
                tempoutf = join(outdir,"temp_"+str(clen)+'.gd')
            remove_cmd = CallGDtools('REMOVE',tempoutf,[tempf],conditions=clist)
            remove_cmd.call()
            ## remove the temporary infile.
            run("rm "+tempf,shell=True,executable='/bin/bash')
            tempf = tempoutf
        
def write_evolved_mutations(evol_pop_clone_labels, inputdir, outf):

    outfh = open(outf,"w")
    outfh.write("Clone,Mutation,Environment,Population\n")

    for path, subdirs, files in walk(inputdir):
        for f in [x for x in files if x.endswith('.gd')]:
            full_f = join(path, f)
    
            infh = open(full_f, 'r', encoding='utf-8')
            name = f.split('_')[0]
            my_row = evol_pop_clone_labels[evol_pop_clone_labels.Name == name].iloc[0]
            env = my_row['Environment']
            population = my_row['PopulationLabel']
            gd = genomediff.GenomeDiff.read(infh)
            muts = gd.mutations
            muttype = ''
            for rec in muts:
                if rec.type == 'SNP':
                    muttype = rec.attributes['snp_type']
                else:
                    muttype = rec.type
                outfh.write(','.join([name,muttype,env,population])+"\n")

def main():

    homedir = expanduser("~")
    projdir = join(homedir,"BoxSync/DM0-evolution")
    resultsdir = join(projdir,"results/genome-analysis")
    breseqoutput_dir = join(projdir,"genomes/consensus")

    pop_clone_label_f = join(projdir,"data/rohan-formatted","populations-and-clones.csv")
    all_pop_clone_labels = pd.read_csv(pop_clone_label_f)
    '''only keep rows for sequenced clones.'''
    all_pop_clone_labels = all_pop_clone_labels[all_pop_clone_labels['Sequenced']==1]

    genomes = [x for x in listdir(breseqoutput_dir) if x.startswith('ZDB') or x.startswith('CZB')]
    genome_paths = [join(breseqoutput_dir,x) for x in genomes]
    genome_path_dict = {x:y for x,y in zip(genomes, genome_paths)}

    ''' assert that all mutations in ancestral strains occur in evolved strains,
    ignoring mutations in these known troublesome regions: 
    rrlA (23S ribosomal RNA)
    rhsA and rhsB (repeat-rich proteins)
    Q48Q mutation in insE-3 or insE-5 (uncertain genomic variation in insE transposon).
    '''    
    contamination_check(genome_path_dict, all_pop_clone_labels,resultsdir)

    ''' Now subtract ancestral mutations from evolved strains, and remove mutations in
    rrlA, rhsA, rhsB, insE-3, insE-5 from the analysis.'''
    subtract_ancestors_and_filter_mutations(genome_path_dict,all_pop_clone_labels,resultsdir)

    ''' tabulate the kinds of mutations in the evolved mutations.'''
    evolved_genomes = [x for x in listdir(breseqoutput_dir) if x.startswith('ZDBp')]
    evol_pop_clone_labels = all_pop_clone_labels[all_pop_clone_labels['Name'].isin(evolved_genomes)]

    mut_table_outf = join(resultsdir,"mutation_types.csv")
    write_evolved_mutations(evol_pop_clone_labels, resultsdir, mut_table_outf)
        
    
main()
