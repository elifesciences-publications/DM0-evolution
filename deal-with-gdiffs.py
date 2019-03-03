'''
deal-with-gdiffs.py by Rohan Maddamsetti.

I first create a directory structure to analyze breseq output for Zach's
DM0 and DM25 Cit+ evolution experiments. Then, I examine the evolved genomes
for possible contamination or cross-contamination. Then, I take the evolved genomes,
with ancestral mutations subtracted out, and write out a table for make-figures.R
to plot the counts and types of mutations in each population.

This script also makes table of the IS insertions in the LTEE and in the MAE,
using the genomes in Jeff Barrick's LTEE-Ecoli GitHub repository.

'''

from os.path import join, basename, exists, isfile, expanduser, dirname
from os import listdir, getcwd, makedirs, walk
from pathlib import Path
import pandas as pd
from subprocess import run
from itertools import product,combinations
from pprint import pprint

## add genomediff parser to path.
import sys
sys.path.append("external/genomediff-python/")
import genomediff

class CallGDtools:

    def __init__(self,calltype,outfile,infiles=[],mut_type=None,conditions=None,overwrite=False,refgenome=None):
        if calltype not in ['INTERSECT','SUBTRACT', 'REMOVE','UNION','ANNOTATE']:
            print('ERROR: NOT IMPLEMENTED.')
            quit()
        self.calltype = calltype
        self.outfile = outfile
        self.infiles = infiles
        self.overwrite = overwrite
        self.refgenome = refgenome
        self.args = ['gdtools', self.calltype, '-o', self.outfile] + self.infiles

        if calltype == 'ANNOTATE': ## only genomediff output format implemented.
            self.args = ['gdtools', self.calltype, '-o', self.outfile] + ['-r', refgenome, '-f', 'GD'] + self.infiles
        elif calltype == 'REMOVE':
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
        return None

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

def write_evolved_mutations(evol_pop_clone_labels, inputdir, outf,polymorphism=False):
    ''' if polymorphism is True, then write out a Frequency field as well.'''
    outfh = open(outf,"w")
    if polymorphism:
        outfh.write("Clone,Mutation,Gene,Position,Frequency,Environment,Population\n")
    else:
        outfh.write("Clone,Mutation,Gene,Position,Environment,Population\n")

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
                pos = str(rec.attributes['position'])
                gene = rec.attributes['gene_name']
                if polymorphism:
                    freq = str(rec.attributes['frequency'])
                if rec.type == 'SNP':
                    muttype = rec.attributes['snp_type']
                else:
                    muttype = rec.type
                if polymorphism:
                    outfh.write(','.join([name, muttype, gene, pos, freq, env, population])+"\n")
                else:
                    outfh.write(','.join([name, muttype, gene, pos, env, population])+"\n")

def make_IS_insertion_tbl(evol_pop_clone_labels,inputdir, outf):

    outfh = open(outf,"w")
    outfh.write("Clone,IS_element,genome_start,genome_end,genes_promoter,genes_inactivated,Environment,Population,GeneName,GenePosition,Generation\n")

    for path, subdirs, files in walk(inputdir):
        for f in [x for x in files if x.endswith('.gd')]:
            full_f = join(path, f)
            infh = open(full_f, 'r', encoding='utf-8')
            name = f.split('_')[0]
            my_row = evol_pop_clone_labels[evol_pop_clone_labels.Name == name].iloc[0]
            env = my_row['Environment']
            population = my_row['PopulationLabel']
            generation = my_row['Generation']
            gd = genomediff.GenomeDiff.read(infh)
            for rec in gd.mutations:
                if rec.type != 'MOB':
                    continue
                try: ## get name of affected promoters.
                    affected_promoter = rec.genes_promoter
                except AttributeError:
                    affected_promoter = 'NA'
                try: ## get the name of inactivated genes.
                    inactivated = rec.genes_inactivated
                except AttributeError:
                    inactivated = 'NA'
                try:
                    gene_pos = rec.gene_position
                except AttributeError:
                    gene_pos = 'NA'
                outfh.write(','.join([name,
                                      rec.repeat_name,
                                      str(rec.position_start),
                                      str(rec.position_end),
                                      affected_promoter,
                                      inactivated,
                                      env,
                                      population,
                                      rec.gene_name,
                                      gene_pos,
                                      str(generation)])+"\n")


def make_LCA_IS_insertion_csv(evol_pop_clone_labels, inputdir, outf, gff):
    outfh = open(outf,"w")
    outfh.write("IS_element,genome_start,genome_end\n")
    ''' parse IS elements in LCA.gff3 annotation file.'''
    for line in open(gff,'r'):
        line = line.strip()
        ldata = line.split('\t')
        if len(ldata) < 3: ## skip DNA sequence
            continue
        if ldata[2] not in ['CDS','repeat_region']:
            continue
        if 'IS' not in line:
            continue
        start_coord = ldata[3]
        end_coord = ldata[4]
        annotation_string = ldata[-1]
        annotations = annotation_string.split(';')
        is_element = ''
        ## two cases: repeat_region annotation produced by gdtools APPLY
        ## and IS element proteins in the reference genome
        if annotations[-1].endswith('repeat region'):
            IS_element_name = annotations[0].split('=')[1]
        elif annotations[-2].startswith('Note'):
            note_text = annotations[-2].split('=')[1].split()
            IS_element_name = [x for x in note_text if x.startswith('IS')].pop()
        else:
            continue
        outfh.write(','.join([IS_element_name, start_coord, end_coord]) + "\n")

def organize_diffs(analysisdir):
    ''' automatically organize files for the dice analysis.'''

    CZB151_paths = list(Path(join(analysisdir,'CZB151')).rglob("*.gd"))
    CZB152_paths = list(Path(join(analysisdir,'CZB152')).rglob("*.gd"))
    CZB154_paths = list(Path(join(analysisdir,'CZB154')).rglob("*.gd"))

    DM0_paths = [x for x in CZB151_paths+CZB152_paths+CZB154_paths if "DM0" in x.parts]
    DM25_paths = [x for x in CZB151_paths+CZB152_paths+CZB154_paths if "DM25" in x.parts]

    def cp_files(comparisontype,treatment,paths):
        for x in paths:
            y = join(analysisdir, comparisontype, treatment, x.name)
            my_args = ['cp',str(x),y]
            my_cmd = ' '.join(my_args)
            makedirs(dirname(y), exist_ok=True)
            run(my_cmd ,executable='/bin/bash',shell=True)

    cp_files("genotype-comparison","CZB151",CZB151_paths)
    cp_files("genotype-comparison","CZB152",CZB152_paths)
    cp_files("genotype-comparison","CZB154",CZB154_paths)
    cp_files("environment-comparison","DM0",DM0_paths)
    cp_files("environment-comparison","DM25",DM25_paths)

def make_LTEE_MAE_IS150_tbl(inputdir, outf):

    outfh = open(outf,"w")
    outfh.write("Clone,IS_element,Position,Environment,Population,GeneName,GenePosition,Generation\n")

    for path, subdirs, files in walk(inputdir):
        for f in [x for x in files if x.endswith('.gd')]:
            if 'Anc' in f: ## skip ancestral strains
                continue
            full_f = join(path, f)
            infh = open(full_f, 'r', encoding='utf-8')
            fname = f.split('.')[0]
            if fname.startswith('JEB'):
                env = 'MAE'
                name = fname
                population = fname
                ## metadata says TIME=14300, but Tenaillon paper says ~13,750 generations.
                generation = '13750'
            else: ## LTEE strain.
                env = 'LTEE'
                population, gen, name = fname.split('_')
                assert gen.endswith('gen')
                generation = gen[:-3]
            gd = genomediff.GenomeDiff.read(infh)
            for rec in gd.mutations:
                if rec.type != 'MOB':
                    continue
                if rec.repeat_name != 'IS150':
                    continue
                try:
                    gene_pos = rec.gene_position
                except AttributeError:
                    gene_pos = 'NA'
                outfh.write(','.join([name,
                                      rec.repeat_name,
                                      str(rec.position),
                                      env,
                                      population,
                                      rec.gene_name,
                                      gene_pos,
                                      generation])+"\n")

def annotate_LTEE_MAE_genomes(LTEE_Ecoli_dir, outdir, refgen):
    ''' run gdtools ANNOTATE on genomes in LTEE-Ecoli directory.'''
    for path, subdirs, files in walk(LTEE_Ecoli_dir):
        for f in [x for x in files if x.endswith('.gd')]:
            full_f = join(path, f)
            out_f = join(outdir,f)
            annotate_cmd = CallGDtools('ANNOTATE',out_f,[full_f],refgenome=refgen)
            annotate_cmd.call()


def main():

    homedir = expanduser("~")
    projdir = join(homedir,"BoxSync/DM0-evolution")
    resultsdir = join(projdir,"results/genome-analysis")
    poly_resultsdir = join(projdir,"results/genome-analysis/polymorphism")
    breseqoutput_dir = join(projdir,"genomes/consensus")
    poly_breseqoutput_dir = join(projdir,"genomes/polymorphism")

    pop_clone_label_f = join(projdir,"data/rohan-formatted","populations-and-clones.csv")
    all_pop_clone_labels = pd.read_csv(pop_clone_label_f)
    '''only keep rows for sequenced clones.'''
    all_pop_clone_labels = all_pop_clone_labels[all_pop_clone_labels['Sequenced']==1]

    genomes = [x for x in listdir(breseqoutput_dir) if x.startswith('ZDB') or x.startswith('CZB')]
    genome_paths = [join(breseqoutput_dir,x) for x in genomes]
    genome_path_dict = {x:y for x,y in zip(genomes, genome_paths)}
    poly_genome_paths = [join(poly_breseqoutput_dir,x) for x in genomes]
    poly_genome_path_dict = {x:y for x,y in zip(genomes, poly_genome_paths)}

    ''' assert that all mutations in ancestral strains occur in evolved strains,
    ignoring mutations in these known troublesome regions:
    rrlA (23S ribosomal RNA)
    rhsA and rhsB (repeat-rich proteins)
    Q48Q mutation in insE-3 or insE-5 (uncertain genomic variation in insE transposon).
    '''
    ##contamination_check(genome_path_dict, all_pop_clone_labels,resultsdir)

    ''' Now subtract ancestral mutations from evolved strains, and remove mutations in
    rrlA, rhsA, rhsB, insE-3, insE-5 from the analysis. Then, organize diffs for the
    Dice analysis. '''
    ##subtract_ancestors_and_filter_mutations(genome_path_dict,all_pop_clone_labels,resultsdir)
    ##organize_diffs(resultsdir)
    ''' do the same thing,using the breseq --polymorphism mode results.'''
    ##subtract_ancestors_and_filter_mutations(poly_genome_path_dict,all_pop_clone_labels,poly_resultsdir)
    ##organize_diffs(poly_resultsdir)

    ''' tabulate the evolved mutations.'''
    evolved_genomes = [x for x in listdir(breseqoutput_dir) if x.startswith('ZDBp')]
    evol_pop_clone_labels = all_pop_clone_labels[all_pop_clone_labels['Name'].isin(evolved_genomes)]
    genomesdir = join(resultsdir,"environment-comparison")
    mut_table_outf = join(resultsdir,"evolved_mutations.csv")
    ##write_evolved_mutations(evol_pop_clone_labels, genomesdir, mut_table_outf)

    poly_genomesdir = join(poly_resultsdir,"environment-comparison")
    poly_mut_table_outf = join(resultsdir,"poly_evolved_mutations.csv")
    ##write_evolved_mutations(evol_pop_clone_labels, poly_genomesdir, poly_mut_table_outf,polymorphism=True)

    ''' make a table of IS-insertions, genome, and start and end locations.
        in the evolved and subtracted clones.'''
    IS_table_outf = join(resultsdir,"IS_insertions.csv")
    make_IS_insertion_tbl(evol_pop_clone_labels, genomesdir, IS_table_outf)

    ''' make a table of IS elements in LCA (including REL606).'''
    LCA_table_outf = join(resultsdir,"LCA_IS_insertions.csv")
    LCAgff = join(projdir,"genomes/curated-diffs/LCA.gff3")
    ##make_LCA_IS_insertion_csv(evol_pop_clone_labels, resultsdir, LCA_table_outf, LCAgff)

    ''' run gdtools ANNOTATE on the LTEE and MAE genomes.'''
    LTEE_MAE_dir = join(projdir,"genomes/LTEE-Ecoli")
    annotated_LTEE_MAE_dir = join(resultsdir,"annotated-LTEE-Ecoli")
    REL606ref = join(projdir,"genomes/REL606.7.gbk")
    if not exists(annotated_LTEE_MAE_dir):
        makedirs(annotated_LTEE_MAE_dir)
    ##annotate_LTEE_MAE_genomes(LTEE_MAE_dir,annotated_LTEE_MAE_dir,REL606ref)


    ''' make a table of IS-insertions, genome, and start and end locations in the
    LTEE and MAE clones. '''

    LTEE_MAE_IS_outf = join(resultsdir, "LTEE_MAE_IS150_insertions.csv")
    make_LTEE_MAE_IS150_tbl(annotated_LTEE_MAE_dir, LTEE_MAE_IS_outf)

main()
