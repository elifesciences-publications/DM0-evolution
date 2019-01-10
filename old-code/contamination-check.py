
from os.path import join, basename, exists
from os import listdir, getcwd, makedirs
import pandas as pd
from subprocess import run
from itertools import product,combinations

class CallGDtools:

    def __init__(self,calltype,outfile,infiles=[],mut_type=None,conditions=None):
        if calltype not in ['INTERSECT','SUBTRACT', 'REMOVE','UNION']:
            print('ERROR: NOT IMPLEMENTED.')
            quit()
        self.calltype = calltype
        self.outfile = outfile
        self.infiles = infiles
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

    def call(self):
        return run(self.args)

    def __str__(self):
        return ' '.join(self.args)

def setup_results_dir_structure(pop_clone_labels,outdir):
    '''
    make a directory tree for the analysis.
    within
    '''
    dirs1 = sorted(set(pop_clone_labels.Founder))
    dirs2 = sorted(set(pop_clone_labels['Ara Status']))
    dirs3 = sorted(set(pop_clone_labels.Environment))
    for x,y,z in product(dirs1,dirs2,dirs3):
        directory = join(outdir,x,y,z)
        if not exists(directory):
            makedirs(directory)

def check_for_secondary_mutations_and_contamination(pathdict,pop_clone_labels,outdir):
    '''
    Q1) Are there secondary mutations in Ara+ revertants?

    For strains grouped by their ancestor:
    function 1) take the intersection of the
    Ara+ strains, and the intersection of the Ara- strains.
    Then subtract X_Ara_Minus_intersection.gd from X_Ara_Plus_intersection.gd.
    This should only contain Ara+ marker there are no secondary mutations.

    Q2) Are all mutations found in ancestral strains found in their descendants?
    (test for outside contamination).

    For strains grouped by their ancestor:
    Take the intersection of all strains. This should be the genotype of the ancestor.

    TODO: USE SEQUENCED ANCESTRAL GENOMES!

    Input: pop_clone_labels := DataFrame of metadata for evolved clones.
           outdir: the subdirectory of this analysis, creatd by setup_results_dir_structure.

    Output: runs gdtools appropriately, and produces the output of those runs.
    '''
    ## make sure that all basenames of the paths are in pop_clone_labels.
    assert set(pop_clone_labels.Name) == set([x for x in pathdict])

    ## calculate the intersections for Ara+ and Ara- separately for each founder.
    the_founders = sorted(set(pop_clone_labels.Founder))
    for founder in the_founders:
        evolved_strains = pop_clone_labels[pop_clone_labels.Founder==founder]
        for arastatus in sorted(set(evolved_strains['Ara Status'])):
            ara_subgroup = evolved_strains[evolved_strains['Ara Status'] == arastatus]
            gdtools_input = [join(pathdict[x],"output/evidence/annotated.gd") for x in ara_subgroup.Name]
            outf = join(outdir,founder,'_'.join([founder,arastatus,'intersection.gd']))
            ara_subgroup_intersection = CallGDtools("INTERSECT",outf,gdtools_input)
            ara_subgroup_intersection.call()
        ## Now, intersect the Ara+ and Ara- intersections to infer the ancestral genome.
        plus_input = join(outdir,founder,'_'.join([founder,'Plus','intersection.gd']))
        minus_input = join(outdir,founder,'_'.join([founder,'Minus','intersection.gd']))
        plus_and_minus = [plus_input,minus_input]
        outf = join(outdir,founder,'_'.join([founder,'inferred.gd']))
        get_ancestor = CallGDtools("INTERSECT",outf,plus_and_minus)
        get_ancestor.call()
    ## finally, infer the last common ancestor of CZB151, 152, 154.
    outf = join(outdir,'LCA_inferred.gd')
    infiles = [join(outdir,x,'_'.join([x,'inferred.gd'])) for x in the_founders]
    get_LCA = CallGDtools("INTERSECT",outf,infiles)
    get_LCA.call()
    ## Now subtract the minus outputs from the plus outputs.
    plus_anc_dict = {'CZB151':'ZDB67','CZB152':'ZDB68','CZB154':'ZDB69'}
    for founder in the_founders:
        outf = join(outdir,founder,'_'.join([plus_anc_dict[founder],'secondary_muts.gd']))
        infiles = [join(outdir,founder,founder+"_Plus_intersection.gd"),
                   join(outdir,founder,founder+"_Minus_intersection.gd")]
        subtract_call = CallGDtools("SUBTRACT",outf,infiles)
        subtract_call.call()
    ## Now subtract the LCA output from the minus output.
    for founder in the_founders:
        outf = join(outdir,founder,founder+'_unique_markers.gd')
        infiles = [join(outdir,founder,founder+"_inferred.gd"),
                   join(outdir,'LCA_inferred.gd')]
        unique_call = CallGDtools("SUBTRACT",outf,infiles)
        unique_call.call()
    ## WARNING!!! There is an apparent bug in that the markers are not truly unique.
    ##            I need to look more closely at this.

def check_for_cross_contamination(pathdict, pop_clone_labels,outdir):
    '''
    Q3) Do any evolved mutations occur in any other evolved backgrounds?
    (test for cross-contamination).

    1): run gd subtract on all genomes, subtracting all mutations found
    in the ancestral clones, including the Ara marker and any secondary mutations.

    2): take subtracted files, and calculate the intersection of each pair.
                any intersection is probable contamination, but could also be a sign
                of strong convergence (unlikely).
    '''
    ## make sure that all basenames of the paths are in pop_clone_labels.
    assert set(pop_clone_labels.Name) == set([x for x in pathdict])

    subtract_files = []
    for clone,arastatus,founder,environment in zip(pop_clone_labels.Name,
                                                   pop_clone_labels['Ara Status'],
                                                   pop_clone_labels.Founder,
                                                   pop_clone_labels.Environment):
        print(clone,arastatus,founder,environment)
        outf = join(outdir,founder,arastatus,environment,clone+"_subtracted.gd")
        print(outf)
        clone_gd = join(pathdict[clone],"output/evidence/annotated.gd")
        most_recent_anc_gd = join(outdir,founder,'_'.join([founder,arastatus,"intersection.gd"]))
        infiles = [clone_gd,most_recent_anc_gd]

        subtract_call = CallGDtools('SUBTRACT',outf,infiles)
        subtract_call.call()
        subtract_files.append(outf)
    ## INTERSECT each pair of subtract_files to check for cross-contamination,
    ## or alternatively, convergence on the molecular level.
    pair_intersect_outdir = join(outdir,'check-cross-contamination')
    if not exists(pair_intersect_outdir):
        makedirs(pair_intersect_outdir)
    for x,y in combinations(subtract_files,2):
        clone1 = basename(x).split('_')[0]
        clone2 = basename(y).split('_')[0]
        my_out = join(pair_intersect_outdir,'_'.join([clone1,clone2,'intersected.gd']))
        CallGDtools('INTERSECT',my_out,[x,y]).call()

    ## what is going on with the possible cross-contamination here?
    ## The safe thing to do is to take the UNION, and the SUBTRACT this
    ## from all of the evolved *.gd files with ancestral mutations already
    ## subtracted, and write output into filtered-annotated-gd dir.
    unionfile = join(outdir,"cross_contamination_union.gd")
    contaminfiles = [join(pair_intersect_outdir,x) for x in listdir(pair_intersect_outdir) if x.endswith('.gd')]
    CallGDtools("UNION",unionfile,contaminfiles).call()
    return unionfile

def filter_cross_contamination(pop_clone_labels,outdir,unionfile):
    '''
    Use gdtools SUBTRACT to get rid of mutations found in independently evolved clones.
    This might get rid of some true positives, but this is the safe thing to do right now.
    '''

    dirs1 = sorted(set(pop_clone_labels.Founder))
    dirs2 = sorted(set(pop_clone_labels['Ara Status']))
    dirs3 = sorted(set(pop_clone_labels.Environment))
    for x,y,z in product(dirs1,dirs2,dirs3):
        directory = join(outdir,x,y,z)
        for i in [f for f in listdir(directory) if f.endswith('gd')]:
            subtract_from_me_path = join(directory,i)
            outf = join(directory,i.split('_')[0] +'.gd')
            infiles = [subtract_from_me_path,unionfile]
            CallGDtools('SUBTRACT',outf,infiles).call()

def subtract_true_ancestors(pathdict,unionfile,analysisdir):
    '''
    hack function to subtract genome diffs of CZB152 and CZB154 from the
    cross-contamination-union file.
    '''
    infiles = [unionfile]
    for anc in pathdict:
        anc_gd = join(pathdict[anc],"output/evidence/annotated.gd")
        infiles.append(anc_gd)
    outf = join(analysisdir,'cross_contamination_minus_ancs.gd')
    CallGDtools('SUBTRACT',outf,infiles).call()


def contamination_check():
    '''
    Note: this is the old main function.
    contamination-check.py by Rohan Maddamsetti.
    
    1) Check for secondary mutations in Ara+ ancestor clones (ZDB67,68,69)
    2) Check for outside contamination (all ancestral mutations found in evolved clones)
    3) Check for cross-contamination (mutations shared between independent evolved clones)
    
    WARNING: I need to debug some weirdness associated with 'unique markers' for the
         different ancestral strains which aren't really unique!
    
    WARNING!!!THERE IS WAY MORE CONVERGENCE THAN EXPECTED!! IS THIS CONTAMINATION??
    
    I discovered that the pykF ancestral mutation is missing from ZDBp898 due to lack of coverage.
    So, lack of coverage issues would explain true positives present in some but not all clones.
    
    TODO: rewrite this code after I have the sequenced genomes of all ancestral strains.
      For now, subtract ancestral 152 and 154 mutations from 'cross_contamination_union.gd' and see
      what mutations are left.
    '''

    projdir = "/Users/Rohandinho/Dropbox (HMS)/DM0-evolution"
    outdir = join(projdir,"results/genome-analysis")
    breseq_outdir = join(projdir,"genomes/consensus")

    pop_clone_label_f = join(projdir,"data/rohan-formatted","populations-and-clones.csv")
    all_pop_clone_labels = pd.read_csv(pop_clone_label_f)

    evol_genomes = [x for x in listdir(breseq_outdir) if x.startswith('ZDBp')]
    ## filter pop_clone_labels by all_genomes.
    evol_pop_clone_labels = all_pop_clone_labels[all_pop_clone_labels['Name'].isin(evol_genomes)]

    anc_genomes = [x for x in listdir(breseq_outdir) if x.startswith('CZB')]
    anc_pop_clone_labels = all_pop_clone_labels[all_pop_clone_labels['Name'].isin(anc_genomes)]

    ## set up file structure for the genome analysis.
    ###setup_results_dir_structure(evol_pop_clone_labels,outdir)

    evol_genome_paths = [join(breseq_outdir,x) for x in evol_genomes]
    evol_path_dict = {x:y for x,y in zip(evol_genomes,evol_genome_paths)}

    ## Answer Q1 and Q2.
    ###check_for_secondary_mutations_and_contamination(pop_clone_labels,outdir)
    ## TODO: subtract ALL ancestral genomes, not just CZB152 and CZB154!
    ###check_for_secondary_mutations_and_contamination(evol_path_dict,evol_pop_clone_labels,outdir)

    ## Answer Q3.
    ###unionfile = check_for_cross_contamination(evol_path_dict,evol_pop_clone_labels,outdir)

    ## filter cross-contamination.
    ###filter_cross_contamination(evol_pop_clone_labels,outdir,unionfile)

    anc_genome_paths = [join(breseq_outdir,x) for x in anc_genomes]
    anc_path_dict = {x:y for x,y in zip(anc_genomes,anc_genome_paths)}


    ## hack to subtract ancestor genomes.
    subtract_true_ancestors(anc_path_dict,unionfile,outdir)
    
    
def main():
    

main()
