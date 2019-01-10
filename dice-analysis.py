'''

dice-analysis.py by Rohan Maddamsetti

1) Generate the input directory structure needed for Dan Deatherage's
TEE.py to run on the Cit+ evolution experiment data.

2) Run TEE.py on the data two ways: DM0 versus DM25, and CZB151 versus CZB152 and CZB154.

It is necessary to have run contamination-check.py before running this program.
In particular, ancestral mutations need to be subtracted from the evolved clones,
and strong parallelism/cross-contamination needs to be removed to have a more
conservative analysis.

'''



'''
Command line used in the TEE analysis:
python TEE.py -n 1000000 -p 150 -dt -s annotated_gd/ -g ../REL606.gbk -e annotated_gd/REL1207.gd -ct 37c

parser = argparse.ArgumentParser(description="Read in gd files, assign mutations to genes, perform statistics where applicable.")
parser.add_argument("-n", "--number", type=int, default=10, help="number of randomizations to perform for statistical models. Set to 0 to skip some tests")
parser.add_argument("-g", "--genbank_ref", help="input genbank file of reference genome, must have blattner numbers as note")
parser.add_argument("-s", "--samples", help="folder containing 'output' folder with gd files to be randomized to detect enrichment")
parser.add_argument("-dt", "--directory_treatments", help="treatment types specified by file architecture", action="store_true")
parser.add_argument("-ct", "--control_treatment", help="specify single treatment to compare all others to for valid and total mutations")
parser.add_argument("--command_line_treatments", help="command line specification of treatments. THIS IS NOT FUNCTIONAL, must specify by directory")
parser.add_argument("-p", "--promoter", type=int, default=0, help="length of promoter to consider as 'valid' mutation default =0")
parser.add_argument("-e", "--excluded_mutations", help="gd file with list of mutations to exclude. Helpful for avoiding using gdtools subtract")
parser.add_argument("--pvalue", type=float, default=0.05, help="p value significance cut off")
parser.add_argument("-pw", "--pairwise", help="perform pairwise dice comparisons among all treatments", action="store_true")
args = parser.parse_args()
'''

from os.path import join ##, basename, exists
##from os import listdir, getcwd, makedirs
from subprocess import run

class CallTEE:
    def __init__(self,gd_dir,control_treatment=''):
        ##self.number = 1000000
        self.number = 10 ## for debugging
        self.promoter_length = 150
        self.genbank_ref = '../genomes/REL606.7.gbk'
        self.control_treatment = control_treatment
        self.samples = gd_dir
        self.exclude = join(gd_dir,'cross_contamination_union.gd')
        ## To get TEE.py working for now, include an excluded file.
        ## do all pairwise comparisons among treatments.
        ## remove control treatment flag for now.

        self.args = ['python', './external/Deatherage-analysis/DiceSimilarity/citrate_dice.py', '-pw', '-n', str(self.number),
                     '-p', str(self.promoter_length), '-dt', '-s', self.samples, '-g', self.genbank_ref, '-e', self.exclude , '-ct', self.control_treatment]

    def call(self):
        return run(self.args)

    def __str__(self):
        return ' '.join(self.args)

def setup_dir_structure():
    '''
    TODO: automatically copy annotated and subtracted
          evolved clone gds into the proper folders to
          run the Dice Similarity analysis. For expediency,
          do by hand at the moment.
'''
    pass

def parse_TEE_output():
    '''
    TODO: parse TEE.py output to make matrix plot. OR turn TEE.py into a module
    '''
    pass

def make_matrix_csv():
    '''
    TODO:
    '''

def main():
    projdir = "/Users/Rohandinho/Dropbox (HMS)/DM0-evolution"
    analysis_dir = join(projdir,"results/genome-analysis/Dice-similarity-analysis")

    do_environment = True
    if do_environment:
        gddir = join(analysis_dir,"environment-comparison")
        ctl_treat = 'DM25'
    else:
        gddir = join(analysis_dir,"genotype-comparison")
        ctl_treat = 'CZB151'

    ##DiceRun = CallTEE(gddir,ctl_treat)
    ##print(DiceRun)
    ##DiceRun.call()

    ## For comparison, run TEE.py on the LTEE to make a matrix for Fig. 1C.
    ## TODO: get rid of FAKE cross_contamination_union.gd for getting TEE.py running.
    ## I do some unnecessary formatting to get TEE.py running at the moment.
    ## NOTE: There's also a bug in the rpy2 calls, but I ignore for now since I only want the matrix.
    ## I should probably refactor citrate_dice.py so that the mutation matrix code runs separately from the statistics.
    ltee_gddir = join(projdir,"genomes/annotated-LTEE-50K-diffs")
    LTEE_DiceRun = CallTEE(ltee_gddir,'Ara+')
    print(LTEE_DiceRun)
    LTEE_DiceRun.call()

main()
