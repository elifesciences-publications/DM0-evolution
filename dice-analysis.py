'''

dice-analysis.py by Rohan Maddamsetti

1) Generate the input directory structure needed for Dan Deatherage's
TEE.py to run on the Cit+ evolution experiment data.

2) Run cirate_dice.py on the data two ways: DM0 versus DM25, and CZB151 versus CZB152 and CZB154.

It is necessary to have run deal-with-diffs.py before running this program.

Command line used in the Daniel Deatherage's TEE analysis:
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
parser.add_argument("-m", "--matrixfile", help='name of file to save matrix of valid mutations in genomes')
parser.add_argument("--dN_only", help='only non-synonymous mutations will be analyzed',action="store_true")
args = parser.parse_args()
'''

from pathlib import Path
from os import makedirs
from os.path import join, expanduser
from subprocess import run
import pandas as pd

class CallDICE:
    def __init__(self, gd_dir, genbank_ref, mutmatrixfile, control_treatment='', excluded=None,dN_only=False):
        self.number = 100 #1000000 ## do one million for final submission.
        self.promoter_length = 150
        self.genbank_ref = genbank_ref
        self.control_treatment = control_treatment
        self.samples = gd_dir
        self.mutmatrixfile = mutmatrixfile
        self.excluded = excluded
        self.dN_only = dN_only
        self.args = ['python', './citrate_dice.py', '-pw', '-n', str(self.number),
                     '-p', str(self.promoter_length), '-dt', '-s', self.samples, '-g', self.genbank_ref, '-ct', self.control_treatment, '--matrixfile', self.mutmatrixfile]
        if self.excluded is not None:
            self.args = self.args + ['-e', self.excluded]
        if self.dN_only:
            self.args = self.args + ['--dN_only']

    def call(self):
        return run(self.args)

    def __str__(self):
        return ' '.join(self.args)

def main():
    homedir = expanduser("~")
    projdir = join(homedir,"BoxSync/DM0-evolution")
    analysis_dir = join(projdir,"results/genome-analysis")

    do_environment = True
    if do_environment:
        gddir = join(analysis_dir,"environment-comparison")
        ctl_treat = 'DM25'
    else: ## No significant difference in targets of selection across parental genotypes.
        gddir = join(analysis_dir,"genotype-comparison")
        ctl_treat = 'CZB151'

    print("********************************** DM0/DM25 Dice analysis **********************************")
    DiceRun = CallDICE(gddir,'../genomes/curated-diffs/LCA.gbk', '../results/DM0-DM25-comparison-mut-matrix.csv', ctl_treat)
    print(DiceRun)
    ##DiceRun.call()

    ## For comparison, run citrate_dice.py on the LTEE to make a matrix for Fig. 1C.
    print("********************************** LTEE Dice analysis for valid mutations **********************************")
    ltee_gddir = join(projdir,"genomes/annotated-LTEE-50K-diffs")
    LTEE_DiceRun = CallDICE(ltee_gddir,'../genomes/REL606.7.gbk','../results/LTEE-mut_matrix.csv','Ara+')
    print(LTEE_DiceRun)
    ##LTEE_DiceRun.call()

    ## run citrate_dice.py on just the 50K Ara-3 clone to get valid mutations for
    ## the Recapitulation Index analysis.
    print("********************************** 50K Ara-3 clone run for valid mutations  **********************************")
    ara_minus3_gddir = join(projdir,"genomes/curated-diffs")
    ara_minus3_DiceRun = CallDICE(ara_minus3_gddir,'../genomes/REL606.7.gbk','../results/ara3_mut_matrix.csv')
    ##ara_minus3_DiceRun.call()

    ## just analyse dN in the DM0/DM25 treatments.
    print("********************************** DM0/DM25 Dice analysis: dN only **********************************")
    dNDiceRun = CallDICE(gddir,'../genomes/curated-diffs/LCA.gbk', '../results/dN-DM0-DM25-comparison-mut-matrix.csv', ctl_treat, dN_only=True)
    print(dNDiceRun)
    dNDiceRun.call()

    ##  dN in the LTEE.
    print("********************************** LTEE Dice analysis for valid dN **********************************")
    dN_LTEE_DiceRun = CallDICE(ltee_gddir,'../genomes/REL606.7.gbk','../results/dN-LTEE-mut_matrix.csv','Ara+', dN_only=True)
    print(dN_LTEE_DiceRun)
    dN_LTEE_DiceRun.call()

    ## run 50K Ara-3 clone to get valid dN for the Recapitulation Index analysis.
    print("********************************** 50K Ara-3 clone for valid dN  **********************************")
    dN_ara_minus3_DiceRun = CallDICE(ara_minus3_gddir,'../genomes/REL606.7.gbk','../results/dN-ara3_mut_matrix.csv', dN_only=True)
    dN_ara_minus3_DiceRun.call()

main()
