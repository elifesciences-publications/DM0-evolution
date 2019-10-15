#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gene centric counts and statistics.
Requires biopython
@author: Rohan Maddamsetti and Daniel Deatherage.

Command line used in Daniel's original TEE analysis:
python TEE.py -n 1000000 -p 150 -dt -s annotated_gd/ -g ../REL606.gbk -e annotated_gd/REL1207.gd -ct 37c

"""
import datetime
import copy
import os
import re
import random
import argparse
import itertools
import numpy
from Bio import SeqIO
from scipy.stats import mannwhitneyu, wilcoxon, kruskal, binom_test, fisher_exact
from collections import Counter
import sys
from pprint import pprint
from tqdm import trange

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


assert args.genbank_ref is not None, "no reference file provided"

if args.directory_treatments is False and args.command_line_treatments is None:
    print("********************************\nWARNING!!!!!!!\nTreatments not specified some functionality will be skipped\n********************************")

print("Number of randomizations: %s" % args.number)
print("Reference file to use: %s" % args.genbank_ref)

GenomeSeq = SeqIO.read(open(args.genbank_ref, "r"), "genbank")


def treatment_read_in():
    temp_dict = {"Total": {}}
    excluded_mutation_list = []
    if args.excluded_mutations is not None:
        with open(args.excluded_mutations) as F:
            for line in F:
                line = line.rstrip().split("\t")
                if re.match("^[A-Z]{3}$", line[0]):  # ignore all non-mutations
                    excluded_mutation_list.append(line)  # store entire line as something to be excluded, note this is stored as list not string

    for treatment in os.listdir(args.samples):  # expects to find directories as treatment labels
        treatment_path = "%s/%s" % (args.samples, treatment)
        if os.path.isdir(treatment_path):  # only look at directories within treatment labels
            temp_dict[treatment] = {}
            for sample in os.listdir(treatment_path):
                sample_path = "%s/%s" % (treatment_path, sample)
                if os.path.isdir(sample_path):  # treatment directories contain directories suggest being run immediately after
                    gd_file = sample_path + "/evidence/annotated.gd"
                    if not os.path.isfile(gd_file):
                        continue
                elif re.findall(".*\.gd", sample) and os.path.isfile(sample_path):  # Program is assuming .gd files are annotated. If not, will produce errors
                    gd_file = sample_path
                else:  # restricted to files, and directories containing annotated gd files
                    continue
                with open(gd_file, "r") as F:
                    mut_list = []
                    for line in F:
                        to_exclude = False
                        line = line.rstrip().split("\t")
                        if re.match("^[A-Z]{3}$", line[0]):
                            for mutation in excluded_mutation_list:  # check sample mutation line against exclusion list
                                if mutation[0] == line[0] and mutation[4] == line[4] and mutation[5] == line[5]:
                                    to_exclude = True
                            if not to_exclude:
                                mut_list.append(line)
                    temp_dict[treatment][sample.replace(".gd", "")] = {"raw_mutations": mut_list}
                    temp_dict["Total"][sample.replace(".gd", "")] = {"raw_mutations": mut_list}
    if args.control_treatment:
        assert args.control_treatment in temp_dict, "control treatment (%s) does not match treatments given (%s)" % (args.control_treatment, list(temp_dict.keys()))
        print("Comparisons to: %s as control treatment will be made." % args. control_treatment)
    return temp_dict


def mutation_gene_assignment(mutation_list, exclusion_list,only_dN=False):
    """Pass in 2 arguments, list of list of mutations (strings split by tabs), and a list of genomic regions to exclude, by default this is regions flagged as repetitive"""
    intergenic_and_multi_gene = []
    valid_mutations = []
    for mutation in mutation_list:
        assigned = False
        mutation[4] = int(mutation[4]) - 1  # pull mutation location and adjust for python starts numbering at 0
        mut_start = int(mutation[4])
        if str(mutation[0]) == "AMP" or str(mutation[0]) == "DEL":  # need to check that amps/del only effect a single gene
            mut_stop = int(mutation[4]) + int(mutation[5])  # add length of amp/del to get stop location
        else:
            mut_stop = int(mutation[4])  # needed for promoter comparison
        #######
        # block to check if mutation within repeat region
        exclude_mut = False  # by default we want to include the mutation
        for region in exclusion_list:  # check each excluded region against the current mutation
            if (region[0] <= mut_start <= region[1]) or (region[0] <= mut_stop <= region[1]):
                exclude_mut = True
                break  # if it overlaps any excluded region stop checking the rest as it doesn't matter if it overlaps multiples
        if exclude_mut:
            x = mutation_line(GenomeSeq.features[1], mutation)  # need placeholder feature to avoid errors
            x[2] = "Repeat_Region"  # need updating if exclusion list to include more than repeat regions
            intergenic_and_multi_gene.append(x)
            continue
        ####### END block
        for feature in GenomeSeq.features:  # loop each mutation through whole genome
            if feature.type == 'gene':
                if feature.location._start.position <= mut_start and mut_stop <= feature.location._end.position:
                    new_line = mutation_line(feature, mutation)
                    if 'snp_type=synonymous' in mutation:
                        intergenic_and_multi_gene.append(new_line)
                    ## if only_dN, then only dN are valid.
                    elif only_dN and 'snp_type=nonsynonymous' not in mutation:
                        intergenic_and_multi_gene.append(new_line)
                    else:
                        valid_mutations.append(new_line)
                    assigned = True
                    break

        ## mutation not assigned to gene; only run if not only_dN
        if not assigned and not only_dN:
            promoter_distance_comparison = []
            ## loop each mutation through whole genome with promoter comparisons
            for feature in GenomeSeq.features:
                if feature.type == 'gene':
                    start_loc = feature.location._start.position
                    end_loc = feature.location._end.position
                    assert start_loc < end_loc, "%s start location greater than gene end location. promoter distance about to be miscalculated" % (feature.qualifiers['gene'][0])
                    if feature.strand == 1:
                        start_loc -= args.promoter  # to prove to yourself that this is correct:print "positive\t%s\t%s\t%s"%(feature.location._start.position, feature.location._end.position, start_loc)
                    elif feature.strand == -1:
                        end_loc += args.promoter  # to prove to yourself that this is correct:print "negative\t%s\t%s\t%s"%(feature.location._start.position, feature.location._end.position, end_loc)
                    if start_loc <= mut_start and mut_stop <= end_loc:
                        promoter_distance_comparison.append(mutation_line(feature, mutation))
                        promoter_distance = []
                        if feature.strand == 1:
                            promoter_distance.extend([mut_start - feature.location._start.position, mut_stop - feature.location._start.position])
                        elif feature.strand == -1:
                            promoter_distance.extend([feature.location._end.position - mut_start, feature.location._end.position - mut_stop])
                        promoter_distance_comparison[-1].append("promoter=" + str(min(promoter_distance)))

            if len(promoter_distance_comparison) == 1:
                valid_mutations.append(promoter_distance_comparison[0])
                assigned = True
            elif len(promoter_distance_comparison) > 1:
                minimal_distances = []
                for possible_nearest_match in promoter_distance_comparison:
                    x = int(possible_nearest_match[-1].replace("promoter=", ''))
                    minimal_distances.append(x)
                y = max(minimal_distances)  # NOTE that distances are promoters and therefore negative distances!
                for possible_nearest_match in promoter_distance_comparison:
                    x = int(possible_nearest_match[-1].replace("promoter=", ''))
                    if x == y:
                        valid_mutations.append(possible_nearest_match)
                        assigned = True
                assert assigned is True, "multiple possible promoter distances compared, but none considered valid"
        if not assigned:
            x = mutation_line(feature, mutation)
            x[2] = "NO_SINGLE_GENE"
            intergenic_and_multi_gene.append(x)
    return {"valid_mutations": valid_mutations, "intergenic_or_multi-gene": intergenic_and_multi_gene}

def gene_name_pull(feature):
    """Attempt to grab gene name, if gene name not found, grab locus_tag instead."""
    try:
        if len(feature.qualifiers['gene']) != 1:
            print(feature.qualifiers['gene'])
            assert False
        return feature.qualifiers['gene'][0]
    except KeyError:
        if len(feature.qualifiers['locus_tag']) != 1:
            print(feature.qualifiers['locus_tag'])
            assert False
        return feature.qualifiers['locus_tag'][0]


def mutation_line(feature, mutation):
    """collapse some of the repeated code used in going from annotated mutation line to concise line with gene names for other analysis"""
    line = [mutation[0], mutation[4], gene_name_pull(feature)]
    if mutation[0] == "SNP":
        assert len(re.findall("snp_type=", str(mutation))) == 1, "'snp_type=' error detected in:\n%s\nIf more than 1 don't know what to assign, if 0 would be reusing old snp_type calls" % (re.findall("snp_type=.*?'", str(mutation)))
        snp_type = re.findall("'snp_type=.*?'", str(mutation))[0].replace("'snp_type=", '').replace("'", "") + ":" + str(GenomeSeq.seq[mutation[4]]) + "_to_" + str(mutation[5])
        line.append(snp_type)
    elif mutation[0] == "MOB":
        line.append(mutation[5])
    elif mutation[0] == "INS" or mutation[0] == "SUB" or mutation[0] == "INV" or mutation[0] == 'CON':
        line.append(len(mutation[5]))
    elif mutation[0] == "AMP" or mutation[0] == "DEL":
        line.append(str(mutation[5]) + "bp")
    else:
        line.append("UNUSED")
    return line


def dice_calc(x, y):
    """calculates dice coefficients 2 lists of mutated gene names (x, y)"""
    x_and_y = [gene for gene in x if gene in y]
    if len(x) <= 0 or len(y) <=0:
        return -1
    dice_value = (2.0 * len(x_and_y)) / (len(x) + len(y))
    return dice_value


def dice_gen(master_dict, treatment_id_dict, excluded_genes, comparison_dictionary):
    """comparison_dictionary always will look at Grand: Mean, Within, and Between. If args.pairwise is used, then pairwise comparisons made as well additional functionality possible in future"""
    dice_dict = copy.deepcopy(comparison_dictionary)
    for combination in itertools.combinations(master_dict["Total"], 2):
        x = [mutation[2] for mutation in master_dict['Total'][combination[0]]['valid_mutations'] if mutation[2] not in excluded_genes]  # pull gene names for 1st of combination pair
        y = [mutation[2] for mutation in master_dict['Total'][combination[1]]['valid_mutations'] if mutation[2] not in excluded_genes]  # pull gene names for 2nd of combination pair
        z = dice_calc(x, y)
        treatment_identification = {"x": ["not_changed_", combination[0]], "y": ["not_changed_", combination[1]]}
        for treatment in treatment_id_dict:  # assign samples to treatment
            if combination[0] in treatment_id_dict[treatment]:
                treatment_identification["x"] = [treatment, combination[0]]
            if combination[1] in treatment_id_dict[treatment]:
                treatment_identification["y"] = [treatment, combination[1]]
        assert treatment_identification["x"][0] != "not_changed_" and treatment_identification["y"][0] != "not_changed_", "The following sample(s) not assigned a treatment %s" % ([x[1] for x in [treatment_identification["x"], treatment_identification["y"]] if x[0] == "not_changed_"])

        if treatment_identification["x"][0] == treatment_identification["y"][0]:
            dice_dict["Grand"]["Within"].append(z)
            dice_dict["Grand"]["Mean"].append(z)
            for treatment_comparison in dice_dict:
                if treatment_identification["x"][0] in treatment_comparison.split(" vs "):
                    dice_dict[treatment_comparison]["Within"].append(z)
                    dice_dict[treatment_comparison]["Mean"].append(z)
            # if z > 0: print combination, "\t", x, "\t", y, "\t", z  # Test to print only dice values when of between treatment groups that show some similarity
        else:
            dice_dict["Grand"]["Between"].append(z)
            dice_dict["Grand"]["Mean"].append(z)
            for treatment_comparison in dice_dict:
                if treatment_identification["x"][0] in treatment_comparison.split(" vs ") and treatment_identification["y"][0] in treatment_comparison.split(" vs "):
                    dice_dict[treatment_comparison]["Between"].append(z)
                    dice_dict[treatment_comparison]["Mean"].append(z)

    assert len(dice_dict["Grand"]["Mean"]) == len(dice_dict["Grand"]["Between"]) + len(dice_dict["Grand"]["Within"]), "separation of between and within not equaling total values"
    for comparison in dice_dict:
        for comparison_type in dice_dict[comparison]:
            # print comparison, "\t", comparison_type, "\t", dice_dict[entry]  # Test to print the raw values of dice calculations
            dice_dict[comparison][comparison_type] = numpy.mean(dice_dict[comparison][comparison_type])  # calculate the average of all the dice values
    return dice_dict


def dice(master_dict, excluded_genes):
    """perform relevant dice calculations with optional list of genes to exclude (ie exclude genes which were individually significant """

    padding = 23
    precision = padding - 3

    if len(excluded_genes) == 0:
        print("\nDice Calculations for all valid mutations")
    else:
        print("\nDice Calculations excluding: %s" % (str(excluded_genes).replace("[", "").replace("]", "")))

    treatment_id_dict = {}  # experimental treatment - [sample list]
    shuffled_treatments = {}  # randomized sample assignment to treatments to be repeated args.number of times
    comparison_dictionary = {"Grand": {"Mean": [], "Within": [], "Between": []}}
    for treatment in [x for x in list(master_dict.keys()) if not x == "Total"]:
        assert treatment != " vs ", " ' vs ' is used as separator for pairwise comparison, but it exists as one of the treatment types. This will cause problems with assignment. Please change vs treatment name to something else"
        treatment_id_dict[treatment] = [sample for sample in master_dict[treatment]]
        shuffled_treatments[treatment] = []  # all treatments must be represented samples to be assigned to each treatment at random later
    if len(treatment_id_dict) < 2:
        print('Warning: less than two treatments supplied. Dice calculations skipped.')
        return -1

    for pairwise_comparison_combination in itertools.combinations(list(treatment_id_dict.keys()), 2):
        comparison_dictionary[str(pairwise_comparison_combination[0] + " vs " + pairwise_comparison_combination[1])] = {"Mean": [], "Within": [], "Between": []}
    z1 = dice_gen(master_dict, treatment_id_dict, excluded_genes, comparison_dictionary)

    ###################
    # random DICE
    random_dice_dict = {}
    for comparison in list(comparison_dictionary.keys()):
        random_dice_dict[comparison] = {"Original_Diff": z1[comparison]['Within'] - z1[comparison]['Between'], "Random_Greater": 0}

    if args.number < 1: ## 0 or fewer comparisons. no point in running.
        return

    ## loop with progress bar.
    for _ in trange(args.number):
        random_sample_list = random.sample(list(master_dict["Total"]), len(master_dict["Total"]))
        random_start = 0
        for treatment in shuffled_treatments:
            random_stop = random_start + len(master_dict[treatment])
            shuffled_treatments[treatment] = random_sample_list[random_start: random_stop]
            random_start += len(master_dict[treatment])
        z2 = dice_gen(master_dict, shuffled_treatments, excluded_genes, comparison_dictionary)
        for comparison in list(comparison_dictionary.keys()):
            if z2[comparison]["Within"] - z2[comparison]["Between"] >= random_dice_dict[comparison]["Original_Diff"]:
                random_dice_dict[comparison]["Random_Greater"] += 1

    print("Calculating dice differences between pairs of treatments.\np-value represents number of times random sample assignment was seen more similar within treatments than between treatments divided by total number of randomizations.\nWith Bonferroni, p value is %f" % (float(args.pvalue / (len(treatment_id_dict) - 1))))  # -1 because treatment not compared to itself
    print(''.join([f'{x:<{padding}.{precision}}' for x in ["Comparison", "Within", "Between", "Mean", "p-value"]]))
    for comparison in sorted(list(comparison_dictionary.keys()), reverse=True):
        about_to_print = [comparison]
        for comparison_type in ["Within", "Between", "Mean"]:
            about_to_print.append(z1[comparison][comparison_type])
        about_to_print.append(random_dice_dict[comparison]["Random_Greater"] / float(args.number))
        print(''.join([f'{x:<{padding}.{precision}}' for x in map(str,about_to_print)]))

    # print "\nOf the %i randomizations to treatment sample assignment, the difference in dice coefficients of within treatment and between treatments was greater than the experimental %i times. This represents a p-value of %f\n" % (args.number, random_dice_more_similar, float(random_dice_more_similar) / args.number)
    ##########


def repeat_id():
    """return List of Lists for repeat regions"""
    big_list = []
    for feature in GenomeSeq.features:
        if feature.type == 'repeat_region':
            big_list.append([feature.location._start.position - 2, feature.location._end.position + 1])  # -2 for 0 space plus 1 for adjacent mutation (ie INS), +1 for adjacent base and INS
    return big_list


def general_mutation_stats(master_dict):
    """
    input: master_dict data structure.
    output: prints out breakdown of mutations over treatments and samples.
    """

    ##set the padding to 23 spaces for now.
    padding = 23
    ## set precision to 20
    precision = padding - 3

    temp_dict = {'Total_Treatment_Mutations': [['Treatment', "Total Mutations", "Valid Mutations", "Average Total", "Average Valid"]], 'Total_Sample_Mutations': [['Sample', 'Treatment', 'Total Mutations', 'Valid Mutations']]}
    for treatment in [x for x in master_dict if not x == "Total"]:
        treatment_print = [treatment, 0, 0]
        for sample in master_dict[treatment]:
            total = len(master_dict[treatment][sample]['total_mutations'])
            valid = len(master_dict[treatment][sample]["valid_mutations"])
            temp_dict['Total_Sample_Mutations'].append([sample] + [treatment] + [total] + [valid])
            treatment_print[1] += total
            treatment_print[2] += valid
        treatment_print.extend([float(treatment_print[1]) / len(master_dict[treatment]), float(treatment_print[2]) / len(master_dict[treatment])])
        temp_dict['Total_Treatment_Mutations'].append(treatment_print)

    print_last = ["Total:", 0, 0]
    for sample in master_dict["Total"]:
        print_last[1] += len(master_dict["Total"][sample]['total_mutations'])
        print_last[2] += len(master_dict["Total"][sample]['valid_mutations'])

    print_last.extend([float(print_last[1]) / len(master_dict["Total"]), float(print_last[2]) / len(master_dict["Total"])])
    for entry in temp_dict["Total_Treatment_Mutations"]:
        print(''.join([f"{x:<{padding}.{precision}}" for x in map(str,entry)]))

    print(''.join([f"{x:<{padding}.{precision}}" for x in map(str,print_last)]))

    if args.control_treatment:
        total_samples = 0
        total_exclude_compare = 0
        valid_exclude_compare = 0

        for treatment in [x for x in master_dict if not x == "Total"]:
            if treatment == args.control_treatment:
                continue
            for sample in master_dict[treatment]:
                total_exclude_compare += len(master_dict[treatment][sample]['total_mutations'])
                valid_exclude_compare += len(master_dict[treatment][sample]["valid_mutations"])
                total_samples += 1

        exclude_compare_print = ["Total excluding " + args.control_treatment + ":"]
        exclude_compare_print.extend([total_exclude_compare, valid_exclude_compare, float(total_exclude_compare) / total_samples, float(valid_exclude_compare) / total_samples])
        print(''.join([f"{x:<{padding}.{precision}}" for x in map(str, exclude_compare_print)]))

    print()  # blank line for better formatting
    for entry in temp_dict['Total_Sample_Mutations']:
        print(''.join([f"{x:<{padding}.{precision}}" for x in map(str, entry)]))

    min_max_dict = {
        "min_total": [min([x[2] for x in temp_dict['Total_Sample_Mutations']][1:])],
        "min_valid": [min([x[3] for x in temp_dict['Total_Sample_Mutations']][1:])],
        "max_total": [max([x[2] for x in temp_dict['Total_Sample_Mutations']][1:])],
        "max_valid": [max([x[3] for x in temp_dict['Total_Sample_Mutations']][1:])]
    }
    for entry in temp_dict['Total_Sample_Mutations']:
        if entry[2] == min_max_dict["min_total"][0]:
            min_max_dict["min_total"].append(entry[0])
        if entry[3] == min_max_dict["min_valid"][0]:
            min_max_dict["min_valid"].append(entry[0])
        if entry[2] == min_max_dict["max_total"][0]:
            min_max_dict["max_total"].append(entry[0])
        if entry[3] == min_max_dict["max_valid"][0]:
            min_max_dict["max_valid"].append(entry[0])

    print("\nMinimum and Maximum mutations detected per sample")
    for entry,val in min_max_dict.items():
     print(''.join([f"{x:<{padding}.{precision}}" for x in map(str,[entry] + val)]))


def mann_whitney_tests(master_dict):
    """Calculate nonparametric Mann_Whitney differences among treatments based on master_dict format"""

    padding = 23
    precision = padding - 3

    print("\nMann Whitney one and two tailed differences:")
    mann_whitney_dict = {'Total': {}, 'Valid': {}}
    for treatment in master_dict:
        if treatment == "Total":
            continue
        xtotal = []
        ytotal = []
        xvalid = []
        yvalid = []
        for sample in master_dict[treatment]:
            xtotal.append(len(master_dict[treatment][sample]['total_mutations']))
            xvalid.append(len(master_dict[treatment][sample]['valid_mutations']))
        for all_other_treatments in master_dict:
            if all_other_treatments == treatment or all_other_treatments == "Total":
                continue
            for sample in master_dict[all_other_treatments]:
                ytotal.append(len(master_dict[all_other_treatments][sample]['total_mutations']))
                yvalid.append(len(master_dict[all_other_treatments][sample]['valid_mutations']))
        mann_whitney_dict['Total'][treatment] = mann_whitney(xtotal, ytotal)
        mann_whitney_dict['Valid'][treatment] = mann_whitney(xvalid, yvalid)
    for entry in mann_whitney_dict:
        print("%s mutations: single test p = %s, Bonferroni corrected p = %s" % (entry, args.pvalue, float(args.pvalue / len(mann_whitney_dict[entry]))))
        print(''.join([f"{x:{padding}.{precision}}" for x in map(str,["Condition", "Less than rest", "Greater than rest", '2 tailed'])]))
        for condition in mann_whitney_dict[entry]:
            printing = [condition, mann_whitney_dict[entry][condition]['less'], mann_whitney_dict[entry][condition]['greater'], mann_whitney_dict[entry][condition]['both']]
            print(''.join([f"{x:{padding}.{precision}}" for x in map(str,printing)]))
        print("")
    if args.control_treatment:
        compare_total = []
        compare_valid = []
        mann_whitney_compare_dict = {'Total': {}, 'Valid': {}}
        for sample in master_dict[args.control_treatment]:
            compare_total.append(len(master_dict[args.control_treatment][sample]['total_mutations']))
            compare_valid.append(len(master_dict[args.control_treatment][sample]['valid_mutations']))
        for treatment in master_dict:
            if treatment == "Total" or treatment == args.control_treatment:
                continue
            treat_total = []
            treat_valid = []
            for sample in master_dict[treatment]:
                treat_total.append(len(master_dict[treatment][sample]['total_mutations']))
                treat_valid.append(len(master_dict[treatment][sample]['valid_mutations']))
            mann_whitney_compare_dict['Total'][treatment] = mann_whitney(compare_total, treat_total)
            mann_whitney_compare_dict['Valid'][treatment] = mann_whitney(compare_valid, treat_valid)
        for entry in mann_whitney_compare_dict:
            print("Comparing %s mutations from each treatment to %s: single test p = %s, Bonferroni corrected p = %s" % (entry, args.control_treatment, args.pvalue, float(args.pvalue / len(mann_whitney_compare_dict[entry]))))
            print(''.join([f"{x:{padding}.{precision}}" for x in map(str,["Comparison", "Control less", "Control greater", '2 tailed'])]))
            for comparison in mann_whitney_compare_dict[entry]:
                comparison_name = args.control_treatment + " Vs " + comparison
                print(''.join([f"{x:{padding}.{precision}}" for x in map(str,[comparison_name, mann_whitney_compare_dict[entry][comparison]['less'], mann_whitney_compare_dict[entry][comparison]['greater'], mann_whitney_compare_dict[entry][comparison]['both']])]))
            print("")


def mann_whitney(x, y):
    """Call scipy.stats Mann-Whitney function."""
    return {"less": mannwhitneyu(x,y,alternative='less').pvalue,
            "greater": mannwhitneyu(x,y,alternative='greater').pvalue,
            "both": mannwhitneyu(x,y,alternative='two-sided').pvalue }


def kruskal_wallis_tests(master_dict):
    """Calculate nonparametric Kruskal-Wallis test values among treatments based on master_dict format"""
    kruskal_test = []
    for treatment in master_dict:
        if treatment == "Total":
            continue
        kruskal_test.append(list(len(master_dict[treatment][sample]['total_mutations']) for sample in master_dict[treatment]))
    if len(kruskal_test) < 2:
        print('Warning: less than two groups. Kruskal-Wallis tests are being skipped.')
        return -1
    print("Total Mutation Kruskal Wallis test pvalue: %f" % (kruskal(*kruskal_test).pvalue))
    kruskal_test = []
    for treatment in master_dict:
        if treatment == "Total":
            continue
        kruskal_test.append(list(len(master_dict[treatment][sample]['valid_mutations']) for sample in master_dict[treatment]))
    print("Valid Mutation Kruskal Wallis test pvalue: %f\n" % (kruskal(*kruskal_test).pvalue))

    if args.control_treatment:
        multiple_experimental_treatments = len(set(master_dict.keys()) - {'Total',args.control_treatment}) > 1
        ## return immediately if there's only one experimental treatment.
        if not multiple_experimental_treatments:
            return
        kruskal_test = []
        for treatment in master_dict:
            if treatment == "Total" or treatment == args.control_treatment:
                continue
            kruskal_test.append(list(len(master_dict[treatment][sample]['total_mutations']) for sample in master_dict[treatment]))
        print("Total Mutation Kruskal Wallis excluding %s treatment test pvalue: %f" % (args.control_treatment, kruskal(*kruskal_test).pvalue))
        kruskal_test = []
        for treatment in master_dict:
            if treatment == "Total" or treatment == args.control_treatment:
                continue
            kruskal_test.append(list(len(master_dict[treatment][sample]['valid_mutations']) for sample in master_dict[treatment]))
        print("Valid Mutation Kruskal Wallis excluding %s treatment test pvalue: %f\n" % (args.control_treatment, kruskal(*kruskal_test).pvalue))

def types_of_mutations(master_dict):
    """print information of number and types of mutations"""
    padding = 23
    precision = padding - 3

    mutation_type_dict = {}
    for treatment in master_dict:
        mutation_type_dict[treatment] = {"SNPs": 0, "Synonymous": 0, "Non-Synonymous": 0, }
        for sample in master_dict[treatment]:
            mutation_type_dict[treatment]['SNPs'] += len(re.findall("SNP", str(master_dict[treatment][sample]['total_mutations'])))
            mutation_type_dict[treatment]['Synonymous'] += len(re.findall("'synonymous", str(master_dict[treatment][sample]['total_mutations'])))
            mutation_type_dict[treatment]['Non-Synonymous'] += len(re.findall("'nonsynonymous", str(master_dict[treatment][sample]['total_mutations'])))
    print("SNP information: ")
    print(''.join([f"{x:{padding}.{precision}}" for x in ["Condition", "SNPs", "Synonymous", "Non-Synonymous"]]))
    for treatment in mutation_type_dict:
        if treatment == "Total":
            continue
        print(''.join([f"{x:{padding}.{precision}}" for x in map(str, [treatment, mutation_type_dict[treatment]["SNPs"], mutation_type_dict[treatment]["Synonymous"], mutation_type_dict[treatment]["Non-Synonymous"]])]))

    print(''.join([f"{x:{padding}.{precision}}" for x in map(str, ["Total", mutation_type_dict["Total"]["SNPs"], mutation_type_dict["Total"]["Synonymous"], mutation_type_dict["Total"]["Non-Synonymous"]])]))

    print("\nSynonymous Mutation information(p value = binomial test if syn greater expected):")
    print(''.join(f'{x:{padding}.{precision}}' for x in ["Treatment", "Sample", "Total Mutations", "Synonymous SNP", "Non-Synonymous SNP"]))

    for treatment in master_dict:
        if treatment == "Total":
            continue
        for sample in master_dict[treatment]:
            if re.findall("'synonymous", str(master_dict[treatment][sample]['total_mutations'])):
                print(''.join([f"{x:{padding}.{precision}}" for x in map(str, [treatment, sample, len(master_dict[treatment][sample]['total_mutations']), str(master_dict[treatment][sample]['total_mutations']).count("'synonymous"), str(master_dict[treatment][sample]['total_mutations']).count("'nonsynonymous")])]))

    mutation_type_dict = {}
    for treatment in master_dict:
        mutation_type_dict[treatment] = { "DEL": [], "AMP": [], "INS": [], "MOB": [], "SUB": [],
                                         "INV": [], "CON":[] }
        for sample in master_dict[treatment]:
            for mutation in master_dict[treatment][sample]['total_mutations']:
                if mutation[0] == "SNP":
                    continue
                elif mutation[0] in ["MOB", "INS", "SUB", "INV", "CON"]:
                    mutation_type_dict[treatment][mutation[0]].append(mutation[3])
                elif mutation[0] in ["AMP", "DEL"]:
                    mutation_type_dict[treatment][mutation[0]].append(int(mutation[3].replace("bp", "")))
                else:
                    print(mutation)
                    assert False, "unused mutation"

    print()
    headings = ["Treatment", "DEL", "AMP", "INS", "SUB", "INV", "CON", "Total MOB"] + list(sorted(set(mutation_type_dict["Total"]["MOB"])))
    print(''.join([f'{x:{10}}' for x in headings]))
    for condition in mutation_type_dict:
        if condition == "Total":
            continue
        body = [condition, len(mutation_type_dict[condition]["DEL"]), len(mutation_type_dict[condition]["AMP"]), len(mutation_type_dict[condition]["INS"]), len(mutation_type_dict[condition]["SUB"]), len(mutation_type_dict[condition]["INV"]), len(mutation_type_dict[condition]["CON"]), len(mutation_type_dict[condition]["MOB"])]
        for mob_type in list(sorted(set(mutation_type_dict["Total"]["MOB"]))):
            body.append(mutation_type_dict[condition]["MOB"].count(mob_type))

        print(''.join([f'{x:{10}}' for x in map(str,body)]))
    body = ["Total", len(mutation_type_dict["Total"]["DEL"]), len(mutation_type_dict["Total"]["AMP"]), len(mutation_type_dict["Total"]["INS"]), len(mutation_type_dict["Total"]["SUB"]), len(mutation_type_dict["Total"]["INV"]), len(mutation_type_dict["Total"]["CON"]), len(mutation_type_dict["Total"]["MOB"])]
    for mob_type in list(sorted(set(mutation_type_dict["Total"]["MOB"]))):
        body.append(mutation_type_dict["Total"]["MOB"].count(mob_type))
    print(''.join([f'{x:{10}}' for x in map(str,body)]))

    print()
    print("Deletion Lengths")
    print(''.join([f"{x:<{10}}" for x in ["Length", "Count"]]))

    deletion_counts = {}
    for deletion in sorted(mutation_type_dict["Total"]["DEL"]):
        if deletion not in deletion_counts:
            deletion_counts[deletion] = 1
        else:
            deletion_counts[deletion] += 1

    for d in sorted(deletion_counts.keys()):
        print(''.join([f"{x:<{10}}" for x in [d,deletion_counts[d]]]))

    print()
    print("Amplication Lengths")
    print(''.join([f"{x:<{10}}" for x in ["Length", "Count"]]))

    amp_counts = {}
    for amplification in sorted(mutation_type_dict["Total"]["AMP"]):
        if amplification not in amp_counts:
            amp_counts[amplification] = 1
        else:
            amp_counts[amplification] += 1

    for a in sorted(amp_counts.keys()):
        print(''.join([f"{x:<{10}}" for x in [a,amp_counts[a]]]))


def run_fishers_exact(fisher_list):
    """2x2 Fisher's exact test from a list of: [success condition 1, fail condition 1, success condition 2, fail condition 2]"""
    assert len(fisher_list) == 4, "fisher list not of correct size"
    y = [ [fisher_list[0],fisher_list[1]],[fisher_list[2],fisher_list[3]] ]
    return fisher_exact(y,alternative='two-sided')[1]

def write_mutation_matrix(master_dict,valid_mutated_genes,outfile):
    outfh = open(outfile,"w")
    print("mutation per sample")
    outfh.write(','.join([x for x in map(str, ['Gene'] + [x for x in master_dict["Total"]])]) + '\n')
    print(''.join([f"{x:<{10}.7}" for x in map(str, ['Gene'] + [x for x in master_dict["Total"]])]))
    for cur_gene in sorted(valid_mutated_genes, key=lambda k: len(valid_mutated_genes[k]), reverse=True):
        printable = [cur_gene]
        for sample in master_dict["Total"]:
            try:
                printable.append(str(master_dict["Total"][sample]["valid_mutations"]).count(cur_gene))
            except KeyError:
                printable.append(0)
        print(''.join([f"{x:<{10}}" for x in map(str,printable)]))
        outfh.write(','.join([x for x in map(str,printable)]) + '\n')


def main():
    repeat_list_of_lists = repeat_id()
    master_dict = None  ## initialize to avoid warning based on ifs

    #Treatment read in
    if args.directory_treatments:  ## In future may have multiple ways to read treatments_and_samples in
        master_dict = treatment_read_in()
        for treatment in master_dict:
            for sample in master_dict[treatment]:
                x = mutation_gene_assignment(master_dict[treatment][sample]['raw_mutations'], repeat_list_of_lists, args.dN_only)
                for mutationlist in x:  ## assign valid_mutations, and intergenic_and_multi_gene mutations to master dict
                    master_dict[treatment][sample][mutationlist] = x[mutationlist]
                k = master_dict[treatment][sample]['valid_mutations'] + master_dict[treatment][sample]['intergenic_or_multi-gene']
                master_dict[treatment][sample]['total_mutations'] = k
    assert master_dict is not None, "master_dict has not been populated. This should never trigger given assertion that args.dictionary_treatments must be supplied"

    valid_mutated_genes = {}
    for treatment in master_dict:
        if treatment == "Total":
            continue
        for sample in master_dict[treatment]:
            for mutation in master_dict[treatment][sample]["valid_mutations"]:
                ## value is a string formatted like so: 'treatment:sample'.
                ## Examples: "Ara+:REL123" or "Ara-:REL456"
                if mutation[2] not in valid_mutated_genes:
                    valid_mutated_genes[mutation[2]] = [treatment + ':' + sample]
                else:
                    valid_mutated_genes[mutation[2]].append(treatment + ':' + sample)

    ''' write out mutation matrix for further analysis.'''
    outfile = args.matrixfile
    write_mutation_matrix(master_dict,valid_mutated_genes,outfile)


    ##  general mutation statistics; print statements within functions
    general_mutation_stats(master_dict)
    types_of_mutations(master_dict)
    mann_whitney_tests(master_dict)
    kruskal_wallis_tests(master_dict)


    dice(master_dict, [])  ## initial dice calculation done without excluding significant genes

    print("\nGene-centric mutation information:\nFisher's exact p-value testing if clustering within treatment greater than expected by chance.\nOnly p values less than %s listed." % args.pvalue)
    print("IMPORTANT NOTE: the numbers reported here are NOT used directly in the Fisher's exact test calculation. The FET is based on number of genomes hit, and not the raw number of mutations.")

    print(''.join([f'{x:{10}}' for x in map(str,["Gene"] + [x for x in master_dict if not x == "Total"] + ["Total", "p-value"])]))
    excluded_genes = []  # list of genes to exclude from repeat dice
    ##  iterate through genes, sorted in descending order based on total number of mutations
    for cur_gene in sorted(valid_mutated_genes, key=lambda k: len(valid_mutated_genes[k]), reverse=True):
        printable = [cur_gene]
        for treatment in [x for x in master_dict if not x == "Total"]:
            try:
                treatments_of_valid_mutations = [x.split(':')[0] for x in valid_mutated_genes[cur_gene]]
                printable.append(treatments_of_valid_mutations.count(treatment))
            except KeyError:
                printable.append(0)
        printable.append(len(valid_mutated_genes[cur_gene]))
        if printable[-1] > 1:
            ## score hits as 1 or 0 (don't count additional hits in a gene).
            hit_genomes = list(set(valid_mutated_genes[cur_gene]))
            treatments_in_hit_genomes = [x.split(':')[0] for x in hit_genomes]
            most_mutated_treatment, most_mutated_count = Counter(treatments_in_hit_genomes).most_common(1)[0]
            total_hit_genomes = len(hit_genomes)
            total_genome_samples = len(master_dict['Total'])
            num_genomes_in_most_mutated_treatment = len(master_dict[most_mutated_treatment])

            s_c1 = most_mutated_count
            f_c1 = num_genomes_in_most_mutated_treatment - most_mutated_count
            s_c2 = total_hit_genomes - most_mutated_count
            ## next line: subtract successes in the remaining treatments from
            ## the number of genomes in the remaining treatments.
            f_c2 = (total_genome_samples - num_genomes_in_most_mutated_treatment) - s_c2
            
            assert f_c1 >= 0, "error 1: negative numbers in f_c1!"
            assert f_c2 >= 0, "error 2: negative numbers in f_c2!"

            '''
            I believe the idea here is the following contingency table:
                                                       DM0   | DM25   ( or whatever treatment is relevant)
            Number of genomes w/ mutation in gene x  | s_c1  | f_c1
            Number of genomes w/o mutation in gene x | s_c2  | f_c2
            '''

            if run_fishers_exact([s_c1, f_c1, s_c2, f_c2]) <= args.pvalue:
                printable.append(run_fishers_exact([s_c1, f_c1, s_c2, f_c2]))
                excluded_genes.append(cur_gene)

            if printable[-1] == 2:
                assert printable[1:-1].count(1) == 2 or printable[1:-1].count(2) == 1, "unclear state"

        print(''.join([f"{x:{10}}" for x in map(str,printable)]))
        ##print("FISHER\'S EXACT TEST COUNTS BASED ON HIT GENOMES: s1, f1, s2, f2: ",s_c1,f_c1,s_c2,f_c2)

    print("Invalid AMP or DEL")
    print(''.join([f"{x:<{23}}" for x in ["Treatment", "Sample", "AMP or DEL", "Size"]]))
    for sample in master_dict["Total"]:
        for mutation in master_dict["Total"][sample]["intergenic_or_multi-gene"]:
            if mutation[0] in ["AMP", "DEL"]:
                sample_treatment = "UNK"
                for treatment in [x for x in master_dict if not x == "Total"]:
                    if sample in master_dict[treatment]:
                        assert sample_treatment is "UNK", "Identical sample name found in multiple treatments. %s" % sample
                        sample_treatment = treatment
                assert sample_treatment is not "UNK", "Sample not identified as belonging to a treatment. %s" % sample
                print(''.join([f"{x:<{23}}" for x in map(str, [sample_treatment, sample, mutation[0], mutation[3]])]))

    print("Valid AMP or DEL")
    print(''.join([f"{x:<{23}}" for x in ["Treatment", "Sample", "AMP or DEL", "Size", "Details"]]))
    for sample in master_dict["Total"]:
        for mutation in master_dict["Total"][sample]["valid_mutations"]:
            if mutation[0] in ["AMP", "DEL"]:
                sample_treatment = "UNK"
                for treatment in [x for x in master_dict if not x == "Total"]:
                    if sample in master_dict[treatment]:
                        assert sample_treatment is "UNK", "Identical sample name found in multiple treatments. %s" % sample
                        sample_treatment = treatment
                assert sample_treatment is not "UNK", "Sample not identified as belonging to a treatment. %s" % sample
                details = ','.join(map(str,mutation))
                print(''.join([f"{x:<{23}}" for x in [sample_treatment, sample, mutation[0], mutation[3], details]]))
    for sample in master_dict["Total"]:
        for mutation in master_dict["Total"][sample]["valid_mutations"]:
            if mutation[2] in ["insJ-5", "insL-2"]:
                print(''.join([f"{x:<{23}}" for x in map(str, [sample] + mutation)]))

main()
