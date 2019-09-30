
rule dice_analysis:
     shell:
	'python dice-analysis.py'

rule prep_growth_data:
     shell:
	'python prep-growth-data.py'

#rule kallisto_index:
#    input:
#        '../results/RNAseq-analysis/LCA_ORFs.fasta'
#    output:
#        '../results/RNAseq-analysis/LCA_ORFs.idx'
#    shell:
#        'kallisto index --make-unique -i {output} {input}'

#rule kallisto_quant:
#    input:
#        index='../results/RNAseq-analysis/LCA_ORFs.idx'
#    output:
#        ''
#    shell:
#        'kallisto quant -i {input.index} -o {output} -b 100 input.readsA1 input.readsA2 input.readsB1 input.readsB2'
