## RNAseq-analysis.R by Rohan Maddamsetti.

## Use sleuth to analyze kallisto quantification of the Kenyon RNAseq data.
## See walkthroughs at:
## https://pachterlab.github.io/sleuth_walkthroughs/pval_agg/analysis.html
## http://achri.blogspot.com/2016/10/using-kallisto-sleuth.html

library(tidyverse)
library(cowplot)
library(sleuth)

home.dir <- path.expand("~")
proj.dir <- file.path(home.dir,"BoxSync/active-projects/DM0-evolution")

pop.clone.labels <- read.csv(
    file.path(proj.dir,
              "data/rohan-formatted/populations-and-clones.csv"),
    stringsAsFactors=FALSE)

metadata <- read.csv(
    file.path(proj.dir,
              "data/rohan-formatted/Kenyon-RNAseq-metadata.csv"),
    header=TRUE,
    stringsAsFactors = FALSE)

## add path names of kallisto output directories to the metadata table.
metadata <- metadata %>%
    mutate(path = file.path(proj.dir,"results","RNAseq-analysis",sample,"abundance.h5"))

## CRITICAL BUG TO FIX: HORRIBLE PARSING PROBLEMS WITH STRINGS IN
## transcript.df. TOP PRIORITY TO FIX!!!

make.transcript.table <- function(s) {
  
  data <- read.table(
    file.path(proj.dir, "results", "RNAseq-analysis", s, "abundance.tsv"),
    header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
      select(target_id, tpm) %>%
          mutate(temp = target_id) %>%
          extract(temp,
                  c("locus_tag", "gene"),
                    ##c("locus_tag", "gene", "note"),
                    ##regex=".+=(.+)\\|.+=(.+)\\|.+=(.+)") %>%
                  ##regex="locus_tag=(.+)\\|gene=(.+)\\|Note=(.+)") %>%
                  regex="^locus_tag=(.+)\\|gene=(.+)\\|") %>%
                      rename(!! s := tpm)
  return(data)
}

## make a df of inferred transcript abundances for each sample.
transcript.df <- map(metadata$sample, make.transcript.table) %>%
  reduce(left_join) %>%
  gather(key='Sample',
         value='tpm',
         -target_id, -locus_tag, -gene) ##%>%
         ##-target_id, -locus_tag, -gene, -note) %>%
  ##separate(Sample, c('Clone','Run'))

## get gene annotation from the kallisto output (target_id already has annotation).
my.annotation <- transcript.df %>%
  select(target_id,locus_tag, gene, note) %>%
  distinct()

## for debugging, write my.annotation to file.
write.csv(transcript.df,file='/Users/Rohandinho/Desktop/transcript_df.csv')
write.csv(my.annotation,file="/Users/Rohandinho/Desktop/annotation.csv")

so <- sleuth_prep(metadata, target_mapping = my.annotation,
                  extra_bootstrap_summary = TRUE) %>%
  sleuth_fit(~ Treatment, 'full') %>%
  sleuth_fit(~ 1, 'reduced') %>%
  sleuth_lrt('reduced','full') %>%
  sleuth_wt("TreatmentEvolved", "full")

lrt.results <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE) %>%
    filter(qval < 0.05) %>% distinct()
## Wald test will return an effect size (regression parameter beta).
wt.results <- sleuth_results(so,'TreatmentEvolved') %>%
    filter(qval < 0.05) %>% distinct()

##sleuth_live(so)

## use this function to plot heatmap.
## plot_transcript_heatmap
## see documentation at:
## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/topics/plot_transcript_heatmap

## see general documentation at:
## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/

## vector of genes that are discussed in the text.

all.genes <- sort(my.annotation$gene)

##discussed.genes <- c('fadA', 'fadB', 'fadH', 'fadR',
##                     'cyoE', 'cyoD','cyoC', 'cyoB', 'cyoA',
##                     'yhiO',
##                     'rpsB')

results.matrix <- sleuth_to_matrix(so, "obs_norm", "tpm")
## change the row names of this matrix for better plotting.
## first turn into a dataframe, 
results.matrix.2 <- data.frame(results.matrix) %>%
rownames_to_column(var = 'target_id') %>%
left_join(my.annotation) %>%
## then filter on genes that are discussed in the text,
####filter(gene %in% discussed.genes) %>%
## and convert back into matrix to plot a heatmap.
####select(-target_id, -locus_tag, -note) %>%
####column_to_rownames(var = 'gene')
