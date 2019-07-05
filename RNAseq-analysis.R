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

make.transcript.table <- function(s) {
  
  data <- read.table(
    file.path(proj.dir, "results", "RNAseq-analysis", s, "abundance.tsv"),
    header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
      select(target_id, tpm) %>%
          mutate(temp = target_id) %>%
            extract(temp,
                    c("locus_tag", "gene", "note"),
                    regex=".+=(.+)\\|.+=(.+)\\|.+=(.+)") %>%
                      rename(!! s := tpm)
  return(data)
}

## make a df of inferred transcript abundances for each sample.
transcript.df <- map(metadata$sample, make.transcript.table) %>%
  reduce(left_join) %>%
  gather(key='Sample',
         value='tpm',
         -target_id, -locus_tag, -gene, -note) %>%
  separate(Sample, c('Clone','Run'))

total.run.transcripts.df <- transcript.df %>%
  group_by(Clone, Run) %>%
  summarize(total.run.transcripts=sum(tpm,na.rm=TRUE))

total.clone.transcripts.df <- total.run.transcripts.df %>%
  summarize(total.transcripts=sum(total.run.transcripts))

total.run.transcripts.df %<>% left_join(total.clone.transcripts.df) %>%
  mutate(weights = total.run.transcripts/total.transcripts)

## weigh the mean by the number of transcripts in each transcriptome.
mean.transcript.df <- transcript.df %>%
  left_join(total.run.transcripts.df) %>%
  group_by(target_id, locus_tag, gene, note, Clone) %>%
  summarize(mean_tpm=weighted.mean(tpm,weights,na.rm=TRUE)) %>%
  ungroup()

## get gene annotation from the kallisto output (target_id already has annotation).
my.annotation <- transcript.df %>%
  select(target_id,locus_tag, gene, note) %>%
  distinct()

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

sleuth_live(so)

## NOTE: I am uncertain whether the following results are
## correct or not. Am I calculating average counts properly??
## check out apparent variability in fadB..
## seems like my calculations are way too naive.
## Hold off on reporting for now, and stick to the sleuth results.

## make a heatmap with clustergrammer.
## filter on RNAs differently expressed between
## ancestors and evolved clones.
clustergrammer.df <- mean.transcript.df %>%
  filter(target_id %in% results.table$target_id) %>%
  mutate(locus = str_c(locus_tag,gene,sep='-')) %>%
  select(-target_id,-locus_tag,-gene,-note) %>%
  drop_na() %>%
  spread(Clone, mean_count) %>%
  column_to_rownames('locus')

## log ratio of evolved with CZB151.
clustergrammer.df2 <- mean.transcript.df %>%
  filter(target_id %in% results.table$target_id) %>%
  mutate(locus = str_c(locus_tag,gene,sep='-')) %>%
  select(-target_id,-locus_tag,-gene,-note) %>%
  spread(Clone, mean_count) %>%
  mutate(ZDBp889 = log2(ZDBp889/CZB151)) %>%
  mutate(ZDBp883 = log2(ZDBp883/CZB151)) %>%
  mutate(ZDBp877 = log2(ZDBp877/CZB151)) %>%
  select(-CZB151, -CZB152) %>%
  filter(is.finite(ZDBp889)) %>%
  filter(is.finite(ZDBp883)) %>%
  filter(is.finite(ZDBp877)) %>%
  drop_na() %>%
  column_to_rownames('locus')
  
## NOTE: add a tab to first row by hand, then save as "tab_clustergrammer.tsv"
## so that not overwritten.
write.table(clustergrammer.df,
            file.path(proj.dir, "results", "sleuth-results", "clustergrammer.tsv"),
            sep='\t', quote= FALSE)

write.table(clustergrammer.df2,
            file.path(proj.dir, "results", "sleuth-results", "clustergrammer2.tsv"),
            sep='\t', quote= FALSE)

clustergrammer.df3 <- mean.transcript.df %>%
  filter(target_id %in% results.table$target_id) %>%
  mutate(locus = str_c(locus_tag,gene,sep='-')) %>%
  select(-target_id,-locus_tag,-gene,-note) %>%
  spread(Clone, mean_count) %>%
  mutate(ZDBp889 = log2(ZDBp889/CZB151)) %>%
  mutate(ZDBp883 = log2(ZDBp883/CZB151)) %>%
  mutate(ZDBp877 = log2(ZDBp877/CZB151)) %>%
  select(-CZB151, -CZB152) %>%
  filter(is.finite(ZDBp889)) %>%
  filter(is.finite(ZDBp883)) %>%
  filter(is.finite(ZDBp877)) %>%
  drop_na() %>%

  filter( ((ZDBp877 > 0) & ((ZDBp883 < 0) | (ZDBp889 < 0))) |
         ((ZDBp877 < 0) & ((ZDBp883 > 0) | (ZDBp889 > 0)))
         )

## ycfR anticorrelates with a lot of genes: rpoB, rps, rpl genes.
## this is a stress response gene that affects outer membrane.
## ycdQ anticorrelates with aceB, pheA, murF, lpdA.
