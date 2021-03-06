## RNAseq-analysis.R by Rohan Maddamsetti.

## Use sleuth to analyze kallisto quantification of the Kenyon RNAseq data.
## See walkthroughs at:
## https://pachterlab.github.io/sleuth_walkthroughs/pval_agg/analysis.html
## http://achri.blogspot.com/2016/10/using-kallisto-sleuth.html

library(tidyverse)
library(cowplot)
library(sleuth)

home.dir <- path.expand("~")
proj.dir <- file.path(home.dir,"Desktop/Research/active-projects/DM0-evolution")

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
                  c("locus_tag", "gene"),
                  regex="^locus_tag=(.+)\\|gene=(.+)$") %>%
                      rename(!! s := tpm)
  return(data)
}

## make a df of inferred transcript abundances for each sample.
transcript.df <- map(metadata$sample, make.transcript.table) %>%
  reduce(left_join) %>%
  gather(key='Sample',
         value='est_counts',
         -target_id, -locus_tag, -gene)

## get gene annotation from the kallisto output (target_id already has annotation).
my.annotation <- transcript.df %>%
  select(target_id,locus_tag, gene) %>%
  distinct()

## for debugging, write my.annotation to file.
##write.csv(transcript.df,file='/Users/Rohandinho/Desktop/transcript_df.csv')
##write.csv(my.annotation,file="/Users/Rohandinho/Desktop/annotation.csv")

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

## write the results to file.
write.csv(lrt.results,file=file.path(proj.dir,"results/sleuth-results/lrt-results.csv"))
write.csv(wt.results,file=file.path(proj.dir,"results/sleuth-results/wald-results.csv"))

## hack the internals of this function to plot heatmap.
## sleuth::plot_transcript_heatmap
## see documentation at:
## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/topics/plot_transcript_heatmap
## see general documentation at:
## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/

all.genes <- sort(my.annotation$gene)
## vector of genes that are discussed in the text.
discussed.genes <- c("fadA", "fadB", "fadD", "fadE", "fadH", "fadR",
                     "cyoA", "cyoB", "cyoC", "cyoD", "cyoE",
                     "sfcA", "dctA", "yhiO",
                     "rpsB", "rpsU", "rpsV", "rpsO", "rpsT",
                     "rplE", "rplJ", "rplN", "rplX",
                     "rpoA", "rpoB", "rpoC", "rpoS",
                     "gyrA", "rho", "sucB", "sucC",
                     "ECB_00826", "ECB_00827",
                     "chpA", "chpR", "recA", "recN",
                     "dinD", "sulA", "umuC", "umuD",
                     "csiE", "sbmC",
                     "glpA", "glpB", "glpC", "glpD", "glpT",
                     ## add genes which are anti-correlated with
                     ## with maeA.
                     "gltS", "dinI", "ECB_03510"
                     )

results.matrix <- sleuth_to_matrix(so, "obs_norm", "est_counts")
## change the row names of this matrix for better plotting.
## first turn into a dataframe, 
results.matrix.2 <- data.frame(results.matrix) %>%
    rownames_to_column(var = 'target_id') %>%
    left_join(my.annotation) %>%
    ## then filter on genes that are discussed in the text,
    filter(gene %in% discussed.genes) %>%
    ## rename sfcA to maeA.
    mutate(gene = replace(gene, gene=="sfcA", "maeA")) %>%
    ## and convert back to plot the heatmap.
    select(-target_id, -locus_tag) %>%
    column_to_rownames(var = 'gene')

## Here I'm hacking some code from the internals of the plot_transcript_heatmap
## function, so that I can use a dataframe rather than a sleuth object to make
## the figure.

## transform counts using: log2(1+est_counts)
trans_mat <- as.matrix(log2(results.matrix.2 + 1))

## default colors
x_axis_angle = 50
color_high <- "#581845"
color_mid <- "#FFC300"
color_low <- "#DAF7A6"
colors <- colorRampPalette(c(color_low, color_mid, color_high))(100)

## make the plot! This is Figure 14.
RNAseq.fig.outf <- file.path(proj.dir,
              "results/figures/Fig14.pdf")

## We want to italicize gene names in the rows.
## See how at: https://www.biostars.org/p/400381/
italicsnames <- lapply(
    rownames(trans_mat),
    function(x) bquote(italic(.(x))))

p <- pheatmap::pheatmap(trans_mat,
                        ##color=colors, ## comment this out if we want red/blue gradient.
                        cluster_cols = TRUE, cluster_rows=TRUE,
                        labels_row = as.expression(italicsnames))

p$gtable$grobs[[3]]$rot <- 360 - x_axis_angle
gridExtra::grid.arrange(p$gtable)

ggsave(RNAseq.fig.outf, p, height=10, width=7)


## run the RShiny web app interface to the sleuth results.
##sleuth_live(so)


## let's do a differential expression analysis to see if anything is anti-correlated
## with sfcA/maeA expression, which is amplified in ZDBp883 and ZDBp889, but not 877.

maeA.metadata <- metadata %>%
filter(Treatment=='Evolved') %>%
select(-Treatment) %>%
mutate(NoAMP = ifelse(Clone %in% c('ZDBp883', 'ZDBp889'),'No','Yes'))

maeA.so <- sleuth_prep(maeA.metadata, target_mapping = my.annotation,
                  extra_bootstrap_summary = TRUE) %>%
  sleuth_fit(~ NoAMP, 'full') %>%
  sleuth_fit(~ 1, 'reduced') %>%
  sleuth_lrt('reduced','full') %>%
  sleuth_wt("NoAMPYes", "full")

maeA.lrt.results <- sleuth_results(maeA.so, 'reduced:full', 'lrt', show_all = TRUE) %>%
    filter(qval < 0.05) %>% distinct()
## Wald test will return an effect size (regression parameter beta).
maeA.wt.results <- sleuth_results(maeA.so,'NoAMPYes') %>%
    filter(qval < 0.05) %>% distinct()

## write maeA AMP lrt results to file.
write.csv(maeA.wt.results,file=file.path("..", "results", "sleuth-results", "maeA-AMP-lrt.csv"))

## Check out the results using the Shiny interface.
##sleuth_live(maeA.so)
