##copy-number-analysis.R by Rohan Maddamsetti.

## 1) use xml2 to get negative binomial fit and dispersion from
## breseq output summary.html. This is H0 null distribution of 1x coverage.

## 2) Find intervals longer than 500bp that reject H0 coverage in genome.
##    at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.

## 3) Do a more rigorous test for each region. take positions in the region separated by more than 500bp,
## and determine the probability that all are independently significant under the null, compared to
## a corrected bonferroni. The 500bp ensures positions cannot be spanned by a single Illumina read.

## 4) Estimate copy number by dividing mean coverage in each region by the mean
##   of the H0 1x coverage distribution.

## 5) return copy number and boundaries for each significant amplification.

## 6) find the genes contained in these regions (and genes chopped at the boundaries).

## 7) do a quick check to see whether polymorphism occurs in the genomes
##    with the largest amplifications (this should be clear from Fig. 3).

## 8) make Fig. 3, showing maeA and citT copy number.

library(xml2)
library(roxygen2)
library(assertthat)
library(IRanges)
library(GenomicRanges)
library(genbankr)

library(data.table)  # faster fread()
library(dplyr)       # consistent data.frame operations
library(purrr)       # consistent & safe list/vector munging
library(tidyr)       # consistent data.frame cleaning
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(DT)          # prettier data.frame output
library(dtplyr)      # dplyr works with data.table now.


#' parse the summary.html breseq output file, and return the mean and dispersion
#' of the negative binomial fit to the read coverage distribution, returned as a
#' data.frame with columns {mean, dispersion}.
#' NOTE: this code has only been tested on the summary file
#' output by breseq 0.30.0. It might fail on earlier or later versions.

#' @export
coverage.nbinom.from.html <- function (breseq.output.dir) {
    summary.html.f <- file.path(breseq.output.dir,"output/summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Reference Sequence Information.
    query <- '//table[./tr/th[contains(text(),"fit dispersion")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    avg <- as.numeric(xml_text(table.data[5]))
    dispersion <- as.numeric(xml_text(table.data[6]))
    print(paste('mean coverage is:',avg))
    print(paste('dispersion is:',dispersion))
    return(data.frame('mean'=avg,'dispersion'=dispersion,'variance'=avg*dispersion))
}


#' Find intervals longer than 500bp that reject H0 coverage in genome.
#' at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.
#' Then do a more rigorous test for each region. take positions in the region separated by more than 500bp,
#' and determine the probability that all are independently significant under the null, compared to
#' a corrected bonferroni. The 500bp ensures positions cannot be spanned by a single Illumina read.
#' Estimate copy number by dividing mean coverage in each region by the mean of the H0 1x coverage distribution.
#' return mean copy number, and boundaries for each region that passes the amplification test.
#' @export
find.amplifications <- function(breseq.output.dir,gnome) { #gnome is not a misspelling.

    ## Use xml2 to get negative binomial fit and dispersion from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir)

    max.read.len <- 500
    genome.length <- 4629812
    alpha <- 0.05
    uncorrected.threshold <- qnbinom(p=alpha,mu=nbinom.fit$mean,size=nbinom.fit$dispersion,lower.tail=FALSE)

    genome.coverage.file <- file.path(breseq.output.dir,"08_mutation_identification/REL606.coverage.tab")

    ## use dtplyr for speed!
    genome.coverage <- tbl_dt(fread(genome.coverage.file)) %>%
        select(position,unique_top_cov,unique_bot_cov) %>% mutate(coverage=unique_top_cov+unique_bot_cov)


    ## find candidate amplifications that pass the uncorrected threshold.
    candidate.amplifications <- filter(genome.coverage,coverage > uncorrected.threshold)
    ## calculate intervals of candidate amplifications.
    boundaries <- candidate.amplifications %>% mutate(left.diff=position - lag(position)) %>%
        mutate(right.diff=lead(position) - position) %>%
        ## corner case: check for the NA values at the endpoints and set them as boundaries.
        mutate(is.right.boundary=is.na(right.diff)|ifelse(right.diff>1,TRUE,FALSE)) %>%
        mutate(is.left.boundary=is.na(left.diff)|ifelse(left.diff>1,TRUE,FALSE)) %>%
        filter(is.left.boundary==TRUE | is.right.boundary==TRUE)

    left.boundaries <- filter(boundaries,is.left.boundary==TRUE) %>% arrange(position)
    right.boundaries <- filter(boundaries,is.right.boundary==TRUE) %>% arrange(position)
    assert_that(nrow(left.boundaries) == nrow(right.boundaries))

    ## helper higher-order function to get min, max, mean coverage of each segment.
    get.segment.coverage <- function(left.bound,right.bound,coverage.table,funcx) {
        seg <- coverage.table %>% filter(position>left.bound) %>% filter(position<right.bound)
        return(funcx(seg$coverage))
    }

    amplified.segments <- data.frame(left.boundary=left.boundaries$position,right.boundary=right.boundaries$position) %>%
        ## filter out intervals less than 500 bp.
        mutate(len=right.boundary-left.boundary) %>% filter(len>max.read.len) %>% mutate(amplication.index=row_number()) %>%
        ## find min, max, and mean coverage of each amplified segment.
        group_by(left.boundary,right.boundary) %>%
        summarise(coverage.min=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,min),
                  coverage.max=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,max),
                  coverage.mean=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,mean)) %>%
        mutate(len=right.boundary-left.boundary) %>%
        mutate(copy.number.min=coverage.min/nbinom.fit$mean,copy.number.max=coverage.max/nbinom.fit$mean,
               copy.number.mean=coverage.mean/nbinom.fit$mean)

    ##print(data.frame(amplified.segments))
    ## divide alpha by the number of tests for the bonferroni correction.
    bonferroni.alpha <- alpha/(genome.length + sum(amplified.segments$len))
    corrected.threshold <- qnbinom(p=bonferroni.alpha,mu=nbinom.fit$mean,size=nbinom.fit$dispersion,lower.tail=FALSE)

    ## This is my test: take the probability of the minimum coverage under H0 to the power of the number of
    ## uncorrelated sites in the amplification (sites more than 500 bp apart). Then see if this is smaller than the
    ## bonferroni corrected p-value for significance..
    significant.amplifications <- amplified.segments %>%
        mutate(pval=(pnbinom(q=coverage.min,
                             mu=nbinom.fit$mean,
                             size=nbinom.fit$dispersion,
                             lower.tail=FALSE))^(len%/%max.read.len)) %>%
        mutate(is.significant=ifelse(pval<bonferroni.alpha,TRUE,FALSE)) %>%
        filter(is.significant==TRUE) %>% mutate(Genome=gnome) %>%
        mutate(bonferroni.corrected.pval=pval*alpha/bonferroni.alpha)

    ##print(data.frame(significant.amplifications))

    return(significant.amplifications)
}

## input: REL606.gbk: file.path of the reference genome,
##    amplifications: data.frame returned by find.amplifications.
annotate.amplifications <- function(amplifications,REL606.gbk) {

    ## create the IRanges object.
    amp.ranges <- IRanges(amplifications$left.boundary,amplifications$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with REL606 genes.
    g.amp.ranges <- GRanges("REL606",ranges=amp.ranges)
    ## and add the data.frame of amplifications as metadata.
    mcols(g.amp.ranges) <- amplifications

    ## find the genes within the amplifications.
    ## The parsing is sloooooow.
    pREL606 <- parseGenBank(REL606.gbk)
    gb <- make_gbrecord(pREL606)
    REL606.genes <- genes(gb)
    ##print(REL606.genes)
    ## find overlaps between annotated genes and amplifications.
    hits <- findOverlaps(REL606.genes,g.amp.ranges,ignore.strand=FALSE)
    ## take the hits, the REL606 annotation, and the amplifications,
    ## and produce a table of genes found in each amplication.

    hits.df <- data.frame(query.index=queryHits(hits),subject.index=subjectHits(hits))

    query.df <- data.frame(query.index=1:length(REL606.genes),
                           gene=REL606.genes$gene,locus_tag=REL606.genes$locus_tag,
                           start=start(ranges(REL606.genes)),end=end(ranges(REL606.genes)))

    subject.df <- bind_cols(data.frame(subject.index=1:length(g.amp.ranges)),data.frame(mcols(g.amp.ranges)))

    amplified.genes.df <- left_join(hits.df,query.df) %>% left_join(subject.df) %>%
        ## if gene is NA, replace with locus_tag. have to change factors to strings!
        mutate(gene = ifelse(is.na(gene),as.character(locus_tag),as.character(gene)))

    return(amplified.genes.df)
}


## see the excellent plots on https://rud.is/projects/facetedheatmaps.html
plot.Fig4A.heatmap <- function(annotated.amps,clone.labels) {

    labeled.annotated.amps <- left_join(annotated.amps,clone.labels,by=c("Genome" = 'Name')) %>%
        select(-query.index,-subject.index,-is.significant,-Sample.Type) %>%
        ## replace 'sfcA' with 'maeA' in the plot.
        mutate(gene = replace(gene, gene == 'sfcA', 'maeA'))
        ## reverse the order genes by start to get axes correct on heatmap.
    labeled.annotated.amps$gene <- with(labeled.annotated.amps, reorder(gene, rev(start)))

    heatmap <- ggplot(labeled.annotated.amps,aes(x=Genome,y=gene,fill=left.boundary,frame=Environment)) +
        geom_tile(color="white",size=0.1) +
        ## draw vertical lines at genes of interest.
        geom_hline(linetype='dashed',yintercept = which(levels(labeled.annotated.amps$gene) %in% c('citT','dctA','maeA'))) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        scale_fill_viridis(name="") +
        coord_equal() +
        facet_wrap(~Environment,ncol=2, scales = "free_x") +
        theme_tufte(base_family='Helvetica') +
        theme(axis.ticks=element_blank()) +
        theme(axis.text.x=element_text(size=5,angle=45,hjust=1)) +
        theme(axis.text.y=element_text(size=5,
                                       ## bold the text for citT, dctA, and maeA
                                       face=ifelse(levels(labeled.annotated.amps$gene) %in% c('citT','dctA','maeA'),"bold.italic","italic"))) +
        theme(panel.border=element_blank()) +
        theme(plot.title=element_text(hjust=0)) +
        theme(strip.text=element_text(hjust=0,size=12)) +
        theme(panel.spacing.x=unit(1, "cm")) +
        theme(panel.spacing.y=unit(0.5, "cm")) +
        theme(legend.title=element_text(size=6)) +
        theme(legend.title.align=1) +
        theme(legend.text=element_text(size=6)) +
        theme(legend.position="bottom") +
        theme(legend.key.size=unit(0.2, "cm")) +
        theme(legend.key.width=unit(1, "cm"))
    return(heatmap)
}

## plot a stacked bar graph to show how big the amplifications are
## in the DM0 and DM25 treatments.
plot.Fig4B.stackedbar <- function(amps,clone.labels) {

amps2 <- amps %>% mutate(total.amp.length = copy.number.mean*len) %>%
    left_join(clone.labels, by=c("Genome" = "Name"))

genome.length <- 4629812

stacked <- ggplot(amps2,aes(x=Genome,y=total.amp.length,fill=left.boundary)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=genome.length/5,linetype="dashed") +
    scale_fill_viridis(name="") +
    coord_equal() +
    facet_wrap(~Environment,ncol=2, scales = "free_x") +
    theme_tufte(base_family='Helvetica') +
    theme(axis.ticks=element_blank()) +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1)) +
        theme(axis.text.y=element_text(size=12)) +
        theme(panel.border=element_blank()) +
        theme(plot.title=element_text(hjust=0)) +
        theme(strip.text=element_text(hjust=0,size=12)) +
        theme(panel.spacing.x=unit(1, "cm")) +
        theme(panel.spacing.y=unit(0.5, "cm")) +
        theme(legend.title=element_text(size=6)) +
        theme(legend.title.align=1) +
        theme(legend.text=element_text(size=6)) +
        theme(legend.position="bottom") +
        theme(legend.key.size=unit(0.2, "cm")) +
        theme(legend.key.width=unit(1, "cm"))

    return(stacked)
}

projdir <- "/Users/Rohandinho/Dropbox (HMS)/DM0-evolution"
outdir <- file.path(projdir,"results")
breseq.output.dir <- file.path(projdir,"genomes/polymorphism")
REL606.gb <- file.path(projdir,"genomes/REL606.7.gbk")
all.genomes <- list.files(breseq.output.dir,pattern='^ZDBp')
all.genome.paths <- sapply(all.genomes, function(x) file.path(breseq.output.dir,x))
genome.input.df <- data.frame(Genome=all.genomes,path=all.genome.paths)

amps <- map2_df(genome.input.df$path,
                genome.input.df$Genome,
                find.amplifications)

write.csv(x=amps,file=file.path(outdir,"amplifications.csv"))

annotated.amps <- amps %>% annotate.amplifications(REL606.gb)
write.csv(x=annotated.amps,file=file.path(outdir,"amplified_genes.csv"))

## check ZDBp895.
## test <- find.amplifications(file.path(breseq.output.dir,"ZDBp895"),"ZDBp895")
## this just fails statistical significance, probably because less data for
## ZDBp895 than the other genomes == less statistical power.

amp.parallelism <- annotated.amps %>% group_by(gene,locus_tag) %>% summarise(count=n()) %>% arrange(desc(count),locus_tag)


## report copy number of maeA, citT, and dctA in a table.
## (I can do this easily by hand).

copy.number.table <- annotated.amps %>%
    filter(gene %in% c('citT','sfcA','dctA')) %>%
    ## change sfcA to maeA
    mutate(Gene=replace(gene,gene=='sfcA','maeA')) %>%
    transform(amplified.segment.length=len) %>%
    arrange(Gene,Genome) %>%
    select(Gene,Genome,copy.number.mean,copy.number.min,copy.number.max, bonferroni.corrected.pval)

write.csv(x=copy.number.table,file=file.path(outdir,"copy_number_table.csv"))

## dctA and maeA amplifications don't occur together. The probability
## of this occurrence is (20/24)*(19/23)*...*(15/19) = 28.7% by chance.

## Make figures.

label.filename <- file.path(projdir,"data/rohan-formatted/populations-and-clones.csv")
clone.labels <- read.csv(label.filename)

## Make a heatmap plot with facet grid on DM0 versus DM25.
## color the matrix based on gene location.

heatmap2 <- plot.Fig4A.heatmap(annotated.amps,clone.labels)
ggsave(heatmap2,filename=file.path(outdir,"figures/Fig4A.pdf"),height=11,width=3)

## Make a stacked bar plot, with facet grid on DM0 versus DM25.
## color the stacked bars based on the left boundary
## of the amplification.

stacked2 <- plot.Fig4B.stackedbar(amps,clone.labels)
ggsave(stacked2,filename=file.path(outdir,"figures/Fig4B.pdf"))
