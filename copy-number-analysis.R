
##copy-number-analysis.R by Rohan Maddamsetti.

## 1) use xml2 to get negative binomial fit and dispersion from
## breseq output summary.html. This is H0 null distribution of 1x coverage.

## 2) Find intervals longer than max.read.len that reject H0 coverage in genome.
##    at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.

## 3) Do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
## and determine the probability that all are independently significant under the null, compared to
## a corrected bonferroni. The max.read.len ensures positions cannot be spanned by a single Illumina read.

## 4) Estimate copy number by dividing mean coverage in each region by the mean
##   of the H0 1x coverage distribution.

## 5) return copy number and boundaries for each significant amplification.

## 6) find the genes contained in these regions (and genes chopped at the boundaries).

## 7) do a quick check to see whether polymorphism occurs in the genomes
##    with the largest amplifications (this should be clear from Fig. 5).

## 8) make Fig. 5, showing maeA and citT copy number.

library(xml2)
library(roxygen2)
library(assertthat)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

library(tidyverse)
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(DT)          # prettier data.frame output
library(data.table)  # faster fread()
library(dtplyr)      # dplyr works with data.table now.
library(cowplot)     # layout figures nicely.

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

#' get the maximum length of a sequencing read from the summary.html breseq
#' output file.
#' @export
max.readlen.from.html <- function (breseq.output.dir) {
    summary.html.f <- file.path(breseq.output.dir,"output/summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Read File Information.
    query <- '//table[./tr/th[contains(text(),"longest")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    readlen.index <- length(table.data) - 1
    max.readlen <- xml_integer(xml_find_all(table.data[readlen.index],".//b//text()"))
    print(paste('max read length is:',max.readlen))
    return(max.readlen)
}


#' Find intervals longer than max.read.len that reject H0 coverage in genome.
#' at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.
#' Then do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
#' and determine the probability that all are independently significant under the null, compared to
#' a corrected bonferroni. max.read.len ensures positions cannot be spanned by a single Illumina read.
#' Estimate copy number by dividing mean coverage in each region by the mean of the H0 1x coverage distribution.
#' return mean copy number, and boundaries for each region that passes the amplification test.
#' @export
find.amplifications <- function(breseq.output.dir,gnome) { #gnome is not a misspelling.

    gnome <- as.character(gnome)
    print(gnome)
    ## Use xml2 to get negative binomial fit and dispersion from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir)

    ## Use xml2 to get max read length from summary.html.
    max.read.len <- max.readlen.from.html(breseq.output.dir)
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
        ## filter out intervals less than max.read.len*2 (400bp).
        mutate(len=right.boundary-left.boundary) %>% filter(len>(2*max.read.len)) %>% mutate(amplication.index=row_number()) %>%
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
    ## uncorrelated sites in the amplification (sites more than max.read.len apart). Then see if this is smaller than the
    ## bonferroni corrected p-value for significance..
    significant.amplifications <- amplified.segments %>%
        mutate(pval=(pnbinom(q=coverage.min,
                             mu=nbinom.fit$mean,
                             size=nbinom.fit$dispersion,
                             lower.tail=FALSE))^(len%/%max.read.len)) %>%
        mutate(is.significant=ifelse(pval<bonferroni.alpha,TRUE,FALSE)) %>%
        filter(is.significant==TRUE) %>% mutate(Genome=as.character(gnome)) %>%
        mutate(bonferroni.corrected.pval=pval*alpha/bonferroni.alpha)

    ##print(data.frame(significant.amplifications))

    return(significant.amplifications)
}

## input: LCA.gbk: file.path of the reference genome,
##    amplifications: data.frame returned by find.amplifications.
annotate.amplifications <- function(amplifications,LCA.gff) {

    ## create the IRanges object.
    amp.ranges <- IRanges(amplifications$left.boundary,amplifications$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with REL606 genes.
    g.amp.ranges <- GRanges("REL606",ranges=amp.ranges)
    ## and add the data.frame of amplifications as metadata.
    mcols(g.amp.ranges) <- amplifications

    ## find the genes within the amplifications.
    pLCA <- import.gff(LCA.gff)
    LCA.Granges <- as(pLCA, "GRanges")

    LCA.genes <- LCA.Granges[LCA.Granges$type == 'gene']
    ## find overlaps between annotated genes and amplifications.
    hits <- findOverlaps(LCA.genes,g.amp.ranges,ignore.strand=FALSE)

    ## take the hits, the LCA annotation, and the amplifications,
    ## and produce a table of genes found in each amplication.

    hits.df <- data.frame(query.index=queryHits(hits),subject.index=subjectHits(hits))

    query.df <- data.frame(query.index=1:length(LCA.genes),
                           gene=LCA.genes$Name,locus_tag=LCA.genes$ID,
                           start=start(ranges(LCA.genes)),end=end(ranges(LCA.genes)))

    subject.df <- bind_cols(data.frame(subject.index=1:length(g.amp.ranges)),data.frame(mcols(g.amp.ranges)))

    amplified.genes.df <- left_join(hits.df,query.df) %>% left_join(subject.df) %>%
        ## if gene is NA, replace with locus_tag. have to change factors to strings!
        mutate(gene = ifelse(is.na(gene),as.character(locus_tag),as.character(gene)))

    return(amplified.genes.df)
}

## see the excellent plots on https://rud.is/projects/facetedheatmaps.html
plot.Fig5A.heatmap <- function(annotated.amps,clone.labels) {

    ## for annotated.amps and clone.labels to play nicely with each other.
    clone.labels$Name <- as.character(clone.labels$Name)

    labeled.annotated.amps <- left_join(annotated.amps,clone.labels,by=c("Genome" = 'Name')) %>%
        select(-query.index,-subject.index,-is.significant,-SampleType, -Population) %>%
        ## replace 'sfcA' with 'maeA' in the plot.
        mutate(gene = replace(gene, gene == 'sfcA', 'maeA')) %>%
        mutate(log.pval=log(bonferroni.corrected.pval)) %>%
        mutate(log2.copy.number.mean=log2(copy.number.mean)) %>%
        transform(Population = PopulationLabel) %>% filter(!(Genome==ParentClone))

    ## order the genes by start to get axes correct on heatmap.
    labeled.annotated.amps$gene <- with(labeled.annotated.amps, reorder(gene, start))
    ## reverse the order of genomes to make axes consistent with Fig5B.
    labeled.annotated.amps$Genome <- factor(labeled.annotated.amps$Genome)
    labeled.annotated.amps$Genome <- factor(labeled.annotated.amps$Genome,
                                            levels=rev(levels(labeled.annotated.amps$Genome)))

    heatmap <- ggplot(labeled.annotated.amps,aes(x=gene,y=Genome,fill=log2.copy.number.mean,frame=Environment)) +
        geom_tile(color="white",size=0.1) +
        ## draw vertical lines at genes of interest.
        geom_vline(linetype='dashed',xintercept = which(levels(labeled.annotated.amps$gene) %in% c('citT','dctA','maeA'))) +
        xlab("Amplified genes") +
        scale_fill_viridis(name="",option="plasma") +
        facet_wrap(~Environment,nrow=2, scales = "free_y") +
        theme_tufte(base_family='Helvetica') +
        theme(axis.ticks=element_line(size=0.1)) +
        theme(axis.text.x=element_text(size=3,angle=90,hjust=1,vjust=0.2,face="italic")) +
        theme(legend.text=element_text(size=6)) +
        guides(fill=FALSE)

    return(heatmap)
}

## plot a stacked bar graph to show how big the amplifications are
## in the DM0 and DM25 treatments.
plot.Fig5B.stackedbar <- function(amps,clone.labels) {

    amps2 <- amps %>% mutate(total.amp.length = copy.number.mean*len) %>%
        mutate(log2.copy.number.mean=log2(copy.number.mean)) %>%
        left_join(clone.labels, by=c("Genome" = "Name")) %>% select(-Population) %>%
        transform(Population=PopulationLabel) %>% filter(!(Genome==ParentClone))

    stacked <- ggplot(amps2,aes(x=Genome,y=total.amp.length,fill=log2.copy.number.mean)) +
        geom_bar(stat="identity") +
        scale_fill_viridis(name=expression("log"[2]*"(copy number)"),option="plasma") +
        facet_grid(~Environment, scales = "free") +
        theme_tufte(base_family='Helvetica') +
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
        theme(legend.key.width=unit(1, "cm")) +
        ylab("Total amplification length")

    return(stacked)
}

home.dir <- path.expand("~")
projdir <- file.path(home.dir,"BoxSync/active-projects/DM0-evolution")
outdir <- file.path(projdir,"results")
breseq.output.dir <- file.path(projdir,"genomes/polymorphism")
LCA.gff3 <- file.path(projdir,"genomes/curated-diffs/LCA.gff3")
all.genomes <- list.files(breseq.output.dir,pattern='^ZDBp|^CZB')
all.genome.paths <- sapply(all.genomes, function(x) file.path(breseq.output.dir,x))
genome.input.df <- data.frame(Genome=all.genomes,path=all.genome.paths)

amps <- map2_df(genome.input.df$path,
                genome.input.df$Genome,
                find.amplifications)

write.csv(x=amps,file=file.path(outdir,"amplifications.csv"))

annotated.amps <- amps %>% annotate.amplifications(LCA.gff3)
write.csv(x=annotated.amps,file=file.path(outdir,"amplified_genes.csv"))

write.csv(x=filter(annotated.amps,Genome=='ZDBp874'),file=file.path(outdir,"ZDBp874_amplified_genes.csv"))

## just report parallel amps in evolved genomes.
amp.parallelism <- annotated.amps %>%
    filter(!(Genome %in% c('CZB151','CZB152','CZB154','ZDB67','ZDB68','ZDB69'))) %>%
    group_by(gene,locus_tag) %>% summarise(count=n()) %>% arrange(locus_tag)
write.csv(x=amp.parallelism,file=file.path(outdir,"amp_parallelism.csv"))

#' report copy number of maeA, citT, and dctA in a table.
#' (I can do this easily by hand).

copy.number.table <- annotated.amps %>%
    filter(gene %in% c('sfcA','dctA')) %>%
    #' change sfcA to maeA
    mutate(Gene=replace(gene,gene=='sfcA','maeA')) %>%
    transform(amplified.segment.length=len) %>%
    arrange(Gene,Genome) %>%
    select(Gene,Genome,amplified.segment.length,copy.number.mean,copy.number.min,copy.number.max, bonferroni.corrected.pval)

write.csv(x=copy.number.table,file=file.path(outdir,"copy_number_table.csv"))


#' Make figures.

label.filename <- file.path(projdir,"data/rohan-formatted/populations-and-clones.csv")
clone.labels <- read.csv(label.filename) %>% mutate(Name=as.character(Name))


#' Make a heatmap plot with facet grid on DM0 versus DM25.
#' color the matrix based on log2(copy number).
Fig5A <- plot.Fig5A.heatmap(annotated.amps,clone.labels)
#' Make a stacked bar plot, with facet grid on DM0 versus DM25.
#' color the stacked bars based on log2(copy number) for consistency.
Fig5B <- plot.Fig5B.stackedbar(amps,clone.labels)

save_plot(file.path(projdir,"results/figures/Fig5A.pdf"),Fig5A,base_aspect_ratio=1.5)
save_plot(file.path(projdir,"results/figures/Fig5B.pdf"),Fig5B,base_aspect_ratio=1.5)

#' Make and save the final Figure 5.
Fig5outf <- file.path(projdir,"results/figures/Fig5.pdf")
Fig5 <- plot_grid(Fig5A, Fig5B, labels = c('A', 'B'), ncol = 1, rel_heights = c(1.1, 1))
save_plot(Fig5outf,Fig5,base_height=7,base_width=10.5)

#' write out a matrix where row is 'maeA-AMP' or 'dctA-AMP'
#' and columns are genome names. This will be used to merge
#' the amplification information with the mutation matrix
#' data in order include these amplifications in the
#' co-occurrence analysis analysis.
write.amp.matrix <- function(annotated.amps,clone.labels, outfile) {
    amp.matrix.df <- left_join(annotated.amps,clone.labels,by=c("Genome" = 'Name')) %>%
    filter(!(Genome %in% ParentClone)) %>%
    mutate(Genome=factor(Genome)) %>%
    mutate(gene = replace(gene, gene == 'sfcA', 'maeA-AMP')) %>%
    mutate(gene = replace(gene, gene == 'dctA', 'dctA-AMP')) %>%
    filter(gene %in% c('maeA-AMP','dctA-AMP')) %>%
    select(Genome,gene)

    maeA.amp.matrix.df <- filter(amp.matrix.df, gene == 'maeA-AMP')
    dctA.amp.matrix.df <- filter(amp.matrix.df, gene == 'dctA-AMP')

    maeA.AMP.binary.vec <- sapply(levels(amp.matrix.df$Genome),function(x) ifelse(x %in% maeA.amp.matrix.df$Genome,1,0))
    dctA.AMP.binary.vec <- sapply(levels(amp.matrix.df$Genome),function(x) ifelse(x %in% dctA.amp.matrix.df$Genome,1,0))

    amp.matrix <- data.frame(rbind(maeA.AMP.binary.vec,dctA.AMP.binary.vec),row.names=NULL)
    amp.matrix$Gene <- c('maeA-AMP','dctA-AMP')
    write.csv(x=amp.matrix,file=outfile,row.names = FALSE)
}

write.amp.matrix(annotated.amps,clone.labels,file.path(outdir,"amp_matrix.csv"))

#' check the genetic background of the maeA and dctA amplifications.
#' strong association with dctA amplifications due to anti-correlation
#' to an ancestral dctA promoter mutation in the CZB151/154 clade.

#' maeA amplifications tend to occur in CZB151/154 rather than 152 background..
#' but Tanush's competition data shows that it's beneficial in both backgrounds?
maeA.dct.amps.df <- left_join(annotated.amps,clone.labels,by=c("Genome" = 'Name')) %>%
    filter(!(Genome %in% ParentClone)) %>%
    mutate(Genome=factor(Genome)) %>%
    mutate(gene = replace(gene, gene == 'sfcA', 'maeA-AMP')) %>%
    mutate(gene = replace(gene, gene == 'dctA', 'dctA-AMP')) %>%
    filter(gene %in% c('maeA-AMP','dctA-AMP')) %>%
    select(Genome,gene,Founder,ParentClone)
