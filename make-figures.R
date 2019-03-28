## make-figures.R by Rohan Maddamsetti.

## go through and remove unneeded dependencies once this script is more polished.
## there's some problem with dplyr::rename due to some of these imports I think.

library(cowplot)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggthemes)
library(viridis)
library(scales)
library(gridExtra)
library(lubridate)
library(assertthat)
library(tidyr)
library(stringr)
library(purrr)
library(tidytext)
library(Matrix)

home.dir <- path.expand("~")
proj.dir <- file.path(home.dir,"BoxSync/DM0-evolution")

pop.clone.labels <- read.csv(
    file.path(proj.dir,
        "data/rohan-formatted/populations-and-clones.csv"))

evolved.mutations <- read.csv(
    file.path(proj.dir,
        "results/genome-analysis/evolved_mutations.csv"))

###############################################
## Figure 1.

## Figure 1A is a schematic of the experimental design, made in Illustrator.
## Figure 1B is the phylogenetic tree constructed using the iTOL webserver.
## Figure 1C: make a stacked bar plot of the different kinds of mutations in each clone.
fig1C.stacked <- ggplot(evolved.mutations,aes(x=Clone,fill=Mutation)) +
    geom_bar() +
    scale_fill_viridis(option="magma",discrete=TRUE) +
    facet_grid(~Environment,scales="free") +
    theme_tufte(base_family='Helvetica') +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=12)) +
    theme(panel.border=element_blank()) +
    theme(plot.title=element_text(hjust=0)) +
    theme(strip.text=element_text(hjust=0,size=12)) +
    theme(panel.spacing.x=unit(1, "cm")) +
    theme(panel.spacing.y=unit(0.5, "cm")) +
    theme(legend.title=element_text(size=8)) +
    theme(legend.title.align=1) +
    theme(legend.text=element_text(size=8))

fig1C.output <- "../results/figures/Fig1C.pdf"
ggsave(fig1C.stacked, file=fig1C.output,width=7,height=7)

##############################################################
## Figure 2: growth curve analysis, showing evolution of
## growth rate and yield in DM0.

## TODO: similar analysis, but with DM25-evolved genomes
## (don't have these data in hand).

## Cite the growthcurver package. https://www.ncbi.nlm.nih.gov/pubmed/27094401
library(growthcurver)

Fig2A.analysis <- function(growth.df) {
    ## The analysis follows what Zack has done:
    ## 1) average blank measurement at every time point.
    ## 2) subtract average blank measurement from each time point.
    ## 3) log2 transform the data.

    ## This next line is just for debugging. keep commented out.
    ##growth.df <- labeled.growth.data

    blanks <- filter(growth.df, Name == 'Blank') %>% group_by(Experiment, Time) %>%
        summarise(blank.time.avg = mean(OD420))
    growth.df <- left_join(growth.df, blanks)
    growth.df <- mutate(growth.df, OD420.minus.blank.avg = OD420 - blank.time.avg)
    ## Add one so that the blank avg maps to zero after the log.
    growth.df <- mutate(growth.df, log2.OD420 = log2(OD420.minus.blank.avg+1))
    ## filter out blank rows.
    growth.df <- filter(growth.df,Name != 'Blank')
    ## Make an hours column with lubridate for plotting.
    growth.df <- mutate(growth.df,Hours = as.numeric(as.duration(hms(Time)))/3600)
    return(growth.df)
}

plot.Fig2A <- function(growth.df) {
    fig2A <- ggplot(growth.df, aes(x=Hours,y=log2.OD420,color=Name)) +
        geom_point(size=0.5) +
        facet_grid(Experiment ~ Founder) +
        theme_classic() +
        #theme_tufte(base_family='Helvetica') +
        #theme(axis.ticks=element_blank()) +
        geom_point(data=filter(growth.df,Name==Founder),size=0.5,color="grey80") +
        guides(color=FALSE) +
        ylab(expression('log'[2]*'(1 + OD420)'))
    return(fig2A)
}

growth.data <- read.csv(file.path(proj.dir,"results/growth-data.csv"))
labeled.growth.data <- left_join(growth.data,pop.clone.labels,by="Name")

final.growth.data <- Fig2A.analysis(labeled.growth.data)
Fig2A <- plot.Fig2A(final.growth.data)
save_plot(file.path(proj.dir,"results/figures/Fig2A.pdf"),Fig2A)
#####################################
## use growthcurver package to fit K and r and derive t0 from growth data.

## This function splits the growth plate data by Well and Experiment,
## Where Experiment is either 'DM0-growth' or 'DM25-growth'.
map.reduce.growth.curve <- function(labeled.plate.data) {

    gdata <- mutate(labeled.plate.data,Hours = as.numeric(as.duration(hms(Time)))/3600)

    ## make a dataframe of all pairs of levels in the Well and Experiment column.
    input.pairs <- expand.grid(levels(gdata$Well),levels(gdata$Experiment))

    ## returns one-row dataframe of growth rate fit and metadata from gdata.
    get.growth.curve.fit <- function(well,experiment) {
        my.data <- gdata %>% filter(Well==well) %>% filter(Experiment==experiment)
        gcfit <- SummarizeGrowth(my.data$Hours,my.data$OD420)
        gcfitv <- gcfit$val
        my.df = data.frame(Name = unique(my.data$Name),
                           Well = well,
                           Experiment = experiment,
                           SampleType = unique(my.data$SampleType),
                           Generation = unique(my.data$Generation),
                           Population = unique(my.data$Population),
                           PopulationLabel = unique(my.data$PopulationLabel),
                           AraStatus = unique(my.data$AraStatus),
                           ParentClone = unique(my.data$ParentClone),
                           Founder = unique(my.data$Founder),
                           Environment = unique(my.data$Environment),
                           Sequenced = unique(my.data$Sequenced),
                           k = gcfitv$k,
                           k_se = gcfitv$k_se,
                           k_p = gcfitv$k_p,
                           n0 = gcfitv$n0,
                           n0_se = gcfitv$n0_se,
                           n0_p = gcfitv$n0_p,
                           r = gcfitv$r,
                           r_se = gcfitv$r_se,
                           r_p = gcfitv$r_p,
                           sigma = gcfitv$sigma,
                           df = gcfitv$df,
                           t_mid = gcfitv$t_mid,
                           t_gen = gcfitv$t_gen,
                           auc_l = gcfitv$auc_l,
                           auc_e = gcfitv$auc_e)
        return(my.df)
    }

    gcfit.df <- bind_rows(map2(input.pairs$Var1,input.pairs$Var2, get.growth.curve.fit))
    return(gcfit.df)
}

growth.curve.fits <- map.reduce.growth.curve(labeled.growth.data)

## lag time t0 (that is, the midpoint of the sigmoid curve) can be
## derived by equating two equivalent parametrizations of the logistic
## growth curve (growthcurver documentation vs. wikipedia).
## A bit of algebra shows that t0 = log((K-N0)/N0)/r.
growth.curve.fit.summary <- growth.curve.fits %>% filter(Name != 'Blank') %>%
    group_by(Experiment,Name) %>% summarize(r.avg = mean(r),
                                            k.avg=mean(k),
                                            n0.avg=mean(n0),
                                            Environment=unique(Environment),
                                            Generation=unique(Generation),
                                            Founder=unique(Founder)) %>%
    mutate(t0=log((k.avg-n0.avg)/n0.avg)/r.avg)


DM0.growth.summary <- filter(growth.curve.fit.summary,Experiment=='DM0-growth') %>%
    rename(DM0.r = r.avg,DM0.k = k.avg,DM0.n0=n0.avg,DM0.t0=t0) %>%
    ungroup(Experiment) %>% select(-Experiment) %>%
    mutate(Generation=as.factor(Generation))

DM25.growth.summary <- filter(growth.curve.fit.summary,Experiment=='DM25-growth') %>%
    rename(DM25.r = r.avg,DM25.k = k.avg, DM25.n0=n0.avg,DM25.t0=t0) %>%
    ungroup(Experiment) %>% select(-Experiment) %>%
    mutate(Generation=as.factor(Generation))

growth.summary <- inner_join(DM0.growth.summary,DM25.growth.summary)

## calculate the log ratio of evolved growth to ancestral growth for both rate and yield,
## and then calculate a confidence intervals around the means, using
## the non-parametric bootstrap.
## This is a better statistical test for an increase in rate.

calc.growth.log.ratios <- function(growth.summary) {

    ## Since the data is small, go ahead and use a for loop
    ## to make vectors corresponding to ancestral R, K, t0
    ## in DM0 and DM25.
    anc.DM0.k <- rep(-1,nrow(growth.summary))
    anc.DM0.r <- rep(-1,nrow(growth.summary))
    anc.DM0.t0 <- rep(-1,nrow(growth.summary))
    anc.DM25.r <- rep(-1,nrow(growth.summary))
    anc.DM25.k <- rep(-1,nrow(growth.summary))
    anc.DM25.t0 <- rep(-1,nrow(growth.summary))
    for (index in 1:nrow(growth.summary)) {
        my.row <- growth.summary[index, ]
        my.anc <- filter(growth.summary,Name==my.row$Founder)
        anc.DM0.k[index] <- my.anc$DM0.k
        anc.DM0.r[index] <- my.anc$DM0.r
        anc.DM0.t0[index] <- my.anc$DM0.t0
        anc.DM25.k[index] <- my.anc$DM25.k
        anc.DM25.r[index] <- my.anc$DM25.r
        anc.DM25.t0[index] <- my.anc$DM25.t0
    }

    growth.summary2 <- growth.summary %>%
        mutate(log.DM0.k.ratio=log(DM0.k/anc.DM0.k)) %>%
        mutate(log.DM0.r.ratio=log(DM0.r/anc.DM0.r)) %>%
        mutate(log.DM0.t0.ratio=log(DM0.t0/anc.DM0.t0)) %>%
        mutate(log.DM25.k.ratio=log(DM25.k/anc.DM25.k)) %>%
        mutate(log.DM25.r.ratio=log(DM25.r/anc.DM25.r)) %>%
        mutate(log.DM25.t0.ratio=log(DM25.t0/anc.DM25.t0))
    return(growth.summary2)
}

final.growth.summary <- calc.growth.log.ratios(growth.summary)
evolved.growth.summary <- filter(final.growth.summary,Name != Founder)

## bootstrap confidence intervals around the mean.
## See http://www.stat.wisc.edu/~larget/stat302/chap3.pdf
## for sample code.
calc.bootstrap.conf.int <- function(vec) {
    B <- 10000
    n <- length(vec)
    boot.samples <- matrix(sample(vec, size = B*n, replace = TRUE),B, n)
    boot.statistics <- apply(boot.samples, 1, mean)
    boot.mean <- mean(boot.statistics)
    boot.se <- sd(boot.statistics)
    boot.confint <- boot.mean + c(-1,1)*2*boot.se
    return(boot.confint)
}

log.DM0.r.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.r.ratio)
log.DM0.k.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.k.ratio)
log.DM0.t0.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.t0.ratio)
log.DM25.r.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.r.ratio)
log.DM25.k.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.k.ratio)
log.DM25.t0.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.t0.ratio)

mean.log.DM0.r.ratio <- mean(evolved.growth.summary$log.DM0.r.ratio)
mean.log.DM0.k.ratio <- mean(evolved.growth.summary$log.DM0.k.ratio)
mean.log.DM0.t0.ratio <- mean(evolved.growth.summary$log.DM0.t0.ratio)
mean.log.DM25.r.ratio <- mean(evolved.growth.summary$log.DM25.r.ratio)
mean.log.DM25.k.ratio <- mean(evolved.growth.summary$log.DM25.k.ratio)
mean.log.DM25.t0.ratio <- mean(evolved.growth.summary$log.DM25.t0.ratio)

## significant increases in growth rate, no significant change in growth yield.
## Make a figure of these confidence intervals (Figure 2D).

bootstrap.results <- data.frame(Parameter=c("DM0 r",
                                            "DM25 r",
                                            "DM0 K",
                                            "DM25 K",
                                            "DM0 t0",
                                            "DM25 t0 "),
                                Estimate = c(mean.log.DM0.r.ratio,
                                             mean.log.DM25.r.ratio,
                                             mean.log.DM0.k.ratio,
                                             mean.log.DM25.k.ratio,
                                             mean.log.DM0.t0.ratio,
                                             mean.log.DM25.t0.ratio),
                                Left = c(log.DM0.r.ratio.conf.int[1],
                                         log.DM25.r.ratio.conf.int[1],
                                         log.DM0.k.ratio.conf.int[1],
                                         log.DM25.k.ratio.conf.int[1],
                                         log.DM0.t0.ratio.conf.int[1],
                                         log.DM25.t0.ratio.conf.int[1]),
                                Right = c(log.DM0.r.ratio.conf.int[2],
                                          log.DM25.r.ratio.conf.int[2],
                                          log.DM0.k.ratio.conf.int[2],
                                          log.DM25.k.ratio.conf.int[2],
                                          log.DM0.t0.ratio.conf.int[2],
                                          log.DM25.t0.ratio.conf.int[2]))

plot.Fig2B <- function (results) {
    the.plot <- ggplot(results,aes(x=Parameter,y=Estimate)) +
        geom_errorbar(aes(ymin=Left,ymax=Right),width=0.1, size=1) +
        geom_line() +
        geom_point(size=2) +
        geom_hline(yintercept=0,linetype='dashed') +
        ylab("log(Evolved/Ancestral)") +
        xlab("Growth parameter") +
        theme_classic()
    return(the.plot)
}

Fig2B <- plot.Fig2B(bootstrap.results)
save_plot(file.path(proj.dir,"results/figures/Fig2B.pdf"),Fig2B)

Fig2C <- ggplot(growth.summary,aes(x=DM25.r,y=DM0.r,color=Founder,shape=Generation)) +
    geom_point() +
    theme_classic() +
    xlab("DM25 growth rate") +
    ylab("DM0 growth rate") +
    scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a")) +
    coord_fixed() +
    guides(color=FALSE,shape=FALSE)

save_plot(file.path(proj.dir,"results/figures/Fig2C.pdf"),Fig2C)

Fig2D <- ggplot(growth.summary,aes(x=DM25.k,y=DM0.k,color=Founder,shape=Generation)) +
    geom_point() +
    theme_classic() +
    xlab("DM25 growth yield") +
    ylab("DM0 growth yield") +
    scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a")) +
    coord_fixed() +
    guides(color=FALSE,shape=FALSE)

save_plot(file.path(proj.dir,"results/figures/Fig2D.pdf"),Fig2D)

####### Make Figure 2 using cowplot.
Fig2outf <- file.path(proj.dir,"results/figures/Fig2.pdf")
Fig2BCD <- plot_grid(Fig2B,Fig2C,Fig2D, labels = c('B','C', 'D'), ncol = 3)
Fig2 <- plot_grid(Fig2A, Fig2BCD, labels = c('A', ''), ncol = 1, rel_heights = c(1.8, 1))
save_plot(Fig2outf,Fig2,base_height=7)

## significant correlation between growth rate in DM25 and DM0:
## Kendall's tau = 0.46, p = 0.0155
cor.test(x=final.growth.summary$DM25.r,y=final.growth.summary$DM0.r,method="kendall")
## no correlation in growth yield in DM25 and DM0
## Kendall's tau = -0.028, p = 0.9226
cor.test(x=final.growth.summary$DM25.k,y=final.growth.summary$DM0.k,method="kendall")

## POSSIBLE TODO IF GIVEN GROWTH RATE DATA FOR CLONES:
## Based on FBA modeling, it seems that maeA side-activity is essential to generate
## NADH during fermentation. Is there a relationship between maeA copy number and
## growth rate?
## Test this hypothesis by comparing growth results to maeA amplifications
## from the copy number analysis.
## These result also have dctA amplifications,
## so do a linear regression on both, with some random effects as well (representing other mutations).
##copy.number.results <- read.csv(file.path(proj.dir,"results/copy_number_table.csv"))

################################################################################
## Figure 3: make a matrix plot of genes with mutations in two or more clones.

## IMPORTANT TODO: refactor code so that I can make figures using just the dN results
## or with all valid mutation in the dice_analysis easily.

#' all valid mutations
raw.matrix <- read.csv(file.path(proj.dir,
        "results/DM0-DM25-comparison-mut-matrix.csv"))

#' just non-synonymous mutations
#' WARNING: THERE IS SOME BUG HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#raw.matrix <- read.csv(file.path(proj.dir,
#        "results/dN-DM0-DM25-comparison-mut-matrix.csv"))

## fix the names of the samples.
names(raw.matrix) <- map_chr(
    names(raw.matrix),
    function (x) str_trunc(x,width=7,side="right",ellipsis=''))

## import rows for the maeA and dctA amplifications.
amp.matrix <- read.csv(file.path(proj.dir,"results/amp_matrix.csv"))
## and merge with the mutation matrix.
merged.with.amps.matrix <- full_join(raw.matrix,amp.matrix)

DM0.DM25.matrix.data <- gather(merged.with.amps.matrix,"Name","mutation.count",2:26) %>%
    left_join(pop.clone.labels) %>%
    select(Gene,Name,mutation.count,Environment,PopulationLabel) %>%
    group_by(Gene) %>% filter(sum(mutation.count)>1)

fig1C.counts <- summarize(DM0.DM25.matrix.data,total.count=sum(mutation.count)) %>%
    arrange(desc(total.count))

## Now: add LTEE mutation matrix as well.
## filter on genes mutated in the DM0 and DM25 matrix data.

#' all valid mutations
##ltee.matrix <- read.csv(
##    file.path(
##        proj.dir,
##"results/LTEE-mut_matrix.csv"))

ltee.matrix <- read.csv(
    file.path(
        proj.dir,
"results/dN-LTEE-mut_matrix.csv"))


## fix the names of the samples.
names(ltee.matrix) <- map_chr(
    names(ltee.matrix),
    function (x) str_trunc(x,width=8,side="left",ellipsis=''))

LTEE.50K.labels <- read.csv(
    file.path(proj.dir,
        "data/rohan-formatted/LTEE-50K-clones.csv"))

ltee.data <- gather(ltee.matrix,"Name","mutation.count",2:13) %>%
    left_join(LTEE.50K.labels) %>%
    group_by(Gene) %>% filter(Gene %in% DM0.DM25.matrix.data$Gene) %>%
    filter(!is.na(Gene))

## keep non-mutators and the Ara-3 50K clone.
non.mutators.and.ara.minus.3 <- ltee.data %>%
    filter(Hypermutator == 0 | PopulationLabel == 'Ara-3') %>%
    select(-Hypermutator)

## Now join LTEE data to the DM0 and DM25 Cit+ data.
## Give the Ara-3 LTEE clone a special label.
matrix.data <- DM0.DM25.matrix.data %>%
    bind_rows(non.mutators.and.ara.minus.3) %>%
    mutate(Name=ifelse(Name=='REL11364','Ara-3: REL11364',Name))

## For the figure, we want to sort the genes by a
## hierarchical clustering. So convert to a matrix data structure,
## do the clustering, and use that ordering for the figure.

## Only use the DM0 and DM25 data for this step.

## note R's strange behavior when casting factors to numbers: have to
## change to character first!
mat.data <- DM0.DM25.matrix.data %>%
    select(-PopulationLabel,-Environment) %>%
    mutate(mutation.count=as.numeric(as.character(mutation.count)))

## turn these data into a matrix of class "dbCMatrix".
mut.matrix <- cast_sparse(mat.data,Gene,Name,mutation.count)

## convert to regular R matrix class
MxG.matrix <- as.matrix(mut.matrix)
## write out to file for analysis in jupyter notebook
## or using the clustergrammer web app:
## http://amp.pharm.mssm.edu/clustergrammer/
write.csv(MxG.matrix,file="../results/MxG_matrix.csv")
write.table(MxG.matrix,file="../results/MxG_matrix.tsv",sep="\t",quote=FALSE)

##and take transpose to turn matrix to Genome x Mutation.
GxM.matrix <- as.matrix(t(mut.matrix))
## write out to file for analysis in jupyter notebook.
write.csv(GxM.matrix,file="../results/GxM_matrix.csv")
write.table(GxM.matrix,file="../results/GxM_matrix.tsv",sep="\t")

## follow this post to use correlation as a distance.
## https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
MxM.scaled.cor.matrix <- cor(scale(GxM.matrix))

MxM.heatmap <- heatmap(MxM.scaled.cor.matrix)
## get the gene order from the heatmap..
sorted.genes <- sapply(MxM.heatmap$rowInd,function(i) rownames(MxM.scaled.cor.matrix)[i])

## now use this gene ordering for the mutation matrix figure.
matrix.data$Gene <- factor(matrix.data$Gene,levels=sorted.genes)
## cast mutation.count into a factor for plotting.
matrix.data$mutation.count <- factor(matrix.data$mutation.count)

matrix.figure <- ggplot(matrix.data,aes(x=Name,
                                        y=Gene,
                                        fill=mutation.count,
                                        frame=Environment)) +
    geom_tile(color="black",size=0.1) +
    ylab("Gene") +
    xlab("Genome") +
    facet_wrap(~Environment,ncol=3, scales = "free_x") +
    theme_tufte(base_family='Helvetica') +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x=element_text(size=10,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=10,hjust=1,face="italic")) +
    theme(legend.text=element_text(size=6)) +
    theme(legend.position="bottom") +
    theme(legend.key.size=unit(0.2, "cm")) +
    theme(legend.key.width=unit(1, "cm")) +
    guides(fill=FALSE) +
    scale_fill_manual(values = c("white", "#ffdf00", "#bebada",
                                 "#fb8072", "#80b1d3", "#fdb462"))

##matrix.outfile <- "../results/figures/Fig3.pdf"
dNmatrix.outfile <- "../results/figures/Fig3B.pdf"
##ggsave(matrix.figure, file=matrix.outfile,width=7,height=15)
ggsave(matrix.figure, file=dNmatrix.outfile,width=7,height=15)

## make a mutation co-occurrence matrix.
## nice to have, but I don't see anything interesting.

MxG.presence <- MxG.matrix
MxG.presence[MxG.presence > 1] <- 1
M.co.occurrence <- MxG.presence %*% t(MxG.presence)
M.co.occurrence.heatmap <- heatmap(M.co.occurrence)

## get the gene order from the heatmap.
M.co.occurrence.sorted.genes <- sapply(M.co.occurrence.heatmap$rowInd,function(i) rownames(M.co.occurrence)[i])

M.co.occurrence.df <- data.frame(row=rownames(M.co.occurrence)[row(M.co.occurrence)],
                                col=colnames(M.co.occurrence)[col(M.co.occurrence)],
                                co.occurrence=c(as.factor(M.co.occurrence)))

M.co.occurrence.df$row <- factor(M.co.occurrence.df$row,levels=M.co.occurrence.sorted.genes)
M.co.occurrence.df$col <- factor(M.co.occurrence.df$col,levels=M.co.occurrence.sorted.genes)

## plot the co-occurrence matrix.
M.co.occurrence.plot <- ggplot(M.co.occurrence.df,aes(row,col,fill=co.occurrence)) +
    geom_tile(color='White') +
    scale_fill_viridis(option="magma",direction=-1) +
    theme(axis.text.x = element_text(angle=90,vjust=-0.1)) +
    guides(fill=FALSE)
M.co.occurrence.plot

ggsave(file.path(proj.dir,"results/figures/co_occurrence.pdf"),
         M.co.occurrence.plot,width=12,height=12)


################################################################################
## analysis of parallel evolution at the same nucleotide.
## discuss numbers and finding in the text (no figure.).

bp.parallel.mutations <- evolved.mutations %>% group_by(Position) %>%
    summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.MOB <- filter(bp.parallel.mutations,Mutation=='MOB')
parallel.DEL <- filter(bp.parallel.mutations,Mutation=='DEL')
parallel.INS <- filter(bp.parallel.mutations,Mutation=='INS')
parallel.dN <- filter(bp.parallel.mutations,Mutation=='nonsynonymous')
parallel.dS <- filter(bp.parallel.mutations,Mutation=='synonymous')

## examine basepair-level parallel evolution in the polymorphism runs.

poly.evolved.mutations <- read.csv(
    file.path(proj.dir,
        "results/genome-analysis/poly_evolved_mutations.csv"))

poly.bp.parallel.mutations <- poly.evolved.mutations %>% filter(Frequency < 1) %>%
    filter(Frequency > 0.1) %>% group_by(Position) %>% summarise(count = n()) %>%
    filter(count>1) %>% inner_join(poly.evolved.mutations)

poly.parallel.MOB <- filter(poly.bp.parallel.mutations,Mutation=='MOB')
poly.parallel.DEL <- filter(poly.bp.parallel.mutations,Mutation=='DEL')
poly.parallel.INS <- filter(poly.bp.parallel.mutations,Mutation=='INS')
poly.parallel.dN <- filter(poly.bp.parallel.mutations,Mutation=='nonsynonymous')
poly.parallel.dS <- filter(poly.bp.parallel.mutations,Mutation=='synonymous')

## no consensus calls in citT, but extensive parallelism in polymorphic calls
## in citT.
citT.mutations <- filter(evolved.mutations,Gene=='citT')
citT.poly.mutations <- filter(poly.evolved.mutations,Gene=='citT')

################################################
## IS element analysis and visualization.

IS.insertions <- read.csv(file.path(proj.dir,
              "results/genome-analysis/IS_insertions.csv")) %>%
    arrange(genome_start)

IS.plot <- ggplot(IS.insertions,aes(x=genome_start,fill=IS_element,frame=Environment)) +
    facet_grid(Environment~.) +
    geom_histogram(bins=400) +
    guides(fill=FALSE) +
    ylab("Count") +
    xlab("Position") +
    theme_tufte(base_family="Helvetica")

save_plot(file.path(proj.dir,"results/figures/Fig4B.pdf"),
         IS.plot,base_aspect_ratio=1.5)


## 81/213 IS insertions recur at the same locations! 38%!
parallel.IS.insertions <- IS.insertions %>%
    group_by(genome_start) %>%
    filter(n()>1)

parallel.IS.summary <- summarize(parallel.IS.insertions,
                                 IS_element=unique(IS_element),
                                 genes_inactivated=unique(genes_inactivated),
                                 genes_promoter=unique(genes_promoter),
                                 count=n()) %>%
    mutate(genes_inactivated=as.character(genes_inactivated)) %>%
    mutate(genes_promoter=as.character(genes_promoter)) %>%
    mutate(annotation=ifelse(is.na(genes_inactivated),
                             genes_promoter,
                             genes_inactivated)) %>%
    mutate(annotation=ifelse(is.na(annotation),'none',annotation))

## plot parallel IS insertions and their annotation.
parallel.IS.plot <- ggplot(parallel.IS.summary,aes(x=genome_start,
                                           y=count,
                                           color=IS_element,
                                           label=annotation)) +
    geom_point() +
    guides(color=FALSE) +
    theme_classic() +
    ylab("Count") +
    xlab("Position") +
    geom_text_repel(fontface = "italic")

save_plot(file.path(proj.dir,"results/figures/Fig4A.pdf"),
         parallel.IS.plot,base_aspect_ratio=1.5)

########
## Plot the rate of increase of IS-elements in the DM0 and DM25 experiments
## in comparison to the rate of increase of IS-elements in Ara-3.

LTEE.MAE.IS.insertions <- read.csv(
    file.path(proj.dir,
              "results/genome-analysis/LTEE_MAE_IS150_insertions.csv")) %>%
    arrange(Position)

LTEE.IS150 <- filter(LTEE.MAE.IS.insertions,Environment=='LTEE')
MAE.IS150 <- filter(LTEE.MAE.IS.insertions,Environment=='MAE')
Ara.minus.3.IS150 <- filter(LTEE.IS150,Population=='Ara-3')

Ara.minus.3.IS150.by.clone <- group_by(Ara.minus.3.IS150,Clone,Generation,Environment,Population) %>%
    summarize(total.count=n())

MAE.IS150.over.time <- group_by(MAE.IS150,Clone,Generation,Environment,Population) %>%
    summarize(total.count=n())

DM0.DM25.over.time <- group_by(IS.insertions,Clone,Generation,Environment,Population) %>%
    summarize(total.count=n())

IS150.rate.df <- rbind(MAE.IS150.over.time,DM0.DM25.over.time,Ara.minus.3.IS150.by.clone)

# The colorblind-friendly palette with black:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

IS150.rate.plot <- ggplot(IS150.rate.df,
                          aes(x=Generation,y=total.count,color=Environment)) +
    theme_classic() +
    scale_colour_manual(values=cbPalette) +
    ylab('IS150 insertions') +
    geom_jitter(width=50) +
    guides(color=FALSE)

save_plot(file.path(proj.dir,"results/figures/Fig4C.pdf"),
         IS150.rate.plot,base_aspect_ratio=1.5)


########
## Combine the IS plots with cowplot to make Figure 4.
Fig4outf <- file.path(proj.dir,"results/figures/Fig4.pdf")
Fig4 <- plot_grid(parallel.IS.plot, IS.plot, IS150.rate.plot,
                  labels = c('A', 'B', 'C'), ncol=1)
save_plot(Fig4outf,Fig4,base_height=10,base_aspect_ratio=0.8)

########
## Conduct the following test for parallel evolution of IS-insertions:
## Assume all IS insertions occur in the set union(A,B), with
## probability = empirical mass distribution over LTEE/MAE/DM0/DM25.
## Set the number of samples to the number of IS150 insertions
## in the DM0 genomes. For 10000 replicates,
## count how often some site gets hit 11 times. I could also repeat this test for the
## chance that a particular site gets hit 11 times (for instance).

## filter out mutations in clones on the same lineage
## by removing duplicates within the same population.
LTEE.MAE.IS150.hit.pos <- LTEE.MAE.IS.insertions %>%
    select(-Clone) %>% distinct() %>%
    group_by(Position,GenePosition) %>%
    summarize(count=n()) %>% ungroup()

## again, filter out duplicates on line of descent
## by removing duplicates within the same population.
DM0.DM25.IS150.hit.pos <- IS.insertions %>%
    filter(IS_element == 'IS150') %>% mutate(Position=genome_start) %>%
    select(Position,Environment,Population,GeneName,GenePosition) %>%
    distinct() %>%
    group_by(Position,GenePosition) %>%
    summarize(count=n()) %>% ungroup()

## use GenePosition to match up IS element insertion positions.
## number of DM0/DM25 IS150 events with sites in LTEE/MAE: 22.
DM0.DM25.pos.intersect <- filter(DM0.DM25.IS150.hit.pos,
                                 GenePosition %in% LTEE.MAE.IS150.hit.pos$GenePosition)
## number of DM0/DM25 IS150 events in total: 124. 22/124 = 0.1774 in LTEE/MAE set.
length(DM0.DM25.IS150.hit.pos$GenePosition)

## make an empirical mass function over GenePosition.
total.hit.pos <- rbind(
    select(DM0.DM25.IS150.hit.pos,-Position),
    select(LTEE.MAE.IS150.hit.pos,-Position)) %>%
    group_by(GenePosition) %>% summarize(count2 = sum(count)) %>%
    arrange(desc(count2)) %>% mutate(prob=count2/sum(count2)) %>%
    mutate(eCDF=cumsum(prob))

## number of draws for each simulation = IS150 in DM0.
DM0.IS150 <- IS.insertions %>% filter(IS_element == 'IS150',Environment=='DM0') %>%
    select(-Clone) %>% distinct()
draws <- nrow(DM0.IS150)

null.parallel.hits <- function(total.hit.pos,empirical.parallel=11,replicates=10000) {

    max.hits <- function(total.hit.pos) {
        my.sample <- sample(x=total.hit.pos$GenePosition,
                            size=draws,replace=TRUE,prob=total.hit.pos$prob)
        my.sample.parallel.hits <- data.table(my.sample)[, .N, keyby = my.sample]
        max.sampled.parallelism <- max(my.sample.parallel.hits$N)
        return(max.sampled.parallelism)
    }

    max.hit.vec <- unlist(map(seq_len(replicates), ~max.hits(total.hit.pos)))
    max.hit.df.col <- data.frame('max.hit'=max.hit.vec)
    past.threshold <- nrow(filter(max.hit.df.col,max.hit>=empirical.parallel))
    return(past.threshold/replicates)
}

## empirical p-value for 11 hits is on the order of 0.0001.
null.parallel.hits(total.hit.pos,replicates=10000)
null.parallel.hits(total.hit.pos,replicates=100000)

## empirical p-value for 8 hits is on the order of 0.01.
null.parallel.hits(total.hit.pos,empirical.parallel=8)

################################################################################
## Figure 8: Fitness and growth effects of plasmid-borne maeA expression.
Fig8.analysis <- function(data, samplesize=6, days.competition=1,rev=FALSE) {
    ## by diluting stationary phase culture 1:100 on day 0, there is another
    ## factor of 100 that multiplies the day 0 plating dilution factor.
    data$D.0 <- 100*data$D.0
    data$Mr <- log2((data$D.1/data$D.0)*data$Red.1*100^days.competition/data$Red.0)/days.competition
    data$Mw <- log2((data$D.1/data$D.0)*data$White.1*100^days.competition/data$White.0)/days.competition
## reverse log ratio when polarity is reversed.
if (rev) data$W <- data$Mr/data$Mw else data$W <- data$Mw/data$Mr

my.mean <- mean(data$W, na.rm=T)
my.sd <- sd(data$W, na.rm=T)
my.confint <- c(my.mean - 1.96*my.sd/sqrt(samplesize), my.mean + 1.96*my.sd/sqrt(samplesize))

print("mean is:")
print(my.mean)

print("confint is:")
print(my.confint)

##return a dataframe of the results for plotting.
left.error <- c(my.confint[1])
right.error <- c(my.confint[2])
results <- data.frame(Fitness=c(my.mean),Left=left.error,Right=right.error)
return(results)
}

fig8.data <- read.csv("../data/rohan-formatted/DM0_Fitness2.csv",header=TRUE)
## these data come from one day competitions that Tanush ran.
res1 <- filter(fig8.data,Red.Pop=='ZDB151_with_maeA') %>% Fig8.analysis()
res2 <- filter(fig8.data,Red.Pop=='ZDB67_with_maeA') %>% Fig8.analysis(rev=TRUE)
res3 <- filter(fig8.data,Red.Pop=='ZDB152_with_maeA') %>% Fig8.analysis()
res4 <- filter(fig8.data,Red.Pop=='ZDB68_with_maeA') %>% Fig8.analysis(rev=TRUE)

## beautiful! All confints overlap with each other, showing no maeA fitness effect
## does not depend on Ara polarity or genetic background.

## let's combine data from the switched polarity competitions (since Ara marker is neutral)
d1 <- filter(fig8.data,Red.Pop=='ZDB151_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d2 <- filter(fig8.data,Red.Pop=='ZDB67_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## now switch column labels.
d2X <- rename(d2,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB151.data <- rbind(d1,d2X)

d3 <- filter(fig8.data,Red.Pop=='ZDB152_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d4 <- filter(fig8.data,Red.Pop=='ZDB68_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## switch column labels.
d4X <- rename(d4,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB152.data <- rbind(d3,d4X)

fres1 <- Fig8.analysis(ZDB151.data,samplesize=6)
fres2 <- Fig8.analysis(ZDB152.data,samplesize=6)
fig8.plot.data <- rbind(fres1,fres2)
#' Correct the strain names here.
fig8.plot.data$Strain <- c('CZB151','CZB152')

## Make Figure 8.
fig8.output <- "../results/figures/Fig8.pdf"

plot.Fig8 <- function (results, output.file) {
    the.plot <- ggplot(results,aes(x=Strain,y=Fitness)) +
        geom_errorbar(aes(ymin=Left,ymax=Right),width=0.1, size=1) +
        geom_line() +
        geom_point(size=2) +
        scale_y_continuous(limits=c(1.0,1.30)) +
        ylab("Fitness of maeA plasmid relative to empty plasmid") +
        theme_classic()
    ggsave(the.plot, file=output.file,width=4,height=4)
}

plot.Fig8(fig8.plot.data,fig8.output)

########################################################

