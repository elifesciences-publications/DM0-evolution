## make-figures.R by Rohan Maddamsetti.

## go through and remove unneeded dependencies once this script is more polished.

## there's some problem with dplyr::rename due to some of these imports I think.

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
library(ggfortify)
library(kernlab)

home.dir <- path.expand("~")
proj.dir <- file.path(home.dir,"BoxSync/DM0-evolution")

pop.clone.labels <- read.csv(
    file.path(proj.dir,
        "data/rohan-formatted/populations-and-clones.csv"))

## Figure 1A is diagram of the DM0 and DM25 experiments made in Illustrator.

evolved.mutations <- read.csv(
    file.path(proj.dir,
        "results/genome-analysis/evolved_mutations.csv"))

## examine basepair-level parallel evolution.

bp.parallel.mutations <- evolved.mutations %>% group_by(Position) %>%
    summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.MOB <- filter(bp.parallel.mutations,Mutation=='MOB')
parallel.DEL <- filter(bp.parallel.mutations,Mutation=='DEL')
parallel.INS <- filter(bp.parallel.mutations,Mutation=='INS')
parallel.dN <- filter(bp.parallel.mutations,Mutation=='nonsynonymous')

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


## Figure 1B: make a stacked bar plot of the different kinds of mutations in each clone.
fig1B.stacked <- ggplot(evolved.mutations,aes(x=Clone,fill=Mutation)) +
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

fig1B.output <- "../results/figures/Fig1B.pdf"
ggsave(fig1B.stacked, file=fig1B.output,width=7,height=7)

## Figure 1C: make a matrix plot of genes with mutations in two or more clones.
fig1C.raw.matrix <- read.csv(file.path(proj.dir,
        "results/DM0-DM25-comparison-mut-matrix.csv"))

## fix the names of the samples.
names(fig1C.raw.matrix) <- map_chr(
    names(fig1C.raw.matrix),
    function (x) str_trunc(x,width=7,side="right",ellipsis=''))

fig1C.DM0.DM25.data <- gather(fig1C.raw.matrix,"Name","mutation.count",2:26) %>%
    left_join(pop.clone.labels) %>%
    select(Gene,Name,mutation.count,Environment,PopulationLabel) %>%
    group_by(Gene) %>% filter(sum(mutation.count)>1)

fig1C.counts <- summarize(fig1C.DM0.DM25.data,total.count=sum(mutation.count)) %>%
    arrange(desc(total.count))

## Now: add LTEE mutation matrix as well.
## filter on genes mutated in the DM0 and DM25 matrix data.
ltee.matrix <- read.csv(
    file.path(
        proj.dir,
        "results/LTEE-mut_matrix.csv"))
## fix the names of the samples.
names(ltee.matrix) <- map_chr(
    names(ltee.matrix),
    function (x) str_trunc(x,width=8,side="left",ellipsis=''))

LTEE.50K.labels <- read.csv(
    file.path(proj.dir,
        "data/rohan-formatted/LTEE-50K-clones.csv"))

ltee.data <- gather(ltee.matrix,"Name","mutation.count",2:13) %>%
    left_join(LTEE.50K.labels) %>%
    group_by(Gene) %>% filter(Gene %in% fig1C.DM0.DM25.data$Gene) %>%
    filter(!is.na(Gene))

## Filter out hypermutators.
non.mutator.ltee.data <- ltee.data %>%
    filter(Hypermutator == 0) %>%
    select(-Hypermutator)

## Now join LTEE data to the DM0 and DM25 Cit+ data.
fig1C.data <- fig1C.DM0.DM25.data %>% bind_rows(non.mutator.ltee.data)

## later I try a different ordering by hierarchical clustering.
## order the genes by number of mutations to get axes correct on heatmap.
fig1C.data$Gene <- factor(fig1C.data$Gene)
fig1C.data$mutation.count <- factor(fig1C.data$mutation.count)

fig1C <-  ggplot(fig1C.data,aes(x=Name,y=Gene,fill=mutation.count,
                                frame=Environment)) +
    geom_tile(color="black",size=0.1) +
    ylab("Gene") +
    facet_wrap(~Environment,ncol=3, scales = "free_x") +
    theme_tufte(base_family='Helvetica') +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x=element_text(size=10,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=10,hjust=1,face="italic")) +
    theme(legend.text=element_text(size=6)) +
    theme(legend.position="bottom") +
    theme(legend.key.size=unit(0.2, "cm")) +
    theme(legend.key.width=unit(1, "cm")) +
    scale_fill_manual(values = c("white", "#ffdf00", "#bebada",
                                 "#fb8072", "#80b1d3", "#fdb462"))

fig1C.output <- "../results/figures/Fig1C.pdf"
ggsave(fig1C, file=fig1C.output,width=7,height=15)

## TODO: add a column corresponding to the valid mutations in Ara-3.
## might be better to make as a separate figure and then assemble together
## using cowplot.

############################################################################
## Recapitulation Index analysis:
## Percent of valid mutations in the DM0 and DM25 treatments that affected genes
## that subsequently acquired a valid mutation in the Ara-3 50K clone.

## PROBLEM: being a hypermutator causes high recapitulation indices in the absence
## of true parallel selection.

A.minus3.50K.df <- read.csv(file.path(proj.dir,"results/ara3_mut_matrix.csv")) %>%
    rename(mutation.count= REL11364_minus_CZB154)


DM0.mutations <- filter(fig1C.DM0.DM25.data,Environment=='DM0')  %>%
    select(-Name,Environment,PopulationLabel)
DM25.mutations <- filter(fig1C.DM0.DM25.data,Environment=='DM25') %>%
    select(-Name,Environment,PopulationLabel)


CalcRecapitulationIndex <- function (treat.muts,future.muts) {
    future.muts <- filter(future.muts,mutation.count > 0)
    recap.mutations <- filter(treat.muts,Gene %in% future.muts$Gene)
    total.recaps <- sum(recap.mutations$mutation.count)
    total.treat.muts <- sum(treat.muts$mutation.count)
    RI <- total.recaps/total.treat.muts
    return(RI)
}

DM0.RI <- CalcRecapitulationIndex(DM0.mutations,A.minus3.50K.df) ## 34.3% for DM0
DM25.RI <- CalcRecapitulationIndex(DM25.mutations,A.minus3.50K.df) ## 32% for DM25

## Then calculate RI for each population separately, excepting Ara-3.
all.but.matrix <- ltee.matrix %>% select(-REL11364)

## initialize the RI vectors.
DM0.RI.vec <- rep(0,ncol(ltee.matrix2)-1)
DM25.RI.vec <- rep(0,ncol(ltee.matrix2)-1)
LTEE.clone.vec <- rep('',ncol(ltee.matrix2)-1)

for (i in 2:12) {
    my.mut.summary <- all.but.matrix[c(1,i)]
    LTEE.clone.vec[i-1] <- names(my.mut.summary)[2]
    names(my.mut.summary)[2] <- 'mutation.count'
    my.mut.summary <- filter(my.mut.summary,mutation.count > 0)
    DM0.RI.vec[i-1] <- CalcRecapitulationIndex(DM0.mutations, my.mut.summary)
    DM25.RI.vec[i-1] <- CalcRecapitulationIndex(DM25.mutations, my.mut.summary)
}

## Now add Ara-3 clone REL11364 to the vectors, and turn into a dataframe.
DM0.RI.col <- c(DM0.RI.vec,DM0.RI)
DM25.RI.col <- c(DM25.RI.vec,DM25.RI)
Name.col <- c(LTEE.clone.vec,'REL11364')
RI.df <- data.frame(Name=Name.col,DM0.RI=DM0.RI.col,DM25.RI=DM25.RI.col)


############################################################################
## let's play with the vectors in fig1C.data.

## note R's strange behavior when casting factors to numbers: have to
## change to character first!
## question: why does PCA result in better clustering with the WRONG encoding?
## my guess is that Linear Discriminant analysis is really what is desired here
## (best view to distiguish classes in genome space).

mat.data <- fig1C.data %>%
    select(-PopulationLabel,-Environment) %>%
    mutate(mutation.count=as.numeric(as.character(mutation.count)))
    ##mutate(mutation.count=as.numeric(mutation.count))

## turn these data into a matrix of class "dbCMatrix".
mut.matrix <- cast_sparse(mat.data,Gene,Name,mutation.count)

## convert to regular R matrix class,
##and take transpose to turn matrix to Genome x Mutation.
GxM.matrix <- as.matrix(t(mut.matrix))
## write out to file for analysis in jupyter notebook.
write.csv(GxM.matrix,file="../results/GxM_matrix.csv")

## cast into a dataframe for annotating the PCA with ggfortify.
GxM.df <- data.frame(GxM.matrix)
GxM.df$Genome <- as.factor(rownames(GxM.matrix))
GxM.df$Environment <- map_chr(
    GxM.df$Genome,
    function(x) unique(filter(fig1C.data,Name==x)$Environment))

## plot a PCA on the GxM matrix.
GxM.PCA <- prcomp(GxM.df[seq(1:ncol(GxM.matrix))],scale=TRUE)

## plot PCA again, using ade4 and factoextra packages.
## see: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/119-pca-in-r-using-ade4-quick-scripts/
GxM.PCA2 <- dudi.pca(GxM.df[seq(1:ncol(GxM.matrix))],scannf = FALSE, nf = 10)
## and visualize.
fviz_eig(GxM.PCA2)

fviz_pca_ind(GxM.PCA2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

fviz_pca_var(GxM.PCA2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

## Do spectral clustering on genomes.
sc <- specc(GxM.matrix, kernel='rbfdot', centers=3)
## add spectral clustering information to dataframe for plotting.
GxM.df$spectral.clustering <- as.factor(sc)

## plot the loadings
## (which mutations distinguish between different environments).
quartz()
autoplot(GxM.PCA,
         data=GxM.df,
         colour='Environment',
         shape='spectral.clustering',
         loadings=TRUE,loadings.label = TRUE) + theme_classic()

autoplot(GxM.PCA,
         data=GxM.df,
         colour='Environment',
         shape='spectral.clustering') + theme_classic()

## do the same for the transpose.
MxG.matrix <- as.matrix(mut.matrix)
## write out to file for analysis in jupyter notebook.
write.csv(MxG.matrix,file="../results/MxG_matrix.csv")
write.table(MxG.matrix,file="../results/MxG_matrix.tsv",sep="\t",quote=FALSE)

## clustergrammer view.
## http://amp.pharm.mssm.edu/clustergrammer/viz_sim_mats/5c59d632c902bb0430741ef0/MxG_matrix.tsv

## cast into a dataframe for annotating the PCA with ggfortify.
MxG.df <- data.frame(MxG.matrix)
MxG.df$Genome <- rownames(MxG.matrix)

quartz() ## for testing.

## plot a PCA on the MxG matrix.
MxG.PCA <- prcomp(MxG.df[seq(1:ncol(MxG.matrix))],scale=TRUE)
autoplot(MxG.PCA, label=TRUE)

## plot PCA again, using ade4 and factoextra packages.
MxG.PCA2 <- dudi.pca(MxG.df[seq(1:ncol(MxG.matrix))],scannf = FALSE, nf = 10)
## and visualize.
fviz_eig(MxG.PCA2)

fviz_pca_ind(MxG.PCA2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

fviz_pca_var(MxG.PCA2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

## make the correlation matrices
GxG.correlation.matrix <- cor(MxG.matrix)
MxM.correlation.matrix <- cor(GxM.matrix)
autoplot(GxG.correlation.matrix)
autoplot(MxM.correlation.matrix)


## https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering

## following this post to turn a correlation into a distance measure.
MxM.scaled.cor.matrix <- cor(scale(GxM.matrix))

MxM.heatmap <- heatmap(MxM.scaled.cor.matrix)
## get the gene order from the heatmap..
sorted.genes <- sapply(MxM.heatmap$rowInd,function(i) rownames(MxM.scaled.cor.matrix)[i])

## does this order improve how Fig 1C looks? Yes it does!
testfig1C.data <- fig1C.data
testfig1C.data$Gene <- factor(fig1C.data$Gene,levels=sorted.genes)

testfig1C <-  ggplot(testfig1C.data,aes(x=Name,y=Gene,fill=mutation.count,
                                frame=Environment)) +
    geom_tile(color="black",size=0.1) +
    ylab("Gene") +
    facet_wrap(~Environment,ncol=3, scales = "free_x") +
    theme_tufte(base_family='Helvetica') +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x=element_text(size=10,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=10,hjust=1,face="italic")) +
    theme(legend.text=element_text(size=6)) +
    theme(legend.position="bottom") +
    theme(legend.key.size=unit(0.2, "cm")) +
    theme(legend.key.width=unit(1, "cm")) +
    scale_fill_manual(values = c("white", "#ffdf00", "#bebada",
                                 "#fb8072", "#80b1d3", "#fdb462"))

testfig1C.output <- "../results/figures/testFig1C.pdf"
ggsave(testfig1C, file=testfig1C.output,width=7,height=15)

## make a mutation co-occurrence matrix.
MxG.presence <- MxG.matrix
MxG.presence[MxG.presence > 1] <- 1
M.co.occurrence <- MxG.presence %*% t(MxG.presence)
M.co.occurrence.heatmap <- heatmap(M.co.occurrence)

## get the gene order from the heatmap..
M.co.occurrence.sorted.genes <- sapply(M.co.occurrence.heatmap$rowInd,function(i) rownames(M.co.occurrence)[i])


M.co.occurrence.df <- data.frame(row=rownames(M.co.occurrence)[row(M.co.occurrence)],
                                col=colnames(M.co.occurrence)[col(M.co.occurrence)],
                                co.occurrence=c(as.factor(M.co.occurrence)))

M.co.occurrence.df$row <- factor(M.co.occurrence.df$row,levels=M.co.occurrence.sorted.genes)
M.co.occurrence.df$col <- factor(M.co.occurrence.df$col,levels=M.co.occurrence.sorted.genes)

M.co.occurrence.plot <- ggplot(M.co.occurrence.df,aes(row,col,fill=co.occurrence)) + geom_tile(color='White') + scale_fill_viridis(option="magma",direction=-1) + theme(axis.text.x = element_text(angle=90,vjust=-0.1))
M.co.occurrence.plot

## look at largest and smallest eigenvalues and eigenvectors.
MxM.eigenresults <- eigen(MxM.correlation.matrix)

autoplot(MxM.eigenresults$vectors)

## do hierarchical clustering on the matrix.
## dist.binary from the ade4 package.
GxM.matrix2 <- GxM.matrix
## turn double hits into single to analyze as binary data.
GxM.matrix2[GxM.matrix2 > 1] <- 1
## method 5 is the Dice distance metric.
GxM.dist.mat <- dist.binary(GxM.matrix2,method=5)
clust.GxM <- hclust(GxM.dist.mat)
new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb")
heatmap(GxM.matrix2,col=new.palette(20))


############################################################################
## Figure 2:
## A) Fitness of evolved clones in DM0 and DM25.
## B) Fitness in DM0 vs. Fitness in DM25 (see if any tradeoff).

## NOTE: This figure still needs work, probably will use clone data instead anyway.

DM0.labels <- read.csv(file.path(proj.dir,
              "data/rohan-formatted/labels-for-DM0-evolved-DM0-fitness.csv"))

DM0.data <- read.csv(file.path(proj.dir,
              "data/rohan-formatted/DM0-evolved-DM0-fitness.csv"))

DM0.data2 <- left_join(DM0.data,DM0.labels) %>%
    mutate(Environment='DM0') %>%
    select(-Set,-Type)

DM25.labels <- read.csv(file.path(proj.dir,
              "data/rohan-formatted/labels-for-DM0-evolved-DM25-fitness.csv"))

DM25.data <- read.csv(file.path(proj.dir,
                                "data/rohan-formatted/DM0-evolved-DM25-fitness.csv"))

DM25.data2 <- left_join(DM25.data,DM25.labels) %>% mutate(Environment='DM25')

founder.labels <- select(pop.clone.labels,Name,Founder) %>% rename(C1=Name)
## put DM0 and DM25 data into a single data frame, and add a Founder column,
## and turn Dilution into a number (Plate Dilution).
Fig2.data <- bind_rows(DM0.data2,DM25.data2) %>%
    left_join(founder.labels) %>%
    mutate(Plate.Dilution = sapply(as.character(Dilution),
                                   function(x) {eval(parse(text=x))}))

## filter out all plates with TMTC colonies, and only analyze one-day competition
## data, and turn Red and White columns (factors) into numeric values
## (since no more TMTC).
## WARNING: as.numeric(Red) is WRONG (gets index of the factor).
## Need to use as.numeric(as.character(Red)) instead.
filtered.Fig2.data <- filter(Fig2.data,Red != 'TMTC') %>%
    filter(White != 'TMTC') %>%
    filter(Day != 'D3') %>% filter(Day != 'D6') %>%
    mutate(Red=as.numeric(as.character(Red))) %>%
    mutate(White=as.numeric(as.character(White)))

## split the data into two data frames for the initial and final day,
## and only include 10^3 dilution D0 and 10^5 dilution D1 data.
Fig2.day0.data <- filter(filtered.Fig2.data,Day=='D0',Dilution=='10^3')
Fig2.day1.data <- filter(filtered.Fig2.data,Day=='D1',Dilution=='10^5')

## change column names of initial.df and final.df before merging them together,
## and drop the Day and Dilution columns.
Fig2.day0.df <- rename(Fig2.day0.data,Red.Pop=C1, White.Pop=C2,
                      Red.0=Red,White.0=White, D.0=Plate.Dilution) %>% select(-Day,-Dilution)
Fig2.day1.df <- rename(Fig2.day1.data,Red.Pop=C1, White.Pop=C2,
                      Red.1=Red,White.1=White, D.1=Plate.Dilution) %>% select(-Day,-Dilution)
## Combine these data frames into a nice format for competition analysis (see Fig7 data).

## This is to help make the Rev column. Reverse if White Pop is a parent clone.
parent.clones <- unique(pop.clone.labels$Parent.Clone)

## Add the Rev and Rep columns.
Fig2.competition.data <- inner_join(Fig2.day0.df,Fig2.day1.df) %>%
    group_by(Environment,Founder,Red.Pop,White.Pop) %>%
    mutate(Rep = row_number()) %>%
    ## The Rev column ensures M_final/M_initial is always evolved/ancestor.
    mutate(Rev=Red.Pop %in% parent.clones)

## by diluting stationary phase culture 1:100 on day 0, there is another
## factor of 100 that multiplies the day 0 plating dilution factor.
Fig2.competition.data2 <- mutate(Fig2.competition.data, D.0=100*D.0) %>%
## always one competitions in these experiments.
    mutate(Mr = log2((D.1/D.0)*Red.1*100^1/Red.0)/1) %>%
    mutate(Mw = log2((D.1/D.0)*White.1*100^1/White.0)/1) %>%
    mutate(W = ifelse(Rev,Mw/Mr,Mr/Mw)) ## reverse log ratio depending on ancestor polarity.
## the data is still grouped_by. so run summarise to get results.

Fig2.results <- summarise(Fig2.competition.data2,
                          Rev=unique(Rev),
                          samplesize=n(),
                          Fitness=mean(W,na.rm=TRUE),
                          W.sd=sd(W,na.rm=TRUE),
                          left.error = Fitness - 1.96*W.sd/sqrt(samplesize),
                          right.error = Fitness + 1.96*W.sd/sqrt(samplesize),
                          Pop = ifelse(Rev,
                                       as.character(White.Pop),
                                       as.character(Red.Pop))) %>% na.omit()


Fig2A.output <- "../results/figures/Fig2A.pdf"

plot.Figure2A <- function (results, output.file) {
    the.plot <- ggplot(results,aes(x=Pop,y=Fitness,color=Rev, label=Pop)) +
        geom_errorbar(aes(ymin=left.error,ymax=right.error),width=0.1, size=1) +
        geom_point(size=1) +
        facet_wrap(~Environment+Founder, scales = "free_x") +
        scale_y_continuous(breaks=seq(0, 2.5, 0.2), limits=c(0.5,2.2)) +
        geom_hline(yintercept=1,linetype="dashed") +
        theme_tufte(base_family='Helvetica') +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        theme(axis.title.x=element_blank()) +
        theme(axis.title.y=element_text(hjust=1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank())

    ggsave(the.plot, file=output.file)
}

plot.Figure2B <- function(results,output.file) {
    #' compare fitness in DM0 to fitness in DM25.

    DM0.fitness.results <- filter(results,Environment=='DM0') %>%
        rename(DM0.Fitness = Fitness) %>%
        rename(DM0.W.sd = W.sd) %>%
        rename(DM0.left.error = left.error) %>%
        rename(DM0.right.error = right.error) %>%
        ungroup(Environment) %>%
        select(-Environment,-samplesize,-Rev)

    DM25.fitness.results <- filter(results,Environment=='DM25') %>%
        rename(DM25.Fitness=Fitness) %>%
        rename(DM25.W.sd=W.sd) %>% ungroup(Environment) %>%
        rename(DM25.left.error = left.error) %>%
        rename(DM25.right.error = right.error) %>%
        select(-Environment,-samplesize,-Rev)

    comparison.results <- inner_join(DM0.fitness.results,DM25.fitness.results)

    the.plot <- ggplot(comparison.results,aes(x=DM25.Fitness,y=DM0.Fitness, label=Pop,color='red')) +
        geom_errorbar(aes(ymin=DM0.left.error,ymax=DM0.right.error), width=0.01, size=1) +
        geom_errorbarh(aes(xmin=DM25.left.error,xmax=DM25.right.error), height=0.01, size=1) +
        geom_point(size=1) +
        guides(colour=FALSE) +
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        theme(axis.title.y=element_text(hjust=1)) +
        geom_abline(slope=1,intercept=0,linetype='dashed') +
        geom_text_repel(size=4,force=50,hjust=-0.5,vjust=-1,color='black')

    ggsave(the.plot, file=output.file)

}

plot.Figure2A(Fig2.results,file.path(home.dir,"/Desktop/Fig2A.pdf"))
plot.Figure2B(Fig2.results,file.path(home.dir,"/Desktop/Fig2B.pdf"))


## print out calculations to compare to Zack's calculations.
## write.csv(Fig2.competition.data2,file.path(home.dir,"/Desktop/Fig2_fitness.csv")_
## write.csv(Fig2.results,file.path(home.dir,"/Desktop/Fig2_fitness_results.csv"))

################################################################################
## Figure 3: Growth curves of evolved populations in DM0 and in DM25.

## Cite the growthcurver package. https://www.ncbi.nlm.nih.gov/pubmed/27094401
library(growthcurver)

Fig3A.analysis <- function(fig3df) {
    ## The analysis follows what Zack has done:
    ## 1) average blank measurement at every time point.
    ## 2) subtract average blank measurement from each time point.
    ## 3) log2 transform the data.

    ## This next line is just for debugging. keep commented out.
    ##fig3df <- labeled.growth.data

    blanks <- filter(fig3df, Name == 'Blank') %>% group_by(Experiment, Time) %>%
        summarise(blank.time.avg = mean(OD420))
    fig3df <- left_join(fig3df, blanks)
    fig3df <- mutate(fig3df, OD420.minus.blank.avg = OD420 - blank.time.avg)
    ## Add one so that the blank avg maps to zero after the log.
    fig3df <- mutate(fig3df, log2.OD420 = log2(OD420.minus.blank.avg+1))
    ## filter out blank rows.
    fig3df <- filter(fig3df,Name != 'Blank')
    ## Make an hours column with lubridate for plotting.
    fig3df <- mutate(fig3df,Hours = as.numeric(as.duration(hms(Time)))/3600)
    return(fig3df)
}

plot.Fig3A <- function(fig3df2) {
    ## This next line is just for debugging. keep commented out.
    ##fig3df2 <- final.fig3.data

    fig3A <- ggplot(fig3df2, aes(x=Hours,y=log2.OD420,color=Name)) + geom_point(size=0.5) +
        facet_grid(Experiment ~ Founder) + theme_classic() #+ ylim(-4,-0.5)
    ggsave(file.path(proj.dir,"results/figures/Fig3A.pdf"),fig3)
}

growth.data <- read.csv(file.path(proj.dir,"results/growth-data-for-fig3.csv"))
labeled.growth.data <- left_join(growth.data,pop.clone.labels,by="Name")

final.fig3.data <- Fig3A.analysis(labeled.growth.data)
plot.Fig3A(final.fig3.data)

#####################################
## use growthcurver package to fit K and r to growth data.

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

growth.curve.fit.summary <- growth.curve.fits %>% filter(Name != 'Blank') %>%
    group_by(Experiment,Name) %>% summarize(r.avg = mean(r),
                                            k.avg=mean(k),
                                            Environment=unique(Environment),
                                            Generation=unique(Generation),
                                            Founder=unique(Founder))

growth.curve.plotA <- ggplot(growth.curve.fit.summary,aes(x=k.avg,y=r.avg,color=Generation,label=Founder)) +
    geom_point() +
    geom_text_repel() +
    facet_wrap(Experiment~Founder) +
    theme_classic()

## examine how rate and yield evolved in the experiment.
ggsave(file.path(proj.dir,"results/figures/Fig3B.pdf"),growth.curve.plotA)

## now directly visualize if there's any tradeoff between DM0 and DM25.
DM0.growth.summary <- filter(growth.curve.fit.summary,Experiment=='DM0-growth') %>%
    rename(DM0.r = r.avg,DM0.k = k.avg) %>% ungroup(Experiment) %>% select(-Experiment) %>%
    mutate(Generation=as.factor(Generation))

DM25.growth.summary <- filter(growth.curve.fit.summary,Experiment=='DM25-growth') %>%
    rename(DM25.r = r.avg,DM25.k = k.avg) %>% ungroup(Experiment) %>% select(-Experiment) %>%
    mutate(Generation=as.factor(Generation))

growth.fit.summary2 <- inner_join(DM0.growth.summary,DM25.growth.summary)

R.comparison.plot <- ggplot(growth.fit.summary2,aes(x=DM25.r,y=DM0.r,color=Founder,shape=Generation)) +
    geom_point() +
    theme_classic()

K.comparison.plot <- ggplot(growth.fit.summary2,aes(x=DM25.k,y=DM0.k,color=Founder,shape=Generation)) +
    geom_point() +
    theme_classic()

ggsave(file.path(proj.dir,"results/figures/Fig3C.pdf"),R.comparison.plot)
ggsave(file.path(proj.dir,"results/figures/Fig3D.pdf"),K.comparison.plot)

################################################################################

## Figure 7: Fitness and growth effects of plasmid-borne maeA expression.
Fig7.analysis <- function(data, samplesize=6, days.competition=1,rev=FALSE) {
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

fig7.data <- read.csv("../data/rohan-formatted/DM0_Fitness2.csv",header=TRUE)
## these data come from one day competitions that Tanush ran.
res1 <- filter(fig7.data,Red.Pop=='ZDB151_with_maeA') %>% Fig7.analysis()
res2 <- filter(fig7.data,Red.Pop=='ZDB67_with_maeA') %>% Fig7.analysis(rev=TRUE)
res3 <- filter(fig7.data,Red.Pop=='ZDB152_with_maeA') %>% Fig7.analysis()
res4 <- filter(fig7.data,Red.Pop=='ZDB68_with_maeA') %>% Fig7.analysis(rev=TRUE)

## beautiful! All confints overlap with each other, showing no maeA fitness effect
## does not depend on Ara polarity or genetic background.

## let's combine data from the switched polarity competitions (since Ara marker is neutral)
d1 <- filter(fig7.data,Red.Pop=='ZDB151_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d2 <- filter(fig7.data,Red.Pop=='ZDB67_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## now switch column labels.
d2X <- rename(d2,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB151.data <- rbind(d1,d2X)

d3 <- filter(fig7.data,Red.Pop=='ZDB152_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d4 <- filter(fig7.data,Red.Pop=='ZDB68_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## switch column labels.
d4X <- rename(d4,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB152.data <- rbind(d3,d4X)

fres1 <- Fig7.analysis(ZDB151.data,samplesize=6)
fres2 <- Fig7.analysis(ZDB152.data,samplesize=6)
fig7.plot.data <- rbind(fres1,fres2)
fig7.plot.data$Strain <- c('ZDB151','ZDB152')

## Make Figure 7.
fig7.output <- "../results/figures/Fig7.pdf"

plot.Figure7 <- function (results, output.file) {
    the.plot <- ggplot(results,aes(x=Strain,y=Fitness)) +
        geom_errorbar(aes(ymin=Left,ymax=Right),width=0.1, size=1) +
        geom_line() +
        geom_point(size=2) +
        scale_y_continuous(limits=c(1.0,1.30)) +
        ylab("Fitness of maeA plasmid relative to empty plasmid") +
        theme_classic()
    ggsave(the.plot, file=output.file,width=4,height=4)
}

plot.Figure7(fig7.plot.data,fig7.output)

########################################################
## IS element analysis and visualization.

IS.insertions <- read.csv(file.path(proj.dir,
                                   "results/genome-analysis/IS_insertions.csv")) %>% arrange(genome_start)


parallel.IS.insertions <- group_by(IS.insertions,genome_start) %>%
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

## 81/213 IS insertions recur at the same locations! 38%!

IS.plot <- ggplot(IS.insertions,aes(x=genome_start,fill=IS_element)) +
    geom_histogram(bins=1000000) +
##    guides(color=FALSE) +
    theme_classic()

ggsave("../results/figures/IS-insertions.pdf",IS.plot)

## plot parallel IS insertions and their annotation.
IS.plot2 <- ggplot(parallel.IS.summary,aes(x=genome_start,
                                           y=count,
                                           color=IS_element,
                                           label=annotation)) +
    geom_point() +
    theme_classic() +
    geom_text_repel()

ggsave("../results/figures/parallel-IS-insertions.pdf",IS.plot2)

## plot genomic distribution of IS elements in LCA (including REL606).
LCA.IS.insertions <- read.csv(file.path(proj.dir,
                                   "results/genome-analysis/LCA_IS_insertions.csv")) %>% arrange(genome_start)

LCA.IS.plot <- ggplot(LCA.IS.insertions,aes(x=genome_start,fill=IS_element)) +
    geom_histogram(bins=1000) +
    theme_classic()

## compare the eCDFs of IS elements in the LCA and in the evolved populations.
## they are not strikingly similar or dissimilar.
IS.insertion.CDF.plot <- ggplot(LCA.IS.insertions,aes(x=genome_start)) +
    stat_ecdf(color='red') + theme_classic() + stat_ecdf(inherit.aes=FALSE,data=IS.insertions,mapping=aes(x=genome_start),color='black')

## p-val: 0.1033. So neither similar nor dissimilar given the data.
ks.test(IS.insertions$genome_start,LCA.IS.insertions$genome_start)
