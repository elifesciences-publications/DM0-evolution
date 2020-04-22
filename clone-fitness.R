## clone-fitness.R by Rohan Maddamsetti.
## make pretty figures for fitness of CZB151/ZDB67-derived clones
## against CZB151 and ZDB67,
## and double-check Zack's fitness calculations.

library(tidyverse)
library(cowplot)

## Fitness calculation code.
## NOTE: this function is specialized for these competitions, and will NOT
## work out of the box on other data.
fitness.analysis <- function(data, samplesize=5, days.competition=3) {
    ## by diluting stationary phase culture 1:100 on day 0, there is another
    ## factor of 100 that multiplies the day 0 plating dilution factor.
    data$D.0 <- 100*data$D.0
    data$Mr <- log((data$D.1/data$D.0)*data$d3.red*100^days.competition/data$d0.red)/days.competition
    data$Mw <- log((data$D.1/data$D.0)*data$d3.white*100^days.competition/data$d0.white)/days.competition

    ## Calculate fitness relative to the reference strain.
    stopifnot(length(unique(data$White.Clone))==1)
    if (unique(data$White.Clone) %in% c('ZDB67','ZDB1262','ZDB1233')) rev <- TRUE else rev <- FALSE
    ## Mw is denominator when have a white reference strain,
    ## Mr is denominator when have a red reference strain.
    if (rev) data$W <- data$Mr/data$Mw else data$W <- data$Mw/data$Mr
    
    my.mean <- mean(data$W, na.rm=T)
    my.sd <- sd(data$W, na.rm=T)
    ## see: https://www.cyclismo.org/tutorial/R/confidence.html#calculating-a-confidence-interval-from-a-t-distribution
    my.t.dist.error <- qt(0.975,df=samplesize-1) * my.sd/sqrt(samplesize)
    
    my.confint <- c(my.mean - my.t.dist.error, my.mean + my.t.dist.error)
        
    ##return a dataframe of the results for plotting.
    left.error <- c(my.confint[1])
    right.error <- c(my.confint[2])

    if (rev) {
        name <- unique(data$Red.Clone)
        reference <- unique(data$White.Clone)
    } else {
        name <- unique(data$White.Clone)
        reference <- unique(data$Red.Clone)
    }
    
    results <- data.frame(Name=name,
                          Reference=reference,
                          Fitness=my.mean,
                          Left=left.error,
                          Right=right.error)
    return(results)
}

plot.clone.fitness <- function(results) {

    ## colorblind-friendly palette.
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    clone.order <- c('Founder: CZB151', 'DM0-1: ZDBp871', 'DM0-2: ZDBp875',
                     'DM0+1: ZDBp889', 'DM0+2: ZDBp892',
                     'DM25-1: ZDBp910', 'DM25-2: ZDBp911',
                     'DM25+1: ZDBp916', 'DM25+2: ZDBp917',
                     'Founder: CZB152', 'DM0-3: ZDBp877', 'DM0-4: ZDBp880',
                     'DM25-3: ZDBp912', 'DM25-4: ZDBp913',
                     'Founder: CZB154', 'DM0-5: ZDBp883', 'DM0-6: ZDBp886',
                     'DM25-5: ZDBp914', 'DM25-6: ZDBp915')                   
    
    results <- results %>%
        ## replace NA values in PopulationLabel to "Founder".
        mutate(PopulationLabel=ifelse(is.na(PopulationLabel),"Founder",PopulationLabel)) %>%
        ## add the PopulationLabel to the Name for plotting.
        mutate(PlotName=paste(PopulationLabel,Name,sep=': ')) %>%
            ## order points from left to right, manually.
        mutate(PlotName=factor(PlotName,levels=clone.order))

    
    the.plot <- ggplot(results,aes(x=PlotName,y=Fitness,color=Founder)) +
        theme_classic() +
        scale_color_manual(values = cbbPalette) +
        geom_errorbar(aes(ymin=Left,ymax=Right),width=0.1, size=1) +
        geom_line() +
        geom_point(size=1.5) +
        theme(axis.text.x = element_text(size=10,angle=90,vjust=0.5)) +
        ylab("Fitness relative to ancestor") +
        geom_hline(yintercept=1,color='gray',linetype='dashed') +
        xlab("Clone") +
        guides(color=FALSE)

    return(the.plot)
}

###################################################################
## Data analysis.

## population and clone metadata.
pop.clone.metadata <- read.csv("../data/rohan-formatted/populations-and-clones.csv", header=TRUE,stringsAsFactors=FALSE)

## These data come from 3-day competitions that Jess Baxter ran.
DM0.fitness.data <- read.csv("../data/rohan-formatted/200314-evolved-clone-fitness-DM0-competitions.csv", header=TRUE,as.is=TRUE)

DM0.groups <- DM0.fitness.data %>%
    group_by(Red.Clone,White.Clone,D.0,D.1) %>%
    ## split by groups of replicates
    group_split()

DM0.results <- DM0.groups %>%
    map_dfr(.f=fitness.analysis) %>%
    ## merge with the metadata.
    left_join(pop.clone.metadata)

easy.comparison.DM0.results <- DM0.results %>% select(Name,Fitness,Left,Right,PopulationLabel,ParentClone)

## now for Jess's 3-day competitions in DM25.
DM25.fitness.data <- read.csv("../data/rohan-formatted/200316-evolved-clone-fitness-DM25-competitions.csv", header=TRUE,as.is=TRUE)

DM25.groups <- DM25.fitness.data %>%
    group_by(Red.Clone,White.Clone,D.0,D.1) %>%
    ## split by groups of replicates
    group_split()

DM25.results <- DM25.groups %>%
    map_dfr(.f=fitness.analysis) %>%
    ## merge with the metadata.
    left_join(pop.clone.metadata)

easy.comparison.DM25.results <- DM25.results %>% select(Name,Fitness,Left,Right,PopulationLabel,ParentClone)

## Make a figure of DM0 fitness for each population.
## Make a figure of DM0 fitness for each population.
Fig3A <- plot.clone.fitness(DM0.results) + ggtitle("Clone fitness measured in DM0 in a three day competition")
Fig3B <- plot.clone.fitness(DM25.results) + ggtitle("Clone fitness measured in DM25 in a three day competition")
## put these figures together to make Figure 3.
Fig3 <- plot_grid(Fig3A,Fig3B,labels=c('A','B'),ncol=1)
ggsave("../results/figures/Fig3.pdf", Fig3, height=8)

write.csv(DM0.results,file="../results/EvolvedCloneFitness-in-DM0.csv")
write.csv(DM25.results,file="../results/EvolvedCloneFitness-in-DM25.csv")
