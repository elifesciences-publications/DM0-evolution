## population-fitness.R by Rohan Maddamsetti.
## make pretty figures for fitness of DM0 project populations
## against CZB151 and ZDB67,
## and double-check Zack's fitness calculations.

library(tidyverse)
library(cowplot)

## Fitness calculation code.
## NOTE: this function is specialized for these competitions, and will NOT
## work out of the box on other data.
fitness.analysis <- function(d, days.competition=1) {
    data <- d ## this is so we modify a copy of the input.
    ## by diluting stationary phase culture 1:100 on day 0, there is another
    ## factor of 100 that multiplies the day 0 plating dilution factor.

    data$D.0 <- 100 * d$D.0 ## don't modify self, esp.  when testing interactively.
    ## now start fitness calculations.
    data$Mr <- log((data$D.1/data$D.0)*data$d1.red*100^days.competition/data$d0.red)/days.competition
    data$Mw <- log((data$D.1/data$D.0)*data$d1.white*100^days.competition/data$d0.white)/days.competition

    ## Calculate fitness relative to the reference strain.
    stopifnot(length(unique(data$White.Pop))==1)
    if (unique(data$White.Pop)=='ZDB67') rev <- TRUE else rev <- FALSE
    ## Mw is denominator when have a white reference strain,
    ## Mr is denominator when have a red reference strain.
    if (rev) data$W <- data$Mr/data$Mw else data$W <- data$Mw/data$Mr

    ## for these competitions, we need to deal with some zero counts.
    ## log(0) = -Inf.
    ## if we're dividing by -Inf, then W = 0.
    ## if we're dividing -Inf on the other hand, 
    ## then W = -Inf.

    ## For the three cases under consideration in these DM0 competitions,
    ## W = 0 --> W = 1, and 
    ## W = -Inf --> W = 0.
    ## This is a case-specific hack, so don't generalize blindly to other data!
    data$W <- sapply(data$W, function(x) ifelse(x == 0, 1, x))
    data$W <- sapply(data$W, function(x) ifelse(x == -Inf, 0, x))
    
    ## This logic comes from Rich's comments in email correspondence:
    ## Like two other cases (ZDBp806 vs CZB151 and ZDBp808 vs CZB151),
    ## the evolved competitor (but not the common competitor) had 0 colonies on Day 1
    ## (but not day 0) in 2 of the 3 reps.  So I would call those cases a relative
    ## fitness of 0.  (It would be infinity if it was reversed,
    ## but in all cases it was the evolved population that failed.) 

    ##In that case, for ZDBp754 I compute mean 0.2585, lower bound 0.0000,
    ## and upper bound 1.3707.

    ##In these 6 cases (3 competing pairs, each with 2 failed cases), I suggest
    ## replacing the common competitor's growth rate with the actual value;
    ## leave the NA for those where the evolved population had 0 colonies on D.1;
    ## and hard code fitness in those cases as 1.  Then you can compute the averages,
    ## etc, as you've done from the 3 reps in every case.
    
    samplesize <- nrow(data)

    my.mean <- mean(data$W, na.rm=FALSE) ## IMPORTANT: return NA if any NA values.
    my.sd <- sd(data$W, na.rm=FALSE)
    ## see: https://www.cyclismo.org/tutorial/R/confidence.html#calculating-a-confidence-interval-from-a-t-distribution
    my.t.dist.error <- qt(0.975,df=samplesize-1) * my.sd/sqrt(samplesize)
    
    my.confint <- c(my.mean - my.t.dist.error, my.mean + my.t.dist.error)
        
    ##return a dataframe of the results for plotting.
    left.error <- c(my.confint[1])
    right.error <- c(my.confint[2])

    if (rev) {
        name <- unique(data$Red.Pop)
        reference <- unique(data$White.Pop)
    } else {
        name <- unique(data$White.Pop)
        reference <- unique(data$Red.Pop)
    }

    ## If mean is negative or undefined, then set to equal 0,
    ## and hard code the error bars to extend from 0 to 4.
    ## Or if one or other limit is within the range,
    ## then hard code it to that value and other to the limit up or down.
    my.mean <- ifelse(is.na(my.mean) | (my.mean <= 0) | (my.mean >= 4), 0, my.mean)
    left.error <- ifelse(is.na(left.error) | (left.error <= 0) | (left.error >= my.mean), 0, left.error)
    right.error <- ifelse(is.na(right.error) | (right.error >= 4), 4, right.error)
    
    results <- data.frame(Treatment=unique(data$Treatment),
                          Name=name,
                          Reference=reference,
                          Fitness=my.mean,
                          Left=left.error,
                          Right=right.error)
    return(results)
}

plot.pop.fitness <- function(results) {
    ## colorblind-friendly palette.
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    pop.order <- c('CZB151', 'DM0-1','DM0-2','DM0+1','DM0+2',
                   'DM25-1','DM25-2','DM25+1','DM25+2',
                   'CZB152','DM0-3','DM0-4','DM25-3','DM25-4',
                   'CZB154','DM0-5','DM0-6','DM25-5','DM25-6',
                   'ZDB68','DM0+3','DM0+4','DM25+3','DM25+4',
                   'ZDB69','DM0+5','DM0+6','DM25+5','DM25+6')
                   
    results <- results %>%
        ## give the ancestor clones a default label.
        mutate(PopulationLabel=ifelse(is.na(PopulationLabel),
                                      Name,PopulationLabel)) %>%
        ## order points from left to right, manually.
        mutate(PopulationLabel=factor(PopulationLabel,levels=pop.order))
    
    the.plot <- ggplot(results,aes(x=PopulationLabel,y=Fitness,color=Founder)) +
        theme_classic() +
        scale_color_manual(values = cbbPalette) +
        geom_errorbar(aes(ymin=Left,ymax=Right),width=0.1, size=1) +
        geom_line() +
        geom_point(size=1.5) +
        theme(axis.text.x = element_text(size=10,angle=90,vjust=0.5)) +
        ylab("Fitness relative to CZB151/ZDB67") +
        xlab("Population") +
        geom_hline(yintercept=1,color='gray',linetype='dashed') +
        geom_hline(yintercept=1.5,color='gray',linetype='dashed') +
        guides(color=FALSE)

    return(the.plot)
}

###################################################################
## Data analysis.

## population and clone metadata.
pop.clone.metadata <- read.csv("../data/rohan-formatted/populations-and-clones.csv", header=TRUE,stringsAsFactors=FALSE)

#### DM0 competitions.
## these data come from one day competitions in DM0 that Zack ran.
DM0.fitness.data <- read.csv("../data/rohan-formatted/200316-evolved-pop-fitness-DM0-competitions.csv", header=TRUE,as.is=TRUE)

DM0.groups <- DM0.fitness.data %>%
    ## drop the intermediate calculations
    select(-Red.M,-White.M,-Fitness) %>%
    group_by(Treatment,Red.Pop,White.Pop,D.0,D.1) %>%
    ## split by groups of replicates
    group_split()

DM0.results <- DM0.groups %>%
    map_dfr(.f=fitness.analysis) %>%
    ## merge with the metadata.
    left_join(pop.clone.metadata)

easy.comparison.DM0.results <- DM0.results %>% select(Name,Fitness,Left,Right,PopulationLabel,ParentClone)

## DM25 competitions
## these data come from one day competitions in DM25 that Zack ran.
DM25.fitness.data <- read.csv("../data/rohan-formatted/200316-evolved-pop-fitness-DM25-competitions.csv", header=TRUE,as.is=TRUE)

DM25.groups <- DM25.fitness.data %>%
    ## drop the intermediate calculations
    select(-Red.M,-White.M,-Fitness) %>%
    group_by(Treatment,Red.Pop,White.Pop,D.0,D.1) %>%
    ## split by groups of replicates
    group_split()

DM25.results <- DM25.groups %>%
    map_dfr(.f=fitness.analysis) %>%
    ## merge with the metadata.
    left_join(pop.clone.metadata)

easy.comparison.DM25.results <- DM25.results %>% select(Name,Fitness,Left,Right,PopulationLabel,ParentClone)

## Make a figure of DM0 fitness for each population.
Fig2A <- plot.pop.fitness(DM0.results) +
    ## remove problematic fitness estimates by setting bounds from 0 to 4.
    ylim(0,4) +
    ggtitle("Population fitness measured in DM0 in a one day competition")
## Make a figure of DM25 fitness for each population.
Fig2B <- plot.pop.fitness(DM25.results) +
    ## set bounds from 0 to 2.
    ylim(0,2) +
    ggtitle("Population fitness measured in DM25 in a one day competition")
## put these figures together to make Figure 2.
Fig2 <- plot_grid(Fig2A,Fig2B,labels=c('A','B'),ncol=1)
ggsave("../results/figures/evolved-pop-fitness.pdf", Fig2, height=7)

write.csv(DM0.results,file="../results/EvolvedPopFitness-in-DM0.csv")
write.csv(DM25.results,file="../results/EvolvedPopFitness-in-DM25.csv")
