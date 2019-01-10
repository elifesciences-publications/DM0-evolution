## make-figures.R by Rohan Maddamsetti.

library(ggplot2)
library(dplyr)
library(ggthemes)
library(viridis)
library(scales)
library(gridExtra)
library(lubridate)
library(assertthat)
library(tidyr)

pop.clone.labels <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/data/rohan-formatted/populations-and-clones.csv")

## Figure 1A is a diagram of the DM0 and DM25 experiments made in Illustrator.

## Figure 1B: make a stacked bar plot of the different kinds of mutations in each clone.

mutation.types <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/results/genome-analysis/mutation_types.csv")

fig1B.stacked <- ggplot(mutation.types,aes(x=Clone,fill=Mutation)) +
    geom_bar() +
    ##scale_fill_viridis(option="magma",discrete=TRUE) +
    coord_equal() +
    facet_wrap(~Environment,ncol=2,scales="free_x") +
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
fig1C.raw.matrix <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/results/genome-analysis/Dice-similarity-analysis/environment-comparison-mutation-matrix.csv")
fig1C.data <- gather(fig1C.raw.matrix,"Name","mutation.count",2:25) %>%
    left_join(pop.clone.labels) %>%
    select(Gene,Name,mutation.count,Environment,Population.Label) %>%
    group_by(Gene) %>% filter(sum(mutation.count)>1)

fig1C.counts <- summarize(fig1C.data,total.count=sum(mutation.count)) %>% arrange(desc(total.count))

## Now: add LTEE mutation matrix as well.
## filter on genes mutated in the DM0 and DM25 matrix data.
fig1C.ltee.matrix <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/results/LTEE-comparison-mutation-matrix.csv")
LTEE.50K.labels <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/data/rohan-formatted/LTEE-50K-clones.csv")

fig1C.ltee.data <- gather(fig1C.ltee.matrix,"Name","mutation.count",2:13) %>%
    left_join(LTEE.50K.labels) %>%
    group_by(Gene) %>% filter(Gene %in% levels(fig1C.data$Gene)) %>%
    filter(!is.na(Gene))

## Now join LTEE data to the DM0 and DM25 Cit+ data.
fig1C.data <- fig1C.data %>% bind_rows(fig1C.ltee.data) %>%
    replace_na(replace=list(Hypermutator=0))

## Possible TODO: label hypermutators in red.

## order the genes by number of mutations to get axes correct on heatmap.
fig1C.data$Gene <- factor(fig1C.data$Gene,levels=fig1C.counts$Gene)
fig1C.data$mutation.count <- factor(fig1C.data$mutation.count)
fig1C <-  ggplot(fig1C.data,aes(x=Population.Label,y=Gene,fill=mutation.count,
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
    scale_fill_manual(values = c("white", "#ffffb3", "#bebada",
                                 "#fb8072", "#80b1d3", "#fdb462"))

fig1C.output <- "../results/figures/Fig1C.pdf"
ggsave(fig1C, file=fig1C.output,width=7,height=7)


############################################################################
## Figure 2: Fitness of evolved clones in DM0 and DM25.
## NOTE: This figure still needs work, probably will use clone data instead anyway.

DM0.labels <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/data/rohan-formatted/labels-for-DM0-evolved-DM0-fitness.csv")
DM0.data <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/data/rohan-formatted/DM0-evolved-DM0-fitness.csv")
DM0.data2 <- left_join(DM0.data,DM0.labels) %>% mutate(Environment='DM0') %>%
    select(-Set,-Type)

DM25.labels <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/data/rohan-formatted/labels-for-DM0-evolved-DM25-fitness.csv")
DM25.data <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/data/rohan-formatted/DM0-evolved-DM25-fitness.csv")
DM25.data2 <- left_join(DM25.data,DM25.labels) %>% mutate(Environment='DM25')

founder.labels <- select(pop.clone.labels,Name,Founder) %>% rename(C1=Name)
## put DM0 and DM25 data into a single data frame, and add a Founder column,
## and turn Dilution into a number (Plate Dilution).
Fig2.data <- bind_rows(DM0.data2,DM25.data2) %>% left_join(founder.labels) %>%
    mutate(Plate.Dilution = sapply(as.character(Dilution),
                                   function(x) {eval(parse(text=x))}))

## filter out all plates with TMTC colonies, and only analyze one-day competition
## data, and turn Red and White columns (factors) into numeric values
## (since no more TMTC).
## WARNING: as.numeric(Red) is WRONG (gets index of the factor).
## Need to use as.numeric(as.character(Red)) instead.
filtered.Fig2.data <- filter(Fig2.data,Red != 'TMTC') %>% filter(White != 'TMTC') %>%
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


Fig2.output <- "../results/figures/Fig2.pdf"

plot.Figure2 <- function (results, output.file) {
    the.plot <- ggplot(results,aes(x=Pop,y=Fitness,color=Rev, label=Pop)) +
        geom_errorbar(aes(ymin=left.error,ymax=right.error),width=0.1, size=1) +
        geom_point(size=1) +
        facet_wrap(~Environment+Founder, scales = "free_x") +
        scale_y_continuous(breaks=seq(0, 2.5, 0.2), limits=c(0.5,2.2)) +
        geom_hline(yintercept=1,linetype="dashed") +
        theme_tufte(base_family='Helvetica') +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        theme(axis.title.x=element_blank()) +
        theme(axis.title.y=element_text(hjust=1))
        theme(axis.ticks = element_blank(), axis.text.x = element_blank())

    ggsave(the.plot, file=output.file)
}

plot.Figure2(Fig2.results,"/Users/Rohandinho/Desktop/test.pdf")


## print out calculations to compare to Zack's calculations.
## write.csv(Fig2.competition.data2,"/Users/Rohandinho/Desktop/Fig2_fitness.csv")
## write.csv(Fig2.results,"/Users/Rohandinho/Desktop/Fig2_fitness_results.csv")

################################################################################
## Figure 3: Growth curves of evolved populations in DM0 and in DM25.

Fig3.analysis <- function(fig3df) {
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

plot.Fig3 <- function(fig3df2) {
    ## This next line is just for debugging. keep commented out.
    ##fig3df2 <- final.fig3.data

    fig3 <- ggplot(fig3df2, aes(x=Hours,y=log2.OD420,color=Name)) + geom_point(size=0.5) +
        facet_grid(Experiment ~ Founder) + theme_classic() #+ ylim(-4,-0.5)
    ggsave("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/results/figures/Fig3.pdf",fig3)
}

growth.data <- read.csv("/Users/Rohandinho/Dropbox (HMS)/DM0-evolution/results/growth-data-for-fig3.csv")
labeled.growth.data <- left_join(growth.data,pop.clone.labels,by="Name")

final.fig3.data <- Fig3.analysis(labeled.growth.data)
plot.Fig3(final.fig3.data)

################################################################################

## Figure 7: Fitness and growth effects of plasmid-borne maeA expression.
Fig7.analysis <- function(data, samplesize=6, days.of.competition=1,rev=FALSE) {
    ## by diluting stationary phase culture 1:100 on day 0, there is another
    ## factor of 100 that multiplies the day 0 plating dilution factor.
    data$D.0 <- 100*data$D.0
    data$Mr <- log2((data$D.1/data$D.0)*data$Red.1*100^days.of.competition/data$Red.0)/days.of.competition
    data$Mw <- log2((data$D.1/data$D.0)*data$White.1*100^days.of.competition/data$White.0)/days.of.competition
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
