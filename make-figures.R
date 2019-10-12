## make-figures.R by Rohan Maddamsetti.

## go through and remove unneeded dependencies once this script is more polished.

library(boot) ## for professional bootstrapping, rather than rolling my own.
library(tidyverse)
library(cowplot)
library(zoo) ## for rollapply function
library(data.table)
library(ggrepel)
library(ggthemes)
library(viridis)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(lubridate)
library(assertthat)
library(growthcurver) ## use growthcurver package to fit r from growth data.
## Cite the growthcurver package. https://www.ncbi.nlm.nih.gov/pubmed/27094401

## IMPORTANT TODO: CHECK HOW the second argument matters in this line of code:
###return(c(mean(x[ind]),var(x[ind])/length(ind)))

calc.bootstrap.conf.int <- function(vec) {
  ## bootstrap confidence intervals around the mean.
  ## Use the boot package to calculate fancy BCA intervals.
  
  mean.boot <- function(x,ind) {
    ## define this type of mean-calculating function to pass to the boot function.
    ## from: https://web.as.uky.edu/statistics/users/pbreheny/621/F12/notes/9-18.pdf
    return(c(mean(x[ind]),var(x[ind])/length(ind)))
  }

  Nbootstraps <- 10000

  ## remove any NA values from vec.
  vec <- vec[!is.na(vec)]

  if (length(vec)) {    
    out <- boot(vec,mean.boot,Nbootstraps)
    ## handle bad inputs in bootstrapping confidence intervals
    ci.result <- tryCatch(boot.ci(out,type="bca"), error= function(c) return(NULL))
    if (is.null(ci.result)) {
      Left <- NA
      Right <- NA
    } else {
      Left <- ci.result$bca[4]
      Right <- ci.result$bca[5]
    }
  } else { ## vec only contained NA values.
    Left <- NA
    Right <- NA
  }   
  final <- c(Left, Right)
  return(final)
}

###################

## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("DM0-evolution","src")))
projdir <- file.path("..")

## For the sake of clarity in exposition, omit the oddball Cit- ZDBp874 clone from
## the analyses. Do keep its growth curve in the supplement, as it was used to calibrate
## those analyses, but omit it from the calculations in the main figures.

pop.clone.labels <- read.csv(
  file.path(projdir,
            "data/rohan-formatted/populations-and-clones.csv"),
  stringsAsFactors=FALSE)

evolved.mutations <- read.csv(
  file.path(projdir,
            "results/genome-analysis/evolved_mutations.csv"),
  stringsAsFactors=FALSE) %>%
mutate(Mutation=as.factor(Mutation)) %>%
## exclude oddball Cit- clone ZDBp874.
filter(Clone != 'ZDBp874')

## Notice parallel mutations at base-pair level in citrate synthase!
filter(evolved.mutations,Gene=='gltA') %>% arrange(Position)
## but remember that fine-tuning was already reported by Erik Quandt.

## any other parallel dN at base-pair level?
parallel.dN <- evolved.mutations %>% filter(Mutation=='nonsynonymous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)

parallel.dN.Table <- filter(evolved.mutations, Position %in% parallel.dN$Position) %>%
arrange(Position)

## 3/5 parallel bp mutations are in citrate synthase, gltA!
## 735765, 735797, 735941.
## These are M172I, A162T, I114F, respectively.
## A162V is reported in Quandt et al. (2015) in eLife as affecting NADH-binding.
## The other two have not been reported before.

## From methods section of Quandt et al. (2015):
## Dimeric models of ligand-free and NADH-bound E. coli citrate synthase (PDB: 1NXE and 1NXG) (Maurus et al., 2003)
## were prepared for analysis by reverting alanine-383 to phenylalanine and processing with the Structure
## Preparation application within MOE (Molecular Operating Environment).

## the other parallel bp mutations are:
## 2209801 is in atoS, 2630053 is in ygaF.

## atoS annotation from Ecocyc:
## AtoS is the sensor histidine kinase of the AtoS/AtoC two-component signal transduction pathway which is best characterised by its induction of the ato operon for metabolism of short chain fatty acids in response to the presence of acetoacetate. AtoS is a homo-dimeric transmembrane protein consisting of an amino-terminal periplasmic sensing domain coupled to a carboxy-terminal cytoplasmic kinase domain [Lioliou05, Filippou08]. AtoS autophosphorylates on a conserved histidine residue by trans-phosphorylation between the monomers of the homodimer [Filippou08]. AtoS tranfers a phosphoryl group to the cytoplasmic response regulator AtoC which controls transcriptional expression of the operon (atoDAEB) involved in short-chain fatty acid catabolism [Jenkins87a, Lioliou04, Lioliou05].

## ygaF annotation from Ecocyc:
## L-2-hydroxyglutarate dehydrogenase (LhgD) is an electron transport chain-coupled dehydrogenase that feeds electrons from the reaction into the membrane quinone pool [Knorr18]. LhgD contains an FAD cofactor which is not covalently attached, and whose reduction potential is relatively high at -25 mV [Kalliri08].
## LhgD was initially thought to be an oxidase, i.e. using molecular oxygen for oxidation of L-2-hydroxyglutarate and producing hydrogen peroxide [Kalliri08].

## LhgD is associated with the cytoplasmic membrane [Zhang07], and its activity is only found in the membrane fraction [Knorr18].

## lhgD is part of an operon whose expression is induced upon carbon starvation and in stationary phase [Marschall98, Becker01, Germer01, Metzner04].

non.MOB.parallel <- evolved.mutations %>% filter(Mutation!='MOB') %>% filter(Mutation!= 'nonsynony
mous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)

## what about other classes of mutations?
Table3 <- filter(evolved.mutations, Position %in% non.MOB.parallel$Position) %>% arrange(Position)

ltee.mutations <- read.csv(
  file.path(projdir,
            "data/rohan-formatted/nature18959-s4.csv"),
  stringsAsFactors=FALSE
  )

## by inspecting evolved.mutations in conjunction with amplified genes, we see that the
## dctA amplification anti-correlates with the promoter mutation.
## Note that the dctA promoter mutation is in CZB151 and CZB154.

dctA.AMPs <- read.csv(file.path(projdir,"results/amplified_genes.csv"),
                      stringsAsFactors=FALSE) %>%
filter(gene=='dctA') %>% dplyr::select(Genome) %>% distinct() %>%
transmute(Name=Genome) %>% left_join(pop.clone.labels)

promoter.mutant <- filter(evolved.mutations, Gene=='dctA/yhjK') %>%
transmute(Name=Clone) %>% left_join(pop.clone.labels) %>%
full_join(filter(pop.clone.labels,Founder %in% c('CZB151','CZB154') & Generation==2500 & SampleType=='Clone' & Sequenced == 1 & Name != 'ZDBp874'))

## out of 24 DM0- and DM25-evolved genomes,
dctA.AMP.promoter.mut.intersection <- inner_join(dctA.AMPs,promoter.mutant)
## 1 has both promoter mutant and dctA AMP.
dctA.AMP.no.promoter.mut <- anti_join(dctA.AMPs,promoter.mutant)
## 5 have dctA AMP but lack the promoter mutation.
promoter.mut.no.dctA.AMP <- anti_join(promoter.mutant,dctA.AMPs)
## 16 have the promoter mutation, and not the dctA AMP.
neither.dctA.AMP.nor.promoter <- filter(pop.clone.labels,Founder == 'CZB152' & Generation==2500 & SampleType=='Clone') %>%
anti_join(dctA.AMPs) %>%
anti_join(promoter.mutant)
## 2 do not have either mutation.

## Fisher's exact test: p = 0.0027.
fisher.test(matrix(c(1,5,16,2),2))

###############################################
## Figure 1.
## Figure 1A is a schematic of the experimental design, made in Illustrator.
## Figure 1B is the phylogenetic tree constructed using the iTOL webserver.
###############################################
## Figure 2.
## Figure 2A: make a stacked bar plot of the different kinds of mutations in each clone.
## Figure 2B: comparison to non-mutator LTEE 5000gen A clones.

plot.evolved.mutations.stackbar <- function(evolved.mutations,DM0=TRUE) {
  if (DM0) {
    muts <- filter(evolved.mutations,Environment=='DM0')
    my.title <- 'DM0 2,500 generations'
  } else {
    muts <- filter(evolved.mutations,Environment=='DM25')
    my.title <- 'DM25 2,500 generations'
  }
  
  fig <- ggplot(muts,aes(x=Clone,fill=Mutation)) +
  geom_bar() +
  ylab("Count") +
  ylim(c(0,30)) +
  scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
  theme_tufte(base_family='Helvetica') +
  theme(axis.ticks=element_blank(),
        axis.text.x=element_text(size=12,angle=45,hjust=1),
        axis.text.y=element_text(size=12),
        panel.border=element_blank(),
        plot.title=element_text(hjust=0, size = 12, face = "bold"),
        strip.text=element_text(hjust=0,size=12),
        panel.spacing.x=unit(1, "cm"),
        panel.spacing.y=unit(0.5, "cm")) +
        ggtitle(my.title) +
        guides(fill=FALSE)
  return(fig)
}

fig2A <- plot.evolved.mutations.stackbar(evolved.mutations,DM0=TRUE)
fig2B <- plot.evolved.mutations.stackbar(evolved.mutations,DM0=FALSE)

ltee.5K.clones <- ltee.mutations %>% filter(Generation==5000) %>%
## take odd rows to get A clones.
filter(row_number() %% 2 == 1) %>%
gather(key='Mutation',value='Count',
       IS.Element.Insertions,Small.Indels,Large.Deletions,
       Large.Duplications,Synonymous.Base.Substitutions,
       Nonsynonymous.Base.Substitutions,Intergenic.Base.Substitutions,
       Pseudogene.Base.Substitutions, Nonsense.Base.Substitutions) %>%
mutate(Mutation=recode(Mutation,
         IS.Element.Insertions='MOB',
         Small.Indels='DEL',
         Large.Deletions='DEL',
         Large.Duplications='INS',
         Synonymous.Base.Substitutions='synonymous',
         Nonsynonymous.Base.Substitutions = 'nonsynonymous',
         Intergenic.Base.Substitutions = 'intergenic',
         Pseudogene.Base.Substitutions = 'pseudogene',
         Nonsense.Base.Substitutions = 'nonsense')) %>%
mutate(Mutation=as.factor(Mutation)) %>%
filter(!(Population %in% c('Ara-2','Ara+3')))

## make sure colors are the same across plots by setting levels.
levels(ltee.5K.clones$Mutation) <- levels(evolved.mutations$Mutation)

fig2C <- ggplot(ltee.5K.clones,aes(x=Population,y=Count,fill=Mutation)) +
geom_bar(stat='identity') +
ylim(c(0,30)) +
xlab('Clone') +
scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +
theme_tufte(base_family='Helvetica') +
theme(axis.ticks=element_blank(),
      axis.text.x=element_text(size=12,angle=45,hjust=1),
      axis.text.y=element_text(size=12),
      panel.border=element_blank(),
      plot.title=element_text(hjust=0, size = 12, face = "bold"),
      strip.text=element_text(hjust=0, size=12),
      panel.spacing.x=unit(1, "cm"),
      panel.spacing.y=unit(0.5, "cm"),
      legend.title=element_text(size=12),
      legend.title.align=1,
      legend.text=element_text(size=12)) +
ggtitle('LTEE 5,000 generations')

## pull the legend from fig2C.
legend <- cowplot::get_legend(fig2C)

## remove the legend from fig2C.
fig2C <- fig2C + guides(fill=FALSE)

fig2 <- plot_grid(fig2A,fig2B,fig2C,legend,
                  labels=c('A','B','C',''),
                  nrow=2)

fig2.output <- "../results/figures/Fig2.pdf"
ggsave(fig2, file=fig2.output,width=7,height=7)

##############################################################
## Figure 3: growth curve analysis, showing evolution of
## growth rate and time lags in DM0, and associated growth
## parameters when grown in DM25.

######################################################################
## let's examine REL606 growth in DM25 to estimate OD420.
## We need these data to provide an empirical basis for the intervals
## chosen to estimate r.glucose and r.citrate in later code.

## IMPORTANT NOTE: Double-check the provenance of these data with Zack.
## This file indicates that this is 96 hour data, but appears to go
## for 166 hours (one week)??
REL606.DM25.growth.data <- read.csv(
  file.path(projdir, "data/rohan-formatted/REL606-DM25-96-hours.csv"),
  stringsAsFactors=FALSE) %>%
mutate(Hours = as.numeric(as.duration(hms(Time)))/3600) %>%
rowwise() %>%
## average the replicate blanks.
mutate(Blank = lift_vd(mean)(DM25.1,DM25.2,DM25.3,DM25.4,DM25.5,DM25.6)) %>%
dplyr::select(-DM25.1,-DM25.2,-DM25.3,-DM25.4,-DM25.5,-DM25.6) %>%
## subtract Blank from well readings.
mutate(REL606.1 = REL606.1 - Blank) %>%
mutate(REL606.2 = REL606.2 - Blank) %>%
mutate(REL606.3 = REL606.3 - Blank) %>%
mutate(REL606.4 = REL606.4 - Blank) %>%
mutate(REL606.5 = REL606.5 - Blank) %>%
## zero out negative values.
mutate_at(vars(REL606.1:REL606.5),function(x)ifelse(x<0,0,x)) %>%
dplyr::select(-Blank) %>%
gather("Replicate","OD420",-Hours,-Time) %>%
mutate(Name='REL606') %>%
mutate(log.OD420=log(OD420)) %>%
## Note: the blanks get messed up several days in. Contamination??
filter(Hours<=24)

plot.REL606.DM25.growth <- function(REL606.df,logscale=FALSE) {
  
  glu.interval.high <- 0.02
  glu.interval.low <- 0.01
  glu.ceiling <- 0.05
  
  if (logscale) {
    fig <- ggplot(REL606.df,
                  aes(x=Hours,
                      y=log.OD420)) +
                      ylab(expression(log~(OD[420]))) +
                      geom_hline(yintercept=log(glu.ceiling),linetype='dashed',color='red') +
                      geom_hline(yintercept=log(glu.interval.high),linetype='dashed',color='black') +
                      geom_hline(yintercept=log(glu.interval.low),linetype='dashed',color='black')

  } else {
    fig <- ggplot(REL606.df,
                  aes(x=Hours,
                      y=OD420)) +
                      ylab(expression(OD[420])) +
                      geom_hline(yintercept=glu.ceiling,linetype='dashed',color='red') +
                      geom_hline(yintercept=glu.interval.high,linetype='dashed',color='black') +
                      geom_hline(yintercept=glu.interval.low,linetype='dashed',color='black')
    
  }
  
  fig <- fig +
  xlab('Time (h)') +
  geom_point(size=0.1, alpha=0.1) +
  theme_classic() +
  guides(color=FALSE) +
  theme(plot.title = element_text(size = 12, face = "bold"))
  
  return(fig)
}

REL606.plot <- plot.REL606.DM25.growth(REL606.DM25.growth.data, logscale=FALSE)
log.REL606.plot <- plot.REL606.DM25.growth(filter(REL606.DM25.growth.data,Hours<24), logscale=TRUE)

S2Fig <- plot_grid(REL606.plot,log.REL606.plot,labels=c('A','B'),ncol=1)

save_plot(file.path(projdir,"results/figures/S2Fig.pdf"),S2Fig,base_height=4,base_width=3)

########################################
prep.growth.df <- function(growth.df) {
  ## The analysis follows what Zack has done:
  ## 1) average blank measurement at every time point.
  ## 2) subtract average blank measurement from each time point.
  ## 3) log transform the data.
  blanks <- filter(growth.df, Name == 'Blank') %>% group_by(Experiment, Time) %>%
  summarise(blank.time.avg = mean(OD420))
  growth.df <- left_join(growth.df, blanks) %>%
  mutate(rawOD420 = OD420) %>%
  mutate(OD420 = rawOD420 - blank.time.avg) %>%
  mutate(log.OD420 = log(OD420)) %>%
  ## filter out blank rows.
  filter(Name != 'Blank') %>%
  ## Make an hours column with lubridate for plotting.
  mutate(Hours = as.numeric(as.duration(hms(Time)))/3600) %>%
  ungroup()
  return(growth.df)
}

plot.single.growthcurve <- function(plot.growth.data, ev.name, fdr, logscale=FALSE) {
  
  if (logscale) {
    fig <- ggplot(plot.growth.data,
                  aes(x=Hours,
                      y=log.OD420,
                      color=Name)) +
                      ylim(c(-6,0)) +
                      ylab(expression(log~(OD[420]))) +
                      geom_hline(yintercept=log(0.05),linetype='dashed',color='red',size=0.1) +
                      geom_hline(yintercept=log(0.1),linetype='dashed',color='red',size=0.1) +
                      ## clever trick to only plot glucose interval on DM25 curves.
                      geom_hline(data=data.frame(yint=log(0.01),Experiment='DM25'), aes(yintercept=yint),linetype='dashed',color='black',size=0.1) +
                      geom_hline(data=data.frame(yint=log(0.02),Experiment='DM25'), aes(yintercept=yint),linetype='dashed',color='black',size=0.1)
    
  } else {
    fig <- ggplot(plot.growth.data,
                  aes(x=Hours,
                      y=OD420,
                      color=Name)) +
                      ylim(c(0,0.6)) +
                      ylab(expression(OD[420]))
  }
  
  ## title for the plot.
  ev.pop <- unique(filter(plot.growth.data,Name==ev.name)$PopulationLabel)
  sample.type <- unique(filter(plot.growth.data,Name==ev.name)$SampleType)
  if (sample.type == 'Clone') {
    title.string <- paste0(ev.pop, ': ', ev.name)
  } else {
    title.string <- paste0(ev.pop, ': population')
  }
  
  ## make Cit- strain stand out.
  if (ev.name == 'ZDBp874') {
    ev.color <- "#E69F00"
  } else {
    ev.color <- "#009E73"
  }
  
  fig <- fig +
  theme_classic() +
  xlab('Time (h)') +
  geom_point(data=filter(plot.growth.data, Name==fdr),
             size=0.1, color="#999999", alpha=0.1) +
  geom_point(data=filter(plot.growth.data,Name==ev.name),
             size=0.1, color=ev.color, alpha=0.1) +
             facet_grid(Experiment ~ Founder) +
             guides(color=FALSE) +
             ggtitle(title.string) +
             theme(plot.title = element_text(size = 12, face = "bold")) +
             theme(strip.text.x = element_blank())
  
  return(fig)
}

plot.growthcurve.figure <- function(growth.df,logscale=FALSE) {
  
  list.of.plots <- vector("list", 12) ## we will have 12 plots.
  i <- 1 ## index for list.
  
  ## recode Experiment column for plotting clarity.
  growth.df <- mutate(growth.df,
                      Experiment=recode_factor(Experiment,
                        `DM0-growth`='DM0',
                        `DM25-growth`='DM25'))
  
  for (fdr in sort(unique(growth.df$Founder))) {
    
    growth.by.founder <- filter(growth.df,Founder == fdr)
    evolved.names <- sort(unique(filter(growth.by.founder, Name != fdr)$Name))
    
    for (ev.name in evolved.names) {
      plot.growth.data <- filter(growth.by.founder, Name %in% c(fdr, ev.name))
      fig <- plot.single.growthcurve(plot.growth.data, ev.name, fdr, logscale)
      ## append fig to the list of plots.
      list.of.plots[[i]] <- fig
      i <- i + 1
    }
  }
  ## use cowplot to plot.
  my.plot <- plot_grid(plotlist=list.of.plots, ncol=4)
  return(my.plot)
}

calc.growth.rates <- function(final.growth.df) {
  ## This is a simpler growth rate calculation code.
  ## use two different intervals for estimating r.glucose and r.citrate,
  ## based on the REL606 growth curves in DM25, the Cit- ZDBp874 growth curve,
  ## and the Cit+ clone growth curves.
  ## This should hopefully provide more reliable estimates.
  
  summarize.well.growth <- function(well.df) {
    
    glu.min.index <- max(min(which(well.df$OD420 > 0.01)), 1)
    glu.max.index <- min(min(which(well.df$OD420 > 0.02)), nrow(well.df))
    
    cit.min.index <- max(min(which(well.df$OD420 > 0.05)), 1)
    cit.max.index <- min(min(which(well.df$OD420 > 0.1)), nrow(well.df))
    ## to compare time lags, measure point when OD420 hits 0.01.
    t.OD.hit.glu.min <- well.df[glu.min.index,]$Hours
    
    if ((glu.max.index <= nrow(well.df)) & (glu.min.index < glu.max.index)) {
      glucose.data <- well.df[glu.min.index:glu.max.index,]
    } else {
      glucose.data <- NULL
    }
    
    if ((cit.max.index <= nrow(well.df)) & (cit.min.index<cit.max.index)) {
      citrate.data <- well.df[cit.min.index:cit.max.index,]
    } else { ## if OD never hits 0.1, then don't measure r.citrate.
      citrate.data <- NULL
    }
    
    ## can't measure growth on glucose if there's no glucose!
    if (unique(well.df$Experiment) == 'DM0-growth'){
      glucose.data <- NULL
    }
    
    get.max.growth.rate <- function(subdata,window=TRUE) {
      ## find max growth rate for glucose or citrate data.
      
      if (is.null(subdata)) return(NA)
      rel.subdata <- dplyr::select(subdata,c('log.OD420','Hours')) %>%
      ##Remove rows that contain nans
      drop_na()
      
      w <- 5 ## require at least 5 datapoints to calculate growth rates.
      if (nrow(rel.subdata) < w) return(NA)
      
      calc.slope <- function(subdf) {
        coef(lm(log.OD420 ~ Hours, data = as.data.frame(subdf)))[2]
      }
      
      if (window) { ## Find the max slope in the interval.
        g.rate <- max(rollapply(zoo(rel.subdata), width = w,
                                FUN = calc.slope,
                                by.column = FALSE, align = "right"))
      }
      else { ## Calculate slope on the whole interval
        g.rate <- calc.slope(rel.subdata)
      }
      return(g.rate)
    }
    
    summarized.df <- data.frame(Name = unique(well.df$Name),
                                Well = unique(well.df$Well),
                                Experiment = unique(well.df$Experiment),
                                SampleType = unique(well.df$SampleType),
                                Generation = unique(well.df$Generation),
                                Environment = unique(well.df$Environment),
                                ParentClone = unique(well.df$ParentClone),
                                Founder = unique(well.df$Founder),
                                r.citrate = get.max.growth.rate(citrate.data),
                                r.glucose = get.max.growth.rate(glucose.data),
                                t.lag = t.OD.hit.glu.min,
                                stringsAsFactors=FALSE)

    return(summarized.df)
  }

  growth.rate.summary <- final.growth.df %>%
  droplevels() %>% ## otherwise runs and dies on empty subsets!
  ## Split on both growth conditions (DM25 or DM0) and Well.
  split(list(.$Experiment,.$Well)) %>%
  map_dfr(.f=summarize.well.growth)
  return(growth.rate.summary)
}

filter.growth.data.on.analysis.domain <- function(growth.rate.df) {
  ## The goal here is to plot the domain of the curves
  ## PER EXPERIMENTAL CONDITIONS AND PER WELL
  ## that is used to calculate growth rates.
  ## We want to plot the data that goes into my rate analysis code.
  ## Also plot intervals used to calculate rates.
  ## This code is EXTREMELY useful for debugging.

  filter.well.data <- function(well.df) {
    ## minimum OD420 > 0.01.
    ## maximum OD420 < 0.1.
    
    ## find first index where OD420 > 0.01.
    min.index <- max(min(which(well.df$OD420 > 0.01)),1)
    ## find first index where OD420 > 0.1.
    max.index <- min(min(which(well.df$OD420 > 0.1)),nrow(well.df))

    ## if OD never hits 0.1, then don't measure.
    if ((max.index > nrow(well.df)) | (min.index >= max.index)) {
      return(NULL)
    }

    return(well.df[min.index:max.index,])
  }

  filtered.df <- growth.rate.df %>%
  droplevels() %>% ## don't map/reduce empty subsets
  split(list(.$Experiment,.$Well)) %>%
  map_dfr(.f=filter.well.data)

  return(filtered.df)
}

map.reduce.growth.curve <- function(final.growth.df) {
  ## This function splits the growth plate data by Well and Experiment,
  ## Where Experiment is either 'DM0-growth' or 'DM25-growth'.

  ## returns one-row dataframe of growth rate fit and metadata from gdata.
  get.growth.curve.fit <- function(well.df) {

    ## since stationary phase is so wacky,
    ## fit data up to 6 hours past OD420 == 0.1,
    ## or to the end of the timecourse, whichever is first.
    cit.max.index <- min(min(which(well.df$OD420 > 0.1)), nrow(well.df))
    t.OD.hit.cit.max <- well.df[cit.max.index,]$Hours
    trim.time <- min(t.OD.hit.cit.max + 6, max(well.df$Hours))

    gcfit <- SummarizeGrowth(well.df$Hours,well.df$OD420,t_trim = trim.time)
    gcfitv <- gcfit$val
    my.df = data.frame(Name = unique(well.df$Name),
      Well = unique(well.df$Well),
      Experiment = unique(well.df$Experiment),
      SampleType = unique(well.df$SampleType),
      Generation = unique(well.df$Generation),
      Population = unique(well.df$Population),
      PopulationLabel = unique(well.df$PopulationLabel),
      AraStatus = unique(well.df$AraStatus),
      ParentClone = unique(well.df$ParentClone),
      Founder = unique(well.df$Founder),
      Environment = unique(well.df$Environment),
      Sequenced = unique(well.df$Sequenced),
      r = gcfitv$r,
      t_mid = gcfitv$t_mid,
      stringsAsFactors=FALSE)
    return(my.df)
  }

  gcfit.df <- final.growth.df %>%
  droplevels() %>% ## so that doesn't run and die on empty subsets
  ## Split on both growth conditions (DM25 or DM0) and Well.
  split(list(.$Experiment,.$Well)) %>%
  map_dfr(.f=get.growth.curve.fit)

  return(gcfit.df)
}

bootstrap.my.growth.confints <- function(growth) {
  ## Calculate and plot growth parameter estimates for each population and clone.
  ## I want to distinguish between noise within estimates versus differences between then.

  bootstrap.my.growth.parameters <- function(df) {
    ## my analysis fits.
    DM0.r.citrate.conf.int <- calc.bootstrap.conf.int(df$DM0.r.citrate)
    DM25.r.citrate.conf.int <- calc.bootstrap.conf.int(df$DM25.r.citrate)
    DM25.r.glucose.conf.int <- calc.bootstrap.conf.int(df$DM25.r.glucose)
    DM0.t.lag.conf.int <- calc.bootstrap.conf.int(df$DM0.t.lag)
    DM25.t.lag.conf.int <- calc.bootstrap.conf.int(df$DM25.t.lag)

    bootstrap.results <- data.frame(Name=rep(unique(df$Name),5),
                                    Founder=rep(unique(df$Founder),5),
                                    Generation=rep(unique(df$Generation),5),
                                    Parameter=c(
                                      "DM0.r.citrate",
                                      "DM25.r.citrate",
                                      "DM25.r.glucose",
                                      "DM0.t.lag",
                                      "DM25.t.lag"),
                                    Estimate = c(
                                      mean(df$DM0.r.citrate, na.rm=TRUE),
                                      mean(df$DM25.r.citrate, na.rm=TRUE),
                                      mean(df$DM25.r.glucose, na.rm=TRUE),
                                      mean(df$DM0.t.lag, na.rm=TRUE),
                                      mean(df$DM25.t.lag, na.rm=TRUE)),
                                    Left = c(
                                      DM0.r.citrate.conf.int[1],
                                      DM25.r.citrate.conf.int[1],
                                      DM25.r.glucose.conf.int[1],
                                      DM0.t.lag.conf.int[1],
                                      DM25.t.lag.conf.int[1]),
                                    Right = c(
                                      DM0.r.citrate.conf.int[2],
                                      DM25.r.citrate.conf.int[2],
                                      DM25.r.glucose.conf.int[2],
                                      DM0.t.lag.conf.int[2],
                                      DM25.t.lag.conf.int[2]),
                                    stringsAsFactors=FALSE)
    return(bootstrap.results)
  }

  ## average values for the same clone or population, over wells of the growth plate.
  confint.df <- growth %>%
  droplevels() %>% ## don't run on empty subsets
  split(.$Name) %>%
  map_dfr(.f=bootstrap.my.growth.parameters)

  return(confint.df)
}

## This function works on estimates from running growthcurver.
bootstrap.growthcurver.confints <- function(growth.fits) {
  ## average values for the same clone or population, over wells of the growth plate.

  bootstrap.growthcurver.parameters <- function(df) {

    ## growthcurver fits.
    DM0.r.conf.int <- calc.bootstrap.conf.int(df$DM0.r)
    DM0.t_mid.conf.int <- calc.bootstrap.conf.int(df$DM0.t_mid)
    DM25.r.conf.int <- calc.bootstrap.conf.int(df$DM25.r)
    DM25.t_mid.conf.int <- calc.bootstrap.conf.int(df$DM25.t_mid)

    bootstrap.results <- data.frame(Name=rep(unique(df$Name),4),
                                    Founder=rep(unique(df$Founder),4),
                                    Generation=rep(unique(df$Generation),4),
                                    Parameter=c("DM0.r",
                                      "DM25.r",
                                      "DM0.t_mid",
                                      "DM25.t_mid"
                                      ),
                                    Estimate = c(mean(df$DM0.r,na.rm=TRUE),
                                      mean(df$DM25.r,na.rm=TRUE),
                                      mean(df$DM0.t_mid, na.rm = TRUE),
                                      mean(df$DM25.t_mid, na.rm = TRUE)
                                      ),
                                    Left = c(DM0.r.conf.int[1],
                                      DM25.r.conf.int[1],
                                      DM0.t_mid.conf.int[1],
                                      DM25.t_mid.conf.int[1]
                                      ),
                                    Right = c(DM0.r.conf.int[2],
                                      DM25.r.conf.int[2],
                                      DM0.t_mid.conf.int[2],
                                      DM25.t_mid.conf.int[2]
                                      ),
                                    stringsAsFactors=FALSE)
    return(bootstrap.results)
  }

  confint.df <- growth.fits %>%
  droplevels() %>% ## don't run on empty subsets
  split(.$Name) %>%
  map_dfr(.f=bootstrap.growthcurver.parameters)

  return(confint.df)
}

plot.growth.confints <- function (plot.df, confints.df, plot.CIs=TRUE) {

  recode.Parameter <- function(df) {
    df <- mutate(df,
                 Parameter=recode(
                   Parameter,
                   DM0.r = "DM0 r",
                   DM25.r = "DM25 r",
                   DM0.r.citrate = "DM0 r citrate",
                   DM25.r.citrate = "DM25 r citrate",
                   DM25.r.glucose = "DM25 r glucose",
                   DM0.t_mid = "DM0 t_mid",
                   DM25.t_mid = "DM25 t_mid",
                   DM0.t.lag = "DM0 t_lag",
                   DM25.t.lag = "DM25 t_lag"
                   ))
    return(df)
  }

  plot.df2 <- recode.Parameter(plot.df) %>% filter(!is.na(Estimate))
  confints.df2 <- recode.Parameter(confints.df) %>% filter(!is.na(Estimate))
  
  list.of.plots <- vector("list",length(unique(plot.df2$Parameter)))
  i <- 1 ## index for list.
  for (param in unique(plot.df2$Parameter)) {
    param.df <- filter(plot.df2, Parameter==param)
    param.confint.df <- filter(confints.df2, Parameter==param)
    fig <- ggplot(param.df, aes(x = Name,y = Estimate,color = Founder)) +
    geom_point(size=0.5) +
    ggtitle(param) +
    ylab("Estimate") +
    xlab("Sample") +
    guides(color=FALSE) +
    scale_color_manual(values = cbbPalette) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(size=8),
          plot.margin = unit(c(0,0,0.25,0.25), "cm"),
          axis.text.x  = element_text(angle=90, vjust=-0.1,size=8)
          )

    ## only add y-axis labels to left-most plot.
    if (i == 1) {
      fig <- fig + theme(axis.title.y = element_text(size=10))
    } else {
      fig <- fig + theme(axis.title.y = element_blank())                  
    }
    
    ## plot CIs.
    if (plot.CIs) {
      fig <- fig +
      geom_errorbar(data = param.confint.df, aes(ymin=Left,ymax=Right), width=1, size=0.5)
    }
    ## set bounds for rate and time parameter plots.
    if (str_detect(param,'r')) {
      fig <- fig + ylim(c(0,2))
    } else {
      fig <- fig + ylim(c(0,25))    
    }
    
    list.of.plots[[i]] <- fig
    i <- i + 1
  }
  ## use cowplot to plot.
  my.plot <- plot_grid(plotlist=list.of.plots, nrow=1)
  return(my.plot)
}

reshape.growth.curve.fits <- function(growth.curve.fits) {
  ## do some gymnastics to plot rates in different media conditions together.
  DM0.growth.curve.fits <- filter(growth.curve.fits,
                                  Experiment=='DM0-growth') %>%
                                  mutate(DM0.r = r) %>%
                                  mutate(DM0.t_mid = t_mid) %>%
                                  dplyr::select(-r,-t_mid,-Experiment)

  DM25.growth.curve.fits <- filter(growth.curve.fits,
                                   Experiment=='DM25-growth') %>%
                                   mutate(DM25.r = r) %>%
                                   mutate(DM25.t_mid = t_mid) %>%
                                   dplyr::select(-r,-t_mid,-Experiment)

  reshaped.growth.curve.fits <- full_join(DM0.growth.curve.fits,
                                          DM25.growth.curve.fits)
  return(reshaped.growth.curve.fits)
}

reshape.my.growth.df <- function(growth) {
  ## do some gymnastics to plot rates in different media conditions together.
  DM0.growth <- filter(growth, Experiment=='DM0-growth') %>%
  mutate(DM0.r.citrate=r.citrate) %>%
  mutate(DM0.t.lag=t.lag) %>%
  dplyr::select(-r.glucose,-r.citrate,-t.lag,-Experiment)

  DM25.growth <- filter(growth, Experiment=='DM25-growth') %>%
  mutate(DM25.r.citrate=r.citrate) %>%
  mutate(DM25.r.glucose=r.glucose) %>%
  mutate(DM25.t.lag=t.lag) %>%
  dplyr::select(-r.glucose,-r.citrate,-t.lag,-Experiment)

  reshaped.growth <- full_join(DM0.growth, DM25.growth)
  return(reshaped.growth)
}

summarize.growthcurver.results <- function(growth.curve.fits) {
  ## lag time t_mid is the midpoint of the fit sigmoid curve.
  ## UNFORTUNATELY: cannot really infer lag time without CFUs for initial inoculum
  ## given Nkrumah's cell death observations.
  ## same goes for cell yield K, can't infer without measuring viable CFUs by plating.
  ## HOWEVER: the lag time estimated here should still be biologically relevant
  ## to experimental conditions, even if the transfer inoculum size is unknown.
  growth.curve.summary <- growth.curve.fits %>%
  group_by(SampleType, Generation, ParentClone,
           Founder, Environment, Dataset, Name) %>%
           summarize(DM0.r = mean(DM0.r, na.rm=TRUE),
                     DM0.t_mid = mean(DM0.t_mid, na.rm=TRUE),
                     DM25.r = mean(DM25.r, na.rm=TRUE),
                     DM25.t_mid = mean(DM25.t_mid, na.rm=TRUE)) %>%
                     ungroup() %>%
                     mutate(Generation=as.factor(Generation))
  return(growth.curve.summary)
}

summarize.growth.results <- function(growth) {
  growth.summary <- growth %>%
  group_by(SampleType, Generation, ParentClone,
           Founder, Environment, Dataset, Name) %>%
           summarise(DM25.r.glucose = mean(DM25.r.glucose, na.rm=TRUE),
                     DM25.r.citrate = mean(DM25.r.citrate, na.rm=TRUE),
                     DM0.r.citrate = mean(DM0.r.citrate, na.rm=TRUE),
                     DM25.t.lag = mean(DM25.t.lag, na.rm=TRUE),
                     DM0.t.lag = mean(DM0.t.lag, na.rm=TRUE)) %>%
                     ungroup() %>%
                     mutate(Generation=as.factor(Generation))
  return(growth.summary)
}

calc.growth.log.ratios <- function(growth.summary) {
  ## This one works on my estimates.
  
  ## calculate the log ratio of evolved growth to ancestral growth for rate,
  ## and time lags,
  ## and then calculate a confidence intervals around the means, using
  ## BCa bootstraps.
  ## This is a better statistical test for an increase in rate.

  ## Since the data is small, go ahead and use a for loop
  ## to make vectors corresponding to ancestral R
  ## in DM0 and DM25.

  ## my methods' estimates
  anc.DM0.r.citrate <- rep(-1,nrow(growth.summary))
  anc.DM0.t.lag <- rep(-1,nrow(growth.summary))

  anc.DM25.r.citrate <- rep(-1,nrow(growth.summary))
  anc.DM25.r.glucose <- rep(-1,nrow(growth.summary))
  anc.DM25.t.lag <- rep(-1,nrow(growth.summary))

  for (index in 1:nrow(growth.summary)) {
    my.row <- growth.summary[index, ]
    my.anc <- filter(growth.summary,Name==my.row$Founder)

    anc.DM0.r.citrate[index] <- my.anc$DM0.r.citrate
    anc.DM0.t.lag[index] <- my.anc$DM0.t.lag
    anc.DM25.r.citrate[index] <- my.anc$DM25.r.citrate
    anc.DM25.r.glucose[index] <- my.anc$DM25.r.glucose
    anc.DM25.t.lag[index] <- my.anc$DM25.t.lag
  }

  growth.summary2 <- growth.summary %>%
  mutate(log.DM0.r.citrate.ratio=log(DM0.r.citrate/anc.DM0.r.citrate)) %>%
  mutate(log.DM0.t.lag.ratio=log(DM0.t.lag/anc.DM0.t.lag)) %>%
  mutate(log.DM25.r.citrate.ratio=log(DM25.r.citrate/anc.DM25.r.citrate)) %>%
  mutate(log.DM25.r.glucose.ratio=log(DM25.r.glucose/anc.DM25.r.glucose)) %>%
  mutate(log.DM25.t.lag.ratio=log(DM25.t.lag/anc.DM25.t.lag))

  return(growth.summary2)
}

calc.growthcurver.log.ratios <- function(growth.summary) {
  ## calculate the log ratio of evolved growth to ancestral growth for rate,
  ## and time lags,
  ## and then calculate a confidence intervals around the means, using
  ## BCa bootstraps.
  ## This is a better statistical test for an increase in rate.

  ## Since the data is small, go ahead and use a for loop
  ## to make vectors corresponding to ancestral R
  ## in DM0 and DM25.

  ## growthcurver estimates
  anc.DM0.r <- rep(-1,nrow(growth.summary))
  anc.DM0.t_mid <- rep(-1,nrow(growth.summary))
  anc.DM25.r <- rep(-1,nrow(growth.summary))
  anc.DM25.t_mid <- rep(-1,nrow(growth.summary))

  for (index in 1:nrow(growth.summary)) {
    my.row <- growth.summary[index, ]
    my.anc <- filter(growth.summary,Name==my.row$Founder)

    anc.DM0.r[index] <- my.anc$DM0.r
    anc.DM0.t_mid[index] <- my.anc$DM0.t_mid
    anc.DM25.r[index] <- my.anc$DM25.r
    anc.DM25.t_mid[index] <- my.anc$DM25.t_mid
  }

  growth.summary2 <- growth.summary %>%
  mutate(log.DM0.r.ratio=log(DM0.r/anc.DM0.r)) %>%
  mutate(log.DM0.t_mid.ratio=log(DM0.t_mid/anc.DM0.t_mid)) %>%
  mutate(log.DM25.r.ratio=log(DM25.r/anc.DM25.r)) %>%
  mutate(log.DM25.t_mid.ratio=log(DM25.t_mid/anc.DM25.t_mid))

  return(growth.summary2)
}

run.growth.ratio.confint.bootstrapping <- function(final.growth.summary) {

  evolved.growth.summary <- filter(final.growth.summary, Name != Founder)
  
  log.DM0.r.citrate.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.r.citrate.ratio)
  mean.log.DM0.r.citrate.ratio <- mean(evolved.growth.summary$log.DM0.r.citrate.ratio,na.rm=TRUE)

  log.DM25.r.citrate.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.r.citrate.ratio)
  mean.log.DM25.r.citrate.ratio <- mean(evolved.growth.summary$log.DM25.r.citrate.ratio,na.rm=TRUE)

  log.DM25.r.glucose.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.r.glucose.ratio)
  mean.log.DM25.r.glucose.ratio <- mean(evolved.growth.summary$log.DM25.r.glucose.ratio,na.rm=TRUE)
  
  log.DM0.t.lag.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.t.lag.ratio)
  mean.log.DM0.t.lag.ratio <- mean(evolved.growth.summary$log.DM0.t.lag.ratio)

  log.DM25.t.lag.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.t.lag.ratio)
  mean.log.DM25.t.lag.ratio <- mean(evolved.growth.summary$log.DM25.t.lag.ratio)
  
  bootstrap.results <- data.frame (
    Parameter = c(
      "log.DM0.r.citrate.ratio",
      "log.DM25.r.citrate.ratio",
      "log.DM25.r.glucose.ratio",
      "log.DM0.t.lag.ratio",
      "log.DM25.t.lag.ratio"),
    Estimate = c(
      mean.log.DM0.r.citrate.ratio,
      mean.log.DM25.r.citrate.ratio,
      mean.log.DM25.r.glucose.ratio,
      mean.log.DM0.t.lag.ratio,
      mean.log.DM25.t.lag.ratio),
    Left = c(
      log.DM0.r.citrate.ratio.conf.int[1],
      log.DM25.r.citrate.ratio.conf.int[1],
      log.DM25.r.glucose.ratio.conf.int[1],
      log.DM0.t.lag.ratio.conf.int[1],
      log.DM25.t.lag.ratio.conf.int[1]),
    Right = c(
      log.DM0.r.citrate.ratio.conf.int[2],
      log.DM25.r.citrate.ratio.conf.int[2],
      log.DM25.r.glucose.ratio.conf.int[2],
      log.DM0.t.lag.ratio.conf.int[2],
      log.DM25.t.lag.ratio.conf.int[2]),
    stringsAsFactors=FALSE)
  
  return(bootstrap.results)
}

run.growthcurver.ratio.confint.bootstrapping <- function(final.growth.summary) {

  evolved.growth.summary <- filter(final.growth.summary, Name != Founder)
  
  log.DM0.r.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.r.ratio)
  mean.log.DM0.r.ratio <- mean(evolved.growth.summary$log.DM0.r.ratio,na.rm=TRUE)

  log.DM0.t_mid.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.t_mid.ratio)
  mean.log.DM0.t_mid.ratio <- mean(evolved.growth.summary$log.DM0.t_mid.ratio)

  log.DM25.r.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.r.ratio)
  mean.log.DM25.r.ratio <- mean(evolved.growth.summary$log.DM25.r.ratio,na.rm=TRUE)
  
  log.DM25.t_mid.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.t_mid.ratio)
  mean.log.DM25.t_mid.ratio <- mean(evolved.growth.summary$log.DM25.t_mid.ratio)
  
  bootstrap.results <- data.frame (
    Parameter = c("log.DM0.r.ratio",
      "log.DM25.r.ratio",
      "log.DM0.t_mid.ratio",
      "log.DM25.t_mid.ratio"),
    Estimate = c(mean.log.DM0.r.ratio,
      mean.log.DM25.r.ratio,
      mean.log.DM0.t_mid.ratio,
      mean.log.DM25.t_mid.ratio),
    Left = c(log.DM0.r.ratio.conf.int[1],
      log.DM25.r.ratio.conf.int[1],
      log.DM0.t_mid.ratio.conf.int[1],
      log.DM25.t_mid.ratio.conf.int[1]),
    Right = c(log.DM0.r.ratio.conf.int[2],
      log.DM25.r.ratio.conf.int[2],
      log.DM0.t_mid.ratio.conf.int[2],
      log.DM25.t_mid.ratio.conf.int[2]),
    stringsAsFactors=FALSE)
  
  return(bootstrap.results)
}

plot.growth.parameters <- function(growth, plot.CIs=TRUE) {    
  growth.CI <- bootstrap.my.growth.confints(growth)
  growth.plot.df <- growth %>%
  gather(key="Parameter",value="Estimate",
         DM0.r.citrate,DM25.r.citrate,DM25.r.glucose,DM0.t.lag,DM25.t.lag)
  growth.plot <- plot.growth.confints(growth.plot.df, growth.CI, plot.CIs)
  return(growth.plot)    
}

plot.growthcurver.parameters <- function(growth, plot.CIs=TRUE) {    
  growth.CI <- bootstrap.growthcurver.confints(growth)
  growth.plot.df <- growth %>%
  gather(key="Parameter",value="Estimate",
         DM0.r,DM25.r,DM0.t_mid,DM25.t_mid)
  growth.plot <- plot.growth.confints(growth.plot.df, growth.CI, plot.CIs)
  return(growth.plot)   
}

plot.parameter.log.ratios <- function (plot.df, confints.df, recode.param.func) {

  plot.df2 <- recode.param.func(plot.df) %>%
  dplyr::select(SampleType, Generation, ParentClone, Founder, Environment,
         Dataset, Name, Parameter, Estimate) %>% drop_na()

  confints.df2 <- recode.param.func(confints.df) %>% drop_na()
  
  the.plot <- ggplot(plot.df2, aes(x=Parameter,y=Estimate)) +
  geom_point(size=1) +
  geom_errorbar(data=confints.df2, aes(ymin=Left,ymax=Right),width=0.5, size=0.5) +
  theme_classic() +
  geom_hline(yintercept=0,linetype='dashed') +
  ylab("log(Evolved/Ancestral)") +
  xlab("Growth parameter") +
  theme(axis.text.x=element_text(size=12))
  
  return(the.plot)
}

plot.growth.parameter.log.ratios <- function (plot.df, confints.df) {

  recode.Growth.Parameter <- function(df) {
    df2 <- mutate(df,
                 Parameter=recode_factor(
                   Parameter,
                   log.DM0.r.citrate.ratio = "DM0 r citrate",              
                   log.DM25.r.citrate.ratio = "DM25 r citrate",
                   log.DM25.r.glucose.ratio = "DM25 r glucose",
                   log.DM0.t.lag.ratio = "DM0 t_lag",
                   log.DM25.t.lag.ratio = "DM25 t_lag")
                 )                   
    return(df2)
  }

  return(plot.parameter.log.ratios(plot.df,confints.df,recode.Growth.Parameter))
}

plot.growthcurver.parameter.log.ratios <- function (plot.df, confints.df) {

  recode.Growthcurver.Parameter <- function(df) {
    df2 <- mutate(df,
                 Parameter=recode_factor(
                   Parameter,
                   log.DM0.r.ratio = "DM0 r",
                   log.DM25.r.ratio = "DM25 r",
                   log.DM0.t_mid.ratio = "DM0 t_mid",
                   log.DM25.t_mid.ratio = "DM25 t_mid")
                 ) 
    return(df2)
  }
  
  return(plot.parameter.log.ratios(plot.df,confints.df,recode.Growthcurver.Parameter))
}

###############################################
## Load the growth data from the plate reader.
## DM0-evolved population data.
DM0.pop.growth.data <- read.csv(
  file.path(projdir,"results/DM0-evolved-pop-growth.csv"),
  stringsAsFactors=FALSE) %>%
left_join(pop.clone.labels, by="Name") %>%
prep.growth.df()

## DM0-evolved clone data.
DM0.clone.growth.data <- read.csv(
  file.path(projdir,"results/DM0-evolved-clone-growth.csv"),
  stringsAsFactors=FALSE) %>%
left_join(pop.clone.labels, by="Name") %>%
prep.growth.df()

## DM25-evolved clone data (only have growth in DM25).
DM25.clone.growth.data <- read.csv(
  file.path(projdir, "results/DM25-evolved-clone-growth.csv"),
  stringsAsFactors=FALSE) %>%
left_join(pop.clone.labels, by="Name") %>%
prep.growth.df()

#######################################
## Examine how growth parameters have changed over evolution.
## Run calc.growth.rates code.
clone.growth <- calc.growth.rates(DM0.clone.growth.data) %>%
mutate(Dataset='CloneGrowth') %>%
reshape.my.growth.df()

pop.growth <- calc.growth.rates(DM0.pop.growth.data) %>%
mutate(Dataset='PopulationGrowth') %>%
reshape.my.growth.df()

DM25.clone.growth <- calc.growth.rates(DM25.clone.growth.data) %>%
mutate(Dataset='CloneGrowth') %>%
reshape.my.growth.df()

## use growthcurver to fit logistic curves.
clone.growth.curve.fits <- map.reduce.growth.curve(DM0.clone.growth.data) %>%
mutate(Dataset='CloneGrowth') %>%
reshape.growth.curve.fits()

pop.growth.curve.fits <- map.reduce.growth.curve(DM0.pop.growth.data)%>%
mutate(Dataset='PopulationGrowth') %>%
reshape.growth.curve.fits()

DM25.clone.growth.curve.fits <- map.reduce.growth.curve(DM25.clone.growth.data) %>%
mutate(Dataset='CloneGrowth') %>%
reshape.growth.curve.fits()

#######################################
## Plot growth curves.

## plot ZDBp874 data to calibrate citrate interval.
ZDBp874.df <- filter(DM0.clone.growth.data,Name %in% c('ZDBp874','CZB151'))
ZDBp874.plot <- plot.single.growthcurve(ZDBp874.df,'ZDBp874','CZB151', logscale=FALSE)

## plot rest of CZB151-descended isolates. Can do good fitness competitions on these.
ZDBp871.df <- filter(DM0.clone.growth.data, Name %in% c('ZDBp871','CZB151'))
ZDBp871.plot <- plot.single.growthcurve(ZDBp871.df,'ZDBp871','CZB151', logscale=FALSE)

ZDBp889.df <- filter(DM0.clone.growth.data, Name %in% c('ZDBp889','CZB151'))
ZDBp889.plot <- plot.single.growthcurve(ZDBp889.df,'ZDBp889','CZB151', logscale=FALSE)

ZDBp892.df <- filter(DM0.clone.growth.data, Name %in% c('ZDBp892','CZB151'))
ZDBp892.plot <- plot.single.growthcurve(ZDBp892.df,'ZDBp892','CZB151', logscale=FALSE)

## SUPER COOL! ZDBp874 grows WORSE than ancestor in DM25!!!!
## ZDBp871 and ZDBp889 also grow worse in DM25 than ancestor!
## This is really interesting!

## S3 Fig: Plot growth curves for DM0-evolved clones.
S3Fig.plot <- plot.growthcurve.figure(filter(DM0.clone.growth.data,Hours<=24), logscale=FALSE)
save_plot(file.path(projdir,"results/figures/S3Fig.pdf"),S3Fig.plot,base_height=7,base_width=11)

## S4 Fig: growth curves for DM0-evolved clones on log-scale.
S4Fig.plot <- plot.growthcurve.figure(DM0.clone.growth.data,logscale=TRUE)
save_plot(file.path(projdir,"results/figures/S4Fig.pdf"),S4Fig.plot,base_height=7,base_width=11)

## S5 Fig: Plot growth curves for DM0-evolved populations.
S5Fig.plot <- plot.growthcurve.figure(DM0.pop.growth.data, logscale=FALSE)
save_plot(file.path(projdir,"results/figures/S5Fig.pdf"), S5Fig.plot,base_height=7,base_width=11)

## S6 Fig: growth curves for DM0-evolved populations on log-scale.
S6Fig.plot <- plot.growthcurve.figure(DM0.pop.growth.data,logscale=TRUE)
save_plot(file.path(projdir,"results/figures/S6Fig.pdf"),S6Fig.plot,base_height=7,base_width=11)

## S7 Figure. Plot growth curves for DM25-evolved clones.
S7Figoutf <- file.path(projdir,"results/figures/S7Fig.pdf")
S7Fig.plot <- plot.growthcurve.figure(DM25.clone.growth.data, logscale=FALSE)
save_plot(S7Figoutf, S7Fig.plot,base_height=7,base_width=11)


## S8 Figure. Plot growth curves for DM25-evolved clones on a log-scale.
S8Fig.plot <- plot.growthcurve.figure(DM25.clone.growth.data, logscale=TRUE)
save_plot(file.path(projdir,"results/figures/S8Fig.pdf"), S8Fig.plot,base_height=7,base_width=11)

## data filtering to plot data included in my growth rate estimation.
## This code is useful for debugging and making sure the rate estimation works right.
filtered.clone.growth.data <- filter.growth.data.on.analysis.domain(DM0.clone.growth.data)
filtered.pop.growth.data <- filter.growth.data.on.analysis.domain(DM0.pop.growth.data)
## filtered log-scale plots, again for debugging purposes.
filtered.log.clone.plot <- plot.growthcurve.figure(filtered.clone.growth.data,logscale=TRUE)
filtered.log.pop.plot <- plot.growthcurve.figure(filtered.pop.growth.data,logscale=TRUE)

################################################################################
## plot my growth estimates and growthcurver estimates as well.
## plot estimates with confidence interval of mean.
## Do the same for the growthcurver fits. Note different confint function needed.

## Filter oddball Cit- ZDBp874 from the clones before estimation.
clone.growth <- filter(clone.growth,Name != 'ZDBp874')

## Figure 3.
Fig3outf <- file.path(projdir,"results/figures/Fig3.pdf")
clone.growth.plot <- plot.growth.parameters(clone.growth)
## See what the figure looks like without CIs.
no.CI.clone.growth.plot <- plot.growth.parameters(clone.growth, plot.CIs=FALSE)

clone.growth.summary <- summarize.growth.results(clone.growth)
final.clone.growth.summary <- calc.growth.log.ratios(clone.growth.summary)
clone.bootstrap.CI <- run.growth.ratio.confint.bootstrapping(final.clone.growth.summary)

Fig3B.plot.df <- final.clone.growth.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
       log.DM25.t.lag.ratio, log.DM0.t.lag.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.

Fig3B <- plot.growth.parameter.log.ratios(Fig3B.plot.df, clone.bootstrap.CI)

Fig3 <- plot_grid(clone.growth.plot,Fig3B,labels=c('A','B'),ncol=1,
                  rel_heights = c(1, 1),rel_widths=c(1,1))
save_plot(Fig3outf,Fig3, base_width=11, base_height=7.5)

## See what the figure looks like without CIs.
no.CI.Fig3outf <- file.path(projdir,"results/figures/noCI.Fig3.pdf")
no.CI.Fig3 <- plot_grid(no.CI.clone.growth.plot,Fig3B,labels=c('A','B'),ncol=1,
                  rel_heights = c(1, 1),rel_widths=c(1,1))
save_plot(no.CI.Fig3outf, no.CI.Fig3, base_width=11, base_height=7.5)


## Supplementary Figure S10. (since growthcurver results are in its own section.)
S10Figoutf <- file.path(projdir,"results/figures/S10Fig.pdf")

## Filter oddball Cit- ZDBp874 from the clones before estimation.
clone.growth.curve.fits <- filter(clone.growth.curve.fits, Name != 'ZDBp874')

clone.growthcurver.plot <- plot.growthcurver.parameters(clone.growth.curve.fits)
clone.growthcurver.summary <- summarize.growthcurver.results(clone.growth.curve.fits)
final.clone.growthcurver.summary <- calc.growthcurver.log.ratios(clone.growthcurver.summary)
clone.growthcurver.bootstrap.CI <- run.growthcurver.ratio.confint.bootstrapping(final.clone.growthcurver.summary)

S10BFig.plot.df <- final.clone.growthcurver.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM0.r.ratio, log.DM25.r.ratio,
       log.DM0.t_mid.ratio, log.DM25.t_mid.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.

S10BFig <- plot.growthcurver.parameter.log.ratios(S10BFig.plot.df, clone.growthcurver.bootstrap.CI)

S10Fig <- plot_grid(clone.growthcurver.plot, S10BFig,labels=c('A','B'), ncol=1,
                   rel_heights = c(1, 1),rel_widths=c(1,1))
save_plot(S10Figoutf,S10Fig, base_width=11, base_height=7.5)

## Figure 4.
Fig4outf <- file.path(projdir,"results/figures/Fig4.pdf")
pop.growth.plot <- plot.growth.parameters(pop.growth)
## See what the figure looks like without CIs.
no.CI.pop.growth.plot <- plot.growth.parameters(pop.growth, plot.CIs=FALSE)

pop.growth.summary <- summarize.growth.results(pop.growth)
final.pop.growth.summary <- calc.growth.log.ratios(pop.growth.summary)

pop.bootstrap.CI <- run.growth.ratio.confint.bootstrapping(final.pop.growth.summary)
Fig4B.plot.df <- final.pop.growth.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
       log.DM25.t.lag.ratio, log.DM0.t.lag.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.
Fig4B <- plot.growth.parameter.log.ratios(Fig4B.plot.df, pop.bootstrap.CI)

Fig4 <- plot_grid(pop.growth.plot,Fig4B,labels=c('A','B'),ncol=1,
                  rel_heights = c(1, 1))
save_plot(Fig4outf,Fig4, base_width=11, base_height=7.5)

## See what the figure looks like without CIs.
no.CI.Fig4outf <- file.path(projdir,"results/figures/noCI.Fig4.pdf")
no.CI.Fig4 <- plot_grid(no.CI.pop.growth.plot,Fig4B,labels=c('A','B'),ncol=1,
                  rel_heights = c(1, 1),rel_widths=c(1,1))
save_plot(no.CI.Fig4outf, no.CI.Fig4, base_width=11, base_height=7.5)


## Supplementary Figure S11 (because in separate section of supplement).
S11Figoutf <- file.path(projdir,"results/figures/S11Fig.pdf")
pop.growthcurver.plot <- plot.growthcurver.parameters(pop.growth.curve.fits)
pop.growthcurver.summary <- summarize.growthcurver.results(pop.growth.curve.fits)
final.pop.growthcurver.summary <- calc.growthcurver.log.ratios(pop.growthcurver.summary)

pop.growthcurver.bootstrap.CI <- run.growthcurver.ratio.confint.bootstrapping(
  final.pop.growthcurver.summary
  )

S11BFig.plot.df <- final.pop.growthcurver.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM0.r.ratio, log.DM25.r.ratio,
       log.DM0.t_mid.ratio, log.DM25.t_mid.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.

S11BFig <- plot.growthcurver.parameter.log.ratios(S11BFig.plot.df, pop.growthcurver.bootstrap.CI)

S11Fig <- plot_grid(pop.growthcurver.plot, S11BFig,labels=c('A','B'), ncol=1,
                   rel_heights = c(1, 1),rel_widths=c(1,1))
save_plot(S11Figoutf,S11Fig, base_width=11, base_height=7.5)


## growth results for DM25-evolved clones.
## Supplementary Figure S9.
S9Figoutf <- file.path(projdir,"results/figures/S9Fig.pdf")
DM25.clone.growth.plot <- plot.growth.parameters(DM25.clone.growth)
DM25.growth.summary <- summarize.growth.results(DM25.clone.growth)
final.DM25.growth.summary <- calc.growth.log.ratios(DM25.growth.summary)
DM25.clone.bootstrap.CI <- run.growth.ratio.confint.bootstrapping(final.DM25.growth.summary)

S9BFig.plot.df <- final.DM25.growth.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
       log.DM25.t.lag.ratio, log.DM0.t.lag.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.
S9BFig <- plot.growth.parameter.log.ratios(S9BFig.plot.df, DM25.clone.bootstrap.CI)

S9Fig <- plot_grid(DM25.clone.growth.plot,S9BFig,labels=c('A','B'),ncol=1,
                    rel_heights = c(1, 1))
save_plot(S9Figoutf,S9Fig, base_width=11, base_height=7.5)


## Supplementary Figure S12. supplementary fig of DM25-evolved clone growth.
S12Figoutf <- file.path(projdir,"results/figures/S12Fig.pdf")
DM25.clone.growthcurver.plot <- plot.growthcurver.parameters(DM25.clone.growth.curve.fits)
DM25.growthcurver.summary <- summarize.growthcurver.results(DM25.clone.growth.curve.fits)
final.DM25.growthcurver.summary <- calc.growthcurver.log.ratios(DM25.growthcurver.summary)
DM25.growthcurver.bootstrap.CI <- run.growthcurver.ratio.confint.bootstrapping(final.DM25.growthcurver.summary)

S12FigB.plot.df <- final.DM25.growthcurver.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.ratio, log.DM25.t_mid.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.

S12FigB <- plot.growthcurver.parameter.log.ratios(S12FigB.plot.df, DM25.growthcurver.bootstrap.CI)

S12Fig <- plot_grid(DM25.clone.growthcurver.plot,
                    S12FigB,
                    labels=c('A','B'),ncol=1,
                    rel_heights = c(1, 1))

save_plot(S12Figoutf,S12Fig, base_width=11, base_height=7.5)

##################### Supplementary figure:
## Compare growthcurver estimates to my estimates,
## in order to show consistent results between the two methods.

growth.summary <- full_join(clone.growth.summary, pop.growth.summary)
growthcurver.summary <- full_join(pop.growthcurver.summary, clone.growthcurver.summary)

growth.estimate.comp.df <- full_join(
  growth.summary,
  growthcurver.summary) %>%
ungroup()

clone.estimates <- filter(growth.estimate.comp.df,Dataset=='CloneGrowth')
pop.estimates <- filter(growth.estimate.comp.df,Dataset=='PopulationGrowth')

DM25.growth.estimate.comp.df <- full_join(
  DM25.growth.summary,
  DM25.growthcurver.summary) %>%
ungroup()

## 
## r = 0.68, 0.70, -0.14, 0.86, 0.37 respectively.
cor(growth.estimate.comp.df$DM0.r, growth.estimate.comp.df$DM0.r.citrate,use="complete.obs")
cor.test(growth.estimate.comp.df$DM0.r, growth.estimate.comp.df$DM0.r.citrate,use="complete.obs",method="pearson")

cor(growth.estimate.comp.df$DM25.r, growth.estimate.comp.df$DM25.r.citrate,use="complete.obs")
cor.test(growth.estimate.comp.df$DM25.r, growth.estimate.comp.df$DM25.r.citrate,use="complete.obs",method="pearson")

cor(growth.estimate.comp.df$DM25.r, growth.estimate.comp.df$DM25.r.glucose,use="complete.obs")
cor.test(growth.estimate.comp.df$DM25.r, growth.estimate.comp.df$DM25.r.glucose,use="complete.obs",method="pearson")

cor(growth.estimate.comp.df$DM0.t_mid, growth.estimate.comp.df$DM0.t.lag,use="complete.obs")
cor.test(growth.estimate.comp.df$DM0.t_mid, growth.estimate.comp.df$DM0.t.lag,use="complete.obs",method="pearson")

cor(growth.estimate.comp.df$DM25.t_mid, growth.estimate.comp.df$DM25.t.lag,use="complete.obs")
cor.test(growth.estimate.comp.df$DM25.t_mid, growth.estimate.comp.df$DM25.t.lag,use="complete.obs",method="pearson")

## r = 0.64, 0.83, -0.04, 0.77, 0.46, respectively.
cor(pop.estimates$DM0.r, pop.estimates$DM0.r.citrate,use="complete.obs")
cor(pop.estimates$DM25.r, pop.estimates$DM25.r.citrate,use="complete.obs")
cor(pop.estimates$DM25.r, pop.estimates$DM25.r.glucose,use="complete.obs")
cor(pop.estimates$DM0.t_mid, pop.estimates$DM0.t.lag,use="complete.obs")
cor(pop.estimates$DM25.t_mid, pop.estimates$DM25.t.lag,use="complete.obs")

## r = 0.48, 0.60, -0.11, 0.90, 0.38, respectively.
cor(clone.estimates$DM0.r, clone.estimates$DM0.r.citrate,use="complete.obs")
cor(clone.estimates$DM25.r, clone.estimates$DM25.r.citrate,use="complete.obs")
cor(clone.estimates$DM25.r, clone.estimates$DM25.r.glucose,use="complete.obs")
cor(clone.estimates$DM0.t_mid, clone.estimates$DM0.t.lag,use="complete.obs")
cor(clone.estimates$DM25.t_mid, clone.estimates$DM25.t.lag,use="complete.obs")

## Note: the following are for the DM0-evolved clones.
param.comp.plot1 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM0.r,y=DM0.r.citrate,color=Dataset)) +
scale_color_manual(values=c('purple','green')) +
geom_point() +
theme_classic() +
guides(color=FALSE) +
xlim(0,2) +
ylim(0,2) +
geom_abline(slope=1,intercept=0,linetype='dotted') +
xlab("DM0 r") +
ylab("DM0 r citrate")

param.comp.plot2 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM25.r,y=DM25.r.citrate,color=Dataset)) +
scale_color_manual(values=c('purple','green')) +
geom_point() +
theme_classic() +
guides(color=FALSE) +
xlim(0,2) +
ylim(0,2) +
geom_abline(slope=1,intercept=0,linetype='dotted') +
xlab("DM25 r") +
ylab("DM25 r citrate")

param.comp.plot3 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM25.r,y=DM25.r.glucose,color=Dataset)) +
scale_color_manual(values=c('purple','green')) +
geom_point() +
theme_classic() +
guides(color=FALSE) +
xlim(0,2) +
ylim(0,2) +
geom_abline(slope=1,intercept=0,linetype='dotted') +
xlab("DM25 r") +
ylab("DM25 r glucose")

param.comp.plot4 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM0.t_mid,y=DM0.t.lag,color=Dataset)) +
scale_color_manual(values=c('purple','green')) +
geom_point() +
theme_classic() +
guides(color=FALSE) +
xlab("DM0 t_mid") +
ylab("DM0 t_lag")

param.comp.plot5 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM25.t_mid,y=DM25.t.lag,color=Dataset)) +
scale_color_manual(values=c('purple','green')) +
geom_point() +
theme_classic() +
guides(color=FALSE) +
xlab("DM25 t_mid") +
ylab("DM25 t_lag")

## Save figure for the Supplement. 
S13Fig <- plot_grid(param.comp.plot1,
                    param.comp.plot2,
                    param.comp.plot3,
                    param.comp.plot4,
                    param.comp.plot5,
                    labels=c('A','B','C','D','E'), ncol=1)
S13Fig.outf <- file.path(projdir, "results/figures/S13Fig.pdf")
save_plot(S13Fig.outf, S13Fig, base_height=10,base_width=4)

####################################
### Supplementary Figure S14:
### no correlation between glucose growth rates and citrate growth rates,
### but  correlation in citrate growth rates in DM0 and DM25.

## no correlation between r citrate and r glucose in DM25 by my estimates.
cor.test(x=clone.estimates$DM25.r.citrate,y=clone.estimates$DM25.r.glucose,method="pearson")
cor.test(x=pop.estimates$DM25.r.citrate,y=pop.estimates$DM25.r.glucose,method="pearson")

## significant correlation between r citrate in DM25 and DM0 by my estimates.
cor.test(x=clone.estimates$DM25.r.citrate,y=clone.estimates$DM0.r.citrate,method="pearson")
cor.test(x=pop.estimates$DM25.r.citrate,y=pop.estimates$DM0.r.citrate,method="pearson")

## correlations between growth rate in DM25 and DM0 by growthcurver:
##  significant for clones and marginially insignificant for populations.
## Pearson's r = 0.74, p = 0.002 for clones; r = 0.50, p = 0.057 for pop.
cor.test(x=clone.estimates$DM25.r,y=clone.estimates$DM0.r,method="pearson")
cor.test(x=pop.estimates$DM25.r,y=pop.estimates$DM0.r,method="pearson")

cit.glucose.cor.plot <- ggplot(growth.summary,
                               aes(x=DM25.r.citrate,y=DM25.r.glucose, color=Founder, shape=Generation)) +
facet_wrap(Founder~Dataset,nrow=1) +
geom_point() +
theme_classic() +
xlab("DM25 r citrate") +
ylab("DM25 r glucose") +
scale_color_manual(values=cbbPalette) +
guides(color=FALSE,shape=FALSE)

## Think more about correlation between growth on citrate in DM0 and growth on citrate in DM25.
cit.cit.cor.plot <- ggplot(growth.summary,
                           aes(x=DM25.r.citrate,y=DM0.r.citrate,color=Founder, shape=Generation)) +
facet_wrap(Founder~Dataset,nrow=1) +
geom_point() +
theme_classic() +
xlab("DM25 r citrate") +
ylab("DM0 r citrate") +
scale_color_manual(values=cbbPalette) +
guides(color=FALSE,shape=FALSE)

## Compare with cit.cit.cor.plot. Very nice! Results look robust.
r.r.cor.plot <- ggplot(growthcurver.summary,
                       aes(x=DM25.r,y=DM0.r,color=Founder, shape=Generation)) +
facet_wrap(Founder~Dataset,nrow=1) +
geom_point() +
theme_classic() +
xlab("DM25 r") +
ylab("DM0 r") +
scale_color_manual(values=cbbPalette) +
guides(color=FALSE,shape=FALSE)

## Save figure for the Supplement. 
S14Fig <- plot_grid(cit.glucose.cor.plot,
                    cit.cit.cor.plot,
                    r.r.cor.plot, labels=c('A','B','C'), ncol=1)
S14Fig.outf <- file.path(projdir,"results/figures/S14Fig.pdf")
save_plot(S14Fig.outf, S14Fig, base_height=8, base_width=12)
######################################################################
####### Figure 5: Nkrumah's cell death results. See CellDeath R script
#######           for analyses and figures.
######################################################################
## Figure 6: make a matrix plot of genes with mutations in two or more clones.

PlotMatrixFigure <- function(raw.matrix.file, amp.matrix.file,
                             ltee.matrix.file, ltee.50k.labels.file,
                             matrix.outfile, pop.clone.labels) {

  raw.matrix <- read.csv(raw.matrix.file, stringsAsFactors=FALSE)

  ## fix the names of the samples.
  names(raw.matrix) <- map_chr(
    names(raw.matrix),
    function (x) str_trunc(x,width=7,side="right",ellipsis=''))

  #' nice helper function
  not_any_na <- function(x) {!any(is.na(x))}

  #' import rows for the maeA and dctA amplifications.
  amp.matrix <- read.csv(amp.matrix.file, stringsAsFactors=FALSE)
  #' and merge with the mutation matrix.
  #' remove columns that contain any NA values (ZDBp875)
  merged.with.amps.matrix <- full_join(raw.matrix,amp.matrix) %>%
  dplyr::select_if(not_any_na)

  DM0.DM25.matrix.data <- gather(merged.with.amps.matrix,"Name","mutation.count",2:ncol(merged.with.amps.matrix)) %>%
  left_join(pop.clone.labels) %>%
  dplyr::select(Gene,Name,mutation.count,Environment,PopulationLabel) %>%
  group_by(Gene) %>% filter(sum(mutation.count)>1)

  #' print out a table of parallelism across the DM0 and DM25 treatments.
  parallel.counts <- summarize(DM0.DM25.matrix.data,total.count=sum(mutation.count)) %>%
  arrange(desc(total.count))
  print('Parallel hits in genes')
  print(filter(parallel.counts,total.count>1),n=Inf)

  #' add LTEE mutation matrix to figure.
  ltee.matrix <- read.csv(ltee.matrix.file, stringsAsFactors=FALSE)

  #' fix the names of the samples.
  names(ltee.matrix) <- map_chr(
    names(ltee.matrix),
    function (x) str_trunc(x,width=8,side="left",ellipsis=''))

  LTEE.50K.labels <- read.csv(ltee.50k.labels.file, stringsAsFactors=FALSE)

  #' filter on genes mutated in the DM0 and DM25 matrix data.
  ltee.data <- gather(ltee.matrix,"Name","mutation.count",2:13) %>%
  left_join(LTEE.50K.labels) %>%
  group_by(Gene) %>% filter(Gene %in% DM0.DM25.matrix.data$Gene) %>%
  filter(!is.na(Gene))

  ## keep non-mutators and the Ara-3 50K clone.
  non.mutators.and.ara.minus.3 <- ltee.data %>%
  filter(Hypermutator == 0 | PopulationLabel == 'Ara-3') %>%
  dplyr::select(-Hypermutator)

  ## Now join LTEE data to the DM0 and DM25 Cit+ data.
  ## Give the Ara-3 LTEE clone a special label.
  matrix.data <- DM0.DM25.matrix.data %>%
  bind_rows(non.mutators.and.ara.minus.3) %>%
  mutate(Name=ifelse(Name=='REL11364','Ara-3: REL11364',Name))

  ## sort genes by number of mutations in each row.
  gene.hit.sort <- group_by(matrix.data,Gene) %>% summarize(hits=sum(mutation.count)) %>%
  arrange(desc(hits))
  matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
  
  ## cast mutation.count into a factor for plotting.
  matrix.data$mutation.count <- factor(matrix.data$mutation.count)
  
  make.matrix.panel <- function(mdata,environment) {
    ggplot(filter(mdata,Environment==environment),
           aes(x=Name,
               y=Gene,
               fill=mutation.count,
               frame=Environment)
           ) +
           geom_tile(color="black",size=0.1) +
           ylab("Gene") +
           xlab("Genome") +
           ggtitle(paste(environment,'genomes')) +
           theme_tufte(base_family='Helvetica') +
           theme(axis.ticks = element_blank(),
                 axis.text.x = element_text(size=10,angle=45,hjust=1),
                 axis.text.y = element_text(size=10,hjust=1,face="italic"),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 ) +
                 guides(fill=FALSE) +
                 scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
                 scale_fill_manual(values = c("white", "#ffdf00", "#bebada",
                                     "#fb8072", "#80b1d3", "#fdb462"))
  }
  
  DM0.matrix.panel <- make.matrix.panel(matrix.data,'DM0')
  DM25.matrix.panel <- make.matrix.panel(matrix.data,'DM25')
  LTEE.matrix.panel <- make.matrix.panel(matrix.data,'LTEE')

  matrix.figure <- plot_grid(DM0.matrix.panel,
                             DM25.matrix.panel,
                             LTEE.matrix.panel,
                             nrow=1,
                             rel_widths=c(13/12,1,7/12),
                             align = 'vh')
                             
  ggsave(matrix.figure, file=matrix.outfile,width=10,height=10)
}

##### Now run this function!
#' on all valid mutations
raw.matrix.f <- file.path(projdir,
                          "results/DM0-DM25-comparison-mut-matrix.csv")

#' on just non-synonymous mutations
dN.raw.matrix.f <- file.path(projdir,
                             "results/dN-DM0-DM25-comparison-mut-matrix.csv")

amp.matrix.f <- file.path(projdir,"results/amp_matrix.csv")
ltee.matrix.f <- file.path(projdir,"results/LTEE-mut_matrix.csv")
dN.ltee.matrix.f <- file.path(projdir,"results/dN-LTEE-mut_matrix.csv")
ltee.50k.labels.f <- file.path(projdir, "data/rohan-formatted/LTEE-50K-clones.csv")
matrix.outf <- "../results/figures/mut_matrix.pdf"
dN.matrix.outf <- "../results/figures/dN_mut_matrix.pdf"

PlotMatrixFigure(raw.matrix.f, amp.matrix.f, ltee.matrix.f, ltee.50k.labels.f, matrix.outf, pop.clone.labels)
PlotMatrixFigure(dN.raw.matrix.f, amp.matrix.f, dN.ltee.matrix.f, ltee.50k.labels.f, dN.matrix.outf, pop.clone.labels)

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
  file.path(projdir,
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
## Figure 7. IS element analysis and visualization.
################################################

IS.palette <- c('#f4a582','#92c5de','#ca0020','black','#0571b0')

IS.insertions <- read.csv(file.path(projdir,
                                    "results/genome-analysis/IS_insertions.csv")) %>%
arrange(genome_start) %>%
## filter oddball Cit- clone ZDBp874.
filter(Clone != 'ZDBp874')

IS.plot <- ggplot(IS.insertions,aes(x=genome_start,fill=IS_element,frame=Environment)) +
facet_grid(Environment~.) +
geom_histogram(bins=200) +
guides(fill=FALSE) +
scale_fill_manual(values=IS.palette) +
ylab("Count") +
xlab("Position") +
theme_tufte(base_family="Helvetica")

## 75/202 IS insertions recur at the same locations! 37.1%!
parallel.IS.insertions <- IS.insertions %>%
group_by(genome_start) %>%
filter(n()>1)

summarize.IS.func <- function(parallel.IS.ins) {
  summarize(parallel.IS.ins,
            IS_element=unique(IS_element),
            genes_inactivated=unique(genes_inactivated),
            genes_promoter=unique(genes_promoter),
            count=n()) %>%
            arrange(desc(count)) %>%
            mutate(genes_inactivated=as.character(genes_inactivated)) %>%
            mutate(genes_promoter=as.character(genes_promoter)) %>%
            mutate(annotation=ifelse(is.na(genes_inactivated),
                     genes_promoter,
                     genes_inactivated)) %>%
                     mutate(annotation=ifelse(is.na(annotation),'none',annotation))
}

parallel.IS.summary <- summarize.IS.func(parallel.IS.insertions)

## plot parallel IS insertions and their annotation.
parallel.IS.plot <- ggplot(parallel.IS.summary,aes(x=genome_start,
                                                   y=count,
                                                   color=IS_element,
                                                   label=annotation)) +
geom_point() +
guides(color=FALSE) +
scale_color_manual(values = IS.palette) +
theme_classic() +
ylab("Count") +
xlab("Position") +
geom_text_repel(fontface = "italic")

## Now just consider IS parallelism in the DM0 environment.
parallel.DM0.IS.insertions <- IS.insertions %>%
filter(Environment == 'DM0') %>%
group_by(genome_start) %>%
filter(n()>1)

parallel.DM0.IS.summary <- summarize.IS.func(parallel.DM0.IS.insertions) 

########
## Plot the rate of increase of IS-elements in the DM0 and DM25 experiments
## in comparison to the rate of increase of IS-elements in Ara-3.

LTEE.MAE.IS.insertions <- read.csv(
  file.path(projdir,
            "results/genome-analysis/LTEE_MAE_IS150_insertions.csv"),
  stringsAsFactors=FALSE) %>%
arrange(Position)

LTEE.IS150 <- filter(LTEE.MAE.IS.insertions,Environment=='LTEE')
## don't bother with MAE comparison, since Jeff Barrick says the rate and spectrum is
## weird due to growth on plates.

Ara.minus.3.IS150 <- filter(LTEE.IS150,Population=='Ara-3')

Ara.minus.3.IS150.by.clone <- group_by(Ara.minus.3.IS150,Clone,Generation,Environment,Population) %>%
summarize(total.count=n())

DM0.DM25.over.time <- group_by(IS.insertions,Clone,Generation,Environment,Population) %>%
summarize(total.count=n())

IS150.rate.df <- rbind(DM0.DM25.over.time,Ara.minus.3.IS150.by.clone)
IS150.rate.plot <- ggplot(IS150.rate.df,
                          aes(x=Generation,
                              y=total.count,color=Environment,shape=Environment)) +
theme_classic() +
scale_color_brewer(palette='YlGnBu',direction=-1) +
ylab(expression(paste('IS',italic('150'),' insertions'))) +
geom_jitter(width=50) +
guides(color=FALSE,shape=FALSE)

########
## Combine the IS plots with cowplot to make Figure 7.
Fig7outf <- file.path(projdir,"results/figures/Fig7.pdf")
Fig7 <- plot_grid(parallel.IS.plot, IS.plot, IS150.rate.plot,
                  labels = c('A', 'B', 'C'), ncol=1)
save_plot(Fig7outf,Fig7,base_height=7,base_aspect_ratio=0.8)

########
## Conduct the following test for parallel evolution of IS-insertions:
## Assume all IS insertions occur in the set union(A,B), with
## probability = empirical mass distribution over LTEE/MAE/DM0/DM25.
## Set the number of samples to the number of IS150 insertions
## in the DM0 genomes. For 10000 replicates,
## count how often some site gets hit 9 times. I could also repeat this test for the
## chance that a particular site gets hit 9 times (for instance).

## filter out mutations in clones on the same lineage
## by removing duplicates within the same population.
LTEE.MAE.IS150.hit.pos <- LTEE.MAE.IS.insertions %>%
dplyr::select(-Clone) %>% distinct() %>%
group_by(Position,GenePosition) %>%
summarize(count=n()) %>% ungroup() %>%
arrange(desc(count))

## again, filter out duplicates on line of descent
## by removing duplicates within the same population.
DM0.DM25.IS150.hit.pos <- IS.insertions %>%
filter(IS_element == 'IS150') %>% mutate(Position=genome_start) %>%
dplyr::select(Position,Environment,Population,GeneName,GenePosition) %>%
distinct() %>%
group_by(Position,GenePosition) %>%
summarize(count=n()) %>% ungroup() %>%
arrange(desc(count))

## use GenePosition to match up IS element insertion positions.
## number of DM0/DM25 IS150 events with sites in LTEE/MAE: 22.
DM0.DM25.pos.intersect <- filter(DM0.DM25.IS150.hit.pos,
                                 GenePosition %in% LTEE.MAE.IS150.hit.pos$GenePosition)
## number of DM0/DM25 IS150 events in total: 124. 22/124 = 0.1774 in LTEE/MAE set.
length(DM0.DM25.IS150.hit.pos$GenePosition)

## make an empirical mass function over GenePosition.
total.hit.pos <- rbind(
  dplyr::select(DM0.DM25.IS150.hit.pos,-Position),
  dplyr::select(LTEE.MAE.IS150.hit.pos,-Position)) %>%
group_by(GenePosition) %>% summarize(count2 = sum(count)) %>%
arrange(desc(count2)) %>% mutate(prob=count2/sum(count2)) %>%
mutate(eCDF=cumsum(prob))

## number of draws for each simulation = IS150 in DM0.
DM0.IS150 <- IS.insertions %>% filter(IS_element == 'IS150',Environment=='DM0') %>%
dplyr::select(-Clone) %>% distinct()
draws <- nrow(DM0.IS150)

null.parallel.hits <- function(total.hit.pos,empirical.parallel=9,replicates=10000) {

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

## empirical p-value for 9 hits is ~ 0.0128.
null.parallel.hits(total.hit.pos,replicates=100000)

## NOTE that this test is really stringent. 9 IS150 insertions in DM0 at position
## 3501576 in the yhiO promoter, 1 in DM0 at position 3501577,
## and ZDBp901 has an insertion annotated as yhiO at position 3501352.
## So, 11/12 have IS150 insertions annotated as yhiO, and 9/12 are at the exact same site.

## empirical p-value for 8 hits is ~ 0.04.
null.parallel.hits(total.hit.pos,empirical.parallel=8)

################################################################################
## Fitness and growth effects of plasmid-borne maeA expression.
maeA.fitness.analysis <- function(data, samplesize=6, days.competition=1,rev=FALSE) {
  ## by diluting stationary phase culture 1:100 on day 0, there is another
  ## factor of 100 that multiplies the day 0 plating dilution factor.
  data$D.0 <- 100*data$D.0
  data$Mr <- log2((data$D.1/data$D.0)*data$Red.1*100^days.competition/data$Red.0)/days.competition
  data$Mw <- log2((data$D.1/data$D.0)*data$White.1*100^days.competition/data$White.0)/days.competition
  ## reverse log ratio when polarity is reversed.
  if (rev) data$W <- data$Mr/data$Mw else data$W <- data$Mw/data$Mr

  my.mean <- mean(data$W, na.rm=T)
  my.sd <- sd(data$W, na.rm=T)
  ## see: https://www.cyclismo.org/tutorial/R/confidence.html#calculating-a-confidence-interval-from-a-t-distribution
  my.t.dist.error <- qt(0.975,df=samplesize-1) * my.sd/sqrt(samplesize)
  
  my.confint <- c(my.mean - my.t.dist.error, my.mean + my.t.dist.error)

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

maeA.fitness.data <- read.csv(file.path(projdir,"data/rohan-formatted/DM0_Fitness2.csv"),
                              header=TRUE,stringsAsFactors=FALSE)
## these data come from one day competitions that Tanush ran.
res1 <- filter(maeA.fitness.data,Red.Pop=='ZDB151_with_maeA') %>% maeA.fitness.analysis()
res2 <- filter(maeA.fitness.data,Red.Pop=='ZDB67_with_maeA') %>% maeA.fitness.analysis(rev=TRUE)
res3 <- filter(maeA.fitness.data,Red.Pop=='ZDB152_with_maeA') %>% maeA.fitness.analysis()
res4 <- filter(maeA.fitness.data,Red.Pop=='ZDB68_with_maeA') %>% maeA.fitness.analysis(rev=TRUE)

## beautiful! All confints overlap with each other, showing no maeA fitness effect
## does not depend on Ara polarity or genetic background.

## let's combine data from the switched polarity competitions (since Ara marker is neutral)
d1 <- maeA.fitness.data %>%
filter(Red.Pop=='ZDB151_with_maeA') %>%
dplyr::select(Red.0,Red.1,White.0,White.1,D.0,D.1)

d2 <- filter(maeA.fitness.data,Red.Pop=='ZDB67_with_maeA') %>% dplyr::select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## now switch column labels.
d2X <- dplyr::rename(d2,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB151.data <- rbind(d1,d2X)

d3 <- filter(maeA.fitness.data,Red.Pop=='ZDB152_with_maeA') %>% dplyr::select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d4 <- filter(maeA.fitness.data,Red.Pop=='ZDB68_with_maeA') %>% dplyr::select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## switch column labels.
d4X <- dplyr::rename(d4,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB152.data <- rbind(d3,d4X)

fres1 <- maeA.fitness.analysis(ZDB151.data,samplesize=6)
fres2 <- maeA.fitness.analysis(ZDB152.data,samplesize=6)
maeA.fitness.plot.data <- rbind(fres1,fres2)
#' Correct the strain names here.
maeA.fitness.plot.data$Strain <- c('CZB151','CZB152')

## Make a figure of the confidence interval.
maeA.fitness.fig <- "../results/figures/maeAFitness.pdf"

plot.maeA.fitness <- function (results, output.file) {
  the.plot <- ggplot(results,aes(x=Strain,y=Fitness)) +
  geom_errorbar(aes(ymin=Left,ymax=Right),width=0.1, size=1) +
  geom_line() +
  geom_point(size=2) +
  scale_y_continuous(limits=c(1.0,1.30)) +
  ylab("Fitness of maeA plasmid relative to empty plasmid") +
  theme_classic()
  ggsave(the.plot, file=output.file,width=4,height=4)
}

plot.maeA.fitness(maeA.fitness.plot.data,maeA.fitness.fig)

########################################################
