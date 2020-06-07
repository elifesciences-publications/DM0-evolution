## make-figures.R by Rohan Maddamsetti.

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

calc.bootstrap.conf.int <- function(vec) {
  ## bootstrap confidence intervals around the mean.
  ## Use the boot package to calculate fancy BCA intervals.
  
  mean.boot <- function(x,ind) {
      ## define this type of mean-calculating function to pass to the boot function.
      ## from: https://web.as.uky.edu/statistics/users/pbreheny/621/F12/notes/9-18.pdf
      ## as it seems this link has been taken down, I downloaded these notes
      ## using the Wayback Machine from the Internet Archive.
      return(c(mean(x[ind]), var(x[ind])/length(x)))
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
mutate(Mutation=as.factor(Mutation))


ltee.mutations <- read.csv(
  file.path(projdir,
            "data/rohan-formatted/nature18959-s4.csv"),
  stringsAsFactors=FALSE
  )

###############################################
## Figure 2.
## Figure 2A: make a stacked bar plot of the different kinds of mutations in each clone.
## Figure 2B: comparison to non-mutator LTEE 5000gen A clones.

plot.fig2.mutations.stackbar <- function(fig.df,panel,leg=FALSE) {
    if (panel == "DM0") {
        muts <- filter(fig.df,Environment=='DM0')
        my.title <- 'DM0 2,500 generations'
    } else if (panel == "DM25") {
        muts <- filter(fig.df,Environment=='DM25')
        my.title <- 'DM25 2,500 generations'
    } else if (panel == "LTEE-5K") {
        muts <- filter(fig.df,Environment=='LTEE-5K')
        my.title <- 'LTEE 5,000 generations'
    } else {
        stopifnot(TRUE == FALSE) ## panic if we get here.
    }
        
    fig <- ggplot(muts,aes(x=Clone, y=Count, fill=Mutation)) +
        geom_bar(stat='identity') +
        ylab("Count") +
        ylim(c(0,30)) +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              panel.border=element_blank(),
              plot.title=element_text(hjust=0, size = 12, face = "bold"),
              strip.text=element_text(hjust=0,size=12),
              panel.spacing.x=unit(1, "cm"),
              panel.spacing.y=unit(0.5, "cm")) +
    ggtitle(my.title)

    if (leg == TRUE) {
        fig <- fig +
            theme(legend.title=element_text(size=12, face="bold"),
                  legend.title.align=0.5,
                  legend.text=element_text(size=12))
    } else {
        fig <- fig + guides(fill=FALSE)
    }
    
    return(fig)
}

fig2.mutation.summary <- evolved.mutations %>%
    ## give nicer names for mutation classes.
    mutate(Mutation=recode(Mutation,
                  MOB = "Mobile element transposition",
                  DEL = "Indel",
                  INS = "Indel",
                  SUB = "Multiple-base substitution",
                  nonsynonymous = "Nonsynonymous",
                  synonymous = "Synonymous",
                  nonsense = "Nonsense",
                  pseudogene = "Pseudogene",
                  intergenic = "Intergenic",
                  )) %>%
    group_by(Clone, Environment, Population, Mutation) %>%
    summarize(Count=n()) %>%
    ungroup() %>%
    data.frame() %>%
    mutate(Mutation=as.character(Mutation))

ltee.5K.clone.mutation.summary <- ltee.mutations %>%
    filter(Generation==5000) %>%
    ## take odd rows to get A clones.
    filter(row_number() %% 2 == 1) %>%
    gather(key='Mutation',value='Count',
           IS.Element.Insertions,Small.Indels,Large.Deletions,
           Large.Duplications,Synonymous.Base.Substitutions,
           Nonsynonymous.Base.Substitutions,Intergenic.Base.Substitutions,
           Pseudogene.Base.Substitutions, Nonsense.Base.Substitutions) %>%
    mutate(Mutation=recode(Mutation,
                           IS.Element.Insertions= "Mobile element transposition",
                           Small.Indels="Indel",
                           Large.Deletions="Indel",
                           Large.Duplications="Large Duplication",
                           Synonymous.Base.Substitutions='Synonymous',
                           Nonsynonymous.Base.Substitutions = 'Nonsynonymous',
                           Intergenic.Base.Substitutions = 'Intergenic',
                           Pseudogene.Base.Substitutions = 'Pseudogene',
                           Nonsense.Base.Substitutions = 'Nonsense')) %>%
    ## to avoid confusion given later CNV analysis, omit Large Duplications here.
    filter(Mutation != "Large Duplication") %>%
    filter(!(Population %in% c('Ara-2','Ara+3'))) %>%
    mutate(Clone = Population) %>% select(-Total.Mutations,Base.Substitutions,
                                          -Total.Deleted.Base.Pairs,
                                          -Total.Inserted.Base.Pairs,
                                          -Estimated.Final.Genome.Size,
                                          -Noncoding.RNA.Base.Substitutions,
                                          -Base.Substitutions,
                                          -Strain.ID,
                                          -Generation) %>%
    mutate(Environment="LTEE-5K")

fig2.df <- rbind(fig2.mutation.summary, ltee.5K.clone.mutation.summary) %>%
    mutate(Mutation=as.factor(Mutation))
    
## make sure colors are the same across plots by setting levels.
fig2A <- plot.fig2.mutations.stackbar(fig2.df,panel="DM0")
fig2B <- plot.fig2.mutations.stackbar(fig2.df,panel="DM25")
fig2C <- plot.fig2.mutations.stackbar(fig2.df,panel="LTEE-5K",leg=TRUE)

## pull the legend from fig2C.
legend <- cowplot::get_legend(fig2C)

## remove the legend from fig1C.
fig2C <- fig2C + guides(fill=FALSE)

fig2 <- plot_grid(fig2A,fig2B,fig2C,legend,
                  labels=c('A','B','C',''),
                  nrow=2)

fig2.output <- "../results/figures/Fig2.pdf"
ggsave(fig2, file=fig2.output,width=7,height=7)

##############################################################
## Figures 3 and 4: Fitness competition results.
## See population-fitness.R and clone-fitness.R scripts.
##############################################################
## Figures 6 and 7: growth curve analysis, showing evolution of
## growth rate and time lags in DM0, and associated growth
## parameters when grown in DM25.

## let's examine REL606 growth in DM25 to estimate OD420.
## We need these data to provide an empirical basis for the intervals
## chosen to estimate r.glucose and r.citrate in later code.

## NOTE: OD is log-transformed using the natural log.
## HOWEVER: log-ratios of growth parameters (and fitness competitions for maeA)
## use log2 for easier interpretability.

## I double-checked the provenance of these data with Zack.
## These data come from a long spec run from back in 2011 for his 2012 Nature paper
## to look at growth phenotypes of Cit+ constructs and Cit+ plasmid transformants,
## with 606, early Cit+, and some evolved Cit- as points of comparison.
## Zack's original file is called "Construct run 4 final.xls"

REL606.DM25.growth.data <- read.csv(
  file.path(projdir, "data/rohan-formatted/REL606-DM25-48-hours.csv"),
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
    filter(Hours<=24) %>%
    ## don't include points before 30 minutes in.
    filter(Hours >= 0.5)

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

Fig5supp1 <- plot_grid(REL606.plot,log.REL606.plot,labels=c('A','B'), ncol=1)
save_plot(file.path(projdir,"results/figures/Fig5-supplement-1.pdf"), Fig5supp1, base_height=4, base_width=3)

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
        ## start using data at the 30 minutes mark.
        filter(Hours >= 0.5) %>%
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
                      geom_hline(yintercept=log(0.05),linetype='dashed',color='red',size=0.4) +
                      geom_hline(yintercept=log(0.1),linetype='dashed',color='red',size=0.4) +
                      ## clever trick to only plot glucose interval on DM25 curves.
                      geom_hline(data=data.frame(yint=log(0.01),Experiment='DM25'), aes(yintercept=yint),linetype='dashed',color='black',size=0.4) +
                      geom_hline(data=data.frame(yint=log(0.02),Experiment='DM25'), aes(yintercept=yint),linetype='dashed',color='black',size=0.4)
    
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
    ev.color <- "#3f007d"
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
    ## My simple growth rate calculation code.
    ## Use two different intervals for estimating r.glucose and r.citrate,
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
        if (unique(well.df$Experiment) == 'DM0-growth') {
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
                                    PopulationLabel = unique(well.df$PopulationLabel),
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
                                        PopulationLabel=rep(unique(df$PopulationLabel),5),
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

plot.growth.confints <- function(plot.df, confints.df, plot.CIs=TRUE, add.pop.name=FALSE) {

    recode.Parameter <- function(df) {
        df <- mutate(df,
                     Parameter=recode(
                         Parameter,
                         DM0.r.citrate = "DM0 citrate growth rate",
                         DM25.r.citrate = "DM25 citrate growth rate",
                         DM25.r.glucose = "DM25 glucose growth rate",
                         DM0.t.lag = "DM0 lag time",
                         DM25.t.lag = "DM25 lag time"
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
        my.title <- bquote(.(param)~(h)^{})
        ## add units of h^-1 to title if param is a rate parameter.
        if (endsWith(param,"rate")) {
            my.title <- bquote(.(param)~(h^-1))
        }
        
        if (add.pop.name) { ## for the clone figure, add the name of its pop.
            param.df <- param.df %>%
                ## replace NA values in PopulationLabel to "Founder".
                mutate(PopulationLabel=ifelse(is.na(PopulationLabel),
                                              "Founder",PopulationLabel)) %>%
                ## add the PopulationLabel to the Name for plotting.
                mutate(PlotName=paste(PopulationLabel,Name,sep=': ')) %>%
                ## order points from left to right, based on order for Name column.
                mutate(PlotName=fct_reorder(PlotName,Name, min))
            
            ## Also have to update the values for param.confint.df.
            param.confint.df <- param.confint.df %>%
                mutate(PopulationLabel=ifelse(is.na(PopulationLabel),
                                              "Founder",PopulationLabel)) %>%
                mutate(PlotName=paste(PopulationLabel,Name,sep=': '))
            
        } else { ## making the population figure.
            param.df <- param.df %>% mutate(PlotName = Name)
            param.confint.df <- param.confint.df %>% mutate(PlotName = Name)
        }
        
        fig <- ggplot(param.df, aes(x = PlotName,y = Estimate,color = Founder)) +
            geom_point(size=0.5) +
            ggtitle(my.title) +
            ylab("Estimate") +
            xlab("Sample") +
            guides(color=FALSE) +
            scale_color_manual(values = cbbPalette) +
            theme_classic() +
            theme(plot.title = element_text(size = 9.5),
                  axis.title.x = element_text(size=10),
                  axis.text.y = element_text(size=8),
                  plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                  axis.text.x  = element_text(angle=90, vjust=0.5,size=8)
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
        ## set bounds for rate and time lag parameter plots.
        if (str_detect(param,'r')) { ## plotting a rate.
            fig <- fig + ylim(c(0,1.5))
        } else { ## plotting time lag.
            if (str_detect(param,'DM25') & str_detect(param,'lag')) { ## for DM25 lag time
                fig <- fig + ylim(c(0,5))
            } else { ## for DM0 lag time.
                fig <- fig + ylim(c(0,25))
            }
        }
        
        list.of.plots[[i]] <- fig
        i <- i + 1
    }
    ## use cowplot to plot.
    my.plot <- plot_grid(plotlist=list.of.plots, nrow=1)
    return(my.plot)
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
  
  ## calculate the log2 ratio of evolved growth to ancestral growth for rate,
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

  for (index in seq_len(nrow(growth.summary))) {
    my.row <- growth.summary[index, ]
    my.anc <- filter(growth.summary,Name==my.row$Founder)

    anc.DM0.r.citrate[index] <- my.anc$DM0.r.citrate
    anc.DM0.t.lag[index] <- my.anc$DM0.t.lag
    anc.DM25.r.citrate[index] <- my.anc$DM25.r.citrate
    anc.DM25.r.glucose[index] <- my.anc$DM25.r.glucose
    anc.DM25.t.lag[index] <- my.anc$DM25.t.lag
  }

  growth.summary2 <- growth.summary %>%
  mutate(log.DM0.r.citrate.ratio=log2(DM0.r.citrate/anc.DM0.r.citrate)) %>%
  mutate(log.DM0.t.lag.ratio=log2(DM0.t.lag/anc.DM0.t.lag)) %>%
  mutate(log.DM25.r.citrate.ratio=log2(DM25.r.citrate/anc.DM25.r.citrate)) %>%
  mutate(log.DM25.r.glucose.ratio=log2(DM25.r.glucose/anc.DM25.r.glucose)) %>%
  mutate(log.DM25.t.lag.ratio=log2(DM25.t.lag/anc.DM25.t.lag))

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
  
  bootstrap.results <- data.frame(
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

plot.growth.parameters <- function(growth, plot.CIs=TRUE, add.pop.name=FALSE) {    
  growth.CI <- bootstrap.my.growth.confints(growth)
  growth.plot.df <- growth %>%
  gather(key="Parameter", value="Estimate",
         DM0.r.citrate,DM25.r.citrate,DM25.r.glucose,DM0.t.lag,DM25.t.lag)
  growth.plot <- plot.growth.confints(growth.plot.df, growth.CI, plot.CIs, add.pop.name)
  return(growth.plot)    
}

plot.parameter.log.ratios <- function(plot.df, confints.df, recode.param.func) {

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

plot.growth.parameter.log.ratios <- function(plot.df, confints.df) {

    recode.Growth.Parameter <- function(df) {
        df2 <- mutate(df,
                      Parameter=recode_factor(
                          Parameter,
                          log.DM0.r.citrate.ratio = "DM0 citrate growth rate",              
                          log.DM25.r.citrate.ratio = "DM25 citrate growth rate",
                          log.DM25.r.glucose.ratio = "DM25 glucose growth rate",
                          log.DM0.t.lag.ratio = "DM0 lag time",
                          log.DM25.t.lag.ratio = "DM25 lag time")
                      )                   
        return(df2)
    }
    
    return(plot.parameter.log.ratios(plot.df,confints.df,recode.Growth.Parameter))
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
    reshape.my.growth.df() %>%
    ## Do some finagling to use Population Labels in Fig. 2 and Supp. Fig. S11.
    ## Set the PopulationLabel of ancestral strains as its Name.
    mutate(PopulationLabel=ifelse(is.na(PopulationLabel),Name,PopulationLabel)) %>%
    ## Then change the Name to PopulationLabel for plotting.
    mutate(Name = PopulationLabel)

DM25.clone.growth <- calc.growth.rates(DM25.clone.growth.data) %>%
    mutate(Dataset='CloneGrowth') %>%
    reshape.my.growth.df()

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

## Fig6-supplement-1: Plot growth curves for DM0-evolved populations.
Fig6supp1 <- plot.growthcurve.figure(DM0.pop.growth.data, logscale=FALSE)
save_plot(file.path(projdir,"results/figures/Fig6-supplement-1.pdf"), Fig6supp1, base_height=7, base_width=11)

## Fig6-supplement-2: growth curves for DM0-evolved populations on log-scale.
Fig6supp2 <- plot.growthcurve.figure(DM0.pop.growth.data,logscale=TRUE)
save_plot(file.path(projdir,"results/figures/Fig6-supplement-2.pdf"), Fig6supp2, base_height=7, base_width=11)


## Fig7-supplement-1: Plot growth curves for DM0-evolved clones.
Fig7supp1 <- plot.growthcurve.figure(filter(DM0.clone.growth.data,Hours<=24), logscale=FALSE)
save_plot(file.path(projdir,"results/figures/Fig7-supplement-1.pdf"), Fig7supp1, base_height=7, base_width=11)

## Fig7-supplement-2: growth curves for DM0-evolved clones on log-scale.
Fig7supp2 <- plot.growthcurve.figure(DM0.clone.growth.data,logscale=TRUE)
save_plot(file.path(projdir,"results/figures/Fig7-supplement-2.pdf"), Fig7supp2, base_height=7, base_width=11)

## Fig8-supplement-1. Plot growth curves for DM25-evolved clones.
Fig8supp1 <- plot.growthcurve.figure(DM25.clone.growth.data, logscale=FALSE)
save_plot(file.path(projdir,"results/figures/Fig8-supplement-1.pdf"), Fig8supp1, base_height=7, base_width=11)

## Fig8-supplement-2. Plot growth curves for DM25-evolved clones on a log-scale.
Fig8supp2 <- plot.growthcurve.figure(DM25.clone.growth.data, logscale=TRUE)
save_plot(file.path(projdir,"results/figures/Fig8-supplement-2.pdf"), Fig8supp2, base_height=7, base_width=11)

## data filtering to plot data included in my growth rate estimation.
## This code is useful for debugging and making sure the rate estimation works right.
filtered.clone.growth.data <- filter.growth.data.on.analysis.domain(DM0.clone.growth.data)
filtered.pop.growth.data <- filter.growth.data.on.analysis.domain(DM0.pop.growth.data)
## filtered log-scale plots, again for debugging purposes.
filtered.log.clone.plot <- plot.growthcurve.figure(filtered.clone.growth.data,logscale=TRUE)
filtered.log.pop.plot <- plot.growthcurve.figure(filtered.pop.growth.data,logscale=TRUE)

################################################################################
## plot my growth estimates.
## plot estimates with confidence interval of mean.

####################### First plot the evolved populations.
## Figure 6: evolved population growth estimates.
Fig6outf <- file.path(projdir,"results/figures/Fig6.pdf")
pop.growth.plot <- plot.growth.parameters(pop.growth)

pop.growth.summary <- summarize.growth.results(pop.growth)
final.pop.growth.summary <- calc.growth.log.ratios(pop.growth.summary)

pop.bootstrap.CI <- run.growth.ratio.confint.bootstrapping(final.pop.growth.summary)
Fig6B.plot.df <- final.pop.growth.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
       log.DM25.t.lag.ratio, log.DM0.t.lag.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.
Fig6B <- plot.growth.parameter.log.ratios(Fig6B.plot.df, pop.bootstrap.CI)

Fig6 <- plot_grid(pop.growth.plot,Fig6B,labels=c('A','B'),ncol=1,
                  rel_heights = c(1, 1))
save_plot(Fig6outf,Fig6, base_width=11, base_height=7.5)

######################  Now plot the clones.
## Filter oddball Cit- ZDBp874 from the clones before estimation.
clone.growth <- filter(clone.growth,Name != 'ZDBp874')

## Figure 7.
Fig7outf <- file.path(projdir,"results/figures/Fig7.pdf")
clone.growth.plot <- plot.growth.parameters(clone.growth, add.pop.name=TRUE)

clone.growth.summary <- summarize.growth.results(clone.growth)
final.clone.growth.summary <- calc.growth.log.ratios(clone.growth.summary)
clone.bootstrap.CI <- run.growth.ratio.confint.bootstrapping(final.clone.growth.summary)

Fig7B.plot.df <- final.clone.growth.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
       log.DM25.t.lag.ratio, log.DM0.t.lag.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.

Fig7B <- plot.growth.parameter.log.ratios(Fig7B.plot.df, clone.bootstrap.CI)

Fig7 <- plot_grid(clone.growth.plot,Fig7B,labels=c('A','B'),ncol=1,
                  rel_heights = c(1, 1),rel_widths=c(1,1))
save_plot(Fig7outf,Fig7, base_width=11, base_height=7.5)

######################
## Figure 8:
## growth results for DM25-evolved clones in DM25.
Fig8outf <- file.path(projdir,"results/figures/Fig8.pdf")
DM25.clone.growth.plot <- plot.growth.parameters(DM25.clone.growth, add.pop.name=TRUE)
DM25.growth.summary <- summarize.growth.results(DM25.clone.growth)
final.DM25.growth.summary <- calc.growth.log.ratios(DM25.growth.summary)
DM25.clone.bootstrap.CI <- run.growth.ratio.confint.bootstrapping(final.DM25.growth.summary)

Fig8B.plot.df <- final.DM25.growth.summary %>%
gather(key="Parameter", value="Estimate",
       log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
       log.DM25.t.lag.ratio, log.DM0.t.lag.ratio) %>%
filter(Name != Founder) ## ancestors always have a ratio of zero.
Fig8B <- plot.growth.parameter.log.ratios(Fig8B.plot.df, DM25.clone.bootstrap.CI)

Fig8 <- plot_grid(DM25.clone.growth.plot, Fig8B, labels=c('A','B'),ncol=1,
                    rel_heights = c(1, 1))
save_plot(Fig8outf, Fig8, base_width=11, base_height=7.5)

####################################
### Figure 9:
### no correlation between glucose growth rates and citrate growth rates,
### but correlation in citrate growth rates in DM0 and DM25.

growth.summary <- full_join(clone.growth.summary, pop.growth.summary) %>%
    mutate(Dataset = case_when(Dataset == 'CloneGrowth' ~ 'Clone Growth',
                           TRUE ~ 'Population Growth'))

## no correlation between r citrate and r glucose in DM25 by my estimates.
cor.test(x=clone.growth.summary$DM25.r.citrate,y=clone.growth.summary$DM25.r.glucose,method="pearson", alternative="two.sided")
cor.test(x=pop.growth.summary$DM25.r.citrate,y=pop.growth.summary$DM25.r.glucose,method="pearson", alternative="two.sided")

## significant correlation between r citrate in DM25 and DM0 by my estimates.
cor.test(x=clone.growth.summary$DM25.r.citrate,y=clone.growth.summary$DM0.r.citrate,method="pearson",, alternative="two.sided")
cor.test(x=pop.growth.summary$DM25.r.citrate,y=pop.growth.summary$DM0.r.citrate,method="pearson", alternative="two.sided")

cit.glucose.cor.plot <- ggplot(growth.summary,
                               aes(x=DM25.r.citrate,y=DM25.r.glucose, color=Founder, shape=Generation)) +
facet_wrap(Founder~Dataset,nrow=1) +
geom_point() +
theme_classic() +
xlab("DM25 citrate growth rate") +
ylab("DM25 glucose growth rate") +
scale_color_manual(values=cbbPalette) +
guides(color=FALSE,shape=FALSE)

## Think more about correlation between growth on citrate in DM0 and growth on citrate in DM25.
cit.cit.cor.plot <- ggplot(growth.summary,
                           aes(x=DM25.r.citrate,y=DM0.r.citrate,color=Founder, shape=Generation)) +
facet_wrap(Founder~Dataset,nrow=1) +
geom_point() +
theme_classic() +
xlab("DM25 citrate growth rate") +
ylab("DM0 citrate growth rate") +
scale_color_manual(values=cbbPalette) +
guides(color=FALSE,shape=FALSE)

## Save Figure 9
Fig9 <- plot_grid(cit.glucose.cor.plot,
                    cit.cit.cor.plot, labels=c('A','B'), ncol=1)
Fig9.outf <- file.path(projdir,"results/figures/Fig9.pdf")
save_plot(Fig9.outf, Fig9, base_height=5, base_width=7.5)

######################################################################
####### Figure 10: Nkrumah's cell death results. See CellDeath.R script
#######           for analyses and figures.
######################################################################
####### Figure 11: make a matrix plot of genes with mutations in two or more clones.
#######
######################################################################
PlotMatrixFigure <- function(raw.matrix.file, amp.matrix.file,
                             ltee.matrix.file, ltee.50k.labels.file,
                             matrix.outfile, pop.clone.labels) {

    raw.matrix <- read.csv(raw.matrix.file, stringsAsFactors=FALSE)
    
    ## fix the names of the samples.
    names(raw.matrix) <- map_chr(
        names(raw.matrix),
        function(x) str_trunc(x,width=7,side="right",ellipsis=''))
    
    #' nice helper function
    not_any_na <- function(x) {!any(is.na(x))}
    
    #' import rows for the maeA and dctA amplifications.
    amp.matrix <- read.csv(amp.matrix.file, stringsAsFactors=FALSE)
    #' and merge with the mutation matrix.
    #' remove columns that contain any NA values (ZDBp875)
    merged.with.amps.matrix <- full_join(raw.matrix,amp.matrix) %>%
        dplyr::select_if(not_any_na)
    
    DM0.DM25.matrix.data <- gather(merged.with.amps.matrix,
                                   "Name",
                                   "mutation.count",
                                   2:ncol(merged.with.amps.matrix)) %>%
        left_join(pop.clone.labels) %>%
        dplyr::select(Gene,Name,mutation.count,Environment,PopulationLabel) %>%
        group_by(Gene) %>% filter(sum(mutation.count)>1)
    
    #' print out a table of parallelism across the DM0 and DM25 treatments.
    parallel.counts <- summarize(DM0.DM25.matrix.data,
                                 total.count=sum(mutation.count)) %>%
        arrange(desc(total.count))
    print('Parallel hits in genes')
    print(filter(parallel.counts,total.count>1),n=Inf)
    
    #' add LTEE mutation matrix to figure.
    ltee.matrix <- read.csv(ltee.matrix.file, stringsAsFactors=FALSE)
    
    #' fix the names of the samples.
    names(ltee.matrix) <- map_chr(
        names(ltee.matrix),
        function(x) str_trunc(x,width=8,side="left",ellipsis=''))
    
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
    
    ## use the Population Names instead of Clone Names for the figure.
    renamed.LTEE.clones <- mutate(non.mutators.and.ara.minus.3,
                                  Name=PopulationLabel)
    
    ## cheap hack to plot empty rows in the LTEE data.
    ## add some zeros for genes with no mutations in renamed.LTEE.clones
    ## that are in the DM0-evolved or DM25-evolved clones.
    ltee.pops <- unique(renamed.LTEE.clones$PopulationLabel)
    get.empty.ltee.pop.df <- function(pop) {
        empty.ltee.genes <- c("yhiO", "ndh", "fadR", "yobG", "ybjR", "ybbO",
                              "rpsV", "maeA-AMP", "dctA-AMP")
        return(data.frame(Gene = empty.ltee.genes,
                   Name = rep(pop, length(empty.ltee.genes)),
                   mutation.count = rep(0,length(empty.ltee.genes)),
                   Environment = rep("LTEE", length(empty.ltee.genes)),
                   PopulationLabel = rep(pop, length(empty.ltee.genes))))
    }
    zero.ltee.muts <- map_dfr(ltee.pops,.f=get.empty.ltee.pop.df)
    
    ## Now join LTEE data to the DM0 and DM25 Cit+ data.
    matrix.data <- bind_rows(DM0.DM25.matrix.data, renamed.LTEE.clones, zero.ltee.muts)
    
    ## sort genes by number of mutations in each row.
    ## this is overwritten by the alternate sorting method that follows,
    ## but keeping this snippet because it's useful to have.
    gene.hit.sort <- group_by(matrix.data,Gene) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(hits))
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    
    ## alternate sorting method: difference in hits between environments.
    DM25mut.count <- filter(matrix.data,Environment=='DM25') %>% group_by(Gene) %>%
        summarize(DM25mut.count=sum(mutation.count))
    
    DM0mut.count <- filter(matrix.data,Environment=='DM0') %>% group_by(Gene) %>%
        summarize(DM0mut.count=sum(mutation.count))
    
    env.mut.sort <- inner_join(DM0mut.count,DM25mut.count) %>%
        mutate(mut.diff=abs(DM0mut.count - DM25mut.count)) %>%
        arrange(desc(mut.diff))
    
    ## now use these calculations to sort genes by the absolute value of the
    ## difference in number of mutations between DM0 and DM25 treatments.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(env.mut.sort$Gene))
    
    
    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)
    
    make.matrix.panel <- function(mdata,environment,leg=FALSE) {
        fig <- ggplot(filter(mdata,Environment==environment),
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
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_manual(name="Mutations",
                              values = c("white", "#ffdf00", "#bebada",
                                         "#fb8072", "#80b1d3", "#fdb462"))

        if (leg == FALSE) {
            fig <- fig + guides(fill=FALSE)
        }
        return(fig)
    }
    
    DM0.matrix.panel <- make.matrix.panel(matrix.data,'DM0')
    DM25.matrix.panel <- make.matrix.panel(matrix.data,'DM25')
    LTEE.matrix.panel <- make.matrix.panel(matrix.data,'LTEE',leg=TRUE)

    ## pull the legend from the LTEE panel.
    legend <- cowplot::get_legend(LTEE.matrix.panel)
    ## remove the legend from the LTEE panel.
    LTEE.matrix.panel <- LTEE.matrix.panel + guides(fill=FALSE)
    
    matrix.figure <- plot_grid(DM0.matrix.panel,
                               DM25.matrix.panel,
                               LTEE.matrix.panel,
                               nrow=1,
                               rel_widths=c(13/12,1,7/12),
                               align = 'vh')
    
    ggsave(matrix.figure, file=matrix.outfile,width=10,height=10)
}

##### Now run this function on all valid mutations to make Figure 5.
raw.matrix.f <- file.path(projdir,
                          "results/DM0-DM25-comparison-mut-matrix.csv")

amp.matrix.f <- file.path(projdir,"results/amp_matrix.csv")
ltee.matrix.f <- file.path(projdir,"results/LTEE-mut_matrix.csv")
ltee.50k.labels.f <- file.path(projdir, "data/rohan-formatted/LTEE-50K-clones.csv")
matrix.outf <- "../results/figures/Fig11.pdf"

PlotMatrixFigure(raw.matrix.f, amp.matrix.f, ltee.matrix.f, ltee.50k.labels.f, matrix.outf, pop.clone.labels)

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

## examine parallel evolution at amino acid level.
parallel.dN <- evolved.mutations %>% filter(Mutation=='nonsynonymous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)
parallel.dN.Table <- filter(evolved.mutations, Position %in% parallel.dN$Position) %>%
arrange(Position)

parallel.dN.Table

## 3/5 parallel bp mutations are in citrate synthase, gltA!
## 735765, 735797, 735941.
## These are M172I, A162T, I114F, respectively.
## A162V is reported in Quandt et al. (2015) in eLife as affecting NADH-binding.
## The other two have not been reported before.
## fine-tuning has been reported by Erik Quandt in that paper.

## From methods section of Quandt et al. (2015):
## Dimeric models of ligand-free and NADH-bound E. coli citrate synthase (PDB: 1NXE and 1NXG) (Maurus et al., 2003)
## were prepared for analysis by reverting alanine-383 to phenylalanine and processing with the Structure
## Preparation application within MOE (Molecular Operating Environment).

## the other parallel bp mutations are:
## 2209801 is in atoS (S351C), 2630053 is in ygaF (I197L).

## IMPORTANT NOTE: the parallel mutation in atoS is missing in the analysis here,
## since it occurs in ZDBp871 and ZDBp874, which is the oddball Cit- clone which
## is excluded from the genome analysis for simplicity.
## This case of parallelism can be seen by
## inspection of the breseq HTML output for these genomes.

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

## There are several mutations associated with the fad regulon-- all of which occur in
## the DM0 treatment.
fad.mutations <- filter(evolved.mutations,!is.na(str_extract(Gene,'fad'))) %>% arrange(Position)
fad.mutations

################################################
## Figure 12. IS element analysis and visualization.
################################################

##IS.palette <- c('#f4a582','#92c5de','#ca0020','black','#76ee00')
IS.palette <- c('#f4a582','#9970ab','#ca0020','black','#76ee00')
IS.insertions <- read.csv(file.path(projdir,
                                    "results/genome-analysis/IS_insertions.csv")) %>%
    arrange(genome_start) %>%
    ## filter oddball Cit- clone ZDBp874.
    filter(Clone != 'ZDBp874') %>%
    mutate(genome_start_MB=genome_start/1000000)

IS.plot <- ggplot(IS.insertions,aes(x=genome_start_MB,fill=IS_element,frame=Environment)) +
    facet_grid(Environment~.) +
    geom_histogram(bins=200) +
    ##guides(fill=FALSE) +
    scale_fill_manual(values=IS.palette, name="IS element",
                      labels = c(expression(paste('IS',italic('1'))),
                                 expression(paste('IS',italic('150'))),
                                 expression(paste('IS',italic('186'))),
                                 expression(paste('IS',italic('3'))),
                                 expression(paste('IS',italic('RSO11'))))
                      ) +
    ylab("Count") +
    xlab("Position (Mb)") +
    theme_classic(base_family="Helvetica")

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
                     mutate(annotation=ifelse(is.na(annotation),'intergenic',annotation))
}

parallel.IS.summary <- summarize.IS.func(parallel.IS.insertions) %>%
    mutate(genome_start_MB=genome_start/1000000)

## plot parallel IS insertions and their annotation.
parallel.IS.plot <- ggplot(parallel.IS.summary,aes(x=genome_start_MB,
                                                   y=count,
                                                   color=IS_element,
                                                   label=annotation)) +
    geom_point() +
    scale_color_manual(values = IS.palette, name = "IS element",
                       labels = c(expression(paste('IS',italic('1'))),
                                  expression(paste('IS',italic('150'))),
                                  expression(paste('IS',italic('186'))),
                                  expression(paste('IS',italic('3'))),
                                  expression(paste('IS',italic('RSO11'))))
                       ) +
theme_classic() +
ylab("Count") +
xlab("Position (Mb)") +
geom_text_repel(fontface = "italic",size=3,show.legend=FALSE)

## Now just consider IS parallelism in the DM0 environment.
parallel.DM0.IS.insertions <- IS.insertions %>%
filter(Environment == 'DM0') %>%
group_by(genome_start) %>%
filter(n()>1)

parallel.DM0.IS.summary <- summarize.IS.func(parallel.DM0.IS.insertions) 

#####################################
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

IS150.insertions <- filter(IS.insertions,IS_element=='IS150')

DM0.DM25.IS150.over.time <- group_by(IS150.insertions,Clone,Generation,Environment,Population) %>%
summarize(total.count=n())

## Calculate Mann-Whitney U-test to compare number of IS150 between DM0 and DM25 treatments.
wilcox.test(total.count ~ Environment, data=DM0.DM25.IS150.over.time)
## significant difference between DM0 and DM25: p = 0.0089.

IS150.rate.df <- rbind(DM0.DM25.IS150.over.time,Ara.minus.3.IS150.by.clone)
IS150.rate.plot <- ggplot(IS150.rate.df,
                          aes(x=Generation,
                              y=total.count,color=Environment,shape=Environment)) +
theme_classic() +
scale_color_manual(values = rev(c("#F1BB7B", "#FD6467", "#5B1A18"))) +
    ylab(expression(paste('IS',italic('150'),' insertions'))) +
    ##guides(color=FALSE,shape=FALSE)
geom_jitter(width=50)


########
## Combine the IS plots with cowplot to make Figure 12.
Fig12outf <- file.path(projdir,"results/figures/Fig12.pdf")
Fig12 <- plot_grid(parallel.IS.plot, IS.plot, IS150.rate.plot,
                  labels = c('A', 'B', 'C'), ncol=1)
save_plot(Fig12outf,Fig12, base_height=7, base_aspect_ratio=1)

## save subpanels to include in talks.
save_plot(file.path(projdir,"results/figures/Fig12A.pdf"),parallel.IS.plot)
save_plot(file.path(projdir,"results/figures/Fig12B.pdf"),IS.plot)
save_plot(file.path(projdir,"results/figures/Fig12C.pdf"),IS150.rate.plot)

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

null.parallel.hits <- function(total.hit.pos,empirical.parallel=9,replicates=100000) {
    
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

## empirical p-value for 9 hits is ~ 0.013.
null.parallel.hits(total.hit.pos)

## NOTE that this test is really stringent. 9 IS150 insertions in DM0 at position
## 3501576 in the yhiO promoter, 1 in DM0 at position 3501577,
## and ZDBp901 has an insertion annotated as yhiO at position 3501352.
## So, 11/12 have IS150 insertions annotated as yhiO, and 9/12 are at the exact same site.

## empirical p-value for 8 hits is ~ 0.04.
null.parallel.hits(total.hit.pos,empirical.parallel=8)

##################
## by inspecting evolved.mutations in conjunction with amplified genes, we see that the
## dctA amplification anti-correlates with the promoter mutation.
## Note that the dctA promoter mutation is in CZB151 and CZB154.

## ALSO: Zack noticed an apparent anti-correlation between the dctA
## amplification and dcuS contingency locus mutation. Calculate all
## three pairwise correlations
## (dctA promoter, dctA amplification, dcuS mutation).

dctA.AMPs <- read.csv(file.path(projdir,"results/amplified_genes.csv"),
                      stringsAsFactors=FALSE) %>%
    filter(gene=='dctA') %>% dplyr::select(Genome) %>% distinct() %>%
    transmute(Name=Genome) %>% left_join(pop.clone.labels)

promoter.mutant <- filter(evolved.mutations, Gene=='dctA/yhjK') %>%
    transmute(Name=Clone) %>% left_join(pop.clone.labels) %>%
    full_join(filter(pop.clone.labels,(Founder %in% c('CZB151','CZB154')) & (Generation==2500) & (SampleType=='Clone') & (Sequenced == 1) & (Name != 'ZDBp874')))

dcuS.mutant <- filter(evolved.mutations, Gene=='dcuS') %>%
    transmute(Name=Clone) %>% left_join(pop.clone.labels) %>%
    full_join(filter(pop.clone.labels,(Founder == 'CZB151') & (Generation==2500) & (SampleType=='Clone') & (Sequenced == 1) & (Name != 'ZDBp874')))

## out of 24 DM0- and DM25-evolved genomes,
## 12 have dcuS mutation. All 12 of these have the dctA promoter mutation.
dcuS.mut.promoter.mut.intersection <- inner_join(dcuS.mutant, promoter.mutant)
## 0 have the dcuS mutation and no dctA promoter mutation.
dcuS.mut.no.promoter.mut <- anti_join(dcuS.mutant, promoter.mutant)
## 5 have the dctA promoter mutation and no dcuS mutation.
promoter.mut.no.dcuS.mut <- anti_join(promoter.mutant, dcuS.mutant)
## 7 have neither the dcuS mutation nor the dctA promoter mutation.
neither.dcuS.mut.nor.promoter <- filter(pop.clone.labels,(Founder %in% c('CZB152','CZB154')) & (Generation==2500) & (SampleType=='Clone') & (Sequenced==1)) %>%
    anti_join(dcuS.mutant) %>%
    anti_join(promoter.mutant)
fisher.test(matrix(c(12,0,5,7),2)) ## p = 0.0046.

## out of 24 DM0- and DM25-evolved genomes,
## 0 have both dcuS mutation and dctA AMP.
dctA.AMP.dcuS.mut.intersection <- inner_join(dctA.AMPs, dcuS.mutant)
## 6 have the dctA amplication and no dcuS mutation.
dctA.AMP.no.dcuS.mut <- anti_join(dctA.AMPs, dcuS.mutant)
## 12 have the dcuS mutation and no dctA amplification.
dcuS.mut.no.dctA.AMP <- anti_join(dcuS.mutant, dctA.AMPs)
## 6 have neither the dcuS mutation nor the dctA amplification.
neither.dctA.AMP.nor.dcuS.mut <- filter(pop.clone.labels,(Founder %in% c('CZB152','CZB154')) & (Generation==2500) & (SampleType=='Clone') & (Sequenced==1)) %>%
    anti_join(dctA.AMPs) %>%
    anti_join(dcuS.mutant)
fisher.test(matrix(c(0,6,12,6),2)) ## p = 0.0137.

## out of 24 DM0- and DM25-evolved genomes,
## 1 has both promoter mutant and dctA AMP.
dctA.AMP.promoter.mut.intersection <- inner_join(dctA.AMPs, promoter.mutant)
## 5 have dctA AMP but lack the promoter mutation.
dctA.AMP.no.promoter.mut <- anti_join(dctA.AMPs, promoter.mutant)
## 16 have the promoter mutation, and not the dctA AMP.
promoter.mut.no.dctA.AMP <- anti_join(promoter.mutant,dctA.AMPs)
## 2 do not have either mutation.
neither.dctA.AMP.nor.promoter <- filter(pop.clone.labels,Founder == 'CZB152' & Generation==2500 & SampleType=='Clone' & Sequenced==1) %>%
    anti_join(dctA.AMPs) %>%
    anti_join(promoter.mutant)

## Fisher's exact test: p = 0.0027.
fisher.test(matrix(c(1,5,16,2),2))
fisher.test(matrix(c(1,16,5,2),2))

## For clarity of exposition in the paper, thinking about the promoter
## mutation in CZB151 and 154, while one of the CZB152-evolved clones ZDBp877)
## also evolved a promoter mutation is complicated to explain (in addition)
## to keeping track here.
## Rich suggested comparing the anticorrelation between CZB151&154 versus 152
## and the dctA amplification instead, so that we have fewer things to keep in mind
## and explain.

## only one of the 16 clones descended from CZB151 and 154 have the dctA AMP.
CZB151.and.154.dctA.AMPs <- filter(dctA.AMPs, Founder %in% c('CZB151', 'CZB154'))
## while 5 out of the 8 clones descended from CZB152 has the dctA AMP.
CZB152.dctA.AMPs <- filter(dctA.AMPs, Founder=='CZB152')
## Fisher's exact test: p = 0.0069.
fisher.test(matrix(c(1,15,5,3),2))
## Of course, one of the CZB152-evolved strains without the dctA AMP
## has the promoter mutation! So further evidence of the the anti-correlation.
## This is the statistical test for this result that is reported in the paper.

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

maeA.fitness.data <- read.csv(file.path(projdir,"data/rohan-formatted/DM0-maeA-plasmid-fitness.csv"),
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

plot.maeA.fitness <- function(results, output.file) {
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
