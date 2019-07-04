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
library(scales)
library(gridExtra)
library(lubridate)
library(assertthat)
library(tidytext)
library(Matrix)
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

    final <- c(Left, Right)
    return(final)
}

## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###################

home.dir <- path.expand("~")
proj.dir <- file.path(home.dir,"BoxSync/active-projects/DM0-evolution")

pop.clone.labels <- read.csv(
    file.path(proj.dir,
              "data/rohan-formatted/populations-and-clones.csv"),
    stringsAsFactors=FALSE)

evolved.mutations <- read.csv(
    file.path(proj.dir,
              "results/genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE)

ltee.mutations <- read.csv(
    file.path(proj.dir,
              "data/rohan-formatted/nature18959-s4.csv"),
    stringsAsFactors=FALSE
)

## by inspecting evolved.mutations in conjunction with amplified genes, we see that the
## dctA amplification anti-correlates with the promoter mutation.
## Note that the dctA promoter mutation is in CZB151 and CZB154.
promoter.mutant <- filter(evolved.mutations, Gene=='dctA/yhjK') %>%
    transmute(Name=Clone) %>% left_join(pop.clone.labels)
dctA.AMPs <- read.csv(file.path(proj.dir,"results/amplified_genes.csv"),
                      stringsAsFactors=FALSE) %>%
    filter(gene=='dctA') %>% select(Genome) %>% distinct() %>%
    transmute(Name=Genome) %>% left_join(pop.clone.labels)
## 5/6 dctA AMP genomes lack the promoter mutation.
## 1/17 with promoter mutation have the AMP.
## 5/7 without promoter mutation have the AMP.
## Fisher's exact test: p = 0.0027.
fisher.test(matrix(c(1,5,16,2),2))

###############################################
## Figure 1.

## Figure 1A is a schematic of the experimental design, made in Illustrator.
## Figure 1B is the phylogenetic tree constructed using the iTOL webserver.
## Figure 1C: make a stacked bar plot of the different kinds of mutations in each clone.
## Figure 1D: comparison to non-mutator LTEE 5000gen and 30000gen A clones.

fig1C <- ggplot(evolved.mutations,aes(x=Clone,fill=Mutation)) +
    geom_bar() +
    ylab("Count") +
    ylim(c(0,80)) +
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
ggsave(fig1C, file=fig1C.output,width=7,height=7)

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
    filter(!(Population %in% c('Ara-2','Ara+3')))

ltee.30K.clones <- ltee.mutations %>% filter(Generation==30000) %>%
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
    filter(!(Population %in% c('Ara-2','Ara+3','Ara-1','Ara-4','Ara+6')))

fig1D.ltee.data <- rbind(ltee.5K.clones,ltee.30K.clones,stringsAsFactors=FALSE) %>%
    mutate(Generation=as.factor(Generation)) %>%
    mutate(Generation=recode_factor(Generation,'5000'='LTEE 5,000 generations','30000'='LTEE 30,000 generations'))

fig1D <- ggplot(fig1D.ltee.data,aes(x=Population,y=Count,fill=Mutation)) +
    geom_bar(stat='identity') +
    scale_fill_viridis(option="magma",discrete=TRUE) +
        facet_grid(~Generation, scales="free") +
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

fig1D.output <- "../results/figures/Fig1D.pdf"
ggsave(fig1D, file=fig1D.output,width=7,height=7)

##############################################################
## Figure 2: growth curve analysis, showing evolution of
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
    file.path(proj.dir, "data/rohan-formatted/REL606-DM25-96-hours.csv"),
    stringsAsFactors=FALSE) %>%
    mutate(Hours = as.numeric(as.duration(hms(Time)))/3600) %>%
    rowwise() %>%
    ## average the replicate blanks.
    mutate(Blank = lift_vd(mean)(DM25.1,DM25.2,DM25.3,DM25.4,DM25.5,DM25.6)) %>%
    select(-DM25.1,-DM25.2,-DM25.3,-DM25.4,-DM25.5,-DM25.6) %>%
    ## subtract Blank from well readings.
    mutate(REL606.1 = REL606.1 - Blank) %>%
    mutate(REL606.2 = REL606.2 - Blank) %>%
    mutate(REL606.3 = REL606.3 - Blank) %>%
    mutate(REL606.4 = REL606.4 - Blank) %>%
    mutate(REL606.5 = REL606.5 - Blank) %>%
    ## zero out negative values.
    mutate_at(vars(REL606.1:REL606.5),function(x)ifelse(x<0,0,x)) %>%
    select(-Blank) %>%
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
            ylab('log(OD420)') +
            geom_hline(yintercept=log(glu.ceiling),linetype='dashed',color='red') +
            geom_hline(yintercept=log(glu.interval.high),linetype='dashed',color='black') +
            geom_hline(yintercept=log(glu.interval.low),linetype='dashed',color='black')

    } else {
        fig <- ggplot(REL606.df,
                      aes(x=Hours,
                          y=OD420)) +
            ylab('OD420') +
            geom_hline(yintercept=glu.ceiling,linetype='dashed',color='red') +
            geom_hline(yintercept=glu.interval.high,linetype='dashed',color='black') +
            geom_hline(yintercept=glu.interval.low,linetype='dashed',color='black')

    }
    title.string <- "REL606 in DM25"

    fig <- fig +
        geom_point(size=0.1, alpha=0.1) +
    theme_classic() +
        guides(color=FALSE) +
        ggtitle(title.string) +
        theme(plot.title = element_text(size = 12, face = "bold"))

    return(fig)
}

REL606.plot <- plot.REL606.DM25.growth(REL606.DM25.growth.data, logscale=FALSE)
log.REL606.plot <- plot.REL606.DM25.growth(filter(REL606.DM25.growth.data,Hours<24), logscale=TRUE)

S1.Fig <- plot_grid(REL606.plot,log.REL606.plot,labels=c('A','B'),ncol=1)

save_plot(file.path(proj.dir,"results/figures/S1Fig.pdf"),S1.Fig,base_height=4,base_width=3)

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
            ylab('log(OD420)') +
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
            ylab('OD420')
    }

    ## title for the plot.
    ev.pop <- unique(filter(plot.growth.data,Name==ev.name)$PopulationLabel)
    sample.type <- unique(filter(plot.growth.data,Name==ev.name)$SampleType)
    if (sample.type == 'Clone') {
        title.string <- paste0(ev.pop,': ', ev.name)
    } else {
        title.string <- ev.pop
    }

    ## make Cit- strain stand out.
    if (ev.name == 'ZDBp874') {
        ev.color <- "#E69F00"
    } else {
        ev.color <- "#009E73"
    }

    fig <- fig +
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
        ## to compare time lags, measure point when OD420 hits 0.05.
        t.OD.hit.cit.min <- well.df[cit.min.index,]$Hours

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
            rel.subdata <- select(subdata,c('log.OD420','Hours')) %>%
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
                                    t.lag = t.OD.hit.cit.min,
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
                                            mean(df$DM0.r.citrate),
                                            mean(df$DM25.r.citrate),
                                            mean(df$DM25.r.glucose),
                                            mean(df$DM0.t.lag),
                                            mean(df$DM25.t.lag)),
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
                                            mean(df$DM0.t_mid),
                                            mean(df$DM25.t_mid)
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

plot.growth.confints <- function (plot.df, confints.df) {

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

    plot.df <- recode.Parameter(plot.df) %>% drop_na()
    confints.df <- recode.Parameter(confints.df) %>% drop_na()

    ggplot(plot.df, aes(x = Name,y = Estimate,color = Founder)) +
        geom_point(size=0.5) +
        ylab("Estimate") +
        xlab("Sequenced Isolate") +
        guides(color=FALSE) +
        scale_color_manual(values = cbbPalette) +
        facet_wrap( Parameter ~ . , scales="free",nrow=1) +
    geom_errorbar(data = confints.df, aes(ymin=Left,ymax=Right), width=1, size=0.5) +
    theme(axis.text.x  = element_text(angle=90, vjust=-0.5, size=5))

}

reshape.growth.curve.fits <- function(growth.curve.fits) {
    ## do some gymnastics to plot rates in different media conditions together.
    DM0.growth.curve.fits <- filter(growth.curve.fits,
                                    Experiment=='DM0-growth') %>%
        mutate(DM0.r = r) %>%
        mutate(DM0.t_mid = t_mid) %>%
        select(-r,-t_mid,-Experiment)

    DM25.growth.curve.fits <- filter(growth.curve.fits,
                                     Experiment=='DM25-growth') %>%
        mutate(DM25.r = r) %>%
        mutate(DM25.t_mid = t_mid) %>%
        select(-r,-t_mid,-Experiment)

    reshaped.growth.curve.fits <- full_join(DM0.growth.curve.fits,
                                            DM25.growth.curve.fits)
    return(reshaped.growth.curve.fits)
}

reshape.my.growth.df <- function(growth) {
    ## do some gymnastics to plot rates in different media conditions together.
    DM0.growth <- filter(growth, Experiment=='DM0-growth') %>%
        mutate(DM0.r.citrate=r.citrate) %>%
        mutate(DM0.t.lag=t.lag) %>%
        select(-r.glucose,-r.citrate,-t.lag,-Experiment)

    DM25.growth <- filter(growth, Experiment=='DM25-growth') %>%
        mutate(DM25.r.citrate=r.citrate) %>%
        mutate(DM25.r.glucose=r.glucose) %>%
        mutate(DM25.t.lag=t.lag) %>%
        select(-r.glucose,-r.citrate,-t.lag,-Experiment)

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

    ## my methods' estimates
    anc.DM0.r.citrate <- rep(-1,nrow(growth.summary))
    anc.DM0.t.lag <- rep(-1,nrow(growth.summary))

    anc.DM25.r.citrate <- rep(-1,nrow(growth.summary))
    anc.DM25.r.glucose <- rep(-1,nrow(growth.summary))
    anc.DM25.t.lag <- rep(-1,nrow(growth.summary))

    for (index in 1:nrow(growth.summary)) {
        my.row <- growth.summary[index, ]
        my.anc <- filter(growth.summary,Name==my.row$Founder)

        anc.DM0.r[index] <- my.anc$DM0.r
        anc.DM0.t_mid[index] <- my.anc$DM0.t_mid
        anc.DM25.r[index] <- my.anc$DM25.r
        anc.DM25.t_mid[index] <- my.anc$DM25.t_mid

        anc.DM0.r.citrate[index] <- my.anc$DM0.r.citrate
        anc.DM0.t.lag[index] <- my.anc$DM0.t.lag
        anc.DM25.r.citrate[index] <- my.anc$DM25.r.citrate
        anc.DM25.r.glucose[index] <- my.anc$DM25.r.glucose
        anc.DM25.t.lag[index] <- my.anc$DM25.t.lag
    }

    growth.summary2 <- growth.summary %>%
        mutate(log.DM0.r.ratio=log(DM0.r/anc.DM0.r)) %>%
        mutate(log.DM0.t_mid.ratio=log(DM0.t_mid/anc.DM0.t_mid)) %>%
        mutate(log.DM25.r.ratio=log(DM25.r/anc.DM25.r)) %>%
        mutate(log.DM25.t_mid.ratio=log(DM25.t_mid/anc.DM25.t_mid)) %>%
        mutate(log.DM0.r.citrate.ratio=log(DM0.r.citrate/anc.DM0.r.citrate)) %>%
        mutate(log.DM0.t.lag.ratio=log(DM0.t.lag/anc.DM0.t.lag)) %>%
        mutate(log.DM25.r.citrate.ratio=log(DM25.r.citrate/anc.DM25.r.citrate)) %>%
        mutate(log.DM25.r.glucose.ratio=log(DM25.r.glucose/anc.DM25.r.glucose)) %>%
        mutate(log.DM25.t.lag.ratio=log(DM25.t.lag/anc.DM25.t.lag))

    return(growth.summary2)
}

run.ratio.confint.bootstrapping <- function(final.growth.summary) {

    evolved.growth.summary <- filter(final.growth.summary, Name != Founder)

    log.DM0.r.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.r.ratio)
    mean.log.DM0.r.ratio <- mean(evolved.growth.summary$log.DM0.r.ratio,na.rm=TRUE)

    log.DM0.t_mid.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM0.t_mid.ratio)
    mean.log.DM0.t_mid.ratio <- mean(evolved.growth.summary$log.DM0.t_mid.ratio)

    log.DM25.r.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.r.ratio)
    mean.log.DM25.r.ratio <- mean(evolved.growth.summary$log.DM25.r.ratio,na.rm=TRUE)

    log.DM25.t_mid.ratio.conf.int <- calc.bootstrap.conf.int(evolved.growth.summary$log.DM25.t_mid.ratio)
    mean.log.DM25.t_mid.ratio <- mean(evolved.growth.summary$log.DM25.t_mid.ratio)

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

    ## Make a figure of these confidence intervals.

    bootstrap.results <- data.frame (
        Parameter = c("log.DM0.r.ratio",
                      "log.DM25.r.ratio",
                      "log.DM0.r.citrate.ratio",
                      "log.DM25.r.citrate.ratio",
                      "log.DM25.r.glucose.ratio",
                      "log.DM0.t_mid.ratio",
                      "log.DM25.t_mid.ratio",
                      "log.DM0.t.lag.ratio",
                      "log.DM25.t.lag.ratio"),
        Estimate = c(mean.log.DM0.r.ratio,
                     mean.log.DM25.r.ratio,
                     mean.log.DM0.r.citrate.ratio,
                     mean.log.DM25.r.citrate.ratio,
                     mean.log.DM25.r.glucose.ratio,
                     mean.log.DM0.t_mid.ratio,
                     mean.log.DM25.t_mid.ratio,
                     mean.log.DM0.t.lag.ratio,
                     mean.log.DM25.t.lag.ratio),
        Left = c(log.DM0.r.ratio.conf.int[1],
                 log.DM25.r.ratio.conf.int[1],
                 log.DM0.r.citrate.ratio.conf.int[1],
                 log.DM25.r.citrate.ratio.conf.int[1],
                 log.DM25.r.glucose.ratio.conf.int[1],
                 log.DM0.t_mid.ratio.conf.int[1],
                 log.DM25.t_mid.ratio.conf.int[1],
                 log.DM0.t.lag.ratio.conf.int[1],
                 log.DM25.t.lag.ratio.conf.int[1]),
        Right = c(log.DM0.r.ratio.conf.int[2],
                  log.DM25.r.ratio.conf.int[2],
                  log.DM0.r.citrate.ratio.conf.int[2],
                  log.DM25.r.citrate.ratio.conf.int[2],
                  log.DM25.r.glucose.ratio.conf.int[2],
                  log.DM0.t_mid.ratio.conf.int[2],
                  log.DM25.t_mid.ratio.conf.int[2],
                  log.DM0.t.lag.ratio.conf.int[2],
                  log.DM25.t.lag.ratio.conf.int[2]),
        stringsAsFactors=FALSE)
    return(bootstrap.results)
}

plot.growth.parameters <- function(growth) {    
    growth.CI <- bootstrap.my.growth.confints(growth)
    growth.plot.df <- growth %>%
        gather(key="Parameter",value="Estimate",
               DM0.r.citrate,DM25.r.citrate,DM25.r.glucose,DM0.t.lag,DM25.t.lag)
    growth.plot <- plot.growth.confints(growth.plot.df, growth.CI)
    return(growth.plot)    
}

plot.growthcurver.parameters <- function(growth) {    
    growth.CI <- bootstrap.growthcurver.confints(growth)
    growth.plot.df <- growth %>%
        gather(key="Parameter",value="Estimate",
               DM0.r,DM25.r,DM0.t_mid,DM25.t_mid)
    growth.plot <- plot.growth.confints(growth.plot.df, growth.CI)
    return(growth.plot)   
}

plot.Fig2D <- function (plot.df, confints.df) {

    recode.Parameter <- function(df) {
        df <- mutate(df,
                     Parameter=recode_factor(
                         Parameter,
                         log.DM0.r.citrate.ratio = "DM0 r citrate",
                         log.DM0.r.ratio = "DM0 r",
                         log.DM0.t.lag.ratio = "DM0 t_lag",
                         log.DM0.t_mid.ratio = "DM0 t_mid",
                         log.DM25.r.citrate.ratio = "DM25 r citrate",
                         log.DM25.r.ratio = "DM25 r",
                         log.DM25.r.glucose.ratio = "DM25 r glucose",
                         log.DM25.t.lag.ratio = "DM25 t_lag",
                         log.DM25.t_mid.ratio = "DM25 t_mid"
                     ))
          return(df)
    }

    plot.df <- recode.Parameter(plot.df)
    confints.df <- recode.Parameter(confints.df) %>% drop_na()
    
    the.plot <- ggplot(plot.df, aes(x=Parameter,y=Estimate)) +
        geom_point(size=1) +
        geom_errorbar(data=confints.df, aes(ymin=Left,ymax=Right),width=0.5, size=0.5) +
        geom_hline(yintercept=0,linetype='dashed') +
        ylab("log(Evolved/Ancestral)") +
        xlab("Growth parameter") +
        theme_classic()
    return(the.plot)
}

###############################################
## Load the growth data from the plate reader.
## DM0-evolved population data.
DM0.pop.growth.data <- read.csv(
  file.path(proj.dir,"results/DM0-evolved-pop-growth.csv"),
  stringsAsFactors=FALSE) %>%
  left_join(pop.clone.labels, by="Name") %>%
  prep.growth.df()

## DM0-evolved clone data.
DM0.clone.growth.data <- read.csv(
  file.path(proj.dir,"results/DM0-evolved-clone-growth.csv"),
  stringsAsFactors=FALSE) %>%
  left_join(pop.clone.labels, by="Name") %>%
  prep.growth.df()

## DM25-evolved clone data (only have growth in DM25).
DM25.clone.growth.data <- read.csv(
  file.path(proj.dir, "results/DM25-evolved-clone-growth.csv"),
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
save_plot(file.path(proj.dir,"results/figures/ZDBp874_growth.pdf"),ZDBp874.plot,base_height=2,base_width=3)

## plot rest of CZB151-descended isolates. Can do good fitness competitions on these.
ZDBp871.df <- filter(DM0.clone.growth.data, Name %in% c('ZDBp871','CZB151'))
ZDBp871.plot <- plot.single.growthcurve(ZDBp871.df,'ZDBp871','CZB151', logscale=FALSE)
save_plot(file.path(proj.dir,"results/figures/ZDBp871_growth.pdf"),ZDBp871.plot,base_height=2,base_width=3)

ZDBp889.df <- filter(DM0.clone.growth.data, Name %in% c('ZDBp889','CZB151'))
ZDBp889.plot <- plot.single.growthcurve(ZDBp889.df,'ZDBp889','CZB151', logscale=FALSE)
save_plot(file.path(proj.dir,"results/figures/ZDBp889_growth.pdf"),ZDBp889.plot,base_height=2,base_width=3)

ZDBp892.df <- filter(DM0.clone.growth.data, Name %in% c('ZDBp892','CZB151'))
ZDBp892.plot <- plot.single.growthcurve(ZDBp892.df,'ZDBp892','CZB151', logscale=FALSE)
save_plot(file.path(proj.dir,"results/figures/ZDBp892_growth.pdf"),ZDBp892.plot,base_height=2,base_width=3)

## SUPER COOL! ZDBp874 grows WORSE than ancestor in DM25!!!!
## ZDBp871 and ZDBp889 also grow worse in DM25 than ancestor!
## This is really interesting!

## Fig. 2A: Plot growth curves for DM0-evolved clones.
Fig2A.plot <- plot.growthcurve.figure(filter(DM0.clone.growth.data,Hours<=24), logscale=FALSE)
save_plot(file.path(proj.dir,"results/figures/Fig2A.pdf"),Fig2A.plot,base_height=7,base_width=11)

## Fig. 2E: Plot growth curves for DM0-evolved populations.
pop.Fig2E.plot <- plot.growthcurve.figure(DM0.pop.growth.data, logscale=FALSE)
save_plot(file.path(proj.dir,"results/figures/pop_Fig2E.pdf"),pop.Fig2E.plot,base_height=7,base_width=11)

## S5 Figure Panel A. Plot growth curves for DM25-evolved clones.
S5FigA.plot <- plot.growthcurve.figure(DM25.clone.growth.data, logscale=FALSE)

## plot log-scale plots.
S3Fig.plot <- plot.growthcurve.figure(DM0.clone.growth.data,logscale=TRUE)
save_plot(file.path(proj.dir,"results/figures/S3Fig.pdf"),S3Fig.plot,base_height=7,base_width=11)
S4Fig.plot <- plot.growthcurve.figure(DM0.pop.growth.data,logscale=TRUE)
save_plot(file.path(proj.dir,"results/figures/S4Fig.pdf"),S4Fig.plot,base_height=7,base_width=11)

S6Fig.plot <- plot.growthcurve.figure(DM25.clone.growth.data, logscale=TRUE)
save_plot(file.path(proj.dir,"results/figures/S6Fig.pdf"), S6Fig.plot,base_height=7,base_width=11)

## data filtering to plot data included in my growth rate estimation.
## This code is useful for debugging and making sure the rate estimation works right.
filtered.clone.growth.data <- filter.growth.data.on.analysis.domain(DM0.clone.growth.data)
filtered.pop.growth.data <- filter.growth.data.on.analysis.domain(DM0.pop.growth.data)
## filtered log-scale plots, again for debugging purposes.
filtered.log.Fig2A.plot <- plot.growthcurve.figure(filtered.clone.growth.data,logscale=TRUE)
filtered.log.Fig2B.plot <- plot.growthcurve.figure(filtered.pop.growth.data,logscale=TRUE)

################################################################################
## plot my growth estimates and growthcurver estimates as well.
## plot estimates with confidence interval of mean.
## Do the same for the growthcurver fits. Note different confint function needed.

clone.growth.plot <- plot.growth.parameters(clone.growth)
clone.growthcurver.plot <- plot.growthcurver.parameters(clone.growth.curve.fits)
clone.growth.parameter.plot <- plot_grid(clone.growth.plot,
                                         clone.growthcurver.plot,ncol=1,
                                         labels=c('B','C'))

pop.growth.plot <- plot.growth.parameters(pop.growth)
pop.growthcurver.plot <- plot.growthcurver.parameters(pop.growth.curve.fits)
pop.growth.parameter.plot <- plot_grid(pop.growth.plot,
                                         pop.growthcurver.plot,ncol=1,
                                         labels=c('F','G'))


DM25.clone.growth.plot <- plot.growth.parameters(DM25.clone.growth)
DM25.clone.growthcurver.plot <- plot.growthcurver.parameters(DM25.clone.growth.curve.fits)
DM25.clone.growth.parameter.plot <- plot_grid(DM25.clone.growth.plot,
                                         DM25.clone.growthcurver.plot,ncol=1,
                                         labels=c('B','C'))

## summarize and take averages over clone and population samples.
clone.growth.summary <- summarize.growth.results(clone.growth)
pop.growth.summary <- summarize.growth.results(pop.growth)
growth.summary <- full_join(clone.growth.summary, pop.growth.summary)

pop.growthcurver.summary <- summarize.growthcurver.results(pop.growth.curve.fits)
clone.growthcurver.summary <- summarize.growthcurver.results(clone.growth.curve.fits) %>%
    ## ZDBp874 is Cit-, so exclude.
    filter(Name != 'ZDBp874')
growthcurver.summary <- full_join(pop.growthcurver.summary, clone.growthcurver.summary)

## growth results for DM25-evolved clones.
DM25.growth.summary <- summarize.growth.results(DM25.clone.growth)
DM25.growthcurver.summary <- summarize.growthcurver.results(DM25.clone.growth.curve.fits)

##################### Supplementary figure:
## Compare growthcurver estimates to my estimates,
## in order to show consistent results between the two methods.
growth.estimate.comp.df <- full_join(
    growth.summary,
    growthcurver.summary) %>%
    ungroup()

DM25.growth.estimate.comp.df <- full_join(
    DM25.growth.summary,
    DM25.growthcurver.summary) %>%
    ungroup()

## r = 0.68, 0.70, -0.14, 0.92, 0.59 respectively.
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


pop.estimates <- filter(growth.estimate.comp.df,Dataset=='PopulationGrowth')
## r = 0.64, 0.83, -0.04, 0.81, 0.92, respectively.
cor(pop.estimates$DM0.r, pop.estimates$DM0.r.citrate,use="complete.obs")
cor(pop.estimates$DM25.r, pop.estimates$DM25.r.citrate,use="complete.obs")
cor(pop.estimates$DM25.r, pop.estimates$DM25.r.glucose,use="complete.obs")
cor(pop.estimates$DM0.t_mid, pop.estimates$DM0.t.lag,use="complete.obs")
cor(pop.estimates$DM25.t_mid, pop.estimates$DM25.t.lag,use="complete.obs")

clone.estimates <- filter(growth.estimate.comp.df,Dataset=='CloneGrowth')
## r = 0.48, 0.60, -0.11, 0.97, 0.68, respectively.
cor(clone.estimates$DM0.r, clone.estimates$DM0.r.citrate,use="complete.obs")
cor(clone.estimates$DM25.r, clone.estimates$DM25.r.citrate,use="complete.obs")
cor(clone.estimates$DM25.r, clone.estimates$DM25.r.glucose,use="complete.obs")
cor(clone.estimates$DM0.t_mid, clone.estimates$DM0.t.lag,use="complete.obs")
cor(clone.estimates$DM25.t_mid, clone.estimates$DM25.t.lag,use="complete.obs")

## Note: the following are for the DM0-evolved clones.
rate.comp.plot1 <- ggplot(growth.estimate.comp.df,
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

rate.comp.plot2 <- ggplot(growth.estimate.comp.df,
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

rate.comp.plot3 <- ggplot(growth.estimate.comp.df,
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

rate.comp.plot4 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM0.t_mid,y=DM0.t.lag,color=Dataset)) +
    scale_color_manual(values=c('purple','green')) +
    geom_point() +
    theme_classic() +
    guides(color=FALSE) +
    geom_abline(slope=1,intercept=0,linetype='dotted') +
    xlab("DM0 t_mid") +
    ylab("DM0 t_lag")

rate.comp.plot5 <- ggplot(growth.estimate.comp.df,
                          aes(x=DM25.t_mid,y=DM25.t.lag,color=Dataset)) +
    scale_color_manual(values=c('purple','green')) +
    geom_point() +
    theme_classic() +
    guides(color=FALSE) +
    geom_abline(slope=1,intercept=0,linetype='dotted') +
    xlab("DM25 t_mid") +
    ylab("DM25 t_lag")

## Save figure for the Supplement. 
S7Fig <- plot_grid(rate.comp.plot1,
                                rate.comp.plot2,
                                rate.comp.plot3,
                                rate.comp.plot4,
                                rate.comp.plot5,
                                labels=c('A','B','C','D','E'), ncol=1)
S7Fig.outf <- file.path(proj.dir, "results/figures/S7Fig.pdf")
save_plot(S7Fig.outf, S7Fig, base_height=10,base_width=4)

####################################
### supplementary figure:
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
S8Fig <- plot_grid(cit.glucose.cor.plot,
                          cit.cit.cor.plot,
                          r.r.cor.plot, labels=c('A','B','C'), ncol=1)
S8Fig.outf <- file.path(proj.dir,"results/figures/S8Fig.pdf")
save_plot(S8Fig.outf, S8Fig, base_height=8, base_width=12)
######################################################################
final.pop.growth.summary <- calc.growth.log.ratios(pop.estimates)
final.clone.growth.summary <- calc.growth.log.ratios(clone.estimates)

## Exclude the Cit- ZDBp874 strain from these calculations.
clone.bootstrap.CI <- run.ratio.confint.bootstrapping(filter(
    final.clone.growth.summary,Name!='ZDBp874'))

clone.Fig2D.plot.df <- final.clone.growth.summary %>%
    gather(key="Parameter", value="Estimate",
           log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
           log.DM25.t.lag.ratio, log.DM0.t.lag.ratio, log.DM0.r.ratio, log.DM0.t_mid.ratio,
           log.DM25.r.ratio, log.DM25.t_mid.ratio) %>%
    filter(Name != Founder) ## ancestors always have a ratio of zero.

Fig2D.clone <- plot.Fig2D(clone.Fig2D.plot.df, clone.bootstrap.CI)
save_plot(file.path(proj.dir,"results/figures/Fig2D_clone.pdf"),Fig2D.clone)

Fig2BCD <- plot_grid(clone.growth.parameter.plot,Fig2D.clone,labels=c('','D'),ncol=1,
                     rel_heights = c(1.8, 1),rel_widths=c(1.5,1))

pop.bootstrap.CI <- run.ratio.confint.bootstrapping(final.pop.growth.summary)
pop.Fig2H.plot.df <- final.pop.growth.summary %>%
    gather(key="Parameter", value="Estimate",
           log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio, log.DM0.r.citrate.ratio,
           log.DM25.t.lag.ratio, log.DM0.t.lag.ratio, log.DM0.r.ratio, log.DM0.t_mid.ratio,
           log.DM25.r.ratio, log.DM25.t_mid.ratio) %>%
    filter(Name != Founder) ## ancestors always have a ratio of zero.
Fig2H.pop <- plot.Fig2D(pop.Fig2H.plot.df, pop.bootstrap.CI)
save_plot(file.path(proj.dir,"results/figures/Fig2D_pop.pdf"),Fig2H.pop)

Fig2FGH <- plot_grid(pop.growth.parameter.plot,Fig2H.pop,labels=c('','H'),ncol=1,
                     rel_heights = c(1.8, 1))

DM25.final.clone.growth.summary <- calc.growth.log.ratios(DM25.growth.estimate.comp.df)
DM25.clone.bootstrap.CI <- run.ratio.confint.bootstrapping(DM25.final.clone.growth.summary)

S5FigD.plot.df <- DM25.final.clone.growth.summary %>%
    gather(key="Parameter", value="Estimate",
           log.DM25.r.glucose.ratio, log.DM25.r.citrate.ratio,
           log.DM25.t.lag.ratio, log.DM25.r.ratio, log.DM25.t_mid.ratio) %>%
    filter(Name != Founder) ## ancestors always have a ratio of zero.

S5FigD <- plot.Fig2D(S5FigD.plot.df, DM25.clone.bootstrap.CI)

S5FigBCD <- plot_grid(DM25.clone.growth.parameter.plot,
                          S5FigD,
                          labels=c('','D'),ncol=1,
                          rel_heights = c(1.8, 1))

####### Make Figure 2 using cowplot.
Fig2outf <- file.path(proj.dir,"results/figures/Fig2.pdf")
Fig2.1 <- plot_grid(Fig2A.plot, Fig2BCD, labels = c('A', ''), nrow = 1, rel_widths = c(1.2, 1))
Fig2.2 <- plot_grid(pop.Fig2E.plot, Fig2FGH, labels = c('E', ''), nrow = 1, rel_widths = c(1.2, 1))
Fig2 <- plot_grid(Fig2.1,Fig2.2,labels=c('',''),nrow=2)
save_plot(Fig2outf,Fig2,base_height=18,base_width=18)

## supplementary fig of DM25-evolved clone growth.
S5Figoutf <- file.path(proj.dir,"results/figures/S5Fig.pdf")
S5Fig <- plot_grid(S5FigA.plot, S5FigBCD, labels = c('A', ''), nrow = 1, rel_widths = c(1.2, 1))
save_plot(S5Figoutf, S5Fig, base_height=9, base_width=18)

####### Figure 3: Nkrumah's cell death results.

### Figure 4 (AFTER cell death results): plot CFU results.
## After looking at the growth curves, let's examine the CFUs in DM25.
calc.CFU.log.ratios <- function(CFU.summary) {
    anc.Counts <- rep(-1,nrow(CFU.summary))
    for (index in 1:nrow(CFU.summary)) {
        my.row <- CFU.summary[index, ]
        my.anc <- filter(CFU.summary,Name==my.row$Founder)
        anc.Counts[index] <- my.anc$Counts
    }
    CFU.summary2 <- CFU.summary %>%
        mutate(log.Count.ratio=log(Counts/anc.Counts))
    return(CFU.summary2)
}

bootstrap.CFU.CI <- function(CFU) {
    evolved.CFU.summary <- filter(CFU, Name != Founder)
    log.CFU.ratio.conf.int <- calc.bootstrap.conf.int(evolved.CFU.summary$log.Count.ratio)
    mean.log.CFU.ratio <- mean(evolved.CFU.summary$log.Count.ratio)

    confint.df <- data.frame (
        Founder = unique(CFU$Founder),
        Environment = unique(CFU$Environment),
        log.Count.ratio = c(mean.log.CFU.ratio),
        Left = c(log.CFU.ratio.conf.int[1]),
        Right = c(log.CFU.ratio.conf.int[2]),
        stringsAsFactors = FALSE
    )
    return(confint.df)
}

CFU.data <- read.csv(file.path(proj.dir,"data/rohan-formatted/DM25-Cell-Counts-6-5-19.csv"),
                     stringsAsFactors=FALSE) %>%
    left_join(pop.clone.labels, by="Name") %>%
    filter(Name != "BLANK") %>% mutate(Generation=as.factor(Generation))

CFU.plotA <- ggplot(CFU.data,aes(x=Name,y=Counts,shape=Environment,color=Generation)) +
    geom_point() +
    facet_grid(.~Founder,scales='free') +
    scale_color_manual(values=cbbPalette) +
    guides(color=FALSE,shape=FALSE) +
    xlab("Sequenced isolate") +
    ylab("Viable cell counts") +
    theme(axis.text.x  = element_text(angle=90, vjust=-0.5, size=5))

CFU.summary <- CFU.data %>%
    filter(Name != 'ZDBp874') %>%
    group_by(Founder,Environment,Generation,Name) %>%
    summarise(Counts=mean(Counts)) %>%
    ungroup() %>%
    calc.CFU.log.ratios()

CFU.bootstrap.CI <- CFU.summary %>%
    droplevels() %>% ## don't run on empty subsets
    split(list(.$Founder, .$Environment)) %>%
    map_dfr(.f=bootstrap.CFU.CI)

CFU.plotB <- ggplot(CFU.summary, aes(x=Environment,y=log.Count.ratio)) +
    geom_point(size=1) +
    geom_errorbar(data=CFU.bootstrap.CI, aes(ymin=Left,ymax=Right),width=0.5, size=0.5) +
    facet_wrap(.~Founder) +
    geom_hline(yintercept=0,linetype='dashed') +
    ylab("log(Evolved/Ancestral)") +
    xlab("Experimental evolution treatment")

CFU.Fig <- plot_grid(CFU.plotA,CFU.plotB,ncol=1,labels=c('A','B'))
CFU.outf <- file.path(proj.dir, "results/figures/CFU_Fig.pdf")
save_plot(CFU.outf, CFU.Fig,base_height=6,base_width=6)

################################################################################
## Figure 5: make a matrix plot of genes with mutations in two or more clones.

PlotMatrixFigure <- function(raw.matrix.file, amp.matrix.file,
                             ltee.matrix.file, ltee.50k.labels.file,
                             matrix.outfile, co.occurrence.outfile) {

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
        select_if(not_any_na)

    DM0.DM25.matrix.data <- gather(merged.with.amps.matrix,"Name","mutation.count",2:ncol(merged.with.amps.matrix)) %>%
        left_join(pop.clone.labels) %>%
        select(Gene,Name,mutation.count,Environment,PopulationLabel) %>%
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
    ## get the gene order from the heatmap,
    sorted.genes <- sapply(MxM.heatmap$rowInd,function(i) rownames(MxM.scaled.cor.matrix)[i])
    ## and use this gene ordering for the mutation matrix figure.
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

    ggsave(matrix.figure, file=matrix.outfile,width=7,height=15)

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
                                     co.occurrence=c(as.factor(M.co.occurrence)),stringsAsFactors=FALSE)

    M.co.occurrence.df$row <- factor(M.co.occurrence.df$row,levels=M.co.occurrence.sorted.genes)
    M.co.occurrence.df$col <- factor(M.co.occurrence.df$col,levels=M.co.occurrence.sorted.genes)

    ## plot the co-occurrence matrix.
    M.co.occurrence.plot <- ggplot(M.co.occurrence.df,aes(row,col,fill=co.occurrence)) +
        geom_tile(color='White') +
        scale_fill_viridis(option="magma",direction=-1) +
        theme(axis.text.x = element_text(angle=90,vjust=-0.1)) +
        guides(fill=FALSE)
    M.co.occurrence.plot

    ggsave(co.occurrence.outfile,
           M.co.occurrence.plot,width=12,height=12)
}

##### Now run this function!
#' on all valid mutations
raw.matrix.f <- file.path(proj.dir,
        "results/DM0-DM25-comparison-mut-matrix.csv")

#' on just non-synonymous mutations
dN.raw.matrix.f <- file.path(proj.dir,
        "results/dN-DM0-DM25-comparison-mut-matrix.csv")

amp.matrix.f <- file.path(proj.dir,"results/amp_matrix.csv")
ltee.matrix.f <- file.path(proj.dir,"results/LTEE-mut_matrix.csv")
dN.ltee.matrix.f <- file.path(proj.dir,"results/dN-LTEE-mut_matrix.csv")
ltee.50k.labels.f <- file.path(proj.dir, "data/rohan-formatted/LTEE-50K-clones.csv")
matrix.outf <- "../results/figures/Fig5.pdf"
dN.matrix.outf <- "../results/figures/Fig5B.pdf"
co.occurrence.outf <- file.path(proj.dir,"results/figures/co_occurrence.pdf")
dN.co.occurrence.outf <- file.path(proj.dir,"results/figures/dN_co_occurrence.pdf")

PlotMatrixFigure(raw.matrix.f, amp.matrix.f, ltee.matrix.f, ltee.50k.labels.f, matrix.outf,co.occurrence.outf)
PlotMatrixFigure(dN.raw.matrix.f, amp.matrix.f, dN.ltee.matrix.f, ltee.50k.labels.f, dN.matrix.outf,dN.co.occurrence.outf)

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
    scale_fill_manual(values=cbbPalette) +
    ylab("Count") +
    xlab("Position") +
    theme_tufte(base_family="Helvetica")

save_plot(file.path(proj.dir,"results/figures/Fig6B.pdf"),
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
    scale_color_manual(values=cbbPalette) +
    theme_classic() +
    ylab("Count") +
    xlab("Position") +
    geom_text_repel(fontface = "italic")

save_plot(file.path(proj.dir,"results/figures/Fig6A.pdf"),
         parallel.IS.plot,base_aspect_ratio=1.5)

########
## Plot the rate of increase of IS-elements in the DM0 and DM25 experiments
## in comparison to the rate of increase of IS-elements in Ara-3.

LTEE.MAE.IS.insertions <- read.csv(
    file.path(proj.dir,
              "results/genome-analysis/LTEE_MAE_IS150_insertions.csv"),
    stringsAsFactors=FALSE) %>%
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

##IS150.rate.df <- rbind(MAE.IS150.over.time,DM0.DM25.over.time,Ara.minus.3.IS150.by.clone)
## don't bother with MAE comparison, since Jeff Barrick says the rate and spectrum is
## weird due to growth on plates.
IS150.rate.df <- rbind(DM0.DM25.over.time,Ara.minus.3.IS150.by.clone)

IS150.rate.plot <- ggplot(IS150.rate.df,
                          aes(x=Generation,y=total.count,color=Environment)) +
    theme_classic() +
    scale_colour_manual(values=c('#1b9e77','#d95f02','#7570b3')) +
    ylab('IS150 insertions') +
    geom_jitter(width=50) +
    guides(color=FALSE)

save_plot(file.path(proj.dir,"results/figures/Fig6C.pdf"),
         IS150.rate.plot,base_aspect_ratio=1.57)


########
## Combine the IS plots with cowplot to make Figure 4.
Fig6outf <- file.path(proj.dir,"results/figures/Fig6.pdf")
Fig6 <- plot_grid(parallel.IS.plot, IS.plot, IS150.rate.plot,
                  labels = c('A', 'B', 'C'), ncol=1)
save_plot(Fig6outf,Fig6,base_height=10,base_aspect_ratio=0.8)

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

maeA.fitness.data <- read.csv(file.path(proj.dir,"data/rohan-formatted/DM0_Fitness2.csv"),
                      header=TRUE,stringsAsFactors=FALSE)
## these data come from one day competitions that Tanush ran.
res1 <- filter(maeA.fitness.data,Red.Pop=='ZDB151_with_maeA') %>% maeA.fitness.analysis()
res2 <- filter(maeA.fitness.data,Red.Pop=='ZDB67_with_maeA') %>% maeA.fitness.analysis(rev=TRUE)
res3 <- filter(maeA.fitness.data,Red.Pop=='ZDB152_with_maeA') %>% maeA.fitness.analysis()
res4 <- filter(maeA.fitness.data,Red.Pop=='ZDB68_with_maeA') %>% maeA.fitness.analysis(rev=TRUE)

## beautiful! All confints overlap with each other, showing no maeA fitness effect
## does not depend on Ara polarity or genetic background.

## let's combine data from the switched polarity competitions (since Ara marker is neutral)
d1 <- filter(maeA.fitness.data,Red.Pop=='ZDB151_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d2 <- filter(maeA.fitness.data,Red.Pop=='ZDB67_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
## now switch column labels.
d2X <- dplyr::rename(d2,Red.0=White.0,Red.1=White.1,White.0=Red.0,White.1=Red.1)
ZDB151.data <- rbind(d1,d2X)

d3 <- filter(maeA.fitness.data,Red.Pop=='ZDB152_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
d4 <- filter(maeA.fitness.data,Red.Pop=='ZDB68_with_maeA') %>% select(Red.0,Red.1,White.0,White.1,D.0,D.1)
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
