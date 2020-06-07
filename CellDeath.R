## CellDeath.R
## Authors: Nkrumah Grant and Rohan Maddamsetti
## Michigan State University 
## 2019-10-11 
## The following script was used to analyze cell death for the paper
## "Genomic and phenotypic evolution of Escherichia coli to a citrate-only
## resource environment".

##Load libraries
library(boot)
library(tidyverse)
library(cowplot)

calc.w.bootstrap.conf.int <- function(vec,w) {
    ## bootstrap weighted confidence intervals around the mean.
    ## Use the boot package to calculate fancy BCA intervals.

    samplewmean <- function(x,ind,j) {
        ## define this type of mean-calculating function to pass to the boot function.
        ## from: https://stackoverflow.com/questions/46231261/bootstrap-weighted-mean-in-r
        ## Also see: https://web.as.uky.edu/statistics/users/pbreheny/621/F12/notes/9-18.pdf
        return(weighted.mean(x[ind],j[ind]))
    }

    Nbootstraps <- 10000
    ## trick is to pass weights to samplewmean using the parameter j
    out <- boot(vec,samplewmean,Nbootstraps,weights=w, j=w)
    ## handle bad inputs in bootstrapping confidence intervals
    ci.result <- tryCatch(boot.ci(out, type="bca"), error= function(c) return(NULL))
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
###################

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("DM0-evolution","src")))
projdir <- file.path("..")

## Load data 
setwd(file.path(projdir, "Nkrumah_MicroscopyData_DM0_project/"))

## Here I create a list of all .txt files as produced by segger output.
seg.dat <- list.files(pattern = "*.txt")

## Here I create a dataframe.
D1 <- as.data.frame(do.call(rbind,Map('cbind', lapply(seg.dat, read.csv), idcol=seg.dat)))

D1$Clone <- as.factor((str_extract(D1$idcol, ("_\\d+_"))))
D1$Treatment <- as.factor((str_extract(D1$idcol,"(?<!not)DM\\d+")))

levels(D1$Clone) <- c("11364", "151", "606", "871", "910") #set levels removing underscores
D1$Clone <- ordered(D1$Clone, levels = c("606", "151", "871", "910", "11364")) #set the order of levels 
levels(D1$Clone) <- c("REL606", "CZB151", "ZDBp871", "ZDBp910", "REL11364")

## Here I use the idcol to create a column called blocks. I retain the idcol for transparency. 
D1$Block <- as.factor(gsub("^.*20190329.*$","1", D1$idcol))
D1$Block <- as.factor(gsub("^.*20190515.*$","2", D1$Block))
D1$Block <- as.factor(gsub("^.*20190516.*$","3", D1$Block))
D1$Block <- as.factor(gsub("^.*20190517.*$","4", D1$Block))
D1$Block <- as.factor(gsub("^.*20190524.*$","5", D1$Block))

## Tidy data:
## 1) Subset the mean fluorescence columns for channels one and two and the Clone and Treatment blocks 
##     - In a future analysis I will not do this 
## 2) Filter dataset to contain only regions with scores greater than 25. 50 is considered a perfect region.
## 3) Create a new column that calculates the ratio of fluorescence on channel's two
##    and one and then the proportion of cells on channel 2. These two relationships
##    basically ask how blue a cell is. 

D3 <- D1[,c(5,9,14,51:53)] %>% ## 1)
    filter(RegionScore > 25 & Block !=0)  %>% ## 2)
    ## Now 3)
    mutate(RelFluor = Fluor2Mean/Fluor1Mean) %>% ## dead/alive
    mutate(PropFluor = Fluor2Mean/(Fluor2Mean + Fluor1Mean)) %>% ## dead/(dead+alive)
    ##4) Assign a "dead" if the mean fluorescence on the dead channel > alive channel
    ## and a "alive" if less than. 
    mutate(RelFluorPrime = as.factor(ifelse(Fluor2Mean>Fluor1Mean,"Dead", "Alive")))

    ##5) Calculate the proportion of dead cells per replicate 
    ##5) Let me adjust this code to also include summary information. 
D4 <- D3 %>% group_by(Clone, Treatment, Block) %>% 
    count(RelFluorPrime) %>% 
    mutate(PropPrime = n / sum(n))

## Here lets generate a couple of plots 
## 1) Lets plot the relationship between blocks of the same clone by environment.
## Will give me a sense of how heterogenous the data is. 

## Here I plot the relative and the propotional fluorescence of each independent replicate. 
## This will give me a sense of the variability in my dataset. 
## Keep in mind that I am using a log scale, thus the intepretations are as followed. 
## interpretations: 
## plot 1: There is greater fluorescence in the death channel when value approaches "0"
## plot 2: The proportion of death is higher as you approach "negative infinity "0"

## Plot 1: Distribution showing relative fluorescences (Dead/Alive)
plot1 <- D3 %>% group_by(Clone, Treatment, Block) %>% 
  ggplot(aes(x=Block, y=log(RelFluor), color=Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Clone + Treatment) +
  labs(y="Log(ratio dead cells)") 

## Plot 2 : Boxplots showing the distribution of proportions (Dead/Dead+Alive) 
plot2 <- D3 %>% group_by(Clone, Treatment, Block) %>% 
    ggplot(aes(x=Block, y=PropFluor, color=Treatment)) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values=c('red3','black')) +
    facet_wrap(~Clone + Treatment) +
    guides(color=FALSE) +
    labs(y="Proportion dead cells")

## Here I plot a 2D density plot to look at whether I can can score cells as dead 
## or alive given some value of fluoresence. 

## Plot 3 : density plots showing seperation 
## log-value(x) > 1 dead. 
plot3 <- D3 %>% group_by(Clone, Treatment, Block) %>% 
  ggplot(aes(x=(log(RelFluor)), y=log(PropFluor), color = RelFluorPrime)) +
  geom_density2d() +
  facet_wrap(~Clone + Treatment) +
  labs(x="Log(proportion dead cells)", y="Log(ratio dead cells)") +
  guides(color=guide_legend(title="Cell viability"))

## Plot 4 : proportion of dead:alive cells on a population level. 

## PI emission is 617 nm when bound.
## SYTO9 emission is 498 nm.
## I used this website to convert to hex.
## https://academo.org/demos/wavelength-to-colour-relationship/
SYTO9.color <- '#00ffa9'
PI.color <- '#ff8200'

cell.death.barplot <- D4 %>%
  group_by(Clone, Treatment) %>% 
  ggplot(aes(x=Treatment, y = n, fill = RelFluorPrime)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(SYTO9.color,PI.color)) +
  theme_classic() +
  facet_wrap(. ~ Clone, nrow=1) +
  guides(fill = FALSE) +
  ylab("Number of cells") +
  labs(fill = "Cell viability") #Changes legend title from "RelFluorPrime" to "Cell viability"
  
## Plot 6: plot of the proportions for independent replicates.  
## First must construct the data frame piece-wise.
D4.2 <- D4 %>%  filter(RelFluorPrime == "Dead") %>% rename(P.Dead = n) %>% rename(PropPrime.Dead = PropPrime)
D4.3 <- D4 %>%  filter(RelFluorPrime == "Alive") %>% rename(P.Alive = n) %>%  rename(PropPrime.Alive = PropPrime)
D4.2$P.Alive <- D4.3$P.Alive
D4.2$PropPrime.Alive <- D4.3$PropPrime.Alive
D4.2 <- D4.2 %>%  
    mutate(N = P.Dead + P.Alive)
D4.2 <- D4.2[,-4]

##Calculate the confidence intervals here. Note: C.I. for proportions are calculated differently. 
P6 <- D4.2 %>% group_by(Clone, Treatment, Block) %>% 
    mutate(B = ((PropPrime.Dead * PropPrime.Alive) / N)) %>% 
    mutate(marg.err = sqrt(B) * 1.96) %>% #1.96 Z score for the 95% C.I. 
    mutate(CIlow = PropPrime.Dead - marg.err, CIHigh = PropPrime.Dead + marg.err) %>%
    mutate(my.text=paste('N =',as.character(N)))

## calculate BCa confidence intervals, and weight proportions
## by the number of cells counted in that block.
w.bootstrap.PropPrime.Dead <- function(df) {
    PropPrime.Dead.confint <- calc.w.bootstrap.conf.int(df$PropPrime.Dead,df$N)
    weighted.bootstrap.results <- data.frame(
        Clone = unique(df$Clone),
        Treatment = unique(df$Treatment),
        PropPrime.Dead = weighted.mean(df$PropPrime.Dead,df$N),
        Left = PropPrime.Dead.confint[1],
        Right = PropPrime.Dead.confint[2],
        stringsAsFactors=FALSE
    )                                        
    return(weighted.bootstrap.results)
}

P6.confint.df <- P6 %>%
    as.data.frame() %>%
    split(list(.$Clone,.$Treatment)) %>%
    keep(function(x) nrow(x) > 0) %>% ## ignore subsets with no rows
    map_dfr(.f = w.bootstrap.PropPrime.Dead)

Fig10B <- P6 %>%
  ggplot(aes(x = Treatment, y = PropPrime.Dead, label = my.text)) +
  geom_point(size=0.5) +
  facet_wrap(. ~ Clone, nrow = 1) +
  theme_classic() +
  ylab("Proportion dead cells") +
  geom_errorbar(data=P6.confint.df,aes(x=Treatment,ymin=Left,ymax=Right), width=0.15, size=0.4,inherit.aes=FALSE) +
  geom_crossbar(data=P6.confint.df,aes(x=Treatment,y=PropPrime.Dead,ymin=PropPrime.Dead,ymax=PropPrime.Dead), width= 0.1, size = 0.4,inherit.aes=FALSE)

labeled.Fig10B <- plot_grid(Fig10B, labels=c('B'),ncol=1)
save_plot(file.path(projdir,"results/figures/nkrumah-figures/Fig10B.pdf"),labeled.Fig10B,base_height=3,base_width=8)
