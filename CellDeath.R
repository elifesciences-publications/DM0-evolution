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
library(ggrepel)
library(devtools)
library(gridExtra)
library(reshape2)
library(grid)
library(ggpubr)
library(PairedData)
library(car)

## IMPORTANT TODO: CHECK HOW the second argument matters in this line of code:
###return(c(mean(x[ind]),var(x[ind])/length(ind)))

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
D1$Block<-  as.factor(gsub("^.*20190329.*$","1",D1$idcol))
D1$Block <- as.factor(gsub("^.*20190515.*$","2",D1$Block))
D1$Block <- as.factor(gsub("^.*20190516.*$","3",D1$Block))
D1$Block <- as.factor(gsub("^.*20190517.*$","4",D1$Block))
D1$Block <- as.factor(gsub("^.*20190524.*$","5",D1$Block))

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

Fig4B <- D4 %>%
  group_by(Clone, Treatment) %>% 
  ggplot(aes(x=Treatment, y = n, fill = RelFluorPrime)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(SYTO9.color,PI.color)) +
  theme_classic() +
  facet_wrap(. ~ Clone, nrow=1) +
  guides(fill = FALSE) +
  ylab("Number of cells") +
  labs(fill = "Cell viability") #Changes legend title from "RelFluorPrime" to "Cell viability"

cell.death.barplot <- D4 %>% group_by(Clone, Treatment) %>% 
  ggplot(aes(x=Treatment, y = n, fill = RelFluorPrime)) +
  geom_bar(stat = "identity") + 
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

Fig4C <- P6 %>%
  ggplot(aes(x = Treatment, y = PropPrime.Dead, label = my.text)) +
  geom_point(size=0.5) +
  facet_wrap(. ~ Clone, nrow = 1) +
  theme_classic() +
  ylab("Proportion dead cells") +
  geom_errorbar(data=P6.confint.df,aes(x=Treatment,ymin=Left,ymax=Right), width=0.15, size=0.4,inherit.aes=FALSE) +
  geom_crossbar(data=P6.confint.df,aes(x=Treatment,y=PropPrime.Dead,ymin=PropPrime.Dead,ymax=PropPrime.Dead), width= 0.1, size = 0.4,inherit.aes=FALSE)

Fig4BC <- plot_grid(Fig4B, Fig4C, labels=c('B','C'),ncol=1)
save_plot(file.path(projdir,"results/figures/nkrumah-figures/Fig4BC.pdf"),Fig4BC,base_height=6,base_width=9)


## Let's look at the the estimated proportion of dead cells to report in the manuscript.
P6.confint.df

## Here I take the P6 dataset and plot the proportion of dead cells for each replicate. MAIN TEXT
## Standard error of the mean = SEM = S/√N = 0.035
## t(α, N-1) = 2.776
## R's way (fun.data) of calculating confidence interval = m +/- (t(α, N-1)*SEM). Need to calculate confidence 
## interval of the proportion. 
old.Fig3C <- P6 %>% 
    ggplot(aes(x= Treatment, y = PropPrime.Dead)) +
    geom_point() +
    facet_wrap(~ Clone) +
    ylab("Proportion dead cells") +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = .15, size = .4) +
    ##line of code defaults to mean_se. This is what I wanted to show. 
    geom_errorbar(stat = "summary", position = "identity", fun.data = mean_ci, width = .050)


## Use the mean proportions to do statistics 
## What is the proportion of death in DM0 and DM25 environments within a clone?
## Is there more death in DM0 or DM25?
## How does death of each of the evolved strains in DM25 compare to that of the ancestor?
## How does death in DM0 and DM25 environments differ between clones?

## geom_errorbar(aes(ymin=PropPrime.Dead - marg.err, ymax=PropPrime.Dead+ marg.err), position = "dodge", colour="black", size=.2, width=.04) 
## the line above is code that shows the dispersion around the mean for each replicate. 
## Allows assessment of homogeneity in variance around mean. 

####################Points to make from visualizing plot P6##########################
## In most cases there is less death in DM25 than in DM0. 
## Most evolved population, 11364 appears to converged on equal proptions of dead cells and similary the common ancestor
## of 871 and 910. 
## 871 displays significantly more death in DM0 although it evolved in this environment, even so it still has a lower
## Proporstion of death when compared to 151.

## Here I calculate the mean and standard deviation of the six samples. Output in console. 
## Shows the values of the means  plotted above in plot 6. 
## P6 %>% 
## group_by(Clone, Treatment) %>% 
## summarise(Mean.Prop.Dead= mean(PropPrime.Dead), StdDev.Prop.Dead = sd(PropPrime.Dead))

## Lets do some statistics! 
 
## One-Way ANOVA test groups: Clone, Treatment. 
## I could filter the data sets by Treatment (will have one for DM0 and the other for DM25) and perform stats on it,
## or I may be able to "group_by" before performing the stats. 

## Here I create a summary dataframe. Not saved as object. 
## P6 %>% group_by(Clone, Treatment) %>% 
##  summarise(
##    count = n(),
##    mean = mean(PropPrime.Dead, na.rm = TRUE),
##    sd = sd(PropPrime.Dead, na.rm = TRUE )
##  )

########Lets do some statistical tests to address some of the questions presented above.##########
##http://www.sthda.com/english/wiki/comparing-means-in-r
##https://www.sheffield.ac.uk/polopoly_fs/1.536445!/file/MASH_ANOVA_in_R.pdf
 
##Is there a significant difference in the proportion of dead cells present in each treatment?
leveneTest(PropPrime.Dead ~ Treatment, data = P6) #p = 0.9499 equal variances can be assumed. 

##histogram(PropPrime.Dead ~ Treatment, data = P6) shows what percent of the total number of cells analyzed
## comes from each treatment. 40% - DM0; 60% - DM25
 
res.aov.trt <- aov(PropPrime.Dead ~ Treatment, data = P6) 
summary(res.aov.trt)  #Significant differences between the treatments. 
##More death in DM0. I can make this interpretation directly.
 
##Plot of the residuals 
res.trt <- res.aov.trt$residuals
hist(res.trt, main="Histogram of standardised residuals for res.trt",xlab="Standardised residuals")
 
P6.DM0 <- P6 %>% filter(Treatment == "DM0")
P6.DM25 <- P6 %>%filter(Treatment == "DM25")
 
## Is there a difference in the proportion of dead cells present when clones are grown in DM0?
res.aov.DM0 <- aov(PropPrime.Dead ~ Clone, data=P6.DM0)
summary(res.aov.DM0) #There is no difference in the proportion of dead cells in DM0. p = 0.638
## The proportion of dead cells is the same across all clones when grown in DM0. 
 
## Is there a difference in the proportion of dead cells present when clones are grown in DM25?
res.aov.DM25 <- aov(PropPrime.Dead ~ Clone, data=P6.DM25)
summary(res.aov.DM25)
## At least one population differs in the proportion of death in DM25 p = 0.00767 which one?
TukeyHSD(res.aov.DM25) ## Tukey honestly significant difference test. Adjusts for multiple comparisons. 
 
## Clones with significantly different mean death in DM25. 
## REL11364 - REL606 -> 0.0257
## REL11364 - ZDBp910 -> 0.028
 
## Paired T.tests for each clone across environments. 

## Filtering dataset by Treatment to allow paired comparisons
## 1) Constructing the dataframes:
## ZDBp910 grew in 1/5 replicates in DM0. Can't run a paired t.test with this clone. 
## REL606 does not grow in DM0. Cant run a paired t.test with this clone. 


P6.DM0.2 <- P6 %>% ungroup() %>% 
    filter(Treatment == "DM0") %>% 
    filter(Clone != "ZDBp910") %>% 
    dplyr::select(Clone, Treatment,PropPrime.Dead)


P6.DM25.2 <- P6 %>% ungroup() %>% 
    filter(Treatment == "DM25") %>% 
    filter(Clone %in% c("CZB151", "ZDBp871", "REL11364")) %>% 
    dplyr::select(Clone, Treatment, PropPrime.Dead)

 
P6.2 <- rbind(P6.DM0.2, P6.DM25.2) 
 
##I now need to create separate dataframes for each of the clones that I will be comparing using t.tests. 
##Grouping variable will be "Treatment."
 
P6.151 <- P6.2 %>% filter(Clone == "CZB151") 
P6.871 <- P6.2 %>%  filter(Clone == "ZDBp871")
P6.11364 <- P6.2 %>% filter(Clone == "REL11364")
 
 
##Shapiro-Wilk Normality test: are the differences normally distributed. If so, I can continue with the t.test,
##if not, I must use a non-parametric alternative to the t.test.

SWNT.151 <- with(P6.151, 
                 PropPrime.Dead[Treatment == "DM0"] - PropPrime.Dead[Treatment == "DM25"])
 
SWNT.871 <- with(P6.871, 
                 PropPrime.Dead[Treatment == "DM0"] - PropPrime.Dead[Treatment == "DM25"])
 
SWNT.11364 <- with(P6.11364, 
                   PropPrime.Dead[Treatment == "DM0"] - PropPrime.Dead[Treatment == "DM25"])
 
 
shapiro.test(SWNT.151)  ##p = 0.1309
shapiro.test(SWNT.871)  ##p = 0.8683
shapiro.test(SWNT.11364) ##p = 0.6543
 
##The results of the Shapiro-Wilk Normality tests indicates that the distribution of the differences in the proportion 
## of dead cells in DM0 and DM25 environments are not significantly different from a normal distribution, 
##thus, we can assume normality and compare mean differences for each clone using the parametric t.test. 

## T.tests 
res.151 <- t.test(PropPrime.Dead ~ Treatment , data = P6.151, paired = TRUE)
res.151 ##p = 0.7753; No differences 

res.871 <- t.test(PropPrime.Dead ~ Treatment, data = P6.871, paired = TRUE) 
res.871 ##p = 0.02585
## Result indicates that there is a differecene in the proporition of dead cells between the two environments.
## Looking at the data plotted from dataframe P6, it is likely that the proporiton of dead cells is higher in DM0
## when compared to DM25. I test this hypothesis in the next line of code. 

res.871.G <- t.test(PropPrime.Dead ~ Treatment, data = P6.871, paired = TRUE, alternative = "greater") 
res.871.G ##Yes, p=0.01293, the propotion of death in DM0 is greater than that in DM25

res.11364 <- t.test(PropPrime.Dead ~ Treatment, data = P6.11364, paired = TRUE)
res.11364 ##p = 0.7651; No differences
