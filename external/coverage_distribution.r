###
##
## NAME
##
## coverage_distribution.r
##
## AUTHORS
##
## Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
## David B. Knoester
##
## LICENSE AND COPYRIGHT
##
## Copyright (c) 2011-20122010 Michigan State University
##
## breseq is free software; you can redistribute it and/or modify it under the terms the 
## GNU General Public License as published by the Free Software Foundation; either 
## version 1, or (at your option) any later version.
##
###

## Arguments:
##   distribution_file=/path/to/input 
##   plot_file=/path/to/output 
##   deletion_propagation_pr_cutoff=float
##   plot_poisson=0 or 1
##   pdf_output=0 or 1

par(family="sans")

plot_poisson = 0;
pdf_output = 1;

this.print.level = 0
#this.print.level = 2

for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
 #   temp = as.numeric(temp) #Im only inputting numbers so I added this to recognize scientific notation
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
  } else {
    assign(ta[[1]][1],TRUE)
    cat("assigned ",ta[[1]][1]," the value of TRUE\n")
  }
}

deletion_propagation_pr_cutoff = as.numeric(deletion_propagation_pr_cutoff);

## initialize values to be filled in
nb_fit_mu = 0
nb_fit_size = 0
m = 0
v = 0
D = 0
deletion_propagation_coverage = -1

#load data
X<-read.table(distribution_file, header=T)

#table might be empty
if (nrow(X) == 0)
{
  #print out statistics
  
  print(nb_fit_size);
  print(nb_fit_mu);
  
  print(m)
  print(v)
  print(D)
  
  print(deletion_propagation_coverage)
  
  q()
}

#create the distribution vector and fit
Y<-rep(X$coverage, X$n)
m<-mean(Y)
v<-var(Y)
D<-v/m

###
## Smooth the distribution with a moving average window of size 5
## so that we can more reliably find it's maximum value
###

ma5 = c(1, 1, 1, 1, 1)/5;

## filtering fails if there are too few points
if (nrow(X) >= 5) {
  X$ma = filter(X$n, ma5)
} else {
	X$ma = X$n
}

i<-0
max_n <- 0;
min_i <- max( trunc(m/4), 1 ); #prevents zero for pathological distributions
max_i <- i;
for (i in min_i:length(X$ma))
{		
  #cat(i, "\n")
	if (!is.na(X$ma[i]) && (X$ma[i] > max_n))
	{
		max_n = X$ma[i];
		max_i = i;
	}
}

##
# Censor data on the right and left of the maximum
##

start_i = max(floor(max_i*0.5), 1);
end_i = min(ceiling(max_i*1.5), length(X$ma));

if (start_i == end_i)
{
  print(nb_fit_size);
  print(nb_fit_mu);
  
  print(m)
  print(v)
  print(D)
  
  print(deletion_propagation_coverage)
  
  q()
}

cat("Fitting from coverage of ", start_i, " to ", end_i, ".\n", sep="")

##
# Coarse grain so that we are only fitting a number of bins that is 1000-2000
#
# The later adjustment for doing the fits this way is to multiply the means
# of the negative binomial and poisson distributions by the binning number.
# (The size parameter of the negative binomial doesn't need to be adjusted.)
##


num_per_bin = trunc((end_i - start_i) / 1000)

if (num_per_bin > 1) 
{
  cat("Coarse-graining for fits\n")
  start_i_for_fits = trunc(start_i/num_per_bin)
  end_i_for_fits = ceiling(end_i/num_per_bin)
  num_bins = end_i - start_i  + 1
  cat("Fitting from coverage in adjusted bins ", start_i_for_fits, " to ", end_i_for_fits, ".\n", sep="")
  cat("Number of bins ", num_bins, ". Each bin has ", num_per_bin, " coverage values.\n", sep="")

  # Create a new vector where we've added together values in bins
  X.for.fits = vector("double", end_i_for_fits)
  for (i in start_i_for_fits:end_i_for_fits)
  {
    for (j in 1:num_per_bin)
    {
      if (i*num_per_bin+j <= length(X$n))
      {
        X.for.fits[i] = X.for.fits[i] + X$n[i*num_per_bin+j]
      }
    }
  }

} else {
  ## AVOID num_per_bin equalling zero!!
  X.for.fits = X$n[1:end_i]
  num_per_bin = 1
  start_i_for_fits = start_i
  end_i_for_fits = end_i
}


##
# Now perform negative binomial fitting to the censored data
##

inner_total<-0;
for (i in start_i_for_fits:end_i_for_fits)
{
	inner_total = inner_total + X.for.fits[i]; 
}
# Yes: it's correct to use X here because we want the overall total total
total_total<-sum(X$n);

## let's preconstruct these for speed
dist = vector("double", end_i_for_fits)

f_nb <- function(par) {

	mu = par[1];
	size = par[2];

  if ((mu <= 0) || (size <= 0))
  {
    return(0);
  }
  
  cat(start_i_for_fits, " ", end_i_for_fits, "\n");
  cat(mu, " ", size, "\n");
  
	dist<-c()
	total <- 0;
	for (i in start_i_for_fits:end_i_for_fits)
	{	
		dist[i] <- dnbinom(i, size=size, mu=mu);
		total <- total + dist[i] 
	}
	#print (mu, size)

 	l <- 0;
	for (i in start_i_for_fits:end_i_for_fits)
	{
		l <- l + ((X.for.fits[i]/inner_total)-(dist[i]/total))^2;
	}
	return(l);
}



## Fit negative binomial 
## - allow fit to fail and set all params to zero/empty if that is the case
nb_fit = NULL
## as.numeric prevents overflow in sums involving integers
mean_estimate = sum((as.numeric(1:end_i_for_fits)*as.numeric(X.for.fits)))/sum(as.numeric(X.for.fits))

nb_fit_mu = -1
nb_fit_size = -1
try_size = 100000
nb_fit = c()

while ( ((nb_fit_mu < 0) || (nb_fit_size < 0) || (nb_fit$code != 1)) && (try_size > 0.001))
{
  try_size = try_size / 10

  ## SIZE ESTIMATE from the censored data can be negative, so try various values instead
  cat("Try Mean: ", mean_estimate, " Size: ", try_size, "\n")

  nb_fit<-nlm(f_nb, c(mean_estimate, try_size), iterlim=1000, print.level=this.print.level )

  nb_fit_mu = nb_fit$estimate[1];
  nb_fit_size = nb_fit$estimate[2];

  cat("Fit Mean: ", nb_fit_mu, " Size: ", nb_fit_size, " Code: ", nb_fit$code, "\n")
}

cat("Final Fit Mean: ", nb_fit_mu, " Size: ", nb_fit_size, " Code: ", nb_fit$code, " Try Size: ", try_size, "\n")

if ((nb_fit_mu < 0) || (nb_fit_size < 0) || (nb_fit$code != 1))
{
  print(0);
  print(0);

  print(m)
  print(v)
  print(D)

  print(deletion_propagation_coverage)

  q()
}


## things can go wrong with fitting and we can still end up with invalid values

fit_nb = c()
if (nb_fit_mu > 0)
{
  end_fract = pnbinom(end_i_for_fits, mu = nb_fit_mu, size=nb_fit_size)
  start_fract = pnbinom(start_i_for_fits, mu = nb_fit_mu, size=nb_fit_size)
  included_fract = end_fract-start_fract;

  ## Adjust so that we are back in full coords before making fit!!
  if (num_per_bin > 1) 
  {
    nb_fit_mu = nb_fit_mu * num_per_bin
  }
  fit_nb = dnbinom(0:max(X$coverage), mu = nb_fit_mu, size=nb_fit_size)*inner_total/included_fract;
}

f_p <- function(par) {

  lambda = par[1];

  if (lambda <= 0)
  {
    return(0);
  }
  
	total <- 0;
	for (i in start_i_for_fits:end_i_for_fits)
	{	
    #cat(i, " ", lambda, "\n");
		dist[i] <- dpois(i, lambda=lambda);
		total <- total + dist[i] 
	}
	#print (total)

 	l <- 0;
	for (i in start_i_for_fits:end_i_for_fits)
	{
		l <- l + ((X.for.fits[i]/inner_total)-(dist[i]/total))^2;
	}
	return(l);
}


## Fit Poisson 
## - allow fit to fail and set all params to zero/empty if that is the case

p_fit = NULL
try(p_fit<-nlm(f_p, c(m), print.level=this.print.level))

fit_p = c()
if (!is.null(p_fit) && (p_fit$estimate[1] > 0))
{
  #print (nb_fit$estimate[1])
  p_fit_lambda = p_fit$estimate[1];
  #print(0:max(X$coverage))

  end_fract = ppois(end_i_for_fits, lambda = p_fit_lambda)
  start_fract = ppois(start_i_for_fits, lambda = p_fit_lambda)
  included_fract = end_fract-start_fract;

  ## Adjust so that we are back in full coords before making fit!!
  if (num_per_bin > 1) 
  {
    p_fit_lambda = p_fit_lambda * num_per_bin
  }
  fit_p<-dpois(0:max(X$coverage), lambda = p_fit_lambda)*inner_total/included_fract;
}


## Graphing
##
## don't graph very high values with very little coverage
i<-max_i
while (i <= length(X$n) && X$n[i]>0.01*max_n)
{		
	i <- i+1;
}
graph_end_i <-i

## Ths leaves enough room to the right of the peak for the legend
graph_end_i = max(floor(2.2 * max_i), graph_end_i);

## graphics settings
my_pch = 21
my_col = "black";
my_col_censored = "red";

if (pdf_output == 0) {
	## units = "px", taa=4, gaa=2 NOT compatible with earlier R versions!
	bitmap(plot_file, height=6, width=7, type = "png16m", res = 72, pointsize=18)
} else {
	pdf(plot_file, height=6, width=7)
}

par(mar=c(5.5,7.5,3,1.5));

max_y = 0
if (plot_poisson) {
	max_y = max(X$n, fit_p, fit_nb)
} else {
	max_y = max(X$n, fit_nb)
}

plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, max_y)*1.05, xlim=c(0, graph_end_i), lwd=1, xaxs="i", yaxs="i", axes=F, las=1, main="Coverage Distribution at Unique-Only Positions", xlab="Coverage depth (reads)", ylab="", cex.lab=1.2, cex.axis=1.2)

mtext(side = 2, text = "Number of reference positions", line = 5.5, cex=1.2)

sciNotation <- function(x, digits = 1) {
    if (length(x) > 1) {
        return(append(sciNotation(x[1]), sciNotation(x[-1])))     
	} 
    if (!x) return(0) 

	exponent <- floor(log10(x)) 
    base <- round(x / 10^exponent, digits)     
	as.expression(substitute(base %*% 10^exponent, list(base = base, exponent = exponent))) 
}

#axis(2, cex.lab=1.2, las=1, cex.axis=1.2, labels=T, at=(0:6)*50000)
axis(2, cex.lab=1.2, las=1, cex.axis=1.2, at = axTicks(2), labels = sciNotation(axTicks(2), 1))
axis(1, cex.lab=1.2, cex.axis=1.2, labels=T)
box()

#graph the coverage as points
fit_data <- subset(X, (coverage>=start_i) & (coverage<=end_i) );
points(fit_data$coverage, fit_data$n, pch=my_pch, col=my_col, bg="white", cex=1.2)

#graph the censored coverage as red points
cat(start_i, " ", end_i, "\n", sep="")

censored_data <- subset(X, (coverage<start_i) | (coverage>end_i) );
points(censored_data$coverage, censored_data$n, pch=my_pch, col=my_col_censored, bg="white", cex=1.2)

#graph the poisson fit IF REQUESTED
if (plot_poisson) {
	lines(0:max(X$coverage), fit_p, lwd=3, lty="22", col="black");
}

#graph the negative binomial fit
lines(0:max(X$coverage), fit_nb, lwd=3, col="black");

if (plot_poisson) {
	legend("topright", c("Coverage distribution", "Censored data", "Negative binomial", "Poisson"), lty=c("blank","blank","solid","22"), lwd=c(1,1,2,2), pch=c(my_pch, my_pch, -1, -1), col=c("black", "red", "black", "black"), bty="n")
} else {
	legend("topright", c("Coverage distribution", "Censored data", "Negative binomial"), lty=c("blank","blank","solid"), lwd=c(1,1,2), pch=c(my_pch, my_pch, -1), col=c("black", "red", "black"), bty="n")
}

dev.off()

## Fit the marginal value that we use for propagating deletions

if (nb_fit_size != 0)
{
  cat(nb_fit_size, " ", nb_fit_mu, "\n")
  deletion_propagation_coverage = qnbinom(deletion_propagation_pr_cutoff, size = nb_fit_size, mu = nb_fit_mu)
}

#print out statistics

print(nb_fit_size);
print(nb_fit_mu);

print(m)
print(v)
print(D)

print(deletion_propagation_coverage)

warnings()
