#The program uses the libraries ape and  geiger they can be installed with install.package("ape") and install.package("geiger") commands
#To produce all figures simply run : Rscript ScriptFiguresDryad.r inputfolder ouputfolder
#The script for each figure can also be run independently, just define input and ouput files and copy paste the script in R.

#Several files are required
#
#MutationTypesThroughTime.txt for figures 1, 3, 4, S2, S3, S4, S5, S6, S9 produced by ComputeMutationThroughTimeDryad.pl
#MutRArray.txt for figure 2 (phylogeny) produced by ComputeMutationThroughTimeDryad.pl
#GenomeComposition.txt and MutationTypesThroughTimeMAE.txt for figures 4, S6, S9 produced by GenomeCompositionComputer.pl and ComputeMutationThroughTimeDryad.pl



args <- commandArgs(trailingOnly=TRUE)
generalinputfile<-args[1]
generaloutputfile<-args[2]

library(ape)
library(geiger)
################################################################################################################################################
######################################################    set colors    ########################################################################
################################################################################################################################################
############ A simple graph to choose the colors for the different populations##################################################################



popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")

rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)
plot(NA,xlim=c(0.5,6.5),ylim=c(-1.5,1.5),xlab="",ylab="")

for ( i in c(0:5))
{
    for (j in c(0:1))
    {
        x00<-i
        y00<-j
              polygon(c(x00,x00,x00+1,x00+1)+0.5,c(y00-1,y00,y00,y00-1)*2,col=popcol2s[(i*2+(j+1))])
    }
    
}




################################################################################################################################################
######################################################    Figure 1    ##########################################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt
outputfile<-generaloutputfile
inputfile<-generalinputfile


#Composition of the genome
#taking into account overlappinf reading frames and start degeneracy and getting rid of non coding genes and pseudogenes
#mutations are sorted into AT->CG AT->GC AT->TA CG->GC CG->TA GC->TA

pdf(paste(outputfile,"Figure1.pdf",sep=""))
nf <- layout(matrix(c(1,3,2,3), 2, 2), height=c(7/10,3/10),respect=TRUE)
#par(mfrow=c(1,2))
par(mar=c(4,5,2,2))
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)



typemuts<-cbind(typemuts,typemuts[,8])
#Compute the total number of mutations
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]


#start the plot
plot(NA,xlim=c(0,50),ylim=c(0,2600),ylab="", xlab="",main="",cex.main=0.8,cex.axis=0.9,las=1,cex.lab=1,yaxt="n",cex.lab=1.2,cex=0.8)
axis(side=2,at=c(0,500,1000,1500,2000,2500),labels=c("0","500","1,000","1,500","2,000","2,500"),las=2,cex=0.9)

title(ylab="Number of mutations", line=3.5, cex.lab=1.2)
title(xlab="Time (thousands of generations)", line=2.9, cex.lab=1.2)

means=c()
mtext("a",3,1,adj=0,at=c(-15),cex=1.5)


#define colors an dtransparency for teh symbols
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

couleur=1

#for each population and each timepoint, compute the mean number of mutations for both sequenced clones of the populations and plot datas
for (pop in pops)
{
    temp<-typemuts[typemuts[,2]==pop,]
    points(temp[,3]/1000,temp[,9],col=popcol2s[couleur],pch=20,cex=1.5)
    means=c()
    count=1
    for (i in c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000))
    {
        means[count]=mean(temp[temp[,3]==i ,9])
        count=count+1;
    }
    
    lines(c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)/1000,means,col=popcols[couleur])
    couleur<-couleur+1
}




#Subsample non mutator clones
means=c()
popcolsnm<-popcols[c(2,4,8,9,10,11)]
rgbcolsnm<-col2rgb(popcolsnm)
popcol2nms<-rgb(t(rgbcolsnm),alpha=150, maxColorValue=255)
pops<-pops[c(2,4,8,9,10,11)]



#second panel

plot(NA,xlim=c(0,50),ylim=c(0,120),ylab="Number of mutations", xlab="",main="",cex.main=0.8,cex.axis=0.9,las=1,cex.lab=1.2)
mtext("b",3,1,adj=0,at=c(-15),cex=1.5)

title(xlab="Time (thousands of generations)", line=2.9, cex.lab=1.2)

pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")

popswitch<-c(25000,   100000,  3000,  100000,  48000,  3000,  8000,   100000,  100000, 100000, 100000,6000)
typemuts<-typemuts[typemuts[,1]!="Aram3_40000_Aram3_40000gen_10979",]

#for each population and each timepoint before mutator swtich, compute the mean number of mutations for both sequenced clones of the populations and plot datas

couleur=1
for (pop in pops)
{
    temp<-typemuts[typemuts[,2]==pop & typemuts$age<popswitch[couleur],]
    points(temp[,3]/1000,temp[,9],col=popcol2s[couleur],pch=20,cex=1.5)
    means=c()
    count=1
    for (i in c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000))
    {
        means[count]=mean(temp[temp[,3]==i ,9])
        count=count+1;
    }
    
    lines(c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)/1000,means,col=popcols[couleur])
    couleur<-couleur+1
}

#plot the color labels
par(mar=c(0,0,2,0))
plot(NA,xlim=c(0,7),ylim=c(0,2),xaxt="n",yaxt="n",frame.plot=FALSE)
yy1<-1.5
yy2<-0.7
delta=1
dep=0.8
text(dep,yy1,"Ara-1",cex=1.2)
text(dep,yy2,"Ara+1",cex=1.2)
text(dep+delta,yy1,"Ara-2",cex=1.2)
text(dep+delta,yy2,"Ara+2",cex=1.2)
text(dep+2*delta,yy1,"Ara-3",cex=1.2)
text(dep+2*delta,yy2,"Ara+3",cex=1.2)
text(dep+3*delta,yy1,"Ara-4",cex=1.2)
text(dep+3*delta,yy2,"Ara+4",cex=1.2)
text(dep+4*delta,yy1,"Ara-5",cex=1.2)
text(dep+4*delta,yy2,"Ara+5",cex=1.2)
text(dep+5*delta,yy1,"Ara-6",cex=1.2)
text(dep+5*delta,yy2,"Ara+6",cex=1.2)
deltabis=0.45
points(c(0,0,1,1,2,2,3,3,4,4,5,5)+dep+deltabis,rep(c(yy1,yy2),6),col=popcol2s,pch=20,cex=2)

dev.off()



################################################################################################################################################
######################################################    Figure 2    ##########################################################################
################################################################################################################################################
#After controlling that, using point mutations, the phylogeny was star like, phylogenies are computed independently for each population and then combined to have clones and population sorted.
#packages ape, geiger and mapplot are required they can be installed with install.package("ape") like command.
#the input file is MutRArray.txt


library(ape)
library(geiger)


outputfile<-generaloutputfile
inputfile<-generalinputfile

MutationTableTotal<-read.table(paste(inputfile,"MutRArray.txt",sep=""))

#sort pop by numbers
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")

#A
treelist<-list()
for (p in pops)
{
    MutationPop<-MutationTableTotal[MutationTableTotal[,1]=="Pop"| MutationTableTotal[,1]=="Ancestor"| MutationTableTotal[,1]==p,]
    MutationPop=t(MutationPop)
    l=length(MutationPop[1,])
    
    #this step filters the mutations that are only present in subsample (useless here)
    MutationPop<-cbind(MutationPop,MutationPop[,3])
    
    MutationPop[,l+1]<-apply(matrix(as.numeric(MutationPop[,c(5:l)]),ncol=(l-4)),1,sum)
    
    
    MutationPopFiltered<-MutationPop# here the filtering is done
    
    MutationPopFiltered[1,]<-MutationPop[1,]
    MutationPopFiltered[2,]<-MutationPop[2,]
    MutationPopFiltered[3,]<-MutationPop[3,]
    MutationPopFiltered[4,]<-MutationPop[4,]
    
    MutationPopFiltered<-MutationPopFiltered[,c(1:l)]
    #isolate the information related to strains
    MutationTransit<-MutationPopFiltered[, c(4:l)]
    #sort according to alphabetic order of strains
    MutationTransit<-MutationTransit[,order(MutationTransit[4,])]
    
    #convert data to a fasta and create teh phylogenny using JC69 model with a fast minimum evolution algorithm
    ll<-length(MutationPopFiltered[,1])
    TreeMatrix<-t(MutationPopFiltered[c(4:ll), c(4:l)])
    TreeMatrix<-replace(TreeMatrix, TreeMatrix=="1","T")
    TreeMatrix<-replace(TreeMatrix, TreeMatrix=="0","A")
    ll<-length(TreeMatrix[,1])
    lll<-length(TreeMatrix[1,])
    AlingmentMatrix<-(TreeMatrix[,c(2:length(TreeMatrix[1,]))])
    dimnames(AlingmentMatrix)[[1]]<-TreeMatrix[,1]
    alntemp<-as.DNAbin(tolower(AlingmentMatrix))
    distmat<-dist.dna(alntemp,model = "JC69")
    
    tree.tmp<-fastme.bal(distmat)
    tree.tmp<-root(tree.tmp,outgroup<-"Ancestor_0_REL606",resolve.root=TRUE)
    
    treelist[[length(treelist)+1]] <-list(pop=p,arbre=tree.tmp)
    
}

#combine the trees
treesum<-treelist[[1]]$arbre

for (i in c(2:12))
treesum<-treesum+treelist[[i]]$arbre

#get rid of the multiple ancestor
i<-2
while(i<length(treesum$tip))
{
    if (treesum$tip[i]=="Ancestor_0_REL606")
    {
        treesum<-drop.tip(treesum,i)
    }
    i<-i+1
}

treesum<-rotate(treesum,266)

#label mutator branches
tree<-treesum #root(treesum,"Ancestor_0_REL606",r=TRUE)
cladep6<-tips(tree,mrca(tree)["Arap6_40000gen_10985","Arap6_10000gen_4535A"])#
cladep3<-tips(tree,mrca(tree)["Arap3_5000gen_2175B","Arap3_50000gen_11345"])#
cladem1<-tips(tree,mrca(tree)["Aram1_40000gen_10938","Aram1_30000gen_10392"])#
#cladem2<-tips(tree,mrca(tree)["Aram2_5000gen_2180B","Aram2_50000gen_11335"])#
cladem2<-tips(tree,mrca(tree)["Aram2_5000gen_2180A","Aram2_50000gen_11335"])#

cladem3<-tips(tree,mrca(tree)["Aram3_40000gen_10979","Aram3_50000gen_11364"])
cladem4<-tips(tree,mrca(tree)["Aram4_10000gen_4539B","Aram4_10000gen_4539A"])#

ancestorp3<-which.edge(tree,mrca(tree)["Arap3_5000gen_2175B","Arap3_50000gen_11345"])
ancestorp6<-which.edge(tree,mrca(tree)["Arap6_40000gen_10985","Arap6_10000gen_4535A"])
ancestorm3<-which.edge(tree,mrca(tree)["Aram3_40000gen_10979","Aram3_50000gen_11364"])
ancestorm1<-which.edge(tree,mrca(tree)["Aram1_40000gen_10938","Aram1_30000gen_10392"])
#ancestorm2<-which.edge(tree,mrca(tree)["Aram2_5000gen_2180B","Aram2_50000gen_11335"])
ancestorm2<-which.edge(tree,mrca(tree)["Aram2_5000gen_2180A","Aram2_50000gen_11335"])
ancestorm4<-which.edge(tree,mrca(tree)["Aram4_10000gen_4539B","Aram4_10000gen_4539A"])

#Modify the branch length of mutators
branchlistmutT<-c(ancestorp6,ancestorm1, which.edge(tree, cladep6),which.edge(tree, cladem1))
branchlistmutS<-c(ancestorp3, ancestorm3,ancestorm2 ,ancestorm4,    which.edge(tree, cladep3),which.edge(tree, cladem2),which.edge(tree, cladem3),which.edge(tree, cladem4))

tree$edge.length[branchlistmutT]<-tree$edge.length[branchlistmutT]/50
tree$edge.length[branchlistmutS]<-tree$edge.length[branchlistmutS]/25

#modify color and width of branches
edgecoul<-rep(1,length(tree$edge.length))
edgewidth<-rep(1,length(tree$edge.length))
edgecoul[branchlistmutT]<-"red"
edgecoul[branchlistmutS]<-"orange"
edgewidth[branchlistmutT]<-2
edgewidth[branchlistmutS]<-1.5


#define the color and symbol for each pop
timestext<-c("_0_","_500gen_","_1000gen_","_1500gen_","_2000gen_","_5000gen_","_10000gen_","_15000gen_","_20000gen_","_30000gen_","_40000gen_","_50000gen_")
times<-c(1,500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)
levels<-c(120,40,50,60,70,80,88,96,102,108,114,120)
popcols<-c("black","dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")

popslabels<-c("REL606","Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")

pops<-c("Ancestor","Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popswitch<-c(25000,100000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

coul<-cbind(tree$tip.label,rep(0,length(tree$tip.label)),rep(0,length(tree$tip.label)),rep(50000,length(tree$tip.label)),rep(0.6,length(tree$tip.label)))
l<-length(coul[,1])
root<-"Ancestor_0_REL606"
for (i in c(1:13))
{
    #GG<-colorRampPalette(c("white", popcol2s[i]),space="Lab")(121)
    for (j in c(1:12))
    {
        #  if (length(coul[grepl(pops[i],coul[,1]) && grepl(times[j],coul[,1]),2])>0)
        #coul[grepl(pops[i],coul[,1]) && grepl(times[j],coul[,1]),2]<-GG[levels[j]]
        for (k in c(1:l))
        {
            if (grepl(pops[i],coul[k,1]) && grepl(timestext[j],coul[k,1]))
            {
                
                coul[k,5]=0.4
                if (root==coul[k,1])
                {
                    print (coul[k,])
                    coul[k,5]<-0.6
                    coul[k,2]<-"black"
                    print (coul[k,])
                }
                coul[k,2]<-popcol2s[i] #
                # Agelimits<-rbind(Agelimits,c(tree$edge[which.edge(tree,coul[k,1]),1],0,times[j]))
                #  Agelimits[k,2]<-0
                # Agelimits[k,3]<-times[j]
            }
        }
        
        
    }
    
}


par(mar=c(0,0,0,0))
#root tree
tree2<-treesum#root(treesum,"Ancestor_0_REL606",r=TRUE)

pdf(paste(outputfile,"Figure2.pdf",sep=""))
par(mar=c(0,1,0,0))

par(mfrow=c(2,1),oma=c(1,1,1,1),mar=c(2,2,2,2))
plot(tree2,type="p",direction="up",show.tip.label=FALSE,edge.col=edgecoul)

#tree2<-multi2di(tree2)
#edgewidth<-rep(1,length(tree2$edge.length))

#plot(tree2,type="p",direction="up",show.tip.label=FALSE,edge.col=edgecoul)

size<-0.8

text(12,0.087,label="Ara-1",cex=size,col=popcols[2])
    text(35,0.012,label="Ara+1",cex=size,col=popcols[3])#207
    
    text(53,0.0760,label="Ara-2",cex=size,col=popcols[4])#260
    #text(255,0.071,label="S",cex=size,col=popcols[3])
    #text(265,0.071,label="L",cex=size,col=popcols[3])
    text(77,0.012,label="Ara+2",cex=size,col=popcols[5])#166
    
    text(102,0.056,label="Ara-3",cex=size,col=popcols[6])
    text(122,0.118,label="Ara+3",cex=size,col=popcols[7])
    
    text(145,0.091,label="Ara-4",cex=size,col=popcols[8])
    text(166,0.012,label="Ara+4",cex=size,col=popcols[9])
    
    text(187,0.012,label="Ara-5",cex=size,col=popcols[10])
    text(209,0.012,label="Ara+5",cex=size,col=popcols[11])
    
    text(233,0.012,label="Ara-6",cex=size,col=popcols[12])
    text(254,0.198,label="Ara+6",cex=size,col=popcols[13])

#text(0,0.197,label="a",cex=1.5)
#text(1,0.008,label="REL606",cex=0.5)
mtext("REL606",at=0,side=2,cex=size,las=1,line=-0.5)
mtext("a",at=0.197,side=2,cex=1.5,las=1,line=1)

coul[,5]<-0.6
tiplabels(pch = 19, col = coul[,2], adj = c(0.5,0.5), cex = as.numeric(coul[,5]))

plot(tree,type="p",direction="up",show.tip.label=FALSE,edge.col=edgecoul,edge.width=edgewidth)

coul[,5]<-0.8
tiplabels(pch = 19, col = coul[,2], adj = c(0.5,0.5), cex = as.numeric(coul[,5]))
size<-0.8
text(12,0.0072/2,label="Ara-1",cex=size,col=popcols[2])
    text(37,0.0072/2,label="Ara+1",cex=size,col=popcols[3])#207
    
    text(53,0.0064/2,label="Ara-2",cex=size,col=popcols[4])#260
    #text(255,0.071,label="S",cex=size,col=popcols[3])
    #text(265,0.071,label="L",cex=size,col=popcols[3])
    text(77,0.0059/2,label="Ara+2",cex=size,col=popcols[5])#166
    
    text(102,0.0068/2,label="Ara-3",cex=size,col=popcols[6])
    text(122,0.0091/2,label="Ara+3",cex=size,col=popcols[7])
    
    text(145,0.0078/2,label="Ara-4",cex=size,col=popcols[8])
    text(166,0.00640/2,label="Ara+4",cex=size,col=popcols[9])
    
    text(186,0.0079/2,label="Ara-5",cex=size,col=popcols[10])
    text(210,0.0068/2,label="Ara+5",cex=size,col=popcols[11])
    
    text(233,0.0062/2,label="Ara-6",cex=size,col=popcols[12])
    text(254,0.0088/2,label="Ara+6",cex=size,col=popcols[13])
    
    
# text(-1,0.0006/2,label="REL606",cex=0.7)
#text(1,0.0089/2,label="b",cex=1.5)
mtext("REL606",at=0,side=2,cex=size,las=1,line=-0.5)
mtext("b",at=0.0089/2,side=2,cex=1.5,las=1,line=1)
dev.off()

################################################################################################################################################
######################################################    Figure 3    ##########################################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt

outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1

timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)

#get rid of IS150 mutator
arap1_and_aram530k_out<-1

pdf(paste(outputfile,"Figure3.pdf",sep=""))
#nf <- layout(matrix(c(1,), 2, 2), height=c(9.9/10,0.1/10),respect=TRUE)
par(mfrow=c(1,1),mar=c(5,5,4,4))

#get data and compute total number of mutations
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)
typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]


#plot
plot(NA,xlim=c(0,50),ylim=c(0,100),ylab="Number of mutations", xlab="Time (thousands of generations)",main="",cex.main=0.8,cex.axis=0.9,las=1,cex.lab=1.2)
means=c()
mtext("b",3,1,adj=0,at=c(-15),cex=1.5)
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popswitch<-c(25000,   100000,  3000,  100000,  38000,  3000,  8000,   100000,  100000, 100000, 100000,6000)

if (arap1_and_aram530k_out)
{ popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}


rgbcols<-col2rgb(popcols)
#rgbcols<-rbind(rgbcols,rep(0.8,length(rgbcols[1,])))
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#plot for each pop and compute the grandmean
grandmeans<-timessampling*0
grandmeanscounts<-grandmeans
dataforfit<-c()
couleur=1
for (pop in pops)
{
    temp<-typemuts[typemuts[,2]==pop & typemuts$age<popswitch[couleur],]
    print (c(temp[,3],temp[,9]))
    points((temp[,3]/1000),temp[,9],col=popcol2s[couleur],pch=20,cex=1.5)
    means=c()
    count=1
    for (i in timessampling[timessampling<popswitch[couleur]])
    {
        means[count]=mean(temp[temp[,3]==i ,9])
        grandmeans[count]<-grandmeans[count]+means[count]
        grandmeanscounts[count]<-grandmeanscounts[count]+1
        dataforfit<-rbind(dataforfit,c(pop,i,means[count]))
        count=count+1;
    }
    
    
    #lines((timessampling[timessampling<popswitch[couleur]]/1000),means[timessampling<popswitch[couleur]],col=popcol2s[couleur],pch=18,cex=1.5)
    couleur<-couleur+1
}
grandmeans<-grandmeans/grandmeanscounts
points(timessampling/1000,grandmeans,col=1,pch=2,cex=1.5)

#convert data for the fitting functions
dataforfit<-as.data.frame(dataforfit)
dimnames(dataforfit)[[2]]<-c("Pop","time","nbmut")
#lines(timessampling/1000,grandmeans,col=rgb(0.5,0.5,0.5,0.5),lwd=4,lty="dashed")
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
dataforfit$time<-as.numeric.factor(dataforfit$time)
dataforfit$nbmut<-as.numeric.factor(dataforfit$nbmut)


dataforfit<-dataforfit[order(dataforfit$time),]

#compute the 3 fits
lfit<-lm(nbmut ~0+time, dataforfit)#linear
sfit<-lm(nbmut ~0+sqrt(time), dataforfit)#square root
lsfit<-lm(nbmut ~0+time+ sqrt(time), dataforfit)#combined

print(extractAIC(lfit))#345.0781
print(extractAIC(sfit))#363.0816
print(extractAIC(lsfit))#275.4998

#draw fits
a<-lfit$coefficients[[1]]
times<-c(0:500)/10
linear_fit<-a*times*1000
lines((times),linear_fit,col=rgb(0.5,0.5,0.5,0.8),lwd=4,lty="dashed")

b<-sfit$coefficients[[1]]
sqrt_fit<-b*sqrt(times*1000)
lines((times),sqrt_fit,col=rgb(0.5,0.5,0.5,0.8),lwd=4)

a_d<-lsfit$coefficients[[1]]
b_d<-lsfit$coefficients[[2]]
combined_fit<-a_d*times*1000+b_d*sqrt(times*1000)
lines((times),combined_fit,col=rgb(0,0,0,1),lwd=4)

#compute p values
var.test(lfit,lsfit)
var.test(sfit,lsfit)

dev.off()


################################################################################################################################################
######################################################    Figure 4 ##############################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt and GenomeComposition.txt
outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1


genomecomp<-read.table(paste(inputfile,"GenomeComposition.txt",sep=""),header=TRUE)
synonymoussitespertypeofmut<-genomecomp[1,2:7]#c(253895,612753,290139,418546,753134,462589)
nonsynonymoussitespertypeofmut<-genomecomp[2,2:7]#c(1636645,1277787,1600401,1616710,1282122,1572667)
intergenicsitespertypeofmut<-genomecomp[3,2:7]#c(284614,284614,284614,201537,201537,201537)
genomicGC<-genomecomp[4,2:7]#c(2175154,2175154,2175154,2236794,2236794,2236794)

#Get data and Compute the total number of mutations
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)
typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]

#use MAE to infer basal rates
typemutsMAE<-read.table(paste(inputfile,"MutationTypesThroughTimeMAE.txt",sep=""),header=TRUE)
refrates<-typemuts[1,]
refrates[5:31]<-colSums(typemutsMAE[5:31])/(13750*15)



#define colors and transparency for the symbols
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#Subsample non mutator clones
means=c()
popcolsnm<-popcols[c(2,4,8,9,10,11)]
#popcolsnm<-popcols[c(4,8,9,10,11)]
rgbcolsnm<-col2rgb(popcolsnm)
popcol2nms<-rgb(t(rgbcolsnm),alpha=150, maxColorValue=255)

popsneutral<-pops[c(2,4,8,9,10,11)]
if (arap1_and_aram530k_out)
{
    popsneutral<-pops[c(4,8,9,10,11)]
   
}
#compute neutral rates used in subsequent figures
NEUTRALRATE<-0
synsrates<-rep(0,6)
couleur=1
for (pop in popsneutral)
{
    temp<-typemuts[typemuts[,2]==pop,]
    
    count=1
    for (i in c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000))
    {
        #means[count]=mean(temp[temp[,3]==i ,9])
        count=count+1;
        if (i==50000)
        {
            print(c(temp[temp[,3]==i,]$syn[1],sum(temp[temp[,3]==i,]$syn)))
            NEUTRALRATE<-NEUTRALRATE+sum(temp[temp[,3]==i,]$syn)
            synsrates<-synsrates+colSums(temp[temp[,3]==i,13:18])
        }
    }
    
    couleur<-couleur+1
}

nbpopsrate<-length(popsneutral)
NEUTRALRATE<-NEUTRALRATE/(50000*2*nbpopsrate)
NSRATE<-sum(synsrates/synonymoussitespertypeofmut*nonsynonymoussitespertypeofmut/(50000*2*nbpopsrate))
INTERRATE<-sum(synsrates/synonymoussitespertypeofmut*intergenicsitespertypeofmut/(50000*2*nbpopsrate))
MUTRATE<-sum(synsrates/synonymoussitespertypeofmut*genomicGC/(50000*2*nbpopsrate))



#set colors and exclure mutators
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popswitch<-c(25000,100000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)

if (arap1_and_aram530k_out)
{
    popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}
timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)
figurelabel<-c("a","b","c","d","e","f","g","h","i","j","k","l")

rgbcols<-col2rgb(popcols)
#rgbcols<-rbind(rgbcols,rep(0.8,length(rgbcols[1,])))
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#get non mutator clones
aranm<-c()
nbpop<-1
for (p in pops)
{
    age<-sort(unique(typemuts[typemuts[,2]==p & typemuts[,3]<popswitch[nbpop],3]))
    
    for (t in age)
    {
        temp<-typemuts[typemuts[,2]==p & typemuts[,3]==t ,]
        labeltemp<-temp[1,1:4]
        labeltemp[2]<-nbpop
        aranm<-rbind(aranm,c(labeltemp,apply(temp[,5:30],2,mean)))
    }
    nbpop<-nbpop+1
}

colnames<-dimnames(aranm)[[2]]
aranm<-as.data.frame(aranm)
new.list<-lapply(1:ncol(aranm), function(x) as.vector(unlist(aranm[,x])))
aranm<-as.data.frame(do.call(cbind, new.list), stringsAsFactors=FALSE)
aranm<-as.data.frame(aranm)
dimnames(aranm)[[2]]<-colnames

l<-length(aranm[1,])
for (i in c(2,3,5:l))
{
    aranm[,i]<-as.numeric(aranm[,i])
}




log2_noinf<-function(x)
{
    y<-x
    y[y<=0]<-0.125
    return(log2(y))
}



pdf(paste(outputfile,"Figure4.pdf",sep=""))
#par (fg="white",col="white",bg="black",col.axis="white", col.main="white",col.lab="white")

age<-sort(unique(aranm$age))
ylimmax<-7
ylimmin<-(-3)

par(mfrow=c(1,3),oma=c(13,2,13,2),mar=c(5,5.5,4,2))
title=c("Observed point mutations relative to neutral expectation","Observed synonymous mutations\nrelative to neutral expectation","Observed nonsynonymous mutations\nrelative to neutral expectation","Observed intergenic point mutations\nrelative to neutral expectation")
coul<-c("purple","red","blue","green")
figurenb<-c("a","a","b","c")

#compute the rates




rate<-cbind(aranm$age,(((aranm$s1+aranm$s2+aranm$s3+aranm$ns1+aranm$ns2+aranm$ns3+aranm$i1+aranm$i2+aranm$i3)+(aranm$s4+aranm$s5+aranm$s6+aranm$ns4+aranm$ns5+aranm$ns6+aranm$i4+aranm$i5+aranm$i6))/(aranm$age))/MUTRATE,(((aranm$s1+aranm$s2+aranm$s3)+(aranm$s4+aranm$s5+aranm$s6))/(aranm$age))/NEUTRALRATE,(((aranm$ns1+aranm$ns2+aranm$ns3)+(aranm$ns4+aranm$ns5+aranm$ns6))/(aranm$age))/NSRATE,(((aranm$i1+aranm$i2+aranm$i3)+(aranm$i4+aranm$i5+aranm$i6))/(aranm$age))/INTERRATE,aranm$pop)

ylimmin<-(-3)
ylimmax<-(5)

for (j in c(2:4))
{
    
    m_rate<-age
    se_rate<-age
    i<-1
    for (a in age)
    {
        m_rate[i]<-mean(rate[rate[,1]==a,j+1])
        l<-length(rate[rate[,1]==a,j+1])
        se_rate[i]<-sd(rate[rate[,1]==a,j+1])/sqrt(l)
        i<-i+1
    }
    plot(rate[,1],log2_noinf(rate[,j+1]),pch=20,ylim=c(ylimmin,ylimmax),main="",ylab=title[j],xlab="",axes=FALSE,col=0,cex.lab=1.2)
    
    title(xlab="Time\n(thousands of generations)", line=3.5, cex.lab=1.2)

    
    mtext(figurenb[j],2,0,at=6,cex=1.5,line=4,las=1)
    polygon(c(age,rev(age)), c(log2_noinf(m_rate + se_rate), rev(log2_noinf(m_rate - se_rate))), col ="grey", border = NA)
    lines(age,log2_noinf(m_rate),lwd=2)
    print(m_rate)
    points(aranm$age,log2_noinf(rate[,j+1]),pch=20,col=popcol2s[rate[,6]],cex=1.5)
    lines(c(0,50000),c(0,0),lty="dashed")
    lines(c(0,50000),c(-2.9,-2.9),lty="dotted",col="grey")
    axis(1, at = c(0,10000,20000,30000,40000,50000), labels = c("0","10","20","30","40","50"))
    axis(2, las=1,at = c(-3:5), labels = c(0,0.25,0.5,1,2,4,8,16,32))
    
}

dev.off()
################################################################################################################################################
######################################################    Figure S1   ##########################################################################
################################################################################################################################################
library(ggplot2)
library(tidyr)
library(dplyr)
outputfile<-generaloutputfile
inputfile<-generalinputfile

plot.width = 7
plot.height = 5
theme_set(theme_bw(base_size = 24))
line_thickness = 0.8
theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
X = read.csv(paste(inputfile,"count.LTEE.final_masked.csv",sep=""))
X$genome_size_final = X$total_bp - X$deleted_bp + X$inserted_bp
X$genome_size_final = as.numeric(X$genome_size_final)
Z = X %>% group_by(population, time) %>% summarize(avg_genome_size_final = mean(genome_size_final))
Z$avg_genome_size_final<-Z$avg_genome_size_final/1000000
p = ggplot(Z, aes(factor(time), avg_genome_size_final))
p + geom_hline(yintercept=X$total_bp[1]/1000000) + geom_boxplot(fill="white") + scale_x_discrete("Time (thousands of generations)", drop=TRUE,labels=c("0.5","1","1.5","2","5","10","15","20","30","40","50"))+scale_y_continuous("Genome Size (MB)")+theme(text=element_text(size=14))

#ggsave(paste(outputfile,"FigureS1.pdf",sep=""), width=8, height=5)

ggsave(paste(outputfile,"FigureS1.eps",sep=""), width=8, height=5,device="eps")

################################################################################################################################################
######################################################    Figure S2   ##########################################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt
outputfile<-generaloutputfile
inputfile<-generalinputfile

#pdf(paste(outputfile,"FigureS2.pdf",sep=""))
cairo_ps(paste(outputfile,"FigureS2.eps",sep=""))
#nf <- layout(matrix(c(1,3,2,4), 2, 2), height=c(9.9/10,0.1/10),respect=TRUE)
par(mfrow=c(1,1),mar=c(5,5,4,4))

#get data
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)

#exclure IS150 mutator
typemuts<-typemuts[typemuts[,1]!="Aram3_40000_Aram3_40000gen_10988",]

#compute total number of mutations
typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]
plot(NA,xlim=c(0,50),ylim=c(0,450),ylab="Number of synonymous mutations", xlab="Time (thousands of generations)",main="",cex.main=0.8,cex.axis=0.9,las=1,cex.lab=1.2)
means=c()
mtext("a",3,1,adj=0,at=c(-15),cex=1.5)
#set colors
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

popswitch<-c(25000,   100000,  3000,  100000,  38000,  3000,  8000,   100000,  100000, 100000, 100000,6000)

#subsample mutator population
popsmut<-pops[popswitch<100000]
popcolsmut<-popcols[popswitch<100000]
rgbcolsmut<-col2rgb(popcolsmut)
popcol2mut<-rgb(t(rgbcolsmut),alpha=150, maxColorValue=255)


couleur=1
for (pop in popsmut)
{
    temp<-typemuts[typemuts[,2]==pop,]
    points(temp[,3]/1000,temp[,10],col=popcol2mut[couleur],pch=20,cex=1.5)
    means=c()
    count=1
    for (i in c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000))
    {
        means[count]=mean(temp[temp[,3]==i ,10])
        count=count+1;
    }
    
    lines(c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)/1000,means,col=popcol2mut[couleur])
    couleur<-couleur+1
}





dev.off()



################################################################################################################################################
######################################################    Figure S3   ##########################################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt
#figure similar to figure 3 but each pop havin gits own fit

outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1


as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#pdf(paste(outputfile,"FigureS3.pdf",sep=""))
cairo_ps(paste(outputfile,"FigureS3.eps",sep=""))

#nf <- layout(matrix(c(1,), 2, 2), height=c(9.9/10,0.1/10),respect=TRUE)
par (mfrow=c(1,1))
opar<-par(no.readonly=TRUE)                               # stockage des paramètre graphiques
par(oma=c(2,2,0,0))                                       # définition des marges, bas, gauche, haut, droit
frame()
mtext("Time (thousands of generations)",side=1,outer=TRUE,cex=1)
mtext("Number of mutations",side=2,outer=TRUE,cex=1)

#par(mar=c(3,3,1,1))
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),4,3,byrow=TRUE),TRUE)               # partage de l'espace graphique en 3/5 2/5

#par(mfrow=c(4,3),mar=c(3,3,1,1),oma=c(2,2,2,2))

#par(mfrow=c(4,3),mar=c(3,3,1,1),oma=c(2,2,2,2))

typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)

typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]
dimnames(typemuts)[[2]][9]<-"nbmut"
means=c()




popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")

popswitch<-c(25000,100000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
if (arap1_and_aram530k_out)
{
    popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}
timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)
figurelabel<-c("a","b","c","d","e","f","g","h","i","j","k","l")

rgbcols<-col2rgb(popcols)
#rgbcols<-rbind(rgbcols,rep(0.8,length(rgbcols[1,])))
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)


couleur=1
for (pop in pops)
{
    par(mfg=c((couleur-1)%/%3+1,(couleur-1)%%3+1))
    par(mar=c(3,3,2,1))
    dataforfit<-c()

    plot(NA,xlim=c(0,50),ylim=c(0,120),ylab="Number of mutations", xlab="Time (thousands of generations)",main=popslabels[couleur],cex.main=0.8,cex.axis=0.5,las=1,cex.lab=0.9)
    mtext(figurelabel[couleur],3,0.5,adj=0,at=c(-5),cex=1)

    temp<-typemuts[typemuts$age< popswitch[couleur] & typemuts[,2]==pop,]
    points(temp[,3]/1000,temp[,9],col=popcol2s[couleur],pch=20,cex=1.5)
    for (i in timessampling[timessampling<popswitch[couleur]])
    {
        dataforfit<-rbind(dataforfit,temp[temp[,3]==i ,c(2,3,9)])
    }
  

#linearFit <- nls(nbmut ~ a*age, dataforfit,start=c(a=1))
#    sqrtFit<-nls(nbmut ~ b*sqrt(age), dataforfit,start=c(b=1))
#    CombienFit<-nls(nbmut ~ a*age+b*sqrt(age), dataforfit,start=c(a=1,b=1))
   
    lfit<-lm(nbmut ~0+age, dataforfit)#linear
    sfit<-lm(nbmut ~0+sqrt(age), dataforfit)#square root
    lsfit<-lm(nbmut ~0+age+ sqrt(age), dataforfit)#combined

    a<-lfit$coefficients[[1]]
    times<-c(0:500)/10
    linear_fit<-a*times*1000
    lines((times),linear_fit,col=rgb(0.5,0.5,0.5,0.8),lwd=1,lty="dashed")

    b<-sfit$coefficients[[1]]
    sqrt_fit<-b*sqrt(times*1000)
    lines((times),sqrt_fit,col=rgb(0.5,0.5,0.5,0.8),lwd=1)

    a_d<-lsfit$coefficients[[1]]
    b_d<-lsfit$coefficients[[2]]
    combined_fit<-a_d*times*1000+b_d*sqrt(times*1000)
    lines((times),combined_fit,col=rgb(0,0,0,1),lwd=1)
    
    #a<-profile(linearFit)$a[[2]][6]
    #times<-c(0:500)/10
    #linear_fit<-a*times*1000
    #lines(times,linear_fit,col=rgb(0.5,0.5,0.5,0.8),lwd=1,lty="dashed")
    
    #b<-profile(sqrtFit)$b[[2]][6]
    #sqrt_fit<-b*sqrt(times*1000)
    #lines(times,sqrt_fit,col=rgb(0.5,0.5,0.5,0.8),lwd=1)
    
    #a_d<-profile(CombienFit)$a[[2]][6,1]
    #b_d<-profile(CombienFit)$b[[2]][6,2]
    #combined_fit<-a_d*times*1000+b_d*sqrt(times*1000)
    #lines(times,combined_fit,col=rgb(0,0,0,1),lwd=1)





couleur<-couleur+1
}

dev.off()

################################################################################################################################################
######################################################    Figure S4               ##############################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt
outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1


timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)

#get rid of IS150 mutator


#get data and compute total number of mutations
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)
typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]



if (arap1_and_aram530k_out)
{ popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}


#plot for each pop and compute the grandmean
grandmeans<-timessampling*0
grandmeanscounts<-grandmeans
dataforfit<-c()
couleur=1
for (pop in pops)
{
    temp<-typemuts[typemuts[,2]==pop & typemuts$age<popswitch[couleur],]
    means=c()
    count=1
    for (i in timessampling[timessampling<popswitch[couleur]])
    {
        means[count]=mean(temp[temp[,3]==i ,9])
        grandmeans[count]<-grandmeans[count]+means[count]
        grandmeanscounts[count]<-grandmeanscounts[count]+1
        dataforfit<-rbind(dataforfit,c(pop,i,means[count]))
        count=count+1;
    }
    
    
    couleur<-couleur+1
}
grandmeans<-grandmeans/grandmeanscounts

#convert data for the fitting functions
dataforfit<-as.data.frame(dataforfit)
dimnames(dataforfit)[[2]]<-c("Pop","time","nbmut")
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
dataforfit$time<-as.numeric.factor(dataforfit$time)
dataforfit$nbmut<-as.numeric.factor(dataforfit$nbmut)





#####


lfit<-lm(nbmut ~0+time, dataforfit)#linear
sfit<-lm(nbmut ~0+sqrt(time), dataforfit)#square root
lsfit<-lm(nbmut ~0+time+ sqrt(time), dataforfit)#combined

#####


dataforfit<-dataforfit[order(dataforfit$time),]

#compute the AIC for all the numbers
lk<-matrix(rep(0,401*401),ncol=401)

lsfittemp<-lsfit
kc<-1
for (aa in c(0:400)/400*0.002)
{
    kl<-1
    for (bb in c(0:400)/400*0.35)
    {
        lk[kc,kl]<-extractAIC(lm(nbmut ~0+offset(aa*time)+ offset(bb*sqrt(time)), dataforfit))[2]
        kl<-kl+1
    }
    print(kc)
    kc<-kc+1
}
#Go back to likelyhood
lk<-lk/(2)



image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
    if(!missing(breaks)){
        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
    }
    if(missing(breaks) & !missing(zlim)){
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    if(missing(breaks) & missing(zlim)){
        zlim <- range(z, na.rm=TRUE)
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    poly <- vector(mode="list", length(col))
    for(i in seq(poly)){
        poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
    }
    xaxt <- ifelse(horiz, "s", "n")
    yaxt <- ifelse(horiz, "n", "s")
    if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
    if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
    if(missing(xlim)) xlim=XLIM
    if(missing(ylim)) ylim=YLIM
    plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i",ylab="",xlab="Log Likelihood",cex.lab=1.2 ,...)
    for(i in seq(poly)){
        if(horiz){
            polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
        }
        if(!horiz){
            polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
        }
    }
}


lkmem<-lk

breaks <- seq(min(lk), max(lk),length.out=100)


#pdf(paste(outputfile,"FigureS4.pdf",sep=""))
setEPS()
postscript(paste(outputfile,"FigureS4.eps",sep=""))

nf<-layout(matrix(c(1,2,3,0,4,0), nrow=2, ncol=3), widths=c(4,4,1), heights=c(4,1))
nf<-layout(matrix(c(1,2),ncol=1), heights=c(4,1))

par(mar=c(4,5,1,2))
image(x=c(0:400)/400*0.002,y=c(0:400)/400*0.35,lk,breaks=breaks,col=rainbow(length(breaks)-1),xlab="Linear coefficient: a", ylab=" Square root coefficient: b",cex.lab=1.2)

#plot the three lk isoclines max 126.6672 + 2, 6 and 10
contour(x=c(0:400)/400*0.002,y=c(0:400)/400*0.35,lk,add=TRUE,nlevels=3,levels=126.6672+c(2,6,10),drawlabels=FALSE)
points(lsfit$coefficients[[1]],lsfit$coefficients[[2]],pch=20,col="black",cex=2)
points(lfit$coefficients[[1]],0,pch=20,col="black",cex=2)
points(0,sfit$coefficients[[1]],pch=20,col="black",cex=2)

par(mar=c(5,5,1,2))
image.scale(lk, col=rainbow(length(breaks)-1), breaks=-breaks, horiz=TRUE)
dev.off()



################################################################################################################################################
######################################################    Figure S5               ##############################################################
################################################################################################################################################
outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1


typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)
typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]

#use MAE to infer basal rates
typemutsMAE<-read.table(paste(inputfile,"MutationTypesThroughTimeMAE.txt",sep=""),header=TRUE)
refrates<-typemuts[1,]
refrates[5:31]<-colSums(typemutsMAE[5:31])/(13750*15)


#set colors and exclure mutators
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popswitch<-c(25000,100000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)


if (arap1_and_aram530k_out)
{
    popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}
timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)
figurelabel<-c("a","b","c","d","e","f","g","h","i","j","k","l")
rgbcols<-col2rgb(popcols)
#rgbcols<-rbind(rgbcols,rep(0.8,length(rgbcols[1,])))
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#get non mutator clones
aranm<-c()
nbpop<-1
for (p in pops)
{
    age<-sort(unique(typemuts[typemuts[,2]==p & typemuts[,3]<popswitch[nbpop],3]))
    
    for (t in age)
    {
        temp<-typemuts[typemuts[,2]==p & typemuts[,3]==t ,]
        labeltemp<-temp[1,1:4]
        labeltemp[2]<-nbpop
        aranm<-rbind(aranm,c(labeltemp,apply(temp[,5:30],2,mean)))
    }
    nbpop<-nbpop+1
}

colnames<-dimnames(aranm)[[2]]
aranm<-as.data.frame(aranm)
new.list<-lapply(1:ncol(aranm), function(x) as.vector(unlist(aranm[,x])))
aranm<-as.data.frame(do.call(cbind, new.list), stringsAsFactors=FALSE)
aranm<-as.data.frame(aranm)
dimnames(aranm)[[2]]<-colnames

l<-length(aranm[1,])
for (i in c(2,3,5:l))
{
    aranm[,i]<-as.numeric(aranm[,i])
}


#pdf(paste(outputfile,"FigureS5.pdf",sep=""))
cairo_ps(paste(outputfile,"FigureS5.eps",sep=""))

#par (fg="white",col="white",bg="black",col.axis="white", col.main="white",col.lab="white")
par(mar=c(5,5,4,4))
age<-sort(unique(aranm$age))
ylimmax<-7
ylimmin<-(-2)

par(mfrow=c(1,1))




m_rate<-age
se_rate<-age
i<-1
for (a in age)
    {
    m_rate[i]<-mean(aranm[aranm$age==a,]$syn)
    i<-i+1
    }

arabb<-aranm
arabb$muts<-aranm$age*(aranm$pop-6.5)*0.03
arabb[arabb$age>4000,]$muts<-(arabb[arabb$age>4000,]$pop-6.5)*150
plot(aranm$age+arabb$muts,aranm$syn,pch=20, ylim=c(0,6),main="",ylab="Number of Synonymous Mutations",xlab="Time (thousands of generations)",axes=FALSE,col=popcol2s[aranm$pop],cex=1.2,cex.lab=1.2)

#plot(aranm$age*(1+0.003*((aranm$pop)-6.5)),aranm$syn,pch=20, ylim=c(0,8),main="",ylab="Number of Synonymous Mutations",xlab="Time (thousands of generations)",axes=FALSE,col=popcol2s[aranm$pop],cex=1.2)
lines(c(0,50000),c(0,0),lty="dashed")
points(age,m_rate,col=1,cex=1.5,pch=2)
#lines(c(0,age),c(0,m_rate),lwd=2,col="black")

lines(c(0,50000),c(0,m_rate[11]),col=rgb(0,0,0,0.5),lwd=3)

axis(1, at = c(0,10000,20000,30000,40000,50000), labels = c("0","10","20","30","40","50"))
axis(2, las=1,at = c(0:5), labels = c(0:5))

boby<-rep(0,11)
k<-1
for (ttime in timessampling)
{mm<-mean(aranm[aranm$age==ttime,]$syn)
vv<-var(aranm[aranm$age==ttime,]$syn)
boby[k]<-vv/mm
k<-k+1
}

#plot(timessampling,boby)
#boby
dev.off()



################################################################################################################################################
######################################################    Figure S6             ##############################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt and GenomeComposition.txt
#as in previous figures we compute rates. Then rather tahn plotting the actual values we plot (n(t)-n(t-1))/(t- (t-1))/ NSrate.
outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1

genomecomp<-read.table(paste(inputfile,"GenomeComposition.txt",sep=""),header=TRUE)
synonymoussitespertypeofmut<-genomecomp[1,2:7]
nonsynonymoussitespertypeofmut<-genomecomp[2,2:7]
intergenicsitespertypeofmut<-genomecomp[3,2:7]
genomicGC<-genomecomp[4,2:7]


#Get data and Compute the total number of mutations
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)
typemuts<-cbind(typemuts,typemuts[,8])
typemuts[,9]<-typemuts[,5]+typemuts[,6]+typemuts[,7]+typemuts[,8]

#use MAE to infer basal rates
typemutsMAE<-read.table(paste(inputfile,"MutationTypesThroughTimeMAE.txt",sep=""),header=TRUE)
refrates<-typemuts[1,]
refrates[5:31]<-colSums(typemutsMAE[5:31])/(13750*15)


#define colors and transparency for the symbols
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#Subsample non mutator clones
means=c()
popcolsnm<-popcols[c(2,4,8,9,10,11)]
rgbcolsnm<-col2rgb(popcolsnm)
popcol2nms<-rgb(t(rgbcolsnm),alpha=150, maxColorValue=255)
popsneutral<-pops[c(2,4,8,9,10,11)]
if (arap1_and_aram530k_out)
{
popsneutral<-pops[c(4,8,9,10,11)]
}
#compute neutral rates used in subsequent figures
NEUTRALRATE<-0
synsrates<-rep(0,6)
couleur=1
for (pop in popsneutral)
{
    temp<-typemuts[typemuts[,2]==pop,]
    
    count=1
    for (i in c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000))
    {
        #means[count]=mean(temp[temp[,3]==i ,9])
        count=count+1;
        if (i==50000)
        {
            print(c(temp[temp[,3]==i,]$syn[1],sum(temp[temp[,3]==i,]$syn)))
            NEUTRALRATE<-NEUTRALRATE+sum(temp[temp[,3]==i,]$syn)
            synsrates<-synsrates+colSums(temp[temp[,3]==i,13:18])
        }
    }
    
    couleur<-couleur+1
}

nbpopsrate<-length(popsneutral)
NEUTRALRATE<-NEUTRALRATE/(50000*2*nbpopsrate)
NSRATE<-sum(synsrates/synonymoussitespertypeofmut*nonsynonymoussitespertypeofmut/(50000*2*nbpopsrate))
INTERRATE<-sum(synsrates/synonymoussitespertypeofmut*intergenicsitespertypeofmut/(50000*2*nbpopsrate))
MUTRATE<-sum(synsrates/synonymoussitespertypeofmut*genomicGC/(50000*2*nbpopsrate))



#set colors and exclure mutators
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popswitch<-c(25000,100000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)


if (arap1_and_aram530k_out)
{
    popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}
timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)
figurelabel<-c("a","b","c","d","e","f","g","h","i","j","k","l")
rgbcols<-col2rgb(popcols)
#rgbcols<-rbind(rgbcols,rep(0.8,length(rgbcols[1,])))
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#get non mutator clones
aranm<-c()
nbpop<-1
for (p in pops)
{
    age<-sort(unique(typemuts[typemuts[,2]==p & typemuts[,3]<popswitch[nbpop],3]))
    
    for (t in age)
    {
        temp<-typemuts[typemuts[,2]==p & typemuts[,3]==t ,]
        labeltemp<-temp[1,1:4]
        labeltemp[2]<-nbpop
        aranm<-rbind(aranm,c(labeltemp,apply(temp[,5:30],2,mean)))
    }
    nbpop<-nbpop+1
}

colnames<-dimnames(aranm)[[2]]
aranm<-as.data.frame(aranm)
new.list<-lapply(1:ncol(aranm), function(x) as.vector(unlist(aranm[,x])))
aranm<-as.data.frame(do.call(cbind, new.list), stringsAsFactors=FALSE)
aranm<-as.data.frame(aranm)
dimnames(aranm)[[2]]<-colnames

l<-length(aranm[1,])
for (i in c(2,3,5:l))
{
    aranm[,i]<-as.numeric(aranm[,i])
}

aranm<-rbind(aranm[1:12,],aranm)
for (p in c(1:12))
aranm[p,]<-c("REL606",p,0,"REL606",rep(0,26))
aranm<-as.data.frame(aranm)


#pdf(paste(outputfile,"FigureS6.pdf",sep=""))
cairo_ps(paste(outputfile,"FigureS6.eps",sep=""))
#par (fg="white",col="white",bg="black",col.axis="white", col.main="white",col.lab="white")
par (mfrow=c(1,1),mar=c(5,6,4,3))
plot(NA,ylim=c(-3,5),xlim=c(0,50000),main="",ylab="Observed nonsynonymous mutations\n relative to neutral expectation",xlab="Time (thousands of generations)",axes=FALSE,col=popcol2s[aranm$pop],cex.lab=1.2,cex=1.2)
axis(1, at = c(0,10000,20000,30000,40000,50000), labels = c("0","10","20","30","40","50"))
axis(2, las=1,at = c(-3:5), labels = c(0,0.25,0.5,1,2,4,8,16,32))
couleur<-1

for (p in c(2:13))
{
    age<-sort(as.numeric(unique(aranm[aranm[,2]==p-1 & as.numeric(aranm[,3])<popswitch[p-1],3])))
    rate_relatif<-age[2:length(age)]

    for (t in c(2:length(age)))
    {
        rate_relatif[t-1]<-((as.numeric(aranm[aranm[,2]==p-1 & aranm[,3]==age[t],]$nonsyn)-as.numeric(aranm[aranm[,2]==p-1 & aranm[,3]==age[t-1],]$nonsyn))/(as.numeric(age[t])-as.numeric(age[t-1])))/NSRATE
    }
    age<-age[2:length(age)]
}

age<-sort(as.numeric(unique(aranm[,3])))
rate_relatif<-age[2:length(age)]*0
se_rate_relatif<-age[2:length(age)]*0
nbcounts<-age[2:length(age)]*0
for (t in c(2:length(age)))
{
    for (p in (1:12))
    {
     if (age[t]<popswitch[p])
     {
         tmp<-((as.numeric(aranm[aranm[,2]==p & aranm[,3]==age[t],]$nonsyn)-as.numeric(aranm[aranm[,2]==p & aranm[,3]==age[t-1],]$nonsyn))/(as.numeric(age[t])-as.numeric(age[t-1])))/NSRATE
         rate_relatif[t-1]<-rate_relatif[t-1]+tmp
         se_rate_relatif[t-1]<-se_rate_relatif[t-1]+tmp^2
         nbcounts[t-1]<-nbcounts[t-1]+1
     }
    }
}
    
for (t in c(2:length(age)))
{
rate_relatif[t-1]<-rate_relatif[t-1]/nbcounts[t-1]
se_rate_relatif[t-1]<-sqrt(se_rate_relatif[t-1]/nbcounts[t-1]-rate_relatif[t-1]*rate_relatif[t-1])/sqrt(nbcounts[t-1])
}

log2_noinf<-function(x)
{
    y<-x
    y[y<=0]<-0.125
    return(log2(y))
}


age<-age[2:length(age)]
polygon(c(age,rev(age)), c(log2_noinf(rate_relatif + se_rate_relatif), rev(log2_noinf(rate_relatif - se_rate_relatif))), col ="grey", border = NA)
lines(age,log2_noinf(rate_relatif),lwd=2)
print(cbind(age,rate_relatif))
for (p in c(2:13))
{
    age<-sort(as.numeric(unique(aranm[aranm[,2]==p-1 & as.numeric(aranm[,3])<popswitch[p-1],3])))
    rate_relatif<-age[2:length(age)]
    
    for (t in c(2:length(age)))
    {
        rate_relatif[t-1]<-((as.numeric(aranm[aranm[,2]==p-1 & aranm[,3]==age[t],]$nonsyn)-as.numeric(aranm[aranm[,2]==p-1 & aranm[,3]==age[t-1],]$nonsyn))/(as.numeric(age[t])-as.numeric(age[t-1])))/NSRATE
    }
    age<-age[2:length(age)]
    print (c(p,rate_relatif))
    points(age*(1+0.003*(p-7.5)),log2_noinf(rate_relatif),pch=20,col=popcol2s[p-1],cex=1.5)
}



lines(c(0,50000),c(0,0),lty="dashed")
lines(c(0,50000),c(-2.9,-2.9),lty="dotted",col="grey")


dev.off()

################################################################################################################################################
######################################################    Figure S7   ##############################################################
################################################################################################################################################
library(ggplot2)
library(tidyr)
library(dplyr)

## load data
## X is the full table

X = read.csv(paste(inputfile,"spectrum_counts.csv",sep=""))

X$population = X$sample

#don't count inversions
X$total = X$total - X$inversion

#large substitutions are all deletions!
X$large_deletion = X$large_deletion + X$large_substitution
X$large_substitution = 0

X$fr_base_substitution = X$base_substitution / X$total
X$fr_small_indel = X$small_indel / X$total
X$fr_large_deletion = X$large_deletion / X$total
X$fr_large_amplification = X$large_amplification / X$total
X$fr_mobile_element_insertion = X$mobile_element_insertion / X$total

## Recategorize into four types of SNPs
## synonymous = synonymous
## nonsynonymous = nonsynymous + nonsense
## intergenic = intergenic
## other = noncoding (RNA) + pseudogene

X$fr_synonymous = X$base_substitution.synonymous / X$total

X$base_substitution.nonsynonymous = X$base_substitution.nonsynonymous + X$base_substitution.nonsense
X$base_substitution.nonsense = 0
X$fr_nonsynonymous = X$base_substitution.nonsynonymous / X$total

X$fr_intergenic = X$base_substitution.intergenic  / X$total

X$base_substitution.other = X$base_substitution.pseudogene + X$base_substitution.noncoding
X$base_substitution.pseudogene = 0
X$base_substitution.noncoding = 0
X$fr_other = X$base_substitution.other  / X$total

X$check_total = X$fr_base_substitution + X$fr_small_indel + X$fr_large_deletion + X$fr_large_amplification + X$fr_mobile_element_insertion

X$check_total = X$fr_synonymous + X$fr_nonsynonymous + X$fr_intergenic + X$fr_other + X$fr_small_indel + X$fr_large_deletion + X$fr_large_amplification + X$fr_mobile_element_insertion

X$time = 0;
X$cumulative = FALSE;
X$population = c("")
for (i in 1:nrow(X)) {
    
    if (grepl(".to.", X$file[i])) {
        
        X$cumulative[i] = TRUE;
        
    } else {
        X$cumulative[i] = FALSE;
    }
    
    file_name = as.character(X$file[i])
    X$time[i] = sub("^.+?\\.(\\d+)gen\\..+", "\\1", file_name , perl=T)
    X$population[i] = sub("^(.+?)\\..+", "\\1", file_name , perl=T)
    
}


X$population = as.factor(X$population)
X$time = as.numeric(X$time)

## save the full table
full_table = X


order_of_pops = c("Ara-5", "Ara-6", "Ara+2", "Ara+4", "Ara+5", "Ara+1", "Ara-2", "Ara-4", "Ara+3", "Ara-3", "Ara+6", "Ara-1")
X = full_table %>% filter(population != "Ara-5-no-alien")
##Reorder based on mutator type
X$population = factor(X$population, levels = order_of_pops)


#plot defaults

plot.width = 5
plot.height = 4
theme_set(theme_bw(base_size = 24))
line_thickness = 0.8
theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))

## Nonmutator spectrum over time

unrolled_counts = full_table %>%
select(population, time, cumulative, base_substitution.synonymous, base_substitution.nonsynonymous, base_substitution.intergenic, base_substitution.other, small_indel, large_deletion, large_amplification, mobile_element_insertion) %>%
gather(mutation.type, count, base_substitution.synonymous:mobile_element_insertion)

instantaneous_counts = subset(unrolled_counts, cumulative==FALSE)

non_mutators = instantaneous_counts %>%
filter( ((population=="Ara+1") & (time<= 2000)) |
((population=="Ara+2") & (time<=50000)) |
((population=="Ara+3") & (time<= 2000)) |
((population=="Ara+4") & (time<=50000)) |
((population=="Ara+5") & (time<=50000)) |
((population=="Ara+6") & (time<= 5000)) |
((population=="Ara-1") & (time<=20000)) |
((population=="Ara-2") & (time<= 2000)) |
((population=="Ara-3") & (time<=30000)) |
((population=="Ara-4") & (time<= 5000)) |
((population=="Ara-5-no-alien") & (time<=50000)) |
((population=="Ara-6") & (time<=50000))
)

non_mutators = non_mutators %>%
group_by(mutation.type, time) %>% summarize(count = sum(count))

non_mutator_totals = non_mutators %>% group_by(time) %>% summarize(total = sum(count))

non_mutators = non_mutators %>% left_join(non_mutator_totals, by="time")

non_mutators$fraction = non_mutators$count / non_mutators$total
non_mutators$time<-non_mutators$time/1000
non_mutators$time = factor(non_mutators$time)
non_mutators$mutation.type = factor(non_mutators$mutation.type, levels = c("base_substitution.synonymous","base_substitution.nonsynonymous","base_substitution.intergenic","base_substitution.other","small_indel","large_deletion","large_amplification","mobile_element_insertion"),labels = c("Synonymous","Nonsynonymous","Intergenic           ","Other            ","Indel (≤50bp)","Deletion (>50bp)","Duplication (>50bp)","IS insertion"))

p = ggplot(arrange(non_mutators,mutation.type), aes(x=time,y=fraction, fill=mutation.type,order=mutation.type))


#k<-p + geom_bar(stat="identity")+guides(col = guide_legend(nrow = 2, byrow = TRUE))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal",legend.title=element_text(size=12,vjust=1,debug=TRUE),legend.text = element_text( size=12))+scale_x_discrete("Time (thousand of generations)", labels=c(0.5,1,1.5,2,5,10,15,20,30,40,50),breaks=c(1:11))+scale_y_continuous("percentage of mutation", breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))
#k+guides(col = guide_legend(nrow = 4, byrow = TRUE))
upperleg<-c("N=","47","69","104","126","144","209","268","304","366","373","448")

cbPalette <- c("#999999", "#D55E00","#0072B2","#E69F00",  "#F0E442",   "#009E73","#56B4E9", "#CC79A7")

#p + geom_bar(stat="identity")+scale_fill_continuous(guide = "legend")
#p+ geom_bar(stat="identity")+scale_fill_discrete(guide = guide_legend(title = "Point mutations"))+scale_fill_manual(values=cbPalette)

#p+geom_bar(stat="identity")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.position="bottom",legend.title=element_text(size=12,vjust=1,debug=TRUE),legend.text = element_text( size=12))+scale_x_discrete("Time (thousand of generations)", labels=c(0.5,1,1.5,2,5,10,15,20,30,40,50),breaks=c(1:11))+scale_y_continuous("percentage of mutation", breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+scale_fill_manual(values=cbPalette)
cairo_ps(paste(outputfile,"FigureS7.eps",sep=""),height=6, width=8,)

p+geom_bar(stat="identity")+theme(plot.title=element_text(size=10,hjust=0),plot.margin = unit(c(4,1,1,1),"line"),axis.text=element_text(size=12),panel.grid.major = element_blank(),panel.grid.minor =element_blank(), axis.title=element_text(size=14), legend.position="bottom", legend.title=element_text(size=10,vjust=1,debug=TRUE),legend.text = element_text( size=10))+scale_x_discrete("Time (thousand of generations)")+scale_y_continuous("Percentage of mutation", breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+scale_fill_manual(values=cbPalette,guide = guide_legend(title = "Base substitutions"))+ggtitle("N=   47           69          104          126          144          209         268         304         366         373        448")
dev.off()
#ggsave(paste(outputfile,"FigureS7.eps",sep=""), height=6, width=8,device="eps")

################################################################################################################################################
######################################################    Figure S8   ##############################################################
################################################################################################################################################


library(ggplot2)
outputfile<-generaloutputfile
inputfile<-generalinputfile


X = read.csv(paste(inputfile,"MAE_fitness.csv",sep=""))
limits<-aes(ymax = CL.95., ymin= CL_95.)
p <- ggplot(X, aes( y= relative_fitness, x= rank))
p + geom_point( shape=15,size=4) + geom_errorbar(limits, width=0.2)+ geom_hline(yintercept=1,linetype = 2,col="red") + scale_x_continuous("Mutation accumulation lineages", labels=c("MAE\nAnc",1:15),breaks=c(1:16))+coord_cartesian(ylim=c(0, 1.25),xlim=c(1,16))+scale_y_continuous("Fitness relative to MAE ancestor",breaks=seq(0,1.25,0.25))+theme(axis.text=element_text(size=12),panel.grid.major = element_blank(),panel.grid.minor =element_blank(),axis.title=element_text(size=14))+geom_text(label=c("","**","","**","","**","","**","","*","","*","**","**","**","**"),x=c(1:16),y=rep(1.22,16),size=9)

#ggsave(paste(outputfile,"FigureS7.pdf",sep=""), width=8, height=5)
ggsave(paste(outputfile,"FigureS8.eps",sep=""), width=8, height=5,device="eps")



################################################################################################################################################
######################################################    Figure S9   ##############################################################
################################################################################################################################################
#requires MutationTypesThroughTime.txt and GenomeComposition.txt
outputfile<-generaloutputfile
inputfile<-generalinputfile
arap1_and_aram530k_out<-1

genomecomp<-read.table(paste(inputfile,"GenomeComposition.txt",sep=""),header=TRUE)
synonymoussitespertypeofmut<-genomecomp[1,2:7]
nonsynonymoussitespertypeofmut<-genomecomp[2,2:7]
intergenicsitespertypeofmut<-genomecomp[3,2:7]
genomicGC<-genomecomp[4,2:7]


#Get data
typemuts<-read.table(paste(inputfile,"MutationTypesThroughTime.txt",sep=""),header=TRUE)

#use MAE to infer basal rates
typemutsMAE<-read.table(paste(inputfile,"MutationTypesThroughTimeMAE.txt",sep=""),header=TRUE)
refrates<-typemuts[1,]
refrates[5:31]<-colSums(typemutsMAE[5:31])/(13750*15)


#define colors and transparency for the symbols
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
rgbcols<-col2rgb(popcols)
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)

#Subsample non mutator clones
means=c()
popcolsnm<-popcols[c(2,4,8,9,10,11)]
rgbcolsnm<-col2rgb(popcolsnm)
popcol2nms<-rgb(t(rgbcolsnm),alpha=150, maxColorValue=255)
pops<-pops[c(2,4,8,9,10,11)]

#compute neutral rates used in subsequent figures
NEUTRALRATE<-0
synsrates<-rep(0,6)
couleur=1
for (pop in pops)
{
    temp<-typemuts[typemuts[,2]==pop,]
    
    count=1
    for (i in c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000))
    {
        #means[count]=mean(temp[temp[,3]==i ,9])
        count=count+1;
        if (i==50000)
        {
            print(c(temp[temp[,3]==i,]$syn[1],sum(temp[temp[,3]==i,]$syn)))
            NEUTRALRATE<-NEUTRALRATE+sum(temp[temp[,3]==i,]$syn)
            synsrates<-synsrates+colSums(temp[temp[,3]==i,13:18])
        }
    }
    
    couleur<-couleur+1
}

nbpopsrate<-length(pops)
NEUTRALRATE<-NEUTRALRATE/(50000*2*nbpopsrate)
NSRATE<-sum(synsrates/synonymoussitespertypeofmut*nonsynonymoussitespertypeofmut/(50000*2*nbpopsrate))
INTERRATE<-sum(synsrates/synonymoussitespertypeofmut*intergenicsitespertypeofmut/(50000*2*nbpopsrate))
MUTRATE<-sum(synsrates/synonymoussitespertypeofmut*genomicGC/(50000*2*nbpopsrate))

log2_noinf<-function(x)
{
    y<-x
    y[y<=0]<-0.25
    return(log2(y))
}

popcols<-c("dark blue","light blue","saddlebrown","lemonchiffon3","red","indianred1","green4","darkolivegreen3","purple", "plum", "darkgoldenrod3","yellow1")
pops<-c("Aram1","Arap1","Aram2","Arap2","Aram3","Arap3","Aram4","Arap4","Aram5","Arap5","Aram6","Arap6")
popslabels<-c("Ara-1","Ara+1","Ara-2","Ara+2","Ara-3","Ara+3","Ara-4","Ara+4","Ara-5","Ara+5","Ara-6","Ara+6")
popswitch<-c(25000,100000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)

arap1_and_aram530k_out<-1
if (arap1_and_aram530k_out)
{
    popswitch<-c(25000,3000,3000,100000,38000,3000,8000,100000,100000,100000,100000,6000)
    typemuts<-typemuts[typemuts[,1]!="Aram5_30000_Aram5_30000gen_10404",]
}
timessampling<-c(500,1000,1500,2000,5000,10000,15000,20000,30000,40000,50000)
figurelabel<-c("a","b","c","d","e","f","g","h","i","j","k","l")
rgbcols<-col2rgb(popcols)
#rgbcols<-rbind(rgbcols,rep(0.8,length(rgbcols[1,])))
popcol2s<-rgb(t(rgbcols),alpha=180, maxColorValue=255)


aranm<-c()
nbpop<-1
for (p in pops)
{
    age<-sort(unique(typemuts[typemuts[,2]==p & typemuts[,3]<popswitch[nbpop],3]))
    
    for (t in age)
    {
        temp<-typemuts[typemuts[,2]==p & typemuts[,3]==t ,]
        labeltemp<-temp[1,1:4]
        labeltemp[2]<-nbpop
        aranm<-rbind(aranm,c(labeltemp,apply(temp[,5:30],2,mean)))
    }
    nbpop<-nbpop+1
}

colnames<-dimnames(aranm)[[2]]
aranm<-as.data.frame(aranm)
new.list<-lapply(1:ncol(aranm), function(x) as.vector(unlist(aranm[,x])))
aranm<-as.data.frame(do.call(cbind, new.list), stringsAsFactors=FALSE)
aranm<-as.data.frame(aranm)
dimnames(aranm)[[2]]<-colnames

l<-length(aranm[1,])
for (i in c(2,3,5:l))
{
aranm[,i]<-as.numeric(aranm[,i])
}




log2_noinf<-function(x)
{
    y<-x
    y[y<=0]<-0.25
    return(log2(y))
}



#pdf(paste(outputfile,"FigureS9.pdf",sep=""))
#setEPS()
#postscript(paste(outputfile,"FigureS9.eps",sep=""))
cairo_ps(paste(outputfile,"FigureS9.eps",sep=""))
#par (fg="white",col="white",bg="black",col.axis="white", col.main="white",col.lab="white")
par (mfrow=c(3,2),mar=c(4,5,2,2),oma=c(2,6,2,2))
age<-sort(unique(aranm$age))
ylimmax<-7
ylimmin<-(-2)

#par(mfrow=c(1,3),oma=c(13,2,13,2))
title=c("Observed nonsynonymous mutations\n relative to neutral expectation","Observed intergenic point mutations\n relative to neutral expectation","Observed IS150 mutations\n relative to neutral expectation","Observed other IS mutations\n relative to neutral expectation","Observed small indel mutations\nrelative to neutral expectation","Observed large indel mutations\n relative to neutral expectation")
coul<-c("purple","red","blue","green")
figurenb<-c("a","b","b","d","e","f")



rate<-cbind(aranm$age,(((aranm$ns1+aranm$ns2+aranm$ns3+aranm$ns4+aranm$ns5+aranm$ns6)/(aranm$age)))/(refrates$nonsyn/refrates$syn*NEUTRALRATE),(((aranm$i1+aranm$i2+aranm$i3+aranm$i4+aranm$i5+aranm$i6)/(aranm$age)))/(refrates$inter/refrates$syn*NEUTRALRATE),(((aranm$is150)/aranm$age)/(((refrates$is150)/refrates$syn)*NEUTRALRATE)),(((aranm$is-aranm$is150)/aranm$age)/(((refrates$is-refrates$is150)/refrates$syn)*NEUTRALRATE)),((aranm$indels/aranm$age)/((refrates$indels/refrates$syn)*NEUTRALRATE)),((aranm$ldels/aranm$age)/((refrates$ldels/refrates$syn)*NEUTRALRATE)),aranm$pop)

ylimmin<-(-2)
ylimmax<-(7)

for (j in c(1:6))
{
    
    m_rate<-age
    se_rate<-age
    i<-1
    for (a in age)
    {
        m_rate[i]<-mean(rate[rate[,1]==a,j+1])
        l<-length(rate[rate[,1]==a,j+1])
        se_rate[i]<-sd(rate[rate[,1]==a,j+1])/sqrt(l)
        i<-i+1
    }
    plot(rate[,1],log2_noinf(rate[,j+1]),pch=20,ylim=c(ylimmin,ylimmax),main=,ylab=title[j],xlab="Time (thousands of generations)",axes=FALSE,col=0,cex.lab=0.9)
    mtext(figurenb[j],3,0,at=c(-22000),cex=1)
    polygon(c(age,rev(age)), c(log2_noinf(m_rate + se_rate), rev(log2_noinf(m_rate - se_rate))), col ="grey", border = NA)
    print(m_rate)
    lines(age,log2_noinf(m_rate),lwd=2)
    points(aranm$age,log2_noinf(rate[,j+1]),pch=20,col=popcol2s[rate[,8]],cex=1.5)
    lines(c(0,50000),c(0,0),lty="dashed")
    lines(c(0,50000),c(-1.9,-1.9),lty="dotted",col="grey")
    axis(1, at = c(0,10000,20000,30000,40000,50000), labels = c("0","10","20","30","40","50"))
    axis(2, las=1,at = c(-2:7), labels = c(0,0.5,1,2,4,8,16,32,64,128))
    
}

dev.off()























