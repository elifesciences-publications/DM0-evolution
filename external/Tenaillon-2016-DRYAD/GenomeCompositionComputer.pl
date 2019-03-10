
#!/usr/bin/perl
use strict;
use warnings;
#MASK and gff3
#this program  computes the type of mutations accessible to be able to compute genome wide Ka or Ks #it takes into account the start codon degenerate code.
#it takes into account the overlapping reading frames
#non coding genes are not taken into account
#masked regions are omitted

#it can be run using the command line perl GenomeCompositionComputer.pl REL606.gff3 REL606.L20.G15.P0.M35.mask.gd
#the first argumlent is the gff3 file of the reference genome
#the second corresponds to the masked regions (repeated regions that are not covered by the sequencing technology)

#it produces the GenomeComposition.txt file which is used for figures.




my $numArgs = $#ARGV + 1;
print "thanks, you gave me $numArgs command-line arguments:\n";


my $gff3file=$ARGV[0];
my $maskfile=$ARGV[1];



my %CodeGenetic=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

our(%CodeGeneticStart)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'M','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'M','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'M','ATC'=>'M','ATT'=>'M','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'M','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

my  (%AAlist)=("F"=>0,"L"=>0,"S"=>0,"Y"=>0, "STOP"=>0,"C"=>0,"W"=>0,"P"=>0,"H"=>0,"Q"=>0,"R"=>0,"I"=>0,"M"=>0,"T"=>0,"N"=>0,"K"=>0,"V"=>0,"A"=>0,"D"=>0,"E"=>0,"G"=>0);

my %MutTypeAction=("1_A"=>"C","1_T"=>"G","2_A"=>"G","2_T"=>"C","3_A"=>"T","3_T"=>"A","4_C"=>"G","4_G"=>"C","5_C"=>"T","5_G"=>"A","6_G"=>"T","6_C"=>"A");
my %Codonlist;
my %CodonlistStart;

my %CodontoAA_mutType;

for my $key (keys %CodeGenetic)
{
    $Codonlist{$key}=0;
    
}
for my $key (keys %CodeGenetic)
{
    $CodonlistStart{$key}=0;
    
}



########################################################################################################################
############################Identify masked regions (regions not covered properly by sequencing technology)#############
########################################################################################################################
my %MaskedBases=();
my @b=();
open (MASKLIST, "<$maskfile");#REL606.L20.G15.P0.M35.mask.gd");
my $line = <MASKLIST>;
while( $line = <MASKLIST>)
{
    chomp $line;
    @b=split(/\t/,$line);
    for (my $i=$b[4]; $i<$b[4]+$b[5]-1;$i++)
    {$MaskedBases{$i}=1;}
}


close MASKLIST;

########################################################################################################################
############################Get genome sequence and create simplified gene position file ###############################
########################################################################################################################

open (GENELIST, "<$gff3file");
open (GENELISTOUT, ">ECBGenesCDS.txt");
$line = <GENELIST>;
$line = <GENELIST>;
$line = <GENELIST>;
while( $line = <GENELIST>)
{
    chomp $line;	# drop \n at end of line
    @b=split(/\t/,$line);
    
    if ($b[0]=~/\#\#/)
    {
        last;
    }
    if ($b[2]eq"CDS" || $b[2]eq"fCDS" || $b[2]eq"tRNA" || $b[2]eq"rRNA")
    {
        my $name;
        $b[8]=~/Name=(.*);Note/;
        $name=$1;
        my $strand="C";
        if ($b[6]eq"+")
        {
            $strand="D"
        }
        print GENELISTOUT "$b[3]\t$b[4]\t$strand\t$name\t$b[2]\n"
        
    }
    
}

my $refsequence="_";
$line = <GENELIST>;
while( $line = <GENELIST>)
{
    chomp $line;
    $refsequence.=$line;
    
}

close GENELISTOUT;

close GENELIST;

my $genomeSize=length $refsequence;
$genomeSize--;

print "$genomeSize\n";
############################################################################################################
#############Attribute to each position a codon, and potentially multiple overlapping frames#################
############################################################################################################
my @GeneName=(1..5000);
my @MutType=(1,2,3,4,5,6);

our %PosToGene=();
our %PosToCodon=();
our %PosToStrand=();
our %PosToPosCodon=();
our %PosToCodonPos=();
our %PosToNbCodon=();
our %BasedTreated=();
open (GENELIST, "<ECBGenesCDS.txt");


my $nb=0;
my $previouspos=1;
my $previousgenenb=$nb-1;
my $overlaptest=0;
my $overlap2=0;
while( $line = <GENELIST>)
{
    chomp $line;	# drop \n at end of line
    @b=split(/\t/,$line);
    
    
    
    for (my $j=$b[0];$j<=$b[1];$j++)
    {
        
        
        
        $PosToGene{$j}=$nb;
        if ($b[4] eq "CDS")
        {
            if ($b[2]eq"D" && $b[4] eq "CDS" )
            {
                my $codon=substr($refsequence, $b[0]+int(($j-$b[0])/3)*3, 3);
                
                if (exists($PosToCodon{$j})==0)
                {
                    $PosToCodon{$j}=$codon;
                    $PosToPosCodon{$j}=($j-$b[0])%3;
                    $PosToCodonPos{$j}=int(($j-$b[0])/3);
                    $PosToStrand{$j}=1;
                    $PosToNbCodon{$j}=1;
                }
                else
                {
                    if (exists($MaskedBases{$j})==0)
                    {
                        $overlap2++;
                    }
                    $PosToCodon{$j}="$PosToCodon{$j},$codon";
                    my $poscod=($j-$b[0])%3;
                    $PosToPosCodon{$j}="$PosToPosCodon{$j},$poscod";
                    $poscod=int(($j-$b[0])/3);
                    $PosToCodonPos{$j}="$PosToCodonPos{$j},$poscod";
                    $PosToStrand{$j}="$PosToStrand{$j},1";
                    $PosToNbCodon{$j}=$PosToNbCodon{$j}+1;
                }
                # print("$j ",substr($refsequence, $j, 1),"   $b[0] $PosToCodon{$j}\t$CodeGenetic{$PosToCodon{$j}}\t\t$PosToPosCodon{$j}\n");
                
            }
            
            if ($b[2]eq"C" && $b[4] eq "CDS" )
            {
                
                my $codon=substr($refsequence, $b[1]-int(($b[1]-$j)/3)*3-2, 3);
                $codon=reverse($codon);
                $codon=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
                if (exists($PosToCodon{$j})==0)
                {
                    $PosToCodon{$j}=$codon;
                    $PosToPosCodon{$j}=($b[1]-$j)%3;
                    $PosToCodonPos{$j}=int(($b[1]-$j)/3);# euclidien division
                    $PosToStrand{$j}=-1;
                    $PosToNbCodon{$j}=1;
                }
                else
                {
                    if (exists($MaskedBases{$j})==0)
                    {
                        $overlap2++;
                    }
                    $PosToCodon{$j}="$PosToCodon{$j},$codon";
                    my $poscod=($b[1]-$j)%3;
                    $PosToPosCodon{$j}="$PosToPosCodon{$j},$poscod";
                    $poscod=int(($b[1]-$j)/3);
                    $PosToCodonPos{$j}="$PosToCodonPos{$j},$poscod";
                    $PosToStrand{$j}="$PosToStrand{$j},-1";
                    $PosToNbCodon{$j}=$PosToNbCodon{$j}+1;
                }
                # my $base=substr($refsequence, $j, 1);
                # $base=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
                #print("$j $b[1] $base $codon\t$CodeGenetic{$codon}\t$PosToPosCodon{$j}\n");
                
            }
        }
        
    }
    
    for (my $j=$previouspos;$j<$b[0];$j++)
    {
        $PosToGene{$j}="$previousgenenb"."_"."$nb";
    }
    $previouspos=$b[1]+1;
    $previousgenenb=$nb;
    $nb++;
}


for (my $j=$previouspos;$j<$genomeSize;$j++)
{
    $PosToGene{$j}="$previousgenenb"."_"."0";
}



close GENELIST;



##########################################################################################################################
#######count the codons of each type and deal here with mutations in the different overlapping reading frames#############
##########################################################################################################################

open (GENELIST, "<ECBGenesCDS.txt");

my $GComited=0;
my $ATomited=0;
my $integenicCGomited=0;
my $integenicATomited=0;
$nb=0;
my $genenb=0;
my $previous=1;
my $GC=0;
my $AT=0;
my $integenicCG=0;
my $integenicAT=0;
my $otherCG=0;
my $otherAT=0;
my $codon;
my $basestotal=0;
my $baseomited=0;
my $baseoverlap=0;

my @Synonymous=(0,0,0,0,0,0);
my @NonSynonymous=(0,0,0,0,0,0);

while( $line = <GENELIST>)
{
    chomp $line;	# drop \n at end of line
    @b=split(/\t/,$line);
    if ($b[4] eq "CDS")
    {
        for (my $j=$b[0];$j<=$b[1];$j++)
        {
           
            if (exists($MaskedBases{$j})==0 )
            {
                
                if (exists($BasedTreated{$j})==0)
                {
                    
                    $basestotal++;
                    
                    my $codonnb=int(($b[1]-$j)/3);
                    my $pos1=$b[1]-$codonnb*3-2;
                    my $pos2=$b[1]-$codonnb*3-1;
                    my $pos3=$b[1]-$codonnb*3;
                    
                    if ($b[2]eq"D")
                    {$codonnb=int(($j-$b[0])/3);
                        $pos1=$b[0]+$codonnb*3;
                        $pos2=$b[0]+$codonnb*3+1;
                        $pos3=$b[0]+$codonnb*3+2;
                    }
                    my $truncatedcodon=0;
                    if (exists($MaskedBases{$pos1}) or exists($MaskedBases{$pos2}) or exists($MaskedBases{$pos3}))
                    {
                        $truncatedcodon=1;
                    }
                    #if one of the base of the codon is covered by two codons
                    if  ($PosToNbCodon{$pos1}>1 or $PosToNbCodon{$pos2}>1 or $PosToNbCodon{$pos3}>1 or $truncatedcodon)
                    {
                        #$baseoverlap++;
                        my @codonslocal=split(/,/,$PosToCodon{$j});
                        my @strandlocal=split(/,/,$PosToStrand{$j});
                        my @codonposlocal=split(/,/,$PosToPosCodon{$j});
                        my @codonnumber=split(/,/,$PosToCodonPos{$j});

                        my $baseref=substr($refsequence,$j,1);
                        
                        foreach my $i (@MutType)
                        {
                            my $mutname="$i"."_$baseref";
                            if (exists $MutTypeAction{$mutname})
                            {
                                my $Syn=0;
                                my $NonSyn=0;
                                for(my $mm=0;$mm<$PosToNbCodon{$j};$mm++)
                                {
                                    #  print ",$codonslocal[$mm]";
                                    my $newbase=$MutTypeAction{$mutname};
                                    
                                    if ($strandlocal[$mm]==-1)
                                    {
                                        $newbase=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
                                    }
                                    my $newcodon=$codonslocal[$mm];
                                    substr $newcodon, $codonposlocal[$mm] ,1, $newbase;
                                    
                                    my $identical=0;
                                    if ($codonnumber[$mm]==0)
                                    {
                                        if ($CodeGeneticStart{$newcodon} eq $CodeGeneticStart{$codonslocal[$mm]})
                                        {$identical=1;}
                                    }
                                    else
                                    {
                                        if ($CodeGenetic{$newcodon} eq $CodeGenetic{$codonslocal[$mm]})
                                        {$identical=1;}
                                    }
                                    if ($identical==1)
                                    {
                                        $Syn++;
                                        
                                    }
                                    else
                                    {
                                        $NonSyn++;
                                    }
                                }
                                if ($NonSyn>0)
                                {
                                    $NonSynonymous[$i-1]++;
                                }
                                else
                                {
                                    $Synonymous[$i-1]++;
                                }
                                
                            }
                        }
                        
                        
                        $BasedTreated{$j}=1;
                    }
                    else
                    {
                        if ($b[2]eq"D")
                        {
                            
                            $codon=substr($refsequence, $b[0]+int(($j-$b[0])/3)*3, 3);
                            $codonnb=int(($j-$b[0])/3);
                            
                            if (int(($j-$b[0])/3)>0)
                            {$Codonlist{$codon}++;}
                            else
                            {$CodonlistStart{$codon}++;}
                            
                        }
                        
                        if ($b[2]eq"C")
                        {
                            
                            $codon=substr($refsequence, $b[1]-int(($b[1]-$j)/3)*3-2, 3);
                            $codonnb=int(($b[1]-$j)/3);
                            # print("$codon\n");
                            $codon=reverse($codon);
                            $codon=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
                            
                            
                            if (int(($b[1]-$j)/3)>0)
                            {$Codonlist{$codon}++;}
                            else
                            {$CodonlistStart{$codon}++;}
                            
                            
                        }
                        
                        
                    }
                    my $base=substr($refsequence,$j,1);
                    if ($base eq "G" || $base eq "C")
                        {
                               $GC++;
                           }
                           else
                           {
                               $AT++;
                           }
                }
            }
            else
            {
                $baseomited++;
                if ($b[2]eq"D" && $b[4] eq "CDS")
                {
                    $codon=substr($refsequence, $b[0]+int(($j-$b[0])/3)*3, 3);
                    
                }
                
                if ($b[2]eq"C" && $b[4] eq "CDS")
                {
                    $codon=substr($refsequence, $b[1]-int(($b[1]-$j)/3)*3-2, 3);
                    
                    # print("$codon\n");
                    $codon=reverse($codon);
                    $codon=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
                    
                }
                for (my $i=0;$i<3;$i++)
                {
                    
                    my $base=substr($codon,$i,1);
                    #  print "$base\t";
                    if ($base eq "G" || $base eq "C")
                    {
                        $GComited++;
                    }
                    else
                    {
                        $ATomited++;
                    }
                    
                }
            }
            
            #}
            
        }
    }
    else
    {
        for (my $j=$b[0];$j<=$b[1];$j++)
        {
            if ($j>=$previous)
            {
                if (exists($MaskedBases{$j})==0)
                {
                    $basestotal++;
                    my $base=substr($refsequence,$j,1);
                    
                    
                    if ($base eq 'G'|| $base eq 'C')
                    {
                        #$GC++;
                        $otherCG++;
                    }
                    else
                    {
                        $otherAT++;
                        #$AT++;
                    }
                }
                
                else
                {
                    $baseomited++;
                    my $base=substr($refsequence,$j,1);
                    if ($base eq'G'|| $base eq'C')
                    {
                        $GComited++;
                        $integenicCGomited++;
                    }
                    else
                    {
                        $ATomited++;
                        $integenicATomited++;
                    }
                }
            }
        }
    }
    
    for (my $j=$previous;$j<$b[0];$j++)
    {
        if (exists($MaskedBases{$j})==0)
        {
            $basestotal++;
            my $base=substr($refsequence,$j,1);
            if ($base eq 'G'|| $base eq 'C')
            {
                $GC++;
                $integenicCG++;
            }
            else
            {
                $integenicAT++;
                $AT++;
            }
        }
        else
        {
            $baseomited++;
            my $base=substr($refsequence,$j,1);
            if ($base eq'G'|| $base eq'C')
            {
                $GComited++;
                $integenicCGomited++;
            }
            else
            {
                $ATomited++;
                $integenicATomited++;
            }
            
        }
        
    }
    
    $previous=$b[1]+1;
    $nb++;
}

print "Synonymous: @Synonymous\n";
print "Non-Synonymous: @NonSynonymous\n";


my $len=length($refsequence);
print " Length seq: $len\n";
for (my $j=$previous;$j<=length($refsequence);$j++)
{
    if (exists($MaskedBases{$j})==0)
    {
        $basestotal++;
        my $base=substr($refsequence,$j,1);
        if ($base eq'G'|| $base eq'C')
        {
            $GC++;
            $integenicCG++;
        }
        else
        {
            $AT++;
            $integenicAT++;
        }
    }
    else
    {
        $baseomited++;
        my $base=substr($refsequence,$j,1);
        if ($base eq'G'|| $base eq'C')
        {
            $GComited++;
            $integenicCGomited++;
        }
        else
        {
            $ATomited++;
            $integenicATomited++;
        }
        
    }
}
close GENELIST;


my $sum=0;
foreach my $i (keys %Codonlist)
{
    #print "$i\t",$Codonlist{$i},"\n";
    $sum+=$Codonlist{$i};
}
print "number of codons in CDS $sum\n";








###################################################################################
###########For each codon find out the impact of all types of mutations ###########
###################################################################################
my %Codon_Mut_S=();
my %Codon_Mut_N=();
my @codonpos=(0,1,2);
print "@MutType\n";

foreach my $i (@MutType)
{
    foreach my $codon (keys %Codonlist)
    {
        my $N=0;
        my $S=0;
        $Codon_Mut_S{"$codon"."_$i"}="_";
        $Codon_Mut_N{"$codon"."_$i"}="_";
        foreach my $j (@codonpos)
        {
            my $newcodon=$codon;
            my $aa=$CodeGenetic{$codon};
            my $newaa=$aa;
            
            my $baseref=substr($codon,$j,1);
            my $mutname="$i"."_$baseref";
            if (exists $MutTypeAction{$mutname})
            {
                substr($newcodon, $j, 1)=$MutTypeAction{$mutname};
                $newaa=$CodeGenetic{$newcodon};
                
                if ($newaa eq $aa)
                {
                    $S++;
                }else
                {
                    $N++;
                }
                
                #here compute the number of AA change
                
                
            }
        }
        $Codon_Mut_S{"$codon"."_$i"}=$S;
        $Codon_Mut_N{"$codon"."_$i"}=$N;
        
    }
    
    
}
###################################################################################
###########For each start codon find out the impact of all types of mutations ###########
###################################################################################
my %Codon_Mut_S_start=();
my %Codon_Mut_N_start=();

foreach my $i (@MutType)
{
    foreach my $codon (keys %CodonlistStart)
    {
        my $N=0;
        my $S=0;
        $Codon_Mut_S_start{"$codon"."_$i"}="_";
        $Codon_Mut_N_start{"$codon"."_$i"}="_";
        foreach my $j (@codonpos)
        {
            my $newcodon=$codon;
            my $aa=$CodeGeneticStart{$codon};
            my $newaa=$aa;
            
            my $baseref=substr($codon,$j,1);
            my $mutname="$i"."_$baseref";
            if (exists $MutTypeAction{$mutname})
            {
                substr($newcodon, $j, 1)=$MutTypeAction{$mutname};
                $newaa=$CodeGeneticStart{$newcodon};
                
                if ($newaa eq $aa)
                {
                    $S++;
                }else
                {
                    $N++;
                }
                
                #here compute the number of AA change
                
                
            }
        }
        $Codon_Mut_S_start{"$codon"."_$i"}=$S;
        $Codon_Mut_N_start{"$codon"."_$i"}=$N;
        
    }
    
    
}





###############################################################################
##############Now deal with the codons found in non overlapping frames #######
###############################################################################


print "Codon\tQuantity\tS1_AT_to_CG\tS2_AT_to_GC\tS3_AT_to_TA\tS4_CG_to_GC\tS5_CG_to_TA\tS6_GC_to_TA\tN1_AT_to_CG\tN2_AT_to_GC\tN3_AT_to_TA\tN4_CG_to_GC\tN5_CG_to_TA\tN6_GC_to_TA\n";



foreach my $i (keys %Codonlist)
{
    print "$i\t$Codonlist{$i}\t";
    foreach my $j (@MutType)
    {
        print $Codon_Mut_S{"$i"."_$j"},"\t";
        $Synonymous[$j-1]+=$Codon_Mut_S{"$i"."_$j"}*$Codonlist{$i}/3;##ajouter les start
    }
    foreach my $j (@MutType)
    {
        print $Codon_Mut_N{"$i"."_$j"},"\t";
        $NonSynonymous[$j-1]+=$Codon_Mut_N{"$i"."_$j"}*$Codonlist{$i}/3;##ajouter les start
        
    }
    print "\n";
    
}
print "\n";print "\n";print "\n";

foreach my $i (keys %CodonlistStart)
{
    print "$i\t$CodonlistStart{$i}\t";
    foreach my $j (@MutType)
    {
        print $Codon_Mut_S_start{"$i"."_$j"},"\t";
        $Synonymous[$j-1]+=$Codon_Mut_S_start{"$i"."_$j"}*$CodonlistStart{$i}/3;##ajouter les start
    }
    foreach my $j (@MutType)
    {
        print $Codon_Mut_N_start{"$i"."_$j"},"\t";
        $NonSynonymous[$j-1]+=$Codon_Mut_N_start{"$i"."_$j"}*$CodonlistStart{$i}/3;##ajouter les start
        
    }
    print "\n";
    
}
print "Synonymous: @Synonymous\n";
print "Non-Synonymous: @NonSynonymous\n";
print "Intergenic: $integenicAT,$integenicAT,$integenicAT,$integenicCG,$integenicCG,$integenicCG\n";
print "Genomic: $AT,$AT,$AT,$GC,$GC,$GC\n";
print "$basestotal\t$baseomited\n";
print "$GComited\t$ATomited\t$integenicCGomited\t$integenicATomited\n";
print "$baseoverlap\t$overlaptest\t$overlap2\n";

open (GENELISTOUT, ">GenomeComposition.txt");
printf GENELISTOUT "Type AT_CG AT_GC AT_TA CG_GC CG_TA GC_TA\n";
printf GENELISTOUT "Synonymous @Synonymous\n";
printf GENELISTOUT "Non-Synonymous @NonSynonymous\n";
printf GENELISTOUT "Intergenic $integenicAT $integenicAT $integenicAT $integenicCG $integenicCG $integenicCG\n";
printf GENELISTOUT "genomicGC $AT $AT $AT $GC $GC $GC\n";
close GENELISTOUT;

