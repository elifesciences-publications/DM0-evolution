
#This file uses several files to work
#the data are stored in either oli.LTEE.final_masked.no_IS_adjacent.tab or oli.MAE.final_masked.no_IS_adjacent.tab if the Mutation Accumalation experiemnt is used
#It also uses REL606.gff3 to have gene positions and genome sequence.
#it is not a generic programs as it uses the specificities of the data such as the sampling time points, which populations and clones are mutators. There are therefore no parameters. 

#The same program was used fo both the Mutation accumulation experiment  (MAE) and the Long term experiemntal evolution. The results of the MAE experiments are been included directly in the code (lines 588, 589) such that the first two lines of the MutationTypesThroughTime.txt are resulst from the MAE.

#the outputs are
#ECBGenesCDS.txt a silplified gene position file
#MutRArray.txt an array like file to be used to produce figure 2 and the phylogeny, it contains presence and absence of all point mutations in all sequenced clones
#TotalRArray.txt an array like file to be used to produce figure 2 and the phylogeny, it contains presence and absence of all  mutations in all sequenced clones
#EventByStrainsSorte.txt a simplified file to store the information on each mutation

#this version only computes mutations that are not within the segments defined by MatchedRegions

use List::Util qw( min max );

######################################################################################
#############################define input files and output files #####################
######################################################################################

my $numArgs = $#ARGV + 1;
my $cmd=$ARGV[0];

open (CMD, "<$cmd") or die "Command file $cmd not found !\n";
#Get genome files

$line = <CMD>;
chomp $line;
@b=split(/\t/,$line);
my $gff3file=$b[1];

#get inputfile
$line = <CMD>;
chomp $line;
@b=split(/\t/,$line);
my $inputfile=$b[1];

#get mutator status and sampling times
$line = <CMD>;
chomp $line;
@b=split(/\t/,$line);
my %MutatorPopulationSwitch=();
my %CloneMutatorStatus=();
my %Timessampled=();
open (CLONE, "<$b[1]") or die "Mutator_info file $b[1] not found !\n";
$line = <CLONE>;
while( $line = <CLONE>)
{
    chomp $line;	# drop \n at end of line
    @b=split(/\t/,$line);
    if (!exists($Timessampled{$b[2]}))
    {
        $Timessampled{$b[2]}=$b[2];
    }
    
    $CloneMutatorStatus{$b[3]}=$b[4];
    if (exists($MutatorPopulationSwitch{$b[0]}))
    {
        if ($b[4]==0 && $b[2]>$MutatorPopulationSwitch{$b[0]})
        {
            $MutatorPopulationSwitch{$b[0]}=$b[2];
        }
    }else
    {
        if ($b[4]==0)
        {
            $MutatorPopulationSwitch{$b[0]}=$b[2];
        }
    }
}
close CLONE;

#foreach my $pp (keys %MutatorPopulationSwitch)
#{
#    print "$pp $MutatorPopulationSwitch{$pp}\n";
#}



@timepoints=sort {$a<=>$b} keys(%Timessampled);
print "@timepoints\n";

#get mutator action
$line = <CMD>;
chomp $line;
@b=split(/\t/,$line);
my $mutatorOnly=$b[1];
print "Mutator to be analysed: $mutatorOnly\n";
#define outputfile names
$line = <CMD>;
chomp $line;
@b=split(/\t/,$line);
my $suffix=$b[1];
if ($suffix eq "none")
{
    $suffix="";
}
$filenameEvent="EventByStrainsSorted".$suffix.".txt";
$MutationSummary= "MutationTypesThroughTime".$suffix.".txt";

close CMD;



my $previouspop;
my $previousage;
my $previousclone;



our(%CodeGenetic)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
our(%CodeGeneticStart)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'M','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'M','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'M','ATC'=>'M','ATT'=>'M','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'M','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');


our (%Degeneracy)=('TCA'=>4,'TCC'=>4,'TCG'=>4,'TCT'=>4,'TTC'=>2,'TTT'=>2,'TTA'=>2,'TTG'=>2,'TAC'=>2,'TAT'=>2,'TAA'=>2,'TAG'=>2,'TGC'=>2,'TGT'=>2,'TGA'=>1,'TGG'=>1,'CTA'=>4,'CTC'=>4,'CTG'=>4,'CTT'=>4,'CCA'=>4,'CCC'=>4,'CCG'=>4,'CCT'=>4,'CAC'=>2,'CAT'=>2,'CAA'=>2,'CAG'=>2,'CGA'=>4,'CGC'=>4,'CGG'=>4,'CGT'=>4,'ATA'=>3,'ATC'=>3,'ATT'=>3,'ATG'=>1,'ACA'=>4,'ACC'=>4,'ACG'=>4,'ACT'=>4,'AAC'=>2,'AAT'=>2,'AAA'=>2,'AAG'=>2,'AGC'=>2,'AGT'=>2,'AGA'=>2,'AGG'=>2,'GTA'=>4,'GTC'=>4,'GTG'=>4,'GTT'=>4,'GCA'=>4,'GCC'=>4,'GCG'=>4,'GCT'=>4,'GAC'=>2,'GAT'=>2,'GAA'=>2,'GAG'=>2,'GGA'=>4,'GGC'=>4,'GGG'=>4,'GGT'=>4);



my %syncodon=();
my %nonsyncodon=();

my %syntotal=();
my %nonsyntotal=();






##############################################################################################################################
#################use the reference genome to get gene positions and genome sequence    #######################################
##############################################################################################################################

open (GENELIST, "<$gff3file") or die "gff3 file $gff3file not found !\n";
open (GENELISTOUT, ">ECBGenesCDS.txt");
$line = <GENELIST>;
$line = <GENELIST>;
$line = <GENELIST>;
$nb=0;
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
        print GENELISTOUT "$b[3]\t$b[4]\t$strand\t$name\t$b[2]\n";
        $nb++;
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

my @GeneName=(1..$nb);
our %PosToGene=();
our %PosToCodon=();
our %PosToStrand=();
our %PosToPosCodon=();
our %PosToCodonPos=();
our %PosToNbCodon=();



my $genomeSize=length $refsequence;
$genomeSize--;

##############################################################################################################################
######### use the simplified gene file to compute for each position in the genome a coding status ############################
##############################################################################################################################

open (GENELIST, "<ECBGenesCDS.txt");

$nb=0;
my $previouspos=1;
my $previousgenenb=$nb-1;
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
                    $PosToCodon{$j}="$PosToCodon{$j},$codon";
                    my $poscod=($b[1]-$j)%3;
                    $PosToPosCodon{$j}="$PosToPosCodon{$j},$poscod";
                    $poscod=int(($b[1]-$j)/3);
                    $PosToCodonPos{$j}="$PosToCodonPos{$j},$poscod";
                    $PosToStrand{$j}="$PosToStrand{$j},-1";
                    $PosToNbCodon{$j}=$PosToNbCodon{$j}+1;
                }
                
                
            }
        }
        else
        {
            if (exists($PosToCodon{$j})==0)
            {
                $PosToCodon{$j}="non_coding";
                $PosToPosCodon{$j}=-1;
                $PosToStrand{$j}=$b[2];
                $PosToCodonPos{$j}=10;
                $PosToNbCodon{$j}=1;
            }
            else
            {
                $PosToCodon{$j}="$PosToCodon{$j},non_coding";
                $PosToPosCodon{$j}="$PosToPosCodon{$j},-1";
                $PosToCodonPos{$j}="$PosToCodonPos{$j},10";
                $PosToStrand{$j}="$PosToStrand{$j},$b[2]";
                $PosToNbCodon{$j}=$PosToNbCodon{$j}+1;
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


##############################################################################################################################
########first read of data to identify all mutations with a unique tag per population ########################################
##############################################################################################################################



my %popevent=();
my %popeventage=();

my %fold4=("A->C"=>0,"A->G"=>0,"A->T"=>0,"C->G"=>0,"C->T"=>0,"G->T"=>0);#,"T->G"=>0,"T->C"=>0,"T->A"=>0,"G->C"=>0,"G->A"=>0,"C->A"=>0);
my %mutationtype=("A->C"=>0,"A->G"=>1,"A->T"=>2,"C->G"=>3,"C->T"=>4,"G->T"=>5,"T->G"=>0,"T->C"=>1,"T->A"=>2,"G->C"=>3,"G->A"=>4,"C->A"=>5);


open STRAINOUT, ">$MutationSummary";
close STRAINOUT;

open STRAINOUT, ">$filenameEvent";
close STRAINOUT;



open STRAINTOEVENT, "<$inputfile" or die "input file $inputfile not found!\n";

while(my $line = <STRAINTOEVENT>)
{
    for $gg (@timepoints)
    {
        my $age=$gg;
        chomp $line;
        @b=split(/\t/,$line);
        
        my $out=0;
        if ($b[5] ne "boubou")
        {  #potential selection for just a mutation Type
            
            
            if ($out==0 && $b[1]==$age)
            {
                
                if (exists($popevent{"$b[0]"."_$b[3]"}))
                {
                    $popevent{"$b[0]"."_$b[3]"}=$popevent{"$b[0]"."_$b[3]"}."_".$CloneMutatorStatus{$b[2]};
                }
                else
                {
                    $popevent{"$b[0]"."_$b[3]"}=$CloneMutatorStatus{$b[2]};
                }
                if (exists($popeventage{"$b[0]"."_$b[3]"}))
                {
                    $popeventage{"$b[0]"."_$b[3]"}=$popeventage{"$b[0]"."_$b[3]"}."_".$age;
                }
                else
                {
                    $popeventage{"$b[0]"."_$b[3]"}=$age;
                }

            }
        }
    }
}
close 	STRAINTOEVENT;
##############################################################################################################################
################################### Attribute to each mutation a coding status and filter#####################################
##############################################################################################################################

open STRAINTOEVENT, "<$inputfile";
open STRAINOUT, ">>$filenameEvent";

my %foundevent=();


while(my $line = <STRAINTOEVENT>)
{
    for $gg (@timepoints)
    {
        my $age=$gg;
        chomp $line;
        @b=split(/\t/,$line);
        
        my $out=0;
        if ($mutatorOnly==1)
        {
            $out=1;#exclusion of all clones
        }
        if ($b[5] ne "boubou"){  #potential selection for just Mutation Type
            
            if ($mutatorOnly==1)
            {
                if ($CloneMutatorStatus{$b[2]}==1  && $b[1]>$MutatorPopulationSwitch{$b[0]})
                {
                $out=0;#inclusion of mutator clones
                }
            }
            if ($mutatorOnly==2)
            {
                if ($CloneMutatorStatus{$b[2]}==1)
                {
                $out=1;#inclusion of mutator clones
                }
            }
            
            
            if ($out==0 && $b[1]==$age)
            {
                #Extract the length of events
                
                my $len;
                if ($b[6]=~/(\d+)_bp/)
                {
                    $len=$1;
                }
                else
                {
                    $len=1;
                    
                }
                # print $len,"\n";
                $b[7]=$len;
                
                # For mutations find the Synonyous or Non Synonymous state
                if ($b[5] eq "Mutation")
                {
                    $b[6]=~/\w-\>(\w)/;
                    
                    $b[7]=$mutationtype{$b[6]};
                    my $newbase=$1;
                    my $type="I";
                    if (exists $PosToCodon{$b[4]})
                    {
                        
                        
                        if ($PosToNbCodon{$b[4]}==1)
                        {
                            
                            if ($PosToCodon{$b[4]} ne "non_coding")
                            {
                                
                                if ($PosToStrand{$b[4]}==-1)
                                {
                                    $newbase=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
                                }
                                my $newcodon=$PosToCodon{$b[4]};
                                substr $newcodon, $PosToPosCodon{$b[4]} ,1, $newbase;
                                #print("newbase: $newbase pos $PosToPosCodon{$b[4]}  codon: $PosToCodon{$b[4]}   AA: $CodeGenetic{$PosToCodon{$b[4]}}   new codon : $newcodon\n");
                                
                                #$PosToCodonPos
                                my $identical=0;
                                if ($PosToCodonPos{$b[4]}==0)
                                {
                                    if ($CodeGeneticStart{$newcodon} eq $CodeGeneticStart{$PosToCodon{$b[4]}})
                                    {$identical=1;}
                                }
                                else
                                {
                                    if ($CodeGenetic{$newcodon} eq $CodeGenetic{$PosToCodon{$b[4]}})
                                    {$identical=1;}
                                }
                                if ($identical==1)#$CodeGenetic{$newcodon} eq $CodeGenetic{$PosToCodon{$b[4]}})
                                {
                                    $type="S";
                                    
                                }
                                else
                                {
                                    $type="N";
                                }
                            }
                            else
                            {
                                $type="O";
                            }
                        }else
                        {
                            #print "AFFFFFF   $b[4],$PosToNbCodon{$b[4]}\n";
                            my $Syn=0;
                            my $NonSyn=0;
                            my $Other=0;
                            $b[6]=~/\w-\>(\w)/;
                            my $newbaseref=$1;
                            #here loop on all codons:
                            my @codonslocal=split(/,/,$PosToCodon{$b[4]});
                            my @strandlocal=split(/,/,$PosToStrand{$b[4]});
                            my @codonposlocal=split(/,/,$PosToPosCodon{$b[4]});
                            my @codonnumber=split(/,/,$PosToCodonPos{$b[4]});
                            
                            for(my $mm=0;$mm<$PosToNbCodon{$b[4]};$mm++)
                            {
                                #  print ",$codonslocal[$mm]";
                                if ($codonslocal[$mm] eq "non_coding")
                                {
                                    $Other++;
                                    next;
                                }
                                my $newbase=$newbaseref;
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
                                # print "non syn $b[4],$Syn,$NonSyn,$Other\n";
                                $type="N";
                            }
                            else
                            {
                                if ($Syn>0)
                                {$type="S";}
                                else
                                {$type="O";}
                                # print "non coding $b[4],$Syn,$NonSyn,$Other\n";
                            }
                            
                        }
                        
                    }
                    $b[6]=$type;
                    
                }
                
                $line=join ("\t",@b);
                
                if (exists($popevent{"$b[0]"."_$b[3]"}))
                {
                    # @a=split("_",$popevent{"$b[0]"."_$b[3]"});
                    
                    #my $min = min @a;
                    #my $max = max @a;
                    
                    @a=split("_",$popeventage{"$b[0]"."_$b[3]"});
                    
                    my $min = min @a;
                    my $max = max @a;
                    
                    
                    
                    if ($mutatorOnly==1)
                    {
                        #get the pop involved:$b[0]
                        if (  $CloneMutatorStatus{$b[2]}==1)
                            {
                                #if ($min>0)
                                #{
                                #$line=~tr/\+\-/pm/;
                                #printf (STRAINOUT "$line\n");
                                #}
                                if ($min>$MutatorPopulationSwitch{$b[0]})
                                {
                                    $line=~tr/\+\-/pm/;
                                    printf (STRAINOUT "$line\n");
                                }
                            }
                    }else
                    {
                        $line=~tr/\+\-/pm/;
                        printf (STRAINOUT "$line\n");
                    }
                    # }
                }
                else
                {
                    $line=~tr/\+\-/pm/;
                    printf (STRAINOUT "$line\n");
                    
                    print "sdgsdlfsdlfqlkhfdl\n";
                    
                    
                }
                $foundevent{"$b[0]"."_$b[3]"}=1;
                
                #  if ($age>=0)
                #{$popevent{"$b[0]"."_$b[3]"}=1;}
            }
            
            
        }
        #print("$b[0]_$b[1]_$b[2]-$b[3]         $b[4],$b[7],$b[5],$b[6],$genelist\n");
        
    }
}
close 	STRAINTOEVENT;
close 	STRAINOUT;



##############################################################################################################################
######################################Count mutations per categories at any given sampling time##############################
##############################################################################################################################

open STRAINTOEVENT, "<$filenameEvent";
open STRAINOUT, ">>$MutationSummary";



#here are the results of the analysis made for the MAE.
printf STRAINOUT "fullclone\tpop\tage\tclone\tmuts\tis\tldels\tindels\tis150\tsyn\tnonsyn\tinter\ts1\ts2\ts3\ts4\ts5\ts6\tns1\tns2\tns3\tns4\tns5\tns6\ti1\ti2\ti3\ti4\ti5\ti6\tothers\n";
#printf STRAINOUT "MAEtotal\tMAE\t14300\tMAE\t199\t41\t28\t41\t14\t39\t117\t39\t1\t4\t1\t3\t17\t13\t9\t13\t8\t5\t34\t48\t3\t14\t3\t0\t8\t11\t2\n";#all muts nofilter
#printf STRAINOUT "MAEtotal\tMAE\t14300\tMAE\t197\t41\t28\t41\t14\t39\t117\t39\t1\t4\t1\t3\t17\t13\t9\t13\t8\t5\t34\t48\t3\t14\t3\t0\t8\t11\t2\n";#IS filtered mutations




for $gg (@timepoints)
{
    my $age=$gg;
    print "$gg\n";
    my $inter;
    my $muts;
    my $is;
    my $is150;
    my $indels;
    my $other;
    my $name;
    my $nameprevious;
    my $first;
    my $largeindels;
    my $syn;
    my $nonsyn;
    my @syns;
    my @nonsyns;
    my @inters;
    my @others;
    $nameprevious="bobo";
    $first=0;
    while(my $line = <STRAINTOEVENT>)
    {
        chomp $line;
        #print $line,"\n";
        @b=split(/\t/,$line);
        my $out=0;
        $name="$b[0]_$b[1]_$b[2]";
        
        if ($name ne $nameprevious)
        {
            #finish analysis of  previous sequence and print
            if ($first>0)
            {
                printf STRAINOUT "$nameprevious\t$previouspop\t$previousage\t$previousclone\t$muts\t$is\t$largeindels\t$indels\t$is150\t$syn\t$nonsyn\t$inter\t@syns\t@nonsyns\t@inters\t$other\n";
                
                
            }
            $first=1;
            # initialise
            $is=$muts=$other=$largeindels=$is150=$syn=$nonsyn=$inter=$indels=0;
            @syns=(0,0,0,0,0,0);
            @nonsyns=(0,0,0,0,0,0);
            @inters=(0,0,0,0,0,0);
            @others=(0,0,0,0,0,0);
        }
        
        
        
        
        
            if ($b[5] eq "Mutation" )
            {
                $muts++;
                if ($b[6] eq "S")
                {
                    $syn++;
                    $syns[$b[7]]++;
                }
                else
                {
                    if ($b[6] eq "N")
                    {
                        $nonsyn++;
                        $nonsyns[$b[7]]++;
                    }
                    else
                    {
                        if ($b[6] eq "O")
                        {
                            $other++;
                            $others[$b[7]]++;
                        }
                        else
                        {
                            $inter++;
                            $inters[$b[7]]++;
                        }
                    }
                }
                #print "$muts\n";
            }
            else
            {
                if ($b[5]=~/IS/)
                {
                    $is++;
                    if ($b[6]=~/IS150/)
                    {
                        $is150++;
                        
                    }
                }
                else
                {
                    if ($b[7]<50)
                    {$indels++;}
                    else
                    {$largeindels++;}
                    
                    #if ($b[5]=~/LargeDel/|| $b[5]=~/Duplication/)
                    #{
                    #    $largeindels++;
                    #}
                    #else
                    #{
                    #    $others++;
                    #}
                }
            }
            
            
            
        
        
        $nameprevious=$name;
        $previouspop=$b[0];
        $previousage=$b[1];
        $previousclone=$b[2];
        
        #print("$b[0]_$b[1]_$b[2]-$b[3]         $b[4],$b[7],$b[5],$b[6],$genelist\n");
        #print("$muts $is $others $largeindels\n");
        
        
    }
    
    printf STRAINOUT "$nameprevious\t$b[0]\t$b[1]\t$b[2]\t$muts\t$is\t$largeindels\t$indels\t$is150\t$syn\t$nonsyn\t$inter\t@syns\t@nonsyns\t@inters\t$other\n";
    
    
    close 	STRAINTOEVENT;
    close 	STRAINOUT;
    
    
    
}



############################################################################
########## Print an array for R to build a phylogeny########################
############################################################################
print "PRINT FOR PHYLOGENY\n";
my %StrainList=();
my %UniqueEventList=();
open STRAINTOEVENT, "<$filenameEvent";
while(my $line = <STRAINTOEVENT>)
{
    chomp $line;
    @b=split(/\t/,$line);
    
    $UniqueEventList{$b[3]}=$line;
    $StrainList{$b[2]}=1;
}
close 	STRAINTOEVENT;


open ALLEVENTARRAY, ">TotalRArray".$suffix.".txt";
open MUTARRAY, ">MutRArray".$suffix.".txt";


#first line gives teh mutation number

printf MUTARRAY "Pop\tAge\tRun\tName\t";
printf ALLEVENTARRAY "Pop\tAge\tRun\tName\t";
my $keys;
foreach $keys (sort keys %UniqueEventList)
{
    printf ALLEVENTARRAY "$keys\t";
    if ($UniqueEventList{$keys}=~/Mutation/)
    {
        printf MUTARRAY "$keys\t";
    }
    
}

#second line  gives the mutation type (SNI, K)
printf MUTARRAY "\nPop\tAge\tRun\tName\t";
printf ALLEVENTARRAY "\nPop\tAge\tRun\tName\t";
foreach $keys (sort keys %UniqueEventList)
{
    
    if ($UniqueEventList{$keys}=~/Mutation/)
    {
        @b=split(/\t/,$UniqueEventList{$keys});
        printf MUTARRAY "$b[6]\t";
        printf ALLEVENTARRAY "$b[6]\t";
        
    }
    else
    {
        printf ALLEVENTARRAY "K\t";
    }
    
}


#third line give 1 for all mutations
printf MUTARRAY "\nPop\tAge\tRun\tName\t";
printf ALLEVENTARRAY "\nPop\tAge\tRun\tName\t";
foreach $keys (sort keys %UniqueEventList)
{
    
    if ($UniqueEventList{$keys}=~/Mutation/)
    {
        printf MUTARRAY "1\t";
    }
    printf ALLEVENTARRAY "1\t";
    
}

#Fourth line give O for the ancestor
printf MUTARRAY "\nAncestor\t0\tREL606\tAncestor_0_REL606\t";
printf ALLEVENTARRAY "\nAncestor\t0\tREL606\tAncestor_0_REL606\t";


foreach $keys (sort keys %UniqueEventList)
{
    if ($UniqueEventList{$keys}=~/Mutation/)
    {
        printf MUTARRAY "0\t";
    }
    
    printf ALLEVENTARRAY "0\t";
    
}

printf MUTARRAY "\n";
printf ALLEVENTARRAY "\n";


#later lines give the presence and absence for each strain
my $strain;
foreach $strain (sort keys %StrainList)
{
    # print("$strain","_gogo","\n");
    $strain=~/(A.*)_(\d+)gen_(.*)/;
    
    printf MUTARRAY "$1\t$2\t$3\t$strain\t";
    printf ALLEVENTARRAY "$1\t$2\t$3\t$strain\t";
    foreach $keys (sort keys %UniqueEventList)
    {
        my  $state=0;
        if ($UniqueEventList{$keys}=~/$strain/)
        {
            $state=1;
        }
        if ($UniqueEventList{$keys}=~/Mutation/)
        {
            printf MUTARRAY "$state\t";
        }
        printf ALLEVENTARRAY "$state\t";
        
    }
    printf MUTARRAY "\n";
    printf ALLEVENTARRAY "\n";
}
close MUTARRAY;
close ALLEVENTARRAY;




