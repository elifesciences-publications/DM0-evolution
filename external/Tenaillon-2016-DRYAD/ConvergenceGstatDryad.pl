


#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/ceil/;


#############################################################################
#####################################read inputs ############################
#############################################################################

my $numArgs = $#ARGV + 1;

if ($numArgs<5)
{
    print "wrong number of arguments: perl ConvergenceGstatDryad.pl REL606.gff3 REL606.L20.G15.P0.M35.mask.gd EventByStrainsSorted.txt 1000 outputfolder\n";
    die;
}
my $gff3=$ARGV[0];
my $maskfile=$ARGV[1];#REL606.L20.G15.P0.M35.mask.gd
my $inputfile=$ARGV[2];#EventByStrainsSorted.txt
my $nbrepetsimuls=$ARGV[3];
my $outputfolder=$ARGV[4];#ConvergenceByTimeDataMutator



my $found;
my @b;
my @c;
my @e;
my $k;
my $keya;
my $keyb;
my $sim;
my ($nb,$sum);
my $listcombined;
my @d;
my %GeneToNb=();
my %NameToB;
my %ColiBSynonymous;
my %NameToBDataSet;
my %BToNameDataSet;
our %Convergenceperlineanddate=();
our %StrainToEvent;
our %StrainToUnit;
our %StrainToGene;
#my %EventToGene;
my %ListFunctionalUnits;
our %StrainEventTemp=();
our %EventToGeneTemp=();
my %GeneToUnit;
my $i;
my $j;
my $UnitName;
our %ConvergenceFromMutationPointofViewEvent;
our %ConvergenceFromMutationPointofViewGene;
our %MutationtoClass;
#our %StrainGene_To_Event;
our(%CodeGenetic)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

sub max {
    my $max = $_[0];
    if ($max>$_[1])
	{
        return($max);
	}
    else
	{
        return($_[1])
	}
}







################################################################################################
#########################Get masked regions#####################################################
################################################################################################


our %PosToGene=();
our %PosToCodon=();
our %PosToStrand=();
our %PosToPosCodon=();
#open (GENELIST, "<ECBgenes.txt");
# the file has the folowing structure
#190	255	D	thrL
#336	2798	D	thrA
#2800	3732	D	thrB
#3733	5019	D	thrC

my %MaskedBases=();
@b=();
open (MASKLIST, "< $maskfile") or die "mask file $maskfile not found!\n";
my $line = <MASKLIST>;
while( $line = <MASKLIST>)
{
    chomp $line;
    @b=split(/\t/,$line);
    for (my $i=$b[4]; $i<$b[4]+$b[5]-1;$i++)
    {$MaskedBases{$i}=1;}
}


close MASKLIST;



################################################################################################
#########################Get genome and create simplified gene file#############################
################################################################################################
$nb=0;
open (GENELIST, "< $gff3") or die "gff3 file $gff3 not found!\n";
;
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
my @GeneLength=(1..$nb);
my @GenePosition=(1..$nb);
my @GeneNonCodingLength=(0) x $nb;


my $genomeSize=length $refsequence;
$genomeSize--;


################################################################################################
#########################Assign a gene to every position in the genome##########################
################################################################################################
open (GENELIST, "<ECBGenesCDS.txt");
$nb=0;
my $totalcoding=0;
my $tatoltast=0;
my $previousend=0;
while( $line = <GENELIST>)
{
    chomp $line;	# drop \n at end of line
    @b=split(/\t/,$line);
    #print "$b[3]\n"; #The gene name is in column 3
    $GeneName[$nb]=$b[3];
    $GenePosition[$nb]=$b[0];
    $GeneLength[$nb]=0;
    for (my $j=$b[0];$j<=$b[1];$j++)
    {
        if (exists($MaskedBases{$j})==0)
        {
            
            if ($b[4]eq "CDS")
            {
                $GeneLength[$nb]++;
                $totalcoding++;
            }
            
        }
    }
    $previousend=$b[1];
    $tatoltast=$tatoltast+$GeneLength[$nb];
    $nb++;
}
print $totalcoding,"\n",$nb,"\n",$tatoltast,"\n";
my $TotalGeneNb=$nb;
close GENELIST;

#print "$refsequence\n";
open (GENELIST, "<ECBGenesCDS.txt");

my %NoncodingtoGenomePos=();

my $nbnoncoding=0;
$nb=0;
my $previouspos=1;
my $previousgenenb=$TotalGeneNb-1;
my $accoutnednbnoncoding=0;

while( $line = <GENELIST>)
{
	chomp $line;	# drop \n at end of line
	#print "$line\n";
	@b=split(/\t/,$line);
	#print "\nvoila $b[3]\n"; The gene name is in column 3
    for (my $j=$b[0];$j<=$b[1];$j++)
    {
        $PosToGene{$j}=$nb;
        
        if ($b[2]eq"D" && $b[4] eq "CDS"&& exists($PosToCodon{$j})==0)
        {
            $PosToCodon{$j}=substr($refsequence, $b[0]+int(($j-$b[0])/3)*3, 3);
            $PosToPosCodon{$j}=($j-$b[0])%3;
            $PosToStrand{$j}=1;
            # print("$j ",substr($refsequence, $j, 1),"   $b[0] $PosToCodon{$j}\t$CodeGenetic{$PosToCodon{$j}}\t\t$PosToPosCodon{$j}\n");
            
        }
        
        if ($b[2]eq"C" && $b[4] eq "CDS"&& exists($PosToCodon{$j})==0)
        {
            my $codon=substr($refsequence, $b[1]-int(($b[1]-$j)/3)*3-2, 3);
            
            # print("$codon\n");
            $codon=reverse($codon);
            $codon=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
            $PosToCodon{$j}=$codon;
            
            $PosToPosCodon{$j}=($b[1]-$j)%3;
            $PosToStrand{$j}=-1;
            my $base=substr($refsequence, $j, 1);
            $base=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
            #print("$j $b[1] $base $codon\t$CodeGenetic{$codon}\t$PosToPosCodon{$j}\n");
            
            
        }
        
        if ($b[4]ne "CDS")
        {
            if (exists($MaskedBases{$j})==0)
            {
                $NoncodingtoGenomePos{$nbnoncoding}=$j;
                $GeneNonCodingLength[$nb]++;
                $nbnoncoding++;
                $accoutnednbnoncoding++;
            }
        }
        
        
        
    }
    
    for (my $j=$previouspos;$j<$b[0];$j++)
    {
        if (!exists($PosToGene{$j}))
        {$PosToGene{$j}="$previousgenenb"."_"."$nb";
       
        
        if (exists($MaskedBases{$j})==0)
        {
        $NoncodingtoGenomePos{$nbnoncoding}=$j;
        $nbnoncoding++;
        $GeneNonCodingLength[$nb]++;
        $GeneNonCodingLength[$previousgenenb]++;
        $accoutnednbnoncoding+=2;
        }
        }
    }
    
    $previouspos=$b[1]+1;
    $previousgenenb=$nb;
	$nb++;
    
}




for (my $j=$previouspos;$j<=$genomeSize;$j++)
{
    if (!exists($PosToGene{$j}))
    {$PosToGene{$j}="$previousgenenb"."_"."0";
    
    if (exists($MaskedBases{$j})==0)
    {
    $NoncodingtoGenomePos{$nbnoncoding}=$j;
    $nbnoncoding++;
    $GeneNonCodingLength[$previousgenenb]++;
    $GeneNonCodingLength[0]++;
    $accoutnednbnoncoding+=2;
    }
    }
}

print "$accoutnednbnoncoding\n";
close GENELIST;





################################################################################################
#########################Compute non coding size################################################
################################################################################################



my $NoncodingSize=0;
for (my $j=0;$j<$TotalGeneNb;$j++)
{
    $NoncodingSize+=$GeneNonCodingLength[$j];
    
}
print "Noncoding $NoncodingSize\n";
$NoncodingSize=$nbnoncoding;




################################################################################################
#########################Find the genes affected and increment the gene_hit table###############
################################################################################################



sub GetGenes3
{
    my ($debut,$length,$genomeSize,$TotalGeneNb,$GeneName_ref,$GeneTargeti_ref,$GeneTargetType_ref,$GeneTargetLength_ref)=@_;
    
    my @GeneName=@$GeneName_ref;
    my @GeneTargeti=@$GeneTargeti_ref;
    my @GeneTargetType=@$GeneTargetType_ref;
    my @GeneTargetLength=@$GeneTargetLength_ref;
    
    my $genelist="";
    my $finaladd="";
    
    
    
    my $firstgenenb=$PosToGene{$debut};
    
    #print("$debut\t$firstgenenb\tggggg\n");
    my @c=split(/_/,$firstgenenb);
    #print("@c\n");
    my $size=@c;
    if ($size>1)
    {
        $genelist="$GeneName[$c[0]]-$c[0],$GeneName[$c[1]]-$c[1]";
        $GeneTargeti[$c[0]]++;
        $GeneTargetType[$c[0]]++;
        $firstgenenb=$c[1];
        $GeneTargeti[$firstgenenb]++;
        $GeneTargetType[$firstgenenb]++;
        
        if ($length>50)
        {
        $GeneTargetLength[$firstgenenb]+=$length;
        }
        
        
    }
    else
    {
        $genelist="$GeneName[$firstgenenb]-$firstgenenb";
        $GeneTargeti[$firstgenenb]++;
        $GeneTargetType[$firstgenenb]++;
        if ($length>50)
        {
            $GeneTargetLength[$firstgenenb]+=$length;
        }
    }
    if ( $PosToGene{$debut} ne $PosToGene{($debut+$length-1) % $genomeSize})
    {
        if ($debut+$length-1<=$genomeSize  )
        {
            
            
            my $lastgenenb;
            @c=split(/_/,$PosToGene{($debut+$length-1)});
            $size=@c;
            $lastgenenb=$c[0];
            if ($size>1)
            {
                $finaladd=",$GeneName[$c[1]]-$c[1]";
                $GeneTargeti[$c[1]]++;
                $GeneTargetType[$c[1]]++;
                if ($length>50)
                {
                    $GeneTargetLength[$c[1]]+=$length;
                }

            }
            
            for (my $u=$firstgenenb+1;$u<=$lastgenenb;$u++)
            {
                $genelist="$genelist,$GeneName[$u]-$u";
                $GeneTargeti[$u]++;
                $GeneTargetType[$u]++;
                if ($length>50)
                {
                    $GeneTargetLength[$u]+=$length;
                }
            }
            $genelist="$genelist$finaladd";
            
            
        }
        else
        {
            
            
            my $lastgenenb;
            @c=split(/_/,$PosToGene{($debut+$length-1) % $genomeSize});
            $size=@c;
            $lastgenenb=$c[0];
            if ($lastgenenb==$TotalGeneNb)
            {
                $lastgenenb=-1;
            }
            if ($size>1)
            {
                $finaladd=",$GeneName[$c[1]]-$c[1]";
                $GeneTargeti[$c[1]]++;
                $GeneTargetType[$c[1]]++;
                if ($length>50)
                {
                    $GeneTargetLength[$c[1]]+=$length;
                }
            }
            for (my $u=$firstgenenb+1;$u<$TotalGeneNb;$u++)
            {
                $genelist="$genelist,$GeneName[$u]-$u";
                $GeneTargeti[$u]++;
                $GeneTargetType[$u]++;
                if ($length>50)
                {
                    $GeneTargetLength[$u]+=$length;
                }
            }
            for (my $u=0;$u<=$lastgenenb;$u++)
            {
                $genelist="$genelist,$GeneName[$u]-$u";
                $GeneTargeti[$u]++;
                $GeneTargetType[$u]++;
                if ($length>50)
                {
                    $GeneTargetLength[$u]+=$length;
                }
            }
            
            $genelist="$genelist$finaladd";
        }
    }
    
    @$GeneTargeti_ref=@GeneTargeti;
    @$GeneTargetType_ref=@GeneTargetType;
    @$GeneTargetLength_ref=@GeneTargetLength;

    
    return($genelist);
    
    
    
}


#################################################################################################
########################## Actions starts here  #################################################
#################################################################################################



############################################################################
#########################Read the input file################################
############################################################################

print "Read EventByStrainsSorted.txt\n";
my %StrainEventRef=();
my $genelist;
my %UniqueEventList=();
my %StrainList=();



open STRAINTOEVENT, ">$outputfolder/GenesGstat_SimulsNS.txt" or die "can not create $outputfolder/GenesGstat_SimulsNS.txt file!\n";
close STRAINTOEVENT;



open STRAINTOEVENT, "<$inputfile" or die "can not find $inputfile!\n";
while(my $line = <STRAINTOEVENT>)
    {
        chomp $line;
        @b=split(/\t/,$line);
        
        $StrainEventRef{"$b[0]_$b[1]_$b[2]§$b[3]"}="$b[4],$b[7],$b[5],$b[6]";
        
        $UniqueEventList{$b[3]}=$line;
        $StrainList{$b[2]}=1;
    }
close 	STRAINTOEVENT;


    my %MutationList;
    my  %MutPosReal;
    my $event;
    my @GeneTargeted;
    my @GeneTargeted_NS;
    my @GeneTargeted_O;
    my @GeneTargeted_I;
    my @GeneTargeted_IS;
    my @GeneTargeted_Sid;
    my @GeneTargeted_Lid;
    my @GeneTargeted_Length;
    my @GeneTargeted_LongDel;
    my @GeneTargeted_LongDupli;
    my @GeneTargeted_S;
    my @GeneExpectedStatsNSyn=(0) x $TotalGeneNb;
    my @GeneGStatsNS=(0) x $TotalGeneNb;
    
    my %UniqueGeneList=();
    my %completeListEvent=();
    

    ############################################################################
    #########Combine mutations into genes and randomise posisitons (i>0)########
    ############################################################################
    print "Start Simulations\n";
    for ($i =0;$i<=$nbrepetsimuls;$i++)
    {
       
        $nb=0;
        %StrainEventTemp=();
        %EventToGeneTemp=();
        %MutationList=();
        %MutPosReal=();
        %UniqueGeneList=();
        %completeListEvent=("NS",0,"I",0,"IS",0,"Indel",0,"Lindels",0,"LengthLIndels",0,"Longdeletions",0,"Longduplication",0,"S",0,"O",0);
        
       

        @GeneTargeted=(0) x $TotalGeneNb;
        @GeneTargeted_NS=(0) x $TotalGeneNb;
        @GeneTargeted_I=(0) x $TotalGeneNb;
        @GeneTargeted_IS=(0) x $TotalGeneNb;
        @GeneTargeted_Lid=(0) x $TotalGeneNb;
        @GeneTargeted_Sid=(0) x $TotalGeneNb;
        @GeneTargeted_Length=(0) x $TotalGeneNb;
        @GeneTargeted_S=(0) x $TotalGeneNb;
        @GeneTargeted_O=(0) x $TotalGeneNb;
        @GeneTargeted_LongDupli=(0) x $TotalGeneNb;
        @GeneTargeted_LongDel=(0) x $TotalGeneNb;

        ####################################
        ####For all the events##############
        ####################################
        foreach $keya (keys %StrainEventRef)
        {
            @b=split(/,/,$StrainEventRef{$keya});
            @c=split(/§/,$keya);
            my @pp=split(/_/,$c[0]);#pop
            my $mem= $b[0];
            #randomise the position for simulations higher than 0
            # Keep the same within population variability: a shared mutaiton within the population has the same new position
           
            #_________________________________________________________________________________________#
            #If the mutations already exists in the same population attritbute the same random position
            #_________________________________________________________________________________________#
            if (exists $MutPosReal{"$pp[0],$b[0],$b[1],$c[1]"})
            {
                $b[0]=$MutPosReal{"$pp[0],$b[0],$b[1],$c[1]"};
            }
            else
            {   #_________________________________________#
                #Counts mutations per type in the sample
                #_________________________________________#
                if ($b[2]eq"Mutation")
                    {
                    if ($b[3]eq"N")
                    {
                        $completeListEvent{"NS"}=$completeListEvent{"NS"}+1;
                    }
                    if ($b[3]eq"I")
                    {
                        $completeListEvent{"I"}=$completeListEvent{"I"}+1;
                    }
                    if ($b[3]eq"S")
                        {
                            $completeListEvent{"S"}=$completeListEvent{"S"}+1;
                        }
                    if ($b[3]eq"O")
                        {
                        $completeListEvent{"O"}=$completeListEvent{"O"}+1;
                        }
                    }
                if ($b[2]eq"IS_Insertion")
                    {
                    $completeListEvent{"IS"}=$completeListEvent{"IS"}+1;
                    }
                
                if ($b[2]=~/Deletion/||  $b[2]eq"Insertion" ||$b[2]=~/Duplication/ ||$b[2]=~/Substitution/)
                    {
                    $b[3]=~/(\d+)_bp/;
                   
                    if ($1>50)
                        {
                        
                        $completeListEvent{"Lindels"}=$completeListEvent{"Lindels"}+1;
                        $completeListEvent{"LengthLIndels"}=$completeListEvent{"LengthLIndels"}+$1;
                        
                        
                        if ($b[2]=~/Deletion/)
                        {$completeListEvent{"Longdeletions"}=$completeListEvent{"Longdeletions"}+1;}
                        else
                            {$completeListEvent{"Longduplication"}=$completeListEvent{"Longduplication"}+1;}
                        
                        }else
                        {
                        $completeListEvent{"Indel"}=$completeListEvent{"Indel"}+1;
                        }
                    }

                
                    
                #________________________________________________________#
                #for simulated dataset draw new position for the mutations
                #________________________________________________________#

                if ($i>0)
                {
                #find a random new position in the genome
                my $rand=ceil(rand()*$genomeSize);
                    while (exists($MaskedBases{$rand}))
                    {
                    $rand=ceil(rand()*$genomeSize)-1;
                    }
                
                #For points mutations mutations are resampled conserving their type
                #for instance, a coding mutation (N or S) is assigned a new possition in the coding genome
                if ($b[3] eq "N" || $b[3] eq "S")
                    {
                        while(!exists $PosToCodon{$rand} || exists($MaskedBases{$rand}))
                        {
                            $rand=ceil(rand()*$genomeSize);
                        }
                    }
                    # an intergenic mutation is assigned a random position among the non-codign part of the genome.
                if ($b[3] eq "I" || $b[3] eq "O")
                    {
                        my $rand=$NoncodingtoGenomePos{ceil(rand()*$nbnoncoding)-1};
                    }
                    
                
                $MutPosReal{"$pp[0],$b[0],$b[1],$c[1]"}=$rand;
                $b[0]=$rand;
                }
                else
                {
                 $MutPosReal{"$pp[0],$b[0],$b[1],$c[1]"}=$b[0];
                }
                
                #________________________________________________________#
                #        For each mutation find the affected genes
                #________________________________________________________#

                
                if ($b[2]eq"Mutation")
                {
                    if ($b[3]eq"N")
                    {
                        
                        $UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_NS,\@GeneTargeted_Length);
                    }
                    if ($b[3]eq"I")
                    {
                        $UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_I,\@GeneTargeted_Length);
                    }
                    if ($b[3]eq"O")
                    {
                        $UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_O,\@GeneTargeted_Length);
                    }

                    if ($b[3]eq"S")
                    {
                        my $firstgenenb=$PosToGene{$b[0]};
                        $GeneTargeted_S[$firstgenenb]=$GeneTargeted_S[$firstgenenb]+1;
                        
                    }
                    
                }
                if ($b[2]eq"IS_Insertion")
                {
                    $UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_IS,\@GeneTargeted_Length);
                }
                
                if ($b[2]=~/Deletion/||  $b[2]eq"Insertion" ||$b[2]=~/Duplication/)
                {
                    $b[3]=~/(\d+)_bp/;
                    
                    if ($1>50)
                    {
                        if ($b[2]=~/Deletion/)
                        {$UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_LongDel,\@GeneTargeted_Length);}
                        else
                        {
                        $UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_LongDupli,\@GeneTargeted_Length);

                        }
                    }else
                    {
                        $UniqueGeneList{"$pp[0],$b[0],$b[1],$c[1]"}= GetGenes3($b[0],$b[1],$genomeSize,$TotalGeneNb,\@GeneName,\@GeneTargeted,\@GeneTargeted_Sid,\@GeneTargeted_Length);
                    }
                }

                

                
            }
                
            #________________________________________________________#
            #    Assign a number to  mutations  (not already found)
            #________________________________________________________#
            
            if (!exists($MutationList{"$b[0],$b[1],$b[2],$c[1]"}))
            
            {
                $MutationList{"$b[0],$b[1],$b[2],$c[1]"}=$nb;
                $event=$nb;
                $EventToGeneTemp{$nb}=$genelist;
                $nb++;
            }
            else
            {
                $event=$MutationList{"$b[0],$b[1],$b[2],$c[1]"};
            }
            
            $StrainEventTemp{"$c[0]§$event"}="$b[0],$b[1],$b[2]\n";
        }
        
        
        
        ##################################################################
        ####For real data, compute Gscore per gene and print##############
        ##################################################################

        if ($i==0)
        {
        foreach $keya (keys %completeListEvent)
        {
            print "$keya\t$completeListEvent{$keya}\n";
        }
            

        @GeneExpectedStatsNSyn=(0) x $TotalGeneNb;
        my $expected=0;
        for ($k=0;$k<$TotalGeneNb;$k++)
        {
            $GeneExpectedStatsNSyn[$k]=$completeListEvent{"NS"}* $GeneLength[$k]/$totalcoding;
        }
            
        @GeneGStatsNS=(0) x $TotalGeneNb;
        open STRAINTOEVENT, ">$outputfolder/GenesGstat.txt";
            
        print STRAINTOEVENT "Gene_Name\tGene_Order\tStart_position\tCoding_Length\tObserved_nonsynonymous_mutation\tExpected_number_of_Non_Synonymous_mutations\tG_score\tSynonymous_mutation\tIntergenic_point_mutation\tPoint_mutation_in_pseudogene_or_noncoding_gene\tIS_insertion\tShort_indel\tLarge_deletion\tLong_duplication\n";
            #print STRAINTOEVENT "Total\t$TotalGeneNb\t0\t$totalcoding\t$NoncodingSize\tnbhits\texpectedNS\t", $completeListEvent{"NS"},"\tGscore(NS)\t", $completeListEvent{"I"},"\t",$completeListEvent{"IS"},"\t",$completeListEvent{"Indel"},"\t",$completeListEvent{"Lindels"},"\t",$completeListEvent{"LengthLIndels"},"\t",$completeListEvent{"Longdeletions"},"\t",$completeListEvent{"Longduplication"},"\t",$completeListEvent{"S"},"\t",$completeListEvent{"O"},"\n";
         
        for ($k=0;$k<$TotalGeneNb;$k++)
            {
                $GeneGStatsNS[$k]=0;
                
                if ($GeneTargeted_NS[$k]>0 and $GeneExpectedStatsNSyn[$k]>0)
                {
                $GeneGStatsNS[$k]=2*$GeneTargeted_NS[$k]*log($GeneTargeted_NS[$k]/$GeneExpectedStatsNSyn[$k]);
                }
                print STRAINTOEVENT "$GeneName[$k]\t$k\t$GenePosition[$k]\t$GeneLength[$k]\t$GeneTargeted_NS[$k]\t$GeneExpectedStatsNSyn[$k]\t$GeneGStatsNS[$k]\t$GeneTargeted_S[$k]\t$GeneTargeted_I[$k]\t$GeneTargeted_O[$k]\t$GeneTargeted_IS[$k]\t$GeneTargeted_Sid[$k]\t$GeneTargeted_LongDel[$k]\t$GeneTargeted_LongDupli[$k]\n";
            }
        close STRAINTOEVENT;

  
            
        open STRAINTOEVENT, ">$outputfolder/GenesObservedHits.txt";
        print STRAINTOEVENT "@GeneName\n";
        print STRAINTOEVENT "@GeneExpectedStatsNSyn\n";
        close STRAINTOEVENT;
    
        }
        
        
        
        ####################################################################################################
        ####For real and simulated datasets store the numbers of Non-synonymous hits per gene and Gscore####
        ####################################################################################################
        

        @GeneGStatsNS=(0) x $TotalGeneNb;
        #compute G afficher Obsersé ôur  chaque gene et premeire ligne les expeceted;
        #un ficher avec les G scores
        open STRAINTOEVENT, ">>$outputfolder/GenesObservedHits.txt";
        print STRAINTOEVENT "@GeneTargeted_NS\n";
        close STRAINTOEVENT;
        
        for ($k=0;$k<$TotalGeneNb;$k++)
        {
            $GeneGStatsNS[$k]=0;
            
            if ($GeneTargeted_NS[$k]>0 and $GeneExpectedStatsNSyn[$k]>0)
            {
            $GeneGStatsNS[$k]=2*$GeneTargeted_NS[$k]*log($GeneTargeted_NS[$k]/$GeneExpectedStatsNSyn[$k]);
            }
        
            
        }
       
        open STRAINTOEVENT, ">>$outputfolder/GenesGstat_SimulsNS.txt";
        print STRAINTOEVENT "@GeneGStatsNS\n";
        close STRAINTOEVENT;
        
        print "$i\n";
    }
    



