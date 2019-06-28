#!/usr/bin/perl
# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard

# This program merges clusters on opposite strands and identifies Duplicated Target Site (TSD) from the supported reads

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use open qw(:std :utf8);
use Utilities;
BEGIN { our $start_run = time(); }

my $in="";
my $ram = "";
my $index ="";
my $rand ="";
my $ncores=4;
my $outprefix="project";
my $algnscore = 20;
my $gap =2;
my $help=0;
my $wd = "";
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'ram=s'            => \$ram,
  'index=s'          => \$index,
  'rand:s'           => \$rand,
  'outprefix=s'      => \$outprefix,
  'algnscore:s'      => \$algnscore,
  'gap=s'            => \$gap,
  'wd:s'             => \$wd,
  'ncores:s'         => \$ncores,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl annotate_breaks.pl -in [FILE_PATH] -ram [FILE_PATH] -outprefix [STRING] -index [TE_MergedIndex] -rand [FILE] -ncores [INT] -algnscore [INT] -gap [INT] -wd [DIR_PATH]\n\n";
  print "This program annotates finalized breaks using the bam file\n\n";
  print "Options:\n\n";
  print "***required:\n";
  print "  -in                    Clusters file in bed format \n";
  print "  -ram                   RAM entry file in sam/bam format \n";
  print "  -index                 TE-consensus index build using BWA [Make sure what you are providing] \n";
  print "  -rand                  random locations this file will be rewritten with addition of a last column\n";
  print "***optional:\n";
  print "  -outprefix             Outputprefix for generating 2 output files [default: project] \n";
  print "  -wd                    Working directory [default: ~]\n";
  print "  -algnscore             alignment score threshold [default: 20] \n";
  print "  -gap                   gap used to merge the breakpoints [default: 2] \n";
  print "  -ncores                Number of threads while reading bam input [default: 4]\n";
  print "  -help|-h               Display usage information.\n\n\n";
  exit 0;
}

#--------------------------------------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------------------------------------
$wd =~ s/\/$//;
my $baseoutprefix =$outprefix;
if($wd ne ""){
  $outprefix = $wd."/".$outprefix;
}
my %chrs =(
    "chr1" => '1',
    "chr2" => '1',
    "chr3" => '1',
    "chr4" => '1',
    "chr5" => '1',
    "chr6" => '1',
    "chr7" => '1',
    "chr8" => '1',
    "chr9" => '1',
    "chr10" => '1',
    "chr11" => '1',
    "chr12" => '1',
    "chr13" => '1',
    "chr14" => '1',
    "chr15" => '1',
    "chr16" => '1',
    "chr17" => '1',
    "chr18" => '1',
    "chr19" => '1',
    "chr20" => '1',
    "chr21" => '1',
    "chr22" => '1',
    "chrX" => '1',
    "chrY" => '1',
);

## Dependancies
my %cov = %{retrieve($outprefix.".coverage.ph")};
print " ERROR: generateRanges.R does not exist in the same directory!\n" unless -e "generateRanges.R";

#--------------------------------------------------------------------------------------------------------------
# I/O
#-------------------------------------------------------------------------------------------------------------
my %transposons;
my %clusters;
my %cliplocs;
my %Locs;
my %randLocs;

my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "[annotate_breaks] START:\t $run_time seconds\n";
print " Command: perl annotate_breaks.pl -in $in -ram $ram -index $index -rand $rand -outprefix -ncores $ncores $outprefix -wd $wd \n";

## get transposons
open (TRANS,"cat $index | grep \">\" | cut -f1 | sort | uniq | sed -e \'s/\>//g\' | ") or die $!;
%transposons = map { chomp; s/\r//; s/\t//; $_ => 1 } <TRANS>;
close(TRANS);
print " transposons being considered in the analyses: " , join(",",keys %transposons),"\n";

## retrieve cluster from the input file
print " reading observed breaks and generating binary\n";
%clusters = %{retrieve($in)};
foreach my $loc (sort keys %clusters){
  my @temp = split("\t",$loc);
  $cliplocs{$temp[0]} .= $temp[1]." ";
} 
while(my($i,$j) = each %cliplocs){
  $j =~ s/ $//;
  $j =~ s/  / /g;
  $j = unique($j," ","false");
  my @temp= split(" ",$j);  
  @temp = sort{$a <=> $b} @temp;
  $cliplocs{$i} = \@temp;
}


## read random location file and add coverage column
print " reading random breaks and generating binary\n";
open FRI, $rand or die "Can't open already genrated random location file $rand\n";
while(<FRI>){
   chomp;
   s/\r//g;
   next if(/^(\#)/);
   next if(/^(\@)/); 
   my @temp = split(/\t/);
   next if(!$chrs{$temp[0]});
   $randLocs{$temp[0]} .= $temp[1]." ";    
   my $valnum=0;
   my $val=0;
   my $ss = int($temp[1]/1000);
   for(my $e=-1;$e<=1;$e++){
     if($cov{$temp[0]}{$ss+$e}){
       $val += $cov{$temp[0]}{$ss+$e};
       $valnum++;
     }
   }
   $val = round($val/(10*$valnum),2) if($valnum>0);   ## readlength 100 vs bin size 1000, hence division by 10
   $Locs{$temp[0]."\t".$temp[1]."\t".$temp[1]."\t".$temp[3]."\t.\t".$temp[5]."\t".$val}++;
}
close(FRI);
while(my($i,$j) = each %randLocs){
  $j =~ s/ $//;
  $j =~ s/  / /g;
  $j = unique($j," ","false");
  my @temp= split(" ",$j);  
  @temp = sort{$a <=> $b} @temp;
  $randLocs{$i} = \@temp;
}
open O, ">",$rand or die "Can't rewrite file $rand\n";
while(my ($i, $j) = each %Locs){
  print O $i,"\n";
}
close(O);


## annotating the clusters using RAM bam file
print  " annotating clusters using bamfile..\n";   
my $clust = count_read_support_from_RAM($ram,\%clusters);
%clusters =%{$clust};
store \%clusters, $outprefix.'.ClustObj_supportRAMcount.ph'; #Save 

$watch_run = time();
$run_time = $watch_run - $start_run;
print  "[annotate_breaks] END:\t $run_time seconds\n";


exit 0;

#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub count_read_support_from_RAM{
  my ($file,$clust) = @_;
  my %clusters = %{$clust};
  print  "  $file\n";
  if($file=~ /.bam/){
    open IN,"samtools view -@ $ncores $file|" or next "Can't open bamfile $file";
  }else{
    open IN,"cat $file|" or next "Can't open file $file";
  }
  
  ## write down supporting reads for each TE element in a separate text file => sort => use supporting read counts for modeling
  my %filehandles;
  foreach my $te (keys %transposons) {
    my $te1 = $te;
    $te1=~ s/\/\S*$//g;
    open( my $fh, ">", $outprefix."_".$te1."_supportingreads.txt") or die "Error: Coudln't open file for writing: $!";
    $filehandles{$te} = $fh;
  }

  while(<IN>){
       next if(/^(\#)/);
       next if(/^(\@)/); 
       chomp;
       s/\r//;  
       my @temp = split(/\t/);
       next if(!exists($cliplocs{$temp[2]})); 
       #(0)1226    (1)145     (2)chr17   (3)81936376        (4)0       (5)20S75M  (6)chr1    (7)245532241       (8)0       (9)GTGTGTGCTGCTTGAGCCCTCCACCCCACCCGGCTAATTTTTGGTATCTTTAGTGGAGACGGGGTTTCACTGTGTTAGCCAGGATGGTCTCGATC (10)//>00000011>B/0>///?////////EE/101B10/A/1B11D1B1F1A10BABFA0000A1HHGBB00FGE1BB11A1FAB11111A>>111 NM:i:3  MD:Z:5A21T6A40  MC:Z:88M7S      AS:i:60 XS:i:60 Yt:Z:MM OP:Z:evi=DE,TP,20,100;side=a;clip=81936376;TE=*;mm=NA;sc=0;str=*;cl=-;end=81936451 
       my $id = $temp[2].$temp[3].$temp[5].$temp[6].$temp[7];
       my ($evi,$side,$clip,$te,$mm,$sc,$te_str,$cl,$end) = $_ =~ m/OP:Z:evi=(\S+?);side=(\S+?);clip=(\d+?);TE=(\S+?);mm=(\S+?);sc=(\d+?);str=(\S+?);cl=(\S+?);end=(\d+)$/;
       my ($class) = $_ =~ m/Yt:Z:(\S+)/;
       my $strand = "+";
       $strand = "-" if($temp[1] & 0x10); 
       my $pnts = &findOverlaps($clip,$cliplocs{$temp[2]}, 2000);
          ## write down events within 2kb of the random or observed breaks
       if($te ne "*" and $pnts ne ""){ 
         my @r = split(";",$cl);
         if(scalar(@r)==0){
          print "   error in reading TE mapping information, $cl\n";
         }
         foreach my $temap (@r){
           my ($ite,$pos) = $temap =~ m/(\w+):(\d+)/;
           $ite =~ s/\r//g; $ite =~ s/\t//g; 
           if($filehandles{$ite}){
             my $fh = $filehandles{$ite};
             print $fh "$temp[2]\t$temp[3]\t$end\t$id\t$strand\t$evi\t$clip\t$temp[4]\t$sc\t$te_str\t$pos\n";  #chr,start,end,id,strand,evi,cliploc,refMapqQ,TEMapScore,TE_strand, TE_pos
           }
         }
       }elsif($te ne "*" and $pnts eq ""){ # check random locations 
         my $alt_pnts = &findOverlaps($clip,$randLocs{$temp[2]}, 2000);
         if($alt_pnts ne ""){
           my @r = split(";",$cl);
           if(scalar(@r)==0){
             print "   error in reading TE mapping information, $cl\n";
           }foreach my $temap (@r){
             my ($ite,$pos) = $temap =~ m/(\w+):(\d+)/;
             $ite =~ s/\r//g; $ite =~ s/\t//g;
             if($filehandles{$ite}){
               my $fh = $filehandles{$ite};
               print $fh "$temp[2]\t$temp[3]\t$end\t$id\t$strand\t$evi\t$clip\t$temp[4]\t$sc\t$te_str\t$pos\n";  #chr,start,end,id,strand,evi,cliploc,refMapqQ,TEMapScore,TE_strand, TE_pos
             }
           }
          }
       }
       next if($pnts eq "");
        
       my @pnts = split(" ",$pnts);
       foreach my $cluster_pos (@pnts){
           my $loc = $temp[2]."\t".$cluster_pos."\t".$side;
           my $dist = abs($cluster_pos-$clip);
           if($clusters{$loc} and $dist <= $gap){
              if($evi=~ m/DE,TP/ and $sc >= $algnscore){ ## reads without any ambiguity about clip sequence
                 $clusters{$loc}{pos1}{$te} .=$id.",";
                 $clusters{$loc}{pos1}{TEClust} .=$cl.";";
              }elsif($evi=~ m/DE,TP/){
                 $clusters{$loc}{pos1}{unmapped} .=$id.",";
              }              
              if($evi=~ m/DE,TP/){
                 ($clusters{$loc}{pos1}{dist2re}) = $evi=~ m/^\S*,\S*,\S*,(\S*)$/;
                 $clusters{$loc}{pos1}{seqstring} = $temp[9].",".$temp[5].",".$strand.",".$temp[3];  #seq, cigar, strand, readstart
                 my $s = Utilities::get_softclipseq($temp[9],$temp[5]);
                 if(exists $clusters{$loc}{pos1}{clipseqs} and length($s)> length($clusters{$loc}{pos1}{clipseqs})){
                    $clusters{$loc}{pos1}{clipseqs}= $s;
                 }elsif(!exists $clusters{$loc}{pos1}{clipseqs}){
                    $clusters{$loc}{pos1}{clipseqs}= $s;
                 }
              }
              
              if( exists($clusters{$loc}{pos1}{startspan}) ){ ## lowest start coord       
                if($clusters{$loc}{pos1}{startspan} > $temp[3]) {
                   $clusters{$loc}{pos1}{startspan} = $temp[3];
                }
              }else{
                $clusters{$loc}{pos1}{startspan} = $temp[3];
              }
              if( exists($clusters{$loc}{pos1}{endspan}) ){ ## highest end coord             
                if($clusters{$loc}{pos1}{endspan} < $end) {
                  $clusters{$loc}{pos1}{endspan} = $end;
                }
              }else{
                $clusters{$loc}{pos1}{endspan} = $end;
              }
           }elsif($clusters{$loc} and $dist < 50 and $evi=~ m/DE,TP/ ){
              $clusters{$loc}{pos1}{fuzzyclipped} .=$id.","; 
           }
           my $alt_side = "a";
           $alt_side = "b" if($side eq "a");
           my $loc1 = $temp[2]."\t".$cluster_pos."\t".$alt_side;
           if($clusters{$loc1} and $cluster_pos>=$temp[3] and $cluster_pos<=$end and $evi =~ /DE,TP/){
              my $l = $temp[2]."_".$clip."_".$side."_".$te."|".$id;   
              $clusters{$loc1}{pos1}{smate} .=$l.",";  ## use this info to tally/write 2nd mate when there is no mapping available for 2nd mate
              $clusters{$loc1}{pos1}{smate_seqstring} = $temp[9].",".$temp[5].",".$strand.",".$temp[3];  #seq, cigar, strand, readstart
              if($te eq "*"){
                $clusters{$loc1}{pos1}{smate_clipseqs} .= Utilities::get_softclipseq($temp[9],$temp[5]).",";  ## capture the supplimentary clip sequences for remapping in the event of PolyA,Unmapped clusters
              }
              if($cl ne "-"){
                 $clusters{$loc1}{pos1}{smateClust} .=$cl.";";
              }
           }
       }
  }
  close(IN);
  foreach my $fh(keys %filehandles){close($fh);}
  
  
  ## generate genomic range objects for the supporting reads
  foreach my $te (keys %transposons){
    my $te1 = $te;
    $te=~ s/\/\S*$//g;
    my $file = $outprefix."_".$te."_supportingreads.txt";
    my $tegr = $wd."/".$baseoutprefix."_".$te;
    print qq[ creating GRange..  ];    ##(1)chr,(2)start,(3)end,(4)id,(5)strand,(6)evi,(7)clip, (8)refMapqQ, (9)TEMapScore, (10) TE_strand
    system( qq[Rscript generateRanges.R $file $tegr] ) == 0 or die qq[ Cound not create GRange Object for $te1\n];
    print qq[ cleaning up..  ];    
    system( qq[rm $file]) == 0 or die qq[ Error in deleting intermetidate files for $te1 \n];            
    print qq[ Done\n];
  }
  return(\%clusters);
}

sub findOverlaps {
  my ($x, $a, $win) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);          # lower, upper end of search interval
  my $i;                               # index of probe
  my $pnts="";

  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
      if(abs($a->[$i] -$x)<=$win){ $pnts .= $a->[$i]." ";}
      if(exists $a->[$i+1] and abs($a->[$i+1] -$x)<=$win){ $pnts .= $a->[$i+1]." ";}
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
      if(abs($a->[$i] -$x)<=$win){ $pnts .= $a->[$i]." ";}
      if(exists $a->[$i+1] and abs($a->[$i+1] -$x)<=$win){ $pnts .= $a->[$i+1]." ";}
    } 
    else {
      if(abs($a->[$i] -$x)<=$win){ $pnts .= $a->[$i]." ";}
      if(exists $a->[$i+1] and abs($a->[$i+1] -$x)<=$win){ $pnts .= $a->[$i+1]." ";}
      $pnts = unique($pnts," ","false");
      return($pnts);
    }
  }
  
  if(exists $a->[$l] and abs($a->[$l] -$x)<=$win){ $pnts .= $a->[$l]." ";}
  if(exists $a->[$l-1] and abs($a->[$l-1] -$x)<=$win){ $pnts .= $a->[$l-1]." ";}
  $pnts = unique($pnts," ","false");
  return($pnts);
}