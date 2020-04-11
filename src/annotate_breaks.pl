#!/usr/bin/perl
# HiTEA
# This program annotates the breaks with all related read features

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use open qw(:std :utf8);
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0);
use src::Utilities;
require "src/vars.pl";
our (%redb);
our (%chrs);
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
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
sub help{
  my $j = shift;
  if($j){
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
   exit 1;
  }
}
$wd =~ s/\/$//;
my $baseoutprefix =$outprefix;
$outprefix = $wd."/".$outprefix if($wd ne "");

if($help or $in eq "" or $ram eq "" or $index eq "" or $rand eq ""){
  print " One or more required inputs are not recognized\n";
  help(1);
}

## Dependancies
my %cov = %{retrieve($outprefix.".coverage.ph")};
print " ERROR: generateRanges.R does not exist in the same directory!\n" unless -e "src/generateRanges.R";

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
%transposons = src::Utilities::get_fasta($index);
print " transposons being considered in the analyses: " , join(",",keys %transposons),"\n";
my @TEELEMENTS= keys %transposons;
push @TEELEMENTS,"unmap";

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
  $j = src::Utilities::unique($j," ","false");
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
   my $ss = int( ($temp[1]/1000)+0.5);
   for(my $e=-1;$e<=1;$e++){
     if($cov{$temp[0]}{$ss+$e}){
       $val += $cov{$temp[0]}{$ss+$e};
       $valnum++;
     }
   }
   $val = src::Utilities::round($val/(10*$valnum),2) if($valnum>0);   ## readlength 100 vs bin size 1000, hence division by 10
   $Locs{$temp[0]."\t".$temp[1]."\t".$temp[1]."\t".$temp[3]."\t.\t".$temp[5]."\t".$val}++;
}
close(FRI);
while(my($i,$j) = each %randLocs){
  $j =~ s/ $//;
  $j =~ s/  / /g;
  $j = src::Utilities::unique($j," ","false");
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
print  " annotating breakpoints using bamfile..\n";   
my $clust = count_read_support_from_RAM($ram,\%clusters);
%clusters =%{$clust};
store \%clusters, $outprefix.'.ClustObj_supportRAMcount.ph'; #Save 

#my $oufile=$outprefix.'.RAM.clusters'; #Save 
#open FOO,">$oufile" or die $!;
#print FOO Dumper %clusters;
#close(FOO);

## generate genomic range objects for the supporting reads
undef %clusters;

foreach my $te (@TEELEMENTS){
  my $te1 = $te;
  $te=~ s/\/\S*$//g;
  my $file = $outprefix."_".$te."_supportingreads.txt";
  my $tegr = $wd."/".$baseoutprefix."_".$te;
  print qq[ creating GRange..  ];    ##(1)chr,(2)start,(3)end,(4)id,(5)strand,(6)evi,(7)clip, (8)refMapqQ, (9)TEMapScore, (10) TE_strand (11) 
  system( qq[Rscript src/generateRanges.R $file $tegr] ) == 0 or die qq[ Cound not create GRange Object for $te1\n];
  print qq[ cleaning up..  ];    
  system( qq[rm $file]) == 0 or die qq[ Error in deleting intermetidate files for $te1 \n];            
  print qq[ Done\n];
}

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
  foreach my $te (@TEELEMENTS){
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

    ## Add deduplication filter on the sorted bam file
    
    my $id = $temp[2].$temp[3].$temp[5].$temp[6].$temp[7];
    my ($evi,$side,$clip,$both,$te,$mm,$sc,$te_str,$end) = $_ =~ m/OP:Z:evi=(\S+?);side=(\S+?);clip=(\d+?);both=(\S+?);TE=(\S+?);mm=(\S+?);sc=(\S+?);str=(\S+?);end=(\d+);$/;
    my ($class) = $_ =~ m/Yt:Z:(\S+)/;
    my $strand = "+";
    $strand = "-" if($temp[1] & 0x10); 
    #print $_,"\n" if(!defined $clip);

    my $pnts = &findOverlaps($clip,$cliplocs{$temp[2]}, 2000);
    my @rr = split(",",$te);
    my @tescores = split(",",$sc);
    
    ## write down events within 2kb of the random or observed breaks
    if($te eq "*" and $pnts ne ""){ 
      my $fh = $filehandles{unmap};
      print $fh "$temp[2]\t$temp[3]\t$end\t$id\t$strand\t$evi\t$clip\t$temp[4]\t0\t*\t*\t*\n";  #chr,start,end,id,strand,evi,cliploc,refMapqQ,TEMapScore,TE_strand, TE_pos, TE
    }elsif($te ne "*" and $pnts ne ""){ 
      if(scalar(@rr)==0){
       print "   error in reading TE mapping information, $te\n";
      }
      my @testrands= split(",",$te_str);
      foreach my $i (0..$#rr){
        my ($ite,$pos) = $rr[$i] =~ m/(\S+):(\d+)/;
        my $iteo = $ite;
        $ite =~ s/~.*$//; $ite =~ s/\r//g; $ite =~ s/\t//g; 
        if($filehandles{$ite}){
          my $fh = $filehandles{$ite};
          print $fh "$temp[2]\t$temp[3]\t$end\t$id\t$strand\t$evi\t$clip\t$temp[4]\t$tescores[$i]\t$testrands[$i]\t$pos\t$iteo\n";  #chr,start,end,id,strand,evi,cliploc,refMapqQ,TEMapScore,TE_strand, TE_pos, TE
        }
      }
    }elsif($te ne "*" and $pnts eq ""){ # check random locations 
      my $alt_pnts = &findOverlaps($clip,$randLocs{$temp[2]}, 2000);
      if($alt_pnts ne ""){
        if(scalar(@rr)==0){
          print "   error in reading TE mapping information, $te\n";
        }
        my @testrands= split(",",$te_str);
        foreach my $i (0..$#rr){
          my ($ite,$pos) = $rr[$i] =~ m/(\S+):(\d+)/;
          my $iteo = $ite;
          $ite =~ s/~.*$//; $ite =~ s/\r//g; $ite =~ s/\t//g;
          if($filehandles{$ite}){
            my $fh = $filehandles{$ite};
            print $fh "$temp[2]\t$temp[3]\t$end\t$id\t$strand\t$evi\t$clip\t$temp[4]\t$tescores[$i]\t$testrands[$i]\t$pos\t$iteo\n";  #chr,start,end,id,strand,evi,cliploc,refMapqQ,TEMapScore,TE_strand, TE_pos, TE
          }
        }
      }
    }
    next if($pnts eq "");
     
    my @pnts = split(" ",$pnts);
    foreach my $cluster_pos (@pnts){
      my $loc = $temp[2]."\t".$cluster_pos."\t".$side;
      my $dist = abs($cluster_pos-$clip);
      ## case-1
      if($clusters{$loc} and $dist <= $gap){
        foreach my $i (0..$#rr){
          if($evi=~ m/DE,TP/ and $tescores[$i] >= $algnscore){ ## reads without any ambiguity about clip sequence
            my ($ite,$pos) = $rr[$i] =~ m/(\S+):(\d+)/;
            $ite =~ s/^\S*?~//;$ite =~ s/\r//g; $ite =~ s/\t//g; 
            $clusters{$loc}{pos1}{$ite}{num} .=$id.",";
            $clusters{$loc}{pos1}{$ite}{TESc} .=$tescores[$i].",";
            $clusters{$loc}{pos1}{$ite}{TEClust} .=$pos.",";
          }
        }
        if($evi=~ m/DE,TP/){
          ($clusters{$loc}{pos1}{dist2re}) = $evi=~ m/^\S*,\S*,\S*,(\S*)$/;
          $clusters{$loc}{pos1}{seqstring} = $temp[9].",".$temp[5].",".$strand.",".$temp[3];  #seq, cigar, strand, readstart
          $clusters{$loc}{pos1}{AllClip} .=$id.";";
          $clusters{$loc}{pos1}{both} .= $both."," if(!defined $clusters{$loc}{pos1}{both} or $both >0 );          
          $clusters{$loc}{pos1}{mapq} = $temp[4] if(!defined $clusters{$loc}{pos1}{mapq} or $clusters{$loc}{pos1}{mapq}<$temp[4]);
          my $s = src::Utilities::get_softclipseq($temp[9],$temp[5]);
          if(exists $clusters{$loc}{pos1}{clipseqs} and length($s)> length($clusters{$loc}{pos1}{clipseqs})){
            $clusters{$loc}{pos1}{clipseqs}= $s;
           }elsif(!exists $clusters{$loc}{pos1}{clipseqs}){
            $clusters{$loc}{pos1}{clipseqs}= $s;
          }
        }
        if( exists($clusters{$loc}{pos1}{startspan}) and $clusters{$loc}{pos1}{startspan} > $temp[3] ){ ## lowest start coord       
          $clusters{$loc}{pos1}{startspan} = $temp[3];
        }else{
          $clusters{$loc}{pos1}{startspan} = $temp[3];
        }
        if( exists($clusters{$loc}{pos1}{endspan}) and $clusters{$loc}{pos1}{endspan} < $end ){ ## highest end coord             
           $clusters{$loc}{pos1}{endspan} = $end;
        }else{
           $clusters{$loc}{pos1}{endspan} = $end;
        }
      }elsif($clusters{$loc} and $dist < 50 and $evi=~ m/DE,TP/ ){
        $clusters{$loc}{pos1}{fuzzyclipped} .=$id.","; 
      }
      ## case-2          
      my $alt_side = "a";
      $alt_side = "b" if($side eq "a");
      my $loc1 = $temp[2]."\t".$cluster_pos."\t".$alt_side;
      if($clusters{$loc1} and $cluster_pos>=$temp[3] and $cluster_pos<=$end and $evi =~ /DE,TP/){
        foreach my $i (0..$#rr){
          my $ite = "*";
          my $pos=0;
          ($ite, $pos) = $rr[$i] =~ m/(\S+):(\d+)/ if($rr[$i] ne "*");
          $ite =~ s/^\S*?~//; $ite =~ s/\r//g; $ite =~ s/\t//g; 
          my $l = $temp[2]."_".$clip."_".$side."_".$ite."|".$id;   
          $clusters{$loc1}{smate}{smate} .=$l.",";  ## use this info to tally/write 2nd mate when there is no mapping available for 2nd mate
          
          if($ite ne "*"){
            $clusters{$loc1}{smate}{$ite}{num} .=$id.",";
            $clusters{$loc1}{smate}{$ite}{TEClust} .=$pos.",";
            $clusters{$loc1}{smate}{$ite}{TESc} .=$tescores[$i].",";
          }
        }
        if($te eq "*"){
          $clusters{$loc1}{smate}{clipseqs} .= src::Utilities::get_softclipseq($temp[9],$temp[5])."|".$clip.",";  ## capture the supplimentary clip sequences for remapping in the event of PolyA,Unmapped clusters
        }
        ($clusters{$loc1}{smate}{dist2re}) = $evi=~ m/^\S*,\S*,\S*,(\S*)$/;          
        $clusters{$loc1}{smate}{AllClip} .=$id.";";
        $clusters{$loc1}{smate}{seqstring} = $temp[9].",".$temp[5].",".$strand.",".$temp[3];  #seq, cigar, strand, readstart  
        $clusters{$loc1}{smate}{both} .= $both."," if(!defined $clusters{$loc1}{smate}{both} or $both >0 );          
        $clusters{$loc1}{smate}{mapq} = $temp[4] if(!defined $clusters{$loc1}{smate}{mapq} or $clusters{$loc1}{smate}{mapq}<$temp[4]);
      }

    }
  }
  close(IN);
  foreach my $fh(keys %filehandles){close($fh);}
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
      $pnts = src::Utilities::unique($pnts," ","false");
      return($pnts);
    }
  }
  
  if(exists $a->[$l] and abs($a->[$l] -$x)<=$win){ $pnts .= $a->[$l]." ";}
  if(exists $a->[$l-1] and abs($a->[$l-1] -$x)<=$win){ $pnts .= $a->[$l-1]." ";}
  $pnts = src::Utilities::unique($pnts," ","false");
  return($pnts);
}