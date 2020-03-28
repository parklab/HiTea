#!/usr/bin/perl

# HiTEA
# This script extacts break positions

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
use utf8;
use open qw(:std :utf8);
use src::Utilities;
require "src/vars.pl";
our (%redb);
our (%chrs);
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
BEGIN { our $start_run = time(); }

my $bam = "";
my $hd ="";
my $wd ="";
my $remap ="F";
my $TE_element="";
my $outprefix="project";
my $help=0;
Getopt::Long::GetOptions(
  'bam=s'            => \$bam,
  'hd:s'             => \$hd,
  'outprefix=s'      => \$outprefix,
  'wd:s'             => \$wd,
  'remap:s'          => \$remap,
  'te:s'             => \$TE_element,
  'help'             => \$help,
  'h'                => \$help,
  ) or die "Incorrect input! Use -h for usage.\n";

sub help{
  my $j = shift;
  if($j){
   print "\nUsage: perl find_breaks.pl -bam [FILE_PATH] -outprefix [STRING] -wd [DIR_PATH] -hd [FILE_PATH] -remap [T/F] -te [STRING] \n\n";
   print "This script generates a bam file (can be used for visualization) and breakpoint list file from mapped sam/bam\n\n";
   print "Options:\n";
   print "required**:\n";
   print "   -bam             sam/bam input file. Read name has information about other mates mapping.  \n";
   print "   -hd              Header file for output bam\n";
   print "optional:\n";
   print "   -wd              Working directory [default: ~]\n";
   print "   -outprefix       Outputprefix for generating the output files \n";
   print "   -remap           Whether to remap the unmapped clipped sequences to PolyMoprhic sequence assembly? (default=F)\n";
   print "   -te              TE for which polymoprhic rempping is being performed\n";
   print "   -help|-h         Display usage information.\n\n\n";
   print "Default outputs:\n";
   print "   (outprefix).clusters.bed\n";
   print "   (outprefix).RAM.bam\n";
   exit 1;
  }
}
$wd =~ s/\/$//;
$outprefix =~ s/\?//g;
$outprefix =~ s/\#//g;
$outprefix =~ s/\'//g;
$outprefix =~ s/\"//g;
$outprefix =~ s/\s*//g;
$outprefix = $wd."/".$outprefix if($wd ne "");
if($help or $bam eq "" or $hd eq ""){
  print " One or more inputs are not recognized**\n";
  help(1);
}
if($remap eq "T" and $TE_element eq ""){
  print " TE element not designated for remapping. Exiting! \n";
  help(1);
}

#--------------------------------------------------------------------------------------------------------------
# main body
#--------------------------------------------------------------------------------------------------------------
my $watch_run=0;
my $run_time=0;
$watch_run = time();
$run_time = $watch_run - our $start_run;
print "[extract_clusters=Start_Time]:\t $run_time seconds\n";

my $out1 = $outprefix.".temp.clusters.bed";
my $out2 = $outprefix.".temp.RAM.bam";
my $out3 =  $outprefix.".temp.remap.fq";
if($TE_element ne "" and $remap eq "T"){
  $out2 = $outprefix.".".$TE_element.".temp.RAM.bam";
  $out3 =  $outprefix.".".$TE_element.".temp.remap.fq.gz";
  open O3, "|gzip -c - > $out3" or die "can't create $out3";
}
open O1, ">> $out1" or die "can't create $out1";
open O2, "| samtools view -bS - > $out2" or die "can't create $out2";
if($hd ne ""){
  open(F, "$hd") or die $!;
  my @header = <F>;
  print O2 @header;
  close F;
  undef(@header);
}

## Get cluster coordinate from Repeat-mapped clipped reads
## Write RAM in a bed format
my %hqual;
my $FIN;
my $badline=0;
if($bam eq "-"){
   print "Input is pipein \n";
   $FIN="STDIN";
   get_breakpoints($FIN);  
}elsif($bam ne "-" and $bam =~ m/.bam$/){
  open($FIN,"samtools view -h $bam | ") or die "no bam: $bam";
  print " Input bam file: $bam\n";
  get_breakpoints($FIN);
}elsif($bam ne "-" and $bam =~ m/.sam$/){
  open($FIN,"$bam") or die "no bam: $bam";
  print " Input sam file: $bam\n";
  get_breakpoints($FIN);
}else{
  print " incorrrect input sam/bam file\n";
  exit 1;
}


## Write clusters into a bed format
my $cls=0;
while( my( $i, $j ) = each %hqual ) {
  $cls++;
  my@temp = split("\t",$i); #chr, start,side
  my $count = scalar(split(" ",$j));
  next if(!exists($chrs{$temp[0]}));  ## only write clusters on chr1-22/X/Y
  if($temp[2] eq "a"){
    print O1 "$temp[0]\t$temp[1]\t$temp[1]\t$count\t$temp[2]\t+\t$j\n";
  }elsif($temp[2] eq "b"){
    print O1 "$temp[0]\t$temp[1]\t$temp[1]\t$count\t$temp[2]\t-\t$j\n";
  }
}
close(O1);
close(O3) if($remap eq "T");
print " Total number of clip read positions in the file: $cls\n";
print " Total badlines in the file: $badline\n";
$watch_run = time();
$run_time = $watch_run - $start_run;
print "[extract_clusters=End_Time]:\t $run_time seconds\n";

exit 0;
#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub get_breakpoints {
 my ($file) = shift;
 my $i=0;
 my $oid ="-";
 my %res;

 while(<$file>){
    next if(/^(\#)/);
    next if(/^(\@)/); 
    chomp;
    s/\r//;  
    next if($_ eq "");
    my @temp = split(/\t/);    
    if($temp[1] ==4 and $temp[0]!~ /DE,TP/){
      next;
    }
    if($temp[1]==4 and $temp[0]=~ m/DE,TP/ and $remap eq "T"){
      print O3 "@","$temp[0]\n$temp[9]\n+\n$temp[10]\n";
      next;
    }elsif($temp[1]==4 and $temp[0]=~ m/DE,TP/ and $remap eq "F"){
      my @sam = split("\x{019}",$temp[0]);
      my $end = src::Utilities::get_read_mapping_span($sam[3],$sam[5]);
      $sam[-1] .= ";TE=*;mm=0;sc=0;str=*;end=".$end.";";
      print O2 join("\t",@sam),"\n";
      next;
    }
    
    my ($mismatch) = "NA";
    my $strand = "*";
    my ($score) = $_ =~ /\tAS:i:(\d*)/;
    if($_=~ /\tNM:i:(\d*)/){ $mismatch=$1;}
    $strand = "+";
    $strand = "-" if($temp[1] & 0x10);

    if($oid eq "-" or ($oid ne "-" and $oid eq $temp[0])){
      $oid =$temp[0] if($oid eq "-");
      if(!exists $res{$temp[2]} or $res{$temp[2]}{score} < $score){
         $res{$temp[2]}{score} = $score;
         $res{$temp[2]}{mm} = $mismatch;
         $res{$temp[2]}{str} = $strand;
         $res{$temp[2]}{start} = $temp[3];
         $res{$temp[2]}{cigar} = $temp[5];
      }
    }

    if($oid ne "-" and $oid ne $temp[0]){
      ## write
      my @sam = split("\x{019}",$oid);
      my $end = src::Utilities::get_read_mapping_span($sam[3],$sam[5]);
      my $read_strand="+";
      $read_strand="-" if($sam[1] & 0x10);       
      my $te;
      my $sc;
      my $mm;
      my $str; 
      my $maxscore;
      my ($cliploc) = $sam[-1] =~ m/;clip=(\d+);/;
      my ($sidedness) = $sam[-1] =~ m/;side=(\S+?);/;     
      my $id = $sam[2].$sam[3].$sam[5].$sam[6].$sam[7].$read_strand;
      foreach my $tel (keys %res){
        if($sidedness eq "a" and $res{$tel}{str} eq "+"){
          $te .= $tel.":".src::Utilities::get_read_mapping_span($res{$tel}{start},$res{$tel}{cigar}).",";   ## end of hte segment that maps on TE
        }elsif($sidedness eq "a" and $res{$tel}{str} eq "-"){
          $te .= $tel.":".$res{$tel}{start}.",";   ## adjust for clip ori 
        }elsif($sidedness eq "b" and $res{$tel}{str} eq "+"){
          $te .= $tel.":".$res{$tel}{start}.",";   ## adjust for clip ori 
        }elsif($sidedness eq "b" and $res{$tel}{str} eq "-"){
          $te .= $tel.":".src::Utilities::get_read_mapping_span($res{$tel}{start},$res{$tel}{cigar}).",";   ## end of hte segment that maps on TE
        }else{
          $te .= $tel.":".src::Utilities::get_clip_coord($res{$tel}{start},$res{$tel}{cigar}).",";   
        }
        $sc .= $res{$tel}{score}.",";
        if(!defined $maxscore or $maxscore<$res{$tel}{score}){
          $maxscore = $res{$tel}{score};
        }
        $mm .= $res{$tel}{mm}.",";
        $str .= $res{$tel}{str}.",";
      }
      $str =~ s/,$//;
      $mm =~ s/,$//;
      $te =~ s/,$//;
      $sc =~ s/,$//;
      $sam[-1] .= ";TE=".$te.";mm=".$mm.";sc=".$sc.";str=".$str.";end=".$end.";";
      print O2 join("\t",@sam),"\n";
      
      if($sam[-1] =~ m/DE,TP/ and $sam[4]>0){ ## only take locations with some TE mapping for cluster generation, that aren;t multimapper
         $hqual{$sam[2]."\t".$cliploc."\t".$sidedness}.=$id."|".$maxscore."|c:".$cliploc." ";
      }
      
      ## reset
      undef %res;
      $oid = $temp[0];  
      $res{$temp[2]}{score} = $score;
      $res{$temp[2]}{mm} = $mismatch;
      $res{$temp[2]}{str} = $strand;
      $res{$temp[2]}{start} = $temp[3];
      $res{$temp[2]}{cigar} = $temp[5];
    }
 }
}
###----------------------------------------------- END OF THE SCRIPT ---------------------------------------------------