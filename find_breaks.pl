#!/usr/bin/perl

# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard
# This script extacts break positions

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
use utf8;
use open qw(:std :utf8);
use Utilities;
#require "subroutines.pl";

BEGIN { our $start_run = time(); }

my $sam = "";
my $bam = "";
my $wd ="";
my $hd ="";
my $outprefix="project";
my $help=0;
Getopt::Long::GetOptions(
  'sam=s'            => \$sam,
  'bam=s'            => \$bam,
  'outprefix=s'      => \$outprefix,
  'wd:s'             => \$wd,
  'hd:s'             => \$hd,
  'help'             => \$help,
  'h'                => \$help,
  ) or die "Incorrect input! Use -h for usage.\n";

if ($help) {
  print "\nUsage: perl find_breaks.pl -bam [FILE_PATH] -outprefix [STRING] -wd [DIR_PATH] -hd [FILE_PATH] \n\n";
  print "This script generates a bam file (can be used for visualization) and breakpoint list file from mapped sam/bam\n\n";
  print "Options:\n";
  print "required**:\n";
  print "   -bam             sam/bam input file. Read name has information about other mates mapping.  \n";
  print "   -hd              Header file for output bam\n";
  print "optional:\n";
  print "   -wd              Working directory [default: ~]\n";
  print "   -outprefix       Outputprefix for generating the output files \n";
  print "   -help|-h         Display usage information.\n\n\n";
  print "Default outputs:\n";
  print "   (outprefix).clusters.bed\n";
  print "   (outprefix).RAM.bam\n";
  exit 0;
}
$wd =~ s/\/$//;
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
my $watch_run=0;
my $run_time=0;

#--------------------------------------------------------------------------------------------------------------
# main body
#--------------------------------------------------------------------------------------------------------------
$watch_run = time();
$run_time = $watch_run - our $start_run;
print "\n[extract_clusters=Start_Time]:\t $run_time seconds\n";
$outprefix =~ s/\?//g;
$outprefix =~ s/\#//g;
$outprefix =~ s/\'//g;
$outprefix =~ s/\"//g;
$outprefix =~ s/\s*//g;


my $out1 = $outprefix.".temp.clusters.bed";
open O1, ">> $out1" or die "can't create $out1";
my $out2 = $outprefix.".temp.RAM.bam";
open O2, "| samtools view -@ 4 -bS - >> $out2" or die "can't create $out2";

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
 
 while(<$file>){
    next if(/^(\#)/);
    next if(/^(\@)/); 
    chomp;
    s/\r//;  
    next if($_ eq "");
    my @temp = split(/\t+/);
    next if($temp[1] >=256);  
    next if($temp[1] ==4 and $temp[0]!~ /DE,TP/);

    my ($score) = $_ =~ /\tAS:i:(\d*)/;
    my $teclust = "-";
    my $strand = "*";
    my ($mismatch) = "NA";
    my $read_strand="+";
    
    my @sam = split("\x{019}",$temp[0]);
    $read_strand="-" if($sam[1] & 0x10);
  
    ## Add map-start on concensus
    if($temp[1] !=4){
      undef($teclust);
      $strand = "+";
      $strand = "-" if($temp[1] & 0x10);
      if($read_strand eq "-" and $strand eq "+"){
        $teclust .= $temp[2].":".$temp[3].";";   ## adjust for clip ori
      }else{
        $teclust .= $temp[2].":".Utilities::get_clip_coord($temp[3],$temp[5]).";";   
      }      
      if($_=~ /XA:Z:(\S*);$/){
        my @altmap = split(";",$1);
        foreach my $a (@altmap){
           my @al = split(",",$a);
           if($al[0] ne $temp[2]){
             my $sstr = "+";
             if($al[1] =~ m/-\d*/){ $sstr = "-";}
             $al[1]=~ s/[+-]//g;
             if($read_strand eq "-" and $sstr eq "+"){
                $teclust .= $al[0].":".$al[1].";"; 
              }else{
                $teclust .= $al[0].":".Utilities::get_clip_coord($al[1],$al[2]).";"; 
             }
           }
        }  
      }
    }
    $teclust =~ s/;$//;
    if($_=~ /\tNM:i:(\d*)/){ $mismatch=$1;}
    $sam[-1] .= ";TE=".$temp[2].";mm=".$mismatch.";sc=".$score.";str=".$strand.";cl=".$teclust;
    my $id = $sam[2].$sam[3].$sam[5].$sam[6].$sam[7];
    my $end = Utilities::get_read_mapping_span($sam[3],$sam[5]);
    $sam[-1] .= ";end=".$end;
    
    my ($cliploc) = $sam[-1] =~ m/;clip=(\d+);/;
    my ($sidedness) = $sam[-1] =~ m/;side=(\S+?);/;     
    if($temp[0]=~ m/DE,TP/ and $temp[2] ne "*" and $sam[4]>0){ ## only take locations with some TE mapping for cluster generation, that aren;t multimapper
       $hqual{$sam[2]."\t".$cliploc."\t".$sidedness}.=$id."|".$score."|c:".$cliploc." ";
    }     
    print O2 join("\t",@sam),"\n";  
 }
}
###----------------------------------------------- END OF THE SCRIPT ---------------------------------------------------