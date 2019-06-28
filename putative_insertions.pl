#!/usr/bin/perl
# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard

# This program identifies putative insertion candidates from the annotated breaks

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use Utilities;
use open qw(:std :utf8);
BEGIN { our $start_run = time(); }

my $in="";
my $index ="";
my $polym ="";
my $repbase ="";
my $bgAnnotations ="";
my $sites ="";
my $motif ="GATC";
my $ncores=4;
my $refMapqQ = 28;
my $TEMapScore = 30;
my $clip = 20;
my $outprefix="project";
my $help=0;
my $wd = "";
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'index=s'          => \$index,
  'polym=s'          => \$polym,
  'repbase:s'        => \$repbase,
  'bgAnnotations:s'  => \$bgAnnotations,
  'outprefix=s'      => \$outprefix,
  'motif=s'          => \$motif,
  'clip=s'           => \$clip,
  'sites=s'          => \$sites,
  'wd:s'             => \$wd,
  'ncores:s'         => \$ncores,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl putative_insertions.pl -in [FILE_PATH] -outprefix [STRING] -polym [PATH] -index [PATH] -repbase [PATH] -bgAnnotations [DIRPATH] -clip [INT] -sites [PATH] -motif [STRING] -ncores [INT] -wd [DIR_PATH]\n\n";
  print "This program annotates finalized breaks using the bam file\n\n";
  print "Options:\n\n";
  print "***required:\n";
  print "  -in                    Clusters file in bed format \n";
  print "  -index                 TE consensus index (generated using BWA)\n";
  print "  -polym                 TE polymorphic example cases index (generated using BWA)\n";
  print "  -repbase               RepBase index (generated using BWA)\n";
  print "  -bgAnnotations         Reference TE annotations\n";
  print "  -sites                 RE motif sites\n";
  print "  -motif                 RE motif (default:GATC)\n";
  print "***optional:\n";
  print "  -clip                  Clip length cutoff [default: 20] \n";
  print "  -outprefix             Outputprefix for generating 2 output files [default: project] \n";
  print "  -wd                    Working directory [default: ~]\n";
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

#--------------------------------------------------------------------------------------------------------------
# I/O
#--------------------------------------------------------------------------------------------------------------
## Dependancies
print " ERROR: Repbase index not provided\n" unless -e $repbase;
print " ERROR: TE-consensus index not provided\n" unless -e $index;
print " ERROR: Restriction motif index is not provided\n" unless -e $sites;
print " ERROR: Background TE annotations are not provided\n" unless -e $bgAnnotations;
print " ERROR: countmodel.R does not exist in the same directory!\n" unless -e "countmodel.R";
print " ERROR: pairLocations.R does not exist in the same directory!\n" unless -e "pairLocations.R";
my %cov;
my $pval_cutoff = 0.05;
my %transposons;
my %polymorphs;
my %crossref;
my %clusters;
my %chromosomes; ## binary
my $olog = $outprefix.".skippedclusters.logs.gz";
open(LOGS,"| gzip -c - > $olog") or die $!;


my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "[putative_insertions] START:\t $run_time seconds\n";
print " Command: perl putative_insertions.pl -in $in -sites $sites -index $index -repbase $repbase -bgAnnotations $bgAnnotations -motif $motif -ncores $ncores -outprefix $outprefix -wd $wd \n";


## get transposons
%transposons = Utilities::get_fasta($index);
print " transposons being considered in the analyses: " , join(",",keys %transposons),"; sizes in bp (",join (",",values %transposons),")\n";
$polym =~ s/\/$//;
foreach my $i (keys %transposons){
  next if($i eq "PolyA");
   my $f = $polym."/".$i.".fa";
   my %r = Utilities::get_fasta($f);
   foreach my $y(keys %r){
     $crossref{$y}=$i;
     $polymorphs{$y}=$r{$y};
   }
}
print " Polymorphic sequences considered in the analyses: " , join(",",keys %polymorphs),"\n";
print " ------------------------------------sizes in bp : " , join(",",values %polymorphs),"\n";
print "  Crossref: ", join(",",keys %crossref), "\n";
print "  -----for: ", join(",",values%crossref),"\n";
## retrieve cluster from the input file
%clusters = %{retrieve($in)};
print " read clusters..\n";

## completing breakpoint object
my $clust = complete_clust_object(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " completed breakpoint object:\t $run_time seconds\n";
print "  -> total breakpoint locations in completed object: ", scalar(keys%clusters),"\n";

## merging breakpoint object
$clust = merge_cluster_entries(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " merged LHS/RHS breakpoints:\t $run_time seconds\n";
print "  -> total putative insertions in completed object: ", scalar(keys%clusters),"\n";
if(scalar(keys%clusters)==0){
  print "  Aborting as there are no clusters identified in the input file \n";
  exit 1;
}

## generate binary for RE sites
%cov = %{retrieve($outprefix.".coverage.ph")};
open FILE, "$sites " or die "Can't open the restriction sites file. Check file path!!";
while (<FILE>) {
  my @locs = split(" ",$_);
  my $key = shift(@locs);
  $chromosomes{$key} = \@locs;
}
close(FILE);
print " generated binary index for RE sites..\n";

## finalize breakpoint object
$clust = finalize_breakpointObject(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " finalized breakpoints with annotations:\t $run_time seconds\n";
print "  -> total putative insertions in completed object: ", scalar(keys%clusters),"\n";

## subfamily annotations of the putative insertions
$clust = subfamily_annotation(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " assessed subfamily annotations of the putative insertions:\t $run_time seconds\n";
print "  -> total putative insertions in completed object: ", scalar(keys%clusters),"\n";


## finding distance to the reference TE copy
my $testring = join(",",keys %transposons);
$testring =~ s/PolyA//;
$testring=~ s/,,/,/;
$testring=~ s/^,//;
$testring=~ s/,$//;
$testring=~ s/^\S*?\///g;
$testring=~ s/,\S*?\//,/g;
my @tes = split(",",Utilities::unique($testring,",","false"));
foreach my $te (@tes){
  next if($te eq "PolyA");
  print " Finding background distances for $te insertions\n";
  if(-e $bgAnnotations."/".$te.".bed"){
    $clust = backgroundTE_overlaps($bgAnnotations."/".$te.".bed",$te,\%clusters);
  }elsif(-e $bgAnnotations."/".$te.".bed.gz"){
    $clust = backgroundTE_overlaps($bgAnnotations."/".$te.".bed.gz",$te,\%clusters);
  }
  %clusters =%{$clust};
}
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " calculated distances to the nearest TE copies:\t $run_time seconds\n";


## background enrichment
$clust = perform_bg_enrichment(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " performed background enrichment:\t $run_time seconds\n";
print "  -> total putative insertions in completed object: ", scalar(keys%clusters),"\n";


store \%clusters, $outprefix.'.ClustObj_PutativeInsertions.ph'; #Save

$watch_run = time();
$run_time = $watch_run - $start_run;
print "[putative_insertions] END:\t $run_time seconds\n";

exit;

#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub complete_clust_object {
  my ($clust) = shift;
  my %clusters = %{$clust};
  my $cfo = $outprefix.".temp.remap.fasta";
  open(FO1,">$cfo") or die "Can't generate $cfo\n";
   

  foreach my $loc (sort keys %clusters){
    my @temp = split("\t",$loc);
    $clusters{$loc}{pos1}{chr} = $temp[0] ;
    $clusters{$loc}{pos1}{cliploc} = $temp[1];
    $clusters{$loc}{pos1}{strand} = "*";
    if(!defined($clusters{$loc}{pos1}{maxscore}) or !exists($clusters{$loc}{pos1}{dist2re}) or !exists($clusters{$loc}{pos1}{seqstring})){
      print qq[ ERROR: Can't find sequence maxscore/dist2re/seqstring for: $loc\n];
      print Dumper $clusters{$loc};
      exit 1;
    }
    
    my @features = ("unmapped","fuzzyclipped",keys %transposons);  
    foreach my $feature (@features){
      if(!exists($clusters{$loc}{pos1}{$feature})){
          $clusters{$loc}{pos1}{$feature} = 0;
        }else{
        $clusters{$loc}{pos1}{$feature} = Utilities::unique($clusters{$loc}{pos1}{$feature}, ",","true"); 
       }
    }
    undef(@features);
    @features = ("clipseqs"); 
    foreach my $feature (@features){
      if(!exists($clusters{$loc}{pos1}{$feature})){
          $clusters{$loc}{pos1}{$feature} = ""; 
        }else{
        $clusters{$loc}{pos1}{$feature} = Utilities::unique($clusters{$loc}{pos1}{$feature},",","false");
       }
    }
    
    ##Unmapped cluster mate
    if(!exists $clusters{$loc}{pos1}{smate}){
       $clusters{$loc}{pos1}{smate} ="-";
       $clusters{$loc}{pos1}{smate_clipseqs} = "-";
    }else{
      $clusters{$loc}{pos1}{smate}=~ s/,$//;
      my @l = split(",",Utilities::unique($clusters{$loc}{pos1}{smate},",","false") );
      my %fq;
      foreach my $i (@l) {
        $i=~ s/\|.*$//;
        $fq{$i}++;
      }
      my @k = keys %fq;
      @k = sort{$fq{$b} <=> $fq{$a}} keys %fq if(scalar @k>1);
      $clusters{$loc}{pos1}{smate} = $k[0]."|".$fq{$k[0]};
    }
    ## identify cluster points within TE assembly
      ## (0)RecrpClustLoc, (1)RecrpClustReads, (2)RecprTE, (3)RecprPolyA
    $clusters{$loc}{pos1}{RecrpClustLoc} = "-";
    $clusters{$loc}{pos1}{RecrpClustReads} = "-";
    $clusters{$loc}{pos1}{RecrpTE} = "-";
    $clusters{$loc}{pos1}{RecrpPolyA} = "FALSE";
    $clusters{$loc}{pos1}{RSLoc} = "-";
    $clusters{$loc}{pos1}{RSReads} = "-";
    $clusters{$loc}{pos1}{RSTE} = "-";
    $clusters{$loc}{pos1}{RSPolyA} = "FALSE";    
    if(exists $clusters{$loc}{pos1}{TEClust} and $clusters{$loc}{pos1}{TEClust} ne "-"){  
      my @array = @{Utilities::get_teMap_clusterfreq($clusters{$loc}{pos1}{TEClust},\%transposons)}; ## TE:NumReads:ClusterLoc
      $clusters{$loc}{pos1}{RecrpClustLoc} = $array[0];
      $clusters{$loc}{pos1}{RecrpClustReads} = $array[1];
      $clusters{$loc}{pos1}{RecrpTE} = $array[2];
      $clusters{$loc}{pos1}{RecrpPolyA} = $array[3];
      $clusters{$loc}{pos1}{TEClust} =$array[2].":".$array[1].":".$array[0];
    }else{
      $clusters{$loc}{pos1}{TEClust} ="-";
    }
    if(exists $clusters{$loc}{pos1}{smateClust} and $clusters{$loc}{pos1}{smateClust} ne "-"){
      my @array = @{Utilities::get_teMap_clusterfreq($clusters{$loc}{pos1}{smateClust},\%transposons)}; ## TE:NumReads:ClusterLoc
      $clusters{$loc}{pos1}{RSLoc} = $array[0];
      $clusters{$loc}{pos1}{RSReads} = $array[1];
      $clusters{$loc}{pos1}{RSTE} = $array[2];
      $clusters{$loc}{pos1}{RSPolyA} = $array[3];      
      $clusters{$loc}{pos1}{smateClust} =$array[2].":".$array[1].":".$array[0];  #TE:NumReads:ClusterLoc
    }else{
      $clusters{$loc}{pos1}{smateClust}="-";
    }
    ## check full length SVA by scanning the repeats in unmapped clip sequences
    #if(exists $clusters{$loc}{pos1}{smate_clipseqs} and $clusters{$loc}{pos1}{RSTE} eq "-" and $clusters{$loc}{pos1}{smate} ne "-"){
    #  my @ttm = $clusters{$loc}{pos1}{smate} =~ /^(\S+)_(\d*)_(\w*)_(\S*)\|(\d*)$/;  ## chr_start_side_te|num
    #  if($ttm[3] eq "*"){
    #    my @array = @{scanSVArepeats($clusters{$loc}{pos1}{smate_clipseqs}) }; ## TE:NumReads:ClusterLoc
    #    if($array[2] eq "SVA"){
    #      $clusters{$loc}{pos1}{RSLoc} = $array[0];
    #      $clusters{$loc}{pos1}{RSReads} = $array[1];
    #      $clusters{$loc}{pos1}{RSTE} = $array[2];
    #      $clusters{$loc}{pos1}{smateClust} =$array[2].":".$array[1].":".$array[0];  #TE:NumReads:ClusterLoc 
    #      $clusters{$loc}{pos1}{smate} = $ttm[0]."_".$ttm[1]."_".$ttm[2]."_".$array[2]."|".$ttm[4];
    #      print " rescued SVA: $loc\n";
    #    }             
    #  }      
    #}
    if(exists $clusters{$loc}{pos1}{smate_clipseqs} and $clusters{$loc}{pos1}{RSTE} eq "-"){
      my $l = $loc;
      $l =~ s/\t/_/g;
      my $cnt=0;
      $clusters{$loc}{pos1}{smate_clipseqs} =~ s/,$//;
      my @t = split(",",$clusters{$loc}{pos1}{smate_clipseqs});
      foreach my $s(@t){
        $cnt++;
        next if($s eq "-" or $s eq "" or $s eq "*");
        print FO1 ">$l","|$cnt\n$s\n";
      }
    }

  }
  close(FO1);
  return(\%clusters);
}


sub merge_cluster_entries{
  my ($clust) = shift;
  my %clusters = %{$clust};
  
  my $out =$outprefix.".temp.clusters.txt";
  my $fname = $outprefix.".temp.clusterTemplates.txt";  
  my %olaps;
  my %new;

  open(FC,">$out") or die $!;
  print FC "chr\tstart\tend\tstrand\tclipcoord\tside\treads\tscore\n";
  foreach my $loc (sort keys %clusters) {
    my @l = split("\t",$loc);
    print FC "$clusters{$loc}{pos1}{chr}\t";
    print FC "$clusters{$loc}{pos1}{startspan}\t";
    print FC "$clusters{$loc}{pos1}{endspan}\t";
    print FC "$clusters{$loc}{pos1}{strand}\t";
    print FC "$clusters{$loc}{pos1}{cliploc}\t";
    print FC "$l[2]\t";
    print FC "$clusters{$loc}{pos1}{clipreads}\t";
    print FC "$clusters{$loc}{pos1}{maxscore}\n";
    $clusters{$loc}{pos1}{clip_sidedness} = $l[2];
  }
  close(FC);
  
  ## Find cluster of clip coordinates based on read mapping span atthe clip location
  if(system("Rscript pairLocations.R $out $fname")==0){
    print "  -> Generated list of breakpoint pairs!\n";
    #system("rm $out");    
  }else{
    print qq[ Error occured while merging the clusters: \n"
              Possible reason (1): GenomicRanges package not installed\n";
                              (2): pairLocations.R script is not present in the working directory\n];
    exit 1;
  }
  
  open (FH,"$fname") or die $!;
  while(<FH>){
    next if(/^(\#)/);
    next if(/^(\@)/); 
    chomp;
    s/\r//;  
    my @temp = split(/\t/);  # id seqnames     coord side  au bu
    my $l = $temp[1]."\t".$temp[2]."\t".$temp[3];
    $olaps{$temp[0]} .= $l.",";
  }
  close(FH);
  print "  -> Total breakpoint pairs in the object: ", scalar(keys %olaps),"\n";

  # Re-generate merged breakpoint object
  #open(FO1,">",$wd."/MyTestClusters_Olaps.txt") or die $!;
  while(my ($l,$j) = each %olaps){
    $j=~ s/,$//;
    $j =Utilities::unique($j,",","false");
    my @temp = split(",",$j);
    my $size = scalar(@temp);
    
    ## only one end softclip at isertion
    if($size ==1){
        %{$new{$j}{pos1}} = %{$clusters{$j}{pos1}};
        %{$new{$j}{pos2}} = %{$clusters{$j}{pos1}};
        delete $clusters{$j};
        $new{$j}{pos1}{bgclusters} = $size;
        $new{$j}{pos2}{unmapped} = 0;
        $new{$j}{pos2}{clipreads} = 0;
        $new{$j}{pos2}{fuzzyclipped} = 0;
        $new{$j}{pos2}{smate} = "-";
        $new{$j}{pos2}{smate_clipseqs} = "-";
        foreach my $te(keys %transposons){
          $new{$j}{pos2}{$te} = 0;
        }        
    }elsif($size == 2){
        %{$new{$temp[0]}{pos1}} = %{$clusters{$temp[0]}{pos1}};
        %{$new{$temp[0]}{pos2}} = %{$clusters{$temp[1]}{pos1}};
        $new{$temp[0]}{pos1}{bgclusters} = $size;
        $new{$temp[0]}{pos1}{smate_clipseqs} ="-";
        $new{$temp[0]}{pos2}{smate_clipseqs} ="-";
    }elsif($size == 3){
      print "  -> Warning: Descrepancy while selecting right and left clip coordinates from the breakpoints: $temp[0]:$temp[1]:$temp[2]\n";
    }
  }
  undef(%clusters);
  %clusters = %new;
  return(\%clusters); 
}

sub finalize_breakpointObject {
   my ($clust) = shift;
   my %clusters = %{$clust};
   my $filtered = 0;
   
   foreach my $loc(sort keys %clusters){
     if(!defined($clusters{$loc}{pos2}{chr}) or !defined($clusters{$loc}{pos2}{cliploc}) ) {
       print " ERROR: LHS/RHS not defined: $loc\n";
       exit 1;
     }
     ## (1) If other breakpoint is unmapped, use the unmapped clip mate for completing the breakpoint annotation
     if($clusters{$loc}{pos1}{cliploc}==$clusters{$loc}{pos2}{cliploc} and $clusters{$loc}{pos1}{clip_sidedness} eq $clusters{$loc}{pos2}{clip_sidedness} and $clusters{$loc}{pos1}{smate} ne "-"){
        my ($coord,$sd,$te,$a) = $clusters{$loc}{pos1}{smate} =~ /^\S+_(\d*)_(\w*)_(\S*)\|(\d*)$/;
        if($a>=1){
          $clusters{$loc}{pos2}{cliploc} = $coord;
          $clusters{$loc}{pos2}{clip_sidedness} = $sd;
          $clusters{$loc}{pos2}{RecrpClustLoc} = "-";
          $clusters{$loc}{pos2}{RecrpClustReads} = "-";
          $clusters{$loc}{pos2}{RecrpTE} = "-";
          $clusters{$loc}{pos2}{RecrpPolyA} = "FALSE";
          if($te ne "*"){
            $clusters{$loc}{pos2}{$te} = $a;
            if($clusters{$loc}{pos1}{smateClust} ne "-"){
              my ($xx,$xy,$xz) = $clusters{$loc}{pos1}{smateClust} =~ m/(\S*):(\S*):(\S*)/; ##TE:NumReads:ClusterLoc
              if($te eq $xx){
                $clusters{$loc}{pos2}{RecrpClustLoc} = $xz;
                $clusters{$loc}{pos2}{RecrpClustReads} = $xy;
                $clusters{$loc}{pos2}{RecrpTE} = $xx;  
                $clusters{$loc}{pos2}{RecrpPolyA} = "TRUE" if($te eq "PolyA");          
                $clusters{$loc}{pos2}{TEClust} = $clusters{$loc}{pos1}{smateClust};
              }
            }
          }
          if(exists $clusters{$loc}{pos1}{smate_seqstring}){
            $clusters{$loc}{pos2}{seqstring} = $clusters{$loc}{pos1}{smate_seqstring};
            my $aa = Utilities::getstrand($clusters{$loc}{pos1}{smate_seqstring});
            if($aa eq "NA" and $te eq "*"){ ## if there is an unmapped polyA tail, use the info to fill slot2
              $clusters{$loc}{pos2}{unmapped} += $a;  
            }elsif($aa ne "NA" and $te eq "*"){
              $clusters{$loc}{pos2}{PolyA} = $a;
              $clusters{$loc}{pos2}{RecrpPolyA} = "TRUE";    
            }
          }
        }
     }
     ## (2) Counts 
     my @features = ("unmapped","fuzzyclipped",keys %transposons);
     foreach my $feature (@features){
        $clusters{$loc}{final}{$feature} = ($clusters{$loc}{pos1}{$feature}+$clusters{$loc}{pos2}{$feature});
     }
     $clusters{$loc}{final}{dist2re} = Utilities::max(($clusters{$loc}{pos1}{dist2re},$clusters{$loc}{pos2}{dist2re}));
     $clusters{$loc}{final}{maxscore} = Utilities::max(($clusters{$loc}{pos1}{maxscore},$clusters{$loc}{pos2}{maxscore}));
   
     ## (3) Locus 
     $clusters{$loc}{final}{chr} = $clusters{$loc}{pos1}{chr};
     $clusters{$loc}{final}{start} = Utilities::min(($clusters{$loc}{pos1}{cliploc},$clusters{$loc}{pos2}{cliploc})) ;
     $clusters{$loc}{final}{end}   = Utilities::max(($clusters{$loc}{pos1}{cliploc},$clusters{$loc}{pos2}{cliploc}));     

     ## (4) Determine inserted TE
     $clusters{$loc}{final}{TE} ="*";
     $clusters{$loc}{final}{TE_mapReads} = 0;
     foreach my $v(keys %transposons){
        if($v ne "PolyA" and $clusters{$loc}{final}{$v}>0){
          if($clusters{$loc}{final}{TE_mapReads} < $clusters{$loc}{final}{$v}){
             $clusters{$loc}{final}{TE} =$v; 
             $clusters{$loc}{final}{TE_mapReads} = $clusters{$loc}{final}{$v};
          }  
        } 
     }
     $clusters{$loc}{final}{TE_mapReads} += $clusters{$loc}{final}{PolyA};       

     ## (5) TE string per side
     $clusters{$loc}{pos1}{TEString}="*";
     $clusters{$loc}{pos2}{TEString}="*";
     my %temp_max;
     foreach my $v(keys %transposons){
       $temp_max{$v} = $clusters{$loc}{pos1}{$v};
     }
     my @maxkeys = sort { $temp_max{$b} <=> $temp_max{$a} } keys %temp_max;
     $clusters{$loc}{pos1}{TEString} = $maxkeys[0] if($temp_max{$maxkeys[0]}>0);
     undef(%temp_max);
     undef(@maxkeys);
     foreach my $v(keys %transposons){
       $temp_max{$v} = $clusters{$loc}{pos2}{$v};
     }
     @maxkeys = sort { $temp_max{$b} <=> $temp_max{$a} } keys %temp_max;
     $clusters{$loc}{pos2}{TEString} = $maxkeys[0] if($temp_max{$maxkeys[0]}>0);
     
     ## (6) strand and tailinfo
     $clusters{$loc}{pos1}{polyAT} = Utilities::getstrand($clusters{$loc}{pos1}{seqstring});
     $clusters{$loc}{pos2}{polyAT} = Utilities::getstrand($clusters{$loc}{pos2}{seqstring});
     $clusters{$loc}{final}{polyAT} = Utilities::getstrand($clusters{$loc}{pos1}{seqstring},$clusters{$loc}{pos2}{seqstring});   
     $clusters{$loc}{final}{strand} = "+" if($clusters{$loc}{final}{polyAT} eq "PolyA") ;
     $clusters{$loc}{final}{strand} = "-" if($clusters{$loc}{final}{polyAT} eq "PolyT") ;
     $clusters{$loc}{final}{strand} = "*" if($clusters{$loc}{final}{polyAT} eq "NA") ;
     if($clusters{$loc}{pos1}{cliploc} == $clusters{$loc}{pos2}{cliploc} and $clusters{$loc}{pos1}{clip_sidedness} eq $clusters{$loc}{pos2}{clip_sidedness}){
        $clusters{$loc}{final}{tailinfo} = $clusters{$loc}{pos1}{polyAT};      
     }else{
        $clusters{$loc}{final}{tailinfo} = $clusters{$loc}{pos1}{polyAT}.",".$clusters{$loc}{pos2}{polyAT}; 
        $clusters{$loc}{final}{tailinfo} =~ s/NA,//;
        $clusters{$loc}{final}{tailinfo} =~ s/,NA//;
     }
     
     ## (7) TSD
     my @l = split(",",($clusters{$loc}{pos1}{seqstring}.",".$clusters{$loc}{final}{start}.",".$clusters{$loc}{final}{end}));
     my $tsd = Utilities::tsd(\@l);
     if($tsd eq "" or $tsd eq "ERROR-TSD"){
         @l = split(",",($clusters{$loc}{pos2}{seqstring}.",".$clusters{$loc}{final}{start}.",".$clusters{$loc}{final}{end}));
         $tsd = Utilities::tsd(\@l);
     }
     if($tsd eq "" or $tsd eq "ERROR-TSD"){
         $clusters{$loc}{final}{TSD} = "*"; 
      }else{
        $clusters{$loc}{final}{TSD} = $tsd;
     }
     
     ## (8) distance to annotated RE motif
     my $dist2remotif="NA";
     if($clusters{$loc}{final}{start} == $clusters{$loc}{final}{end}){
         my $a = &binarysearch($clusters{$loc}{final}{start}, $chromosomes{$clusters{$loc}{final}{chr}});
         if($a <= $clusters{$loc}{final}{start} and ($a+length($motif)) >= $clusters{$loc}{final}{start} ) {
           $dist2remotif =0;
         }elsif($a < $clusters{$loc}{final}{start}){
           $dist2remotif = abs($a-$clusters{$loc}{final}{start}) - length($motif);
         }else{
          $dist2remotif= abs($a-$clusters{$loc}{final}{start});
         } 
     }else{
         my $x = 0;
         my $a = &binarysearch($clusters{$loc}{final}{start}, $chromosomes{$clusters{$loc}{final}{chr}});
         if($a <=$clusters{$loc}{final}{start} and ($a+length($motif)) >= $clusters{$loc}{final}{start}){
           $x=0;
         }elsif($a < $clusters{$loc}{final}{start}){
            $x = abs($a-$clusters{$loc}{final}{start}) - length($motif);
         }else{
            $x = abs($a-$clusters{$loc}{final}{start});
         }
         my $y = 0;
         my $b = &binarysearch($clusters{$loc}{final}{end}, $chromosomes{$clusters{$loc}{final}{chr}});
         if($b <=$clusters{$loc}{final}{end} and ($b+length($motif)) >= $clusters{$loc}{final}{end}){
           $y=0;
         }elsif($a < $clusters{$loc}{final}{end}){
            $y = abs($a-$clusters{$loc}{final}{end}) - length($motif);
         }else{
            $y = abs($a-$clusters{$loc}{final}{end});
         }
         $dist2remotif= Utilities::min( ($x,$y) );
     }
     $clusters{$loc}{final}{dist2remotif} = $dist2remotif;
     
     ## (9) Purity of the cluster
     my $purity_ratio = 0;
     my $tsp  = $clusters{$loc}{final}{TE};
     $clusters{$loc}{final}{AllClip}=0;
     foreach my $v(keys %transposons){
         $clusters{$loc}{final}{AllClip} += $clusters{$loc}{final}{$v};
     }
     $clusters{$loc}{final}{AllClip} +=$clusters{$loc}{final}{unmapped};
     if($tsp ne "*" and $clusters{$loc}{final}{AllClip} > 0){
       $purity_ratio =  Utilities::round( ($clusters{$loc}{final}{$tsp}+$clusters{$loc}{final}{PolyA})/($clusters{$loc}{final}{AllClip}),2);
     }
     $clusters{$loc}{final}{purity} = $purity_ratio;     
     
     ## (10) Insertion size
     $clusters{$loc}{final}{insert_size}="-";
     if($tsp ne "*" and $clusters{$loc}{final}{start} != $clusters{$loc}{final}{end}){
        if($clusters{$loc}{pos1}{RecrpPolyA} eq "TRUE" and $clusters{$loc}{pos2}{RecrpTE} eq $tsp){
           $clusters{$loc}{final}{insert_size}= abs($transposons{$tsp} - $clusters{$loc}{pos2}{RecrpClustLoc} );
        }elsif($clusters{$loc}{pos2}{RecrpPolyA} eq "TRUE" and $clusters{$loc}{pos1}{RecrpTE} eq $tsp){
          $clusters{$loc}{final}{insert_size}= abs($transposons{$tsp} - $clusters{$loc}{pos1}{RecrpClustLoc} );
        }elsif($clusters{$loc}{pos2}{RecrpPolyA} ne "TRUE" and $clusters{$loc}{pos1}{RecrpPolyA} ne "TRUE"){
          if($clusters{$loc}{pos1}{RecrpTE} eq $tsp and $clusters{$loc}{pos2}{RecrpTE} eq $tsp){
            $clusters{$loc}{final}{insert_size}= abs($clusters{$loc}{pos1}{RecrpClustLoc} - $clusters{$loc}{pos2}{RecrpClustLoc});
          }
        }elsif($clusters{$loc}{final}{polyAT} !~ m/NA/){
          my $m = Utilities::min(($clusters{$loc}{pos1}{RecrpClustLoc},$clusters{$loc}{pos2}{RecrpClustLoc}));
          $clusters{$loc}{final}{insert_size}= abs($transposons{$tsp} - $m ) if (defined $m);
        }  
     }

     ## (11) Reciprocal clusters percent. Include PolyA reads as well in the analysis
     $clusters{$loc}{pos1}{RecrpClustPercent}=0;
     $clusters{$loc}{pos2}{RecrpClustPercent}=0;
     if($tsp ne "*"){
       if($clusters{$loc}{pos1}{cliploc}==$clusters{$loc}{pos2}{cliploc} and $clusters{$loc}{pos1}{clip_sidedness} eq $clusters{$loc}{pos2}{clip_sidedness}){
         if($clusters{$loc}{pos1}{RecrpTE} eq $tsp){
            $clusters{$loc}{pos1}{RecrpClustPercent} = Utilities::round( ($clusters{$loc}{pos1}{RecrpClustReads}*100)/$clusters{$loc}{pos1}{$tsp},2);
         } 
       }else{
         if($clusters{$loc}{pos1}{RecrpTE} eq $tsp){
           if($clusters{$loc}{pos1}{$tsp}>0){
             $clusters{$loc}{pos1}{RecrpClustPercent} = Utilities::round( ($clusters{$loc}{pos1}{RecrpClustReads}*100)/$clusters{$loc}{pos1}{$tsp},2);
           }else{
            print Dumper $clusters{$loc};
           }
         }
         if($clusters{$loc}{pos2}{RecrpTE} eq $tsp){
           if($clusters{$loc}{pos2}{$tsp}>0){          
            $clusters{$loc}{pos2}{RecrpClustPercent} = Utilities::round( ($clusters{$loc}{pos2}{RecrpClustReads}*100)/$clusters{$loc}{pos2}{$tsp},2);
           }else{
            print Dumper $clusters{$loc};
           }
         }
       }
     }
     $clusters{$loc}{final}{RecrpClustPercent}= Utilities::max(($clusters{$loc}{pos1}{RecrpClustPercent},$clusters{$loc}{pos2}{RecrpClustPercent}));

     ## a) remove cluster that have no TE mapping reads 
     if(!defined $tsp or $tsp eq "*"){
       print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\tNo TE-mapping clip reads (filter step)\tMapReads:-;-\t$clusters{$loc}{pos1}{PolyA};$clusters{$loc}{pos2}{PolyA} ";
       print LOGS "$clusters{$loc}{pos1}{Alu},$clusters{$loc}{pos1}{L1},$clusters{$loc}{pos1}{SVA},$clusters{$loc}{pos1}{HERVK};";
       print LOGS "$clusters{$loc}{pos2}{Alu},$clusters{$loc}{pos2}{L1},$clusters{$loc}{pos2}{SVA},$clusters{$loc}{pos2}{HERVK};\n";
       $filtered++;
       delete($clusters{$loc});
       next;
     }   
   }
   print " breakpoints where clipped part failed to align on non PolyA TE <2 reads mapping to TE assembly : $filtered \n";
   return(\%clusters);

}

sub binarysearch {
  my ($x, $a) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);          # lower, upper end of search interval
  my $i;                               # index of probe
  my $q1=0;
  my $q2=0;


  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
    } 
    else {
      # return $i+1; #original 
      $q1 = abs($a->[$i] - $x);
      $q2 = $q1;  
      if($a->[$i+1]){
        $q2 = abs($a->[$i+1] - $x);
      }
      if($q1 > $q2){
        return($a->[$i+1]); # found
      }else{
        return($a->[$i]); # found
      }
      #if($upper==1){
      #  return($q2);
      #}else{
      #  return($q1);
      #}
      
    }
  }
  $q2 = abs($a->[$l-1] - $x);
  $q1 = $q2;
  if($a->[$l]){
    $q1 = abs($a->[$l] - $x);  
  }
  if($q1 < $q2){
    return($a->[$l]); # found
  }else{
    return($a->[$l-1]); # found
  }
  #if($upper==1){
  #  return($q1);
  #}else{
  #  return($q2);
  #}
  #return($l); #original         
}


sub subfamily_annotation{
   my ($clust) = shift;
   my %clusters = %{$clust};
   ## (2) : Identify subfamily annotation (a subfamily fasta index library is needed)
   my $cfo = $outprefix.".temp.consensus.fasta";
   open(FO1,">$cfo") or die "Can't generate $cfo\n";
   
   foreach my $loc (keys %clusters){
     ## longest reads for subfamily annotation
     if($clusters{$loc}{pos1}{cliploc} < $clusters{$loc}{pos2}{cliploc}){
        $clusters{$loc}{final}{RC_consensus} = Utilities::get_longestseq($clusters{$loc}{pos1}{clipseqs});  
        $clusters{$loc}{final}{LC_consensus} = Utilities::get_longestseq($clusters{$loc}{pos2}{clipseqs});  
     }else{
        $clusters{$loc}{final}{RC_consensus} = Utilities::get_longestseq($clusters{$loc}{pos2}{clipseqs});  
        $clusters{$loc}{final}{LC_consensus} = Utilities::get_longestseq($clusters{$loc}{pos1}{clipseqs});  
     }
     my $l = $loc;
     $l=~ s/\t/_/g;
     print FO1 ">$l","_RC\n$clusters{$loc}{final}{RC_consensus}\n";
     print FO1 ">$l","_LC\n$clusters{$loc}{final}{LC_consensus}\n";   
     #$clusters{$loc}{final}{transduction} =;   
   }
   close(FO1);

   if(defined $repbase and $repbase ne ""){
     open(IN,"bwa mem -v 1 $repbase $cfo|") or die $!;
     my $qw = 0.4; ## arbitrary threshold - fraction of longest softclipped sequence mapping to reference 
     while (<IN>){
       next if(/^(\#)/);
       next if(/^(\@)/); 
       s/\n//; 
       s/\r//;  
       my @temp = split(/\t+/);
       next if($temp[1]>=256);

       my ($chr,$start,$strand,$side) = $temp[0]=~ m/(\S*)_(\S*)_(\S*)_(\S*)/;
       my ($score) = $_ =~ /AS:i:(\d*)/;
       my ($mismatch) = "NA";
       if($_=~ /NM:i:(\d*)/){
         $mismatch=$1;
       }
       my $loc = $chr."\t".$start."\t".$strand;
       my $matches=0;
       if($temp[5]=~m/(\d*)M/){ $matches+=$1;}
       if($temp[5]=~m/(\d*)I/){ $matches+=$1;}
       if($temp[5]=~m/(\d*)D/){ $matches+=$1;}
       my $rat = Utilities::round($score/length($temp[9]),2);
       $clusters{$loc}{final}{subfamilyMapInfo} .= $temp[2]."|".$side."|".$matches."|".length($temp[9])."|".$score."|".$mismatch."|".$rat.",";
  
       if($clusters{$loc}{final}{TE} eq "*"){
          $clusters{$loc}{final}{subfamily}=$temp[2];   
       }elsif(defined $clusters{$loc}{final}{subfamily}){
         if($rat > 0.5 and $rat>$qw and $temp[2] =~ m/$clusters{$loc}{final}{TE}/){
           $clusters{$loc}{final}{subfamily} = $temp[2]; 
         }
       }elsif($rat > 0.5 and $temp[2] =~ m/$clusters{$loc}{final}{TE}/ ){
          $clusters{$loc}{final}{subfamily} = $temp[2]; 
          $qw = $rat;
       }else{
           $clusters{$loc}{final}{subfamily} = "*"; 
           $qw = 0.4;
       }
     }
     close(IN);
   }else{
     foreach my $loc(sort keys %clusters){
       $clusters{$loc}{final}{subfamily} = "*";
       $clusters{$loc}{final}{subfamilyMapInfo} = "*";
     }
   }
   return(\%clusters);
}

sub backgroundTE_overlaps{
  my ($file,$te,$clust) =@_;
  my %clusters =%{$clust};
    
  my $cnt=0;
  my $fo = $outprefix.".temp.".$te.".mergedclusters.bed";
  open(FO,"| sort -k1,1 -k2,2n - > $fo") or die $!;
  foreach my $loc (sort keys %clusters){
     my $tsp = $clusters{$loc}{final}{TE};
     next if($tsp eq "*" or $tsp eq "PolyA");
     $tsp =~ s/^\S*\///g;
     if($tsp eq $te){  
        $cnt++;
        print FO "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
        print FO ".\t.\t$clusters{$loc}{final}{strand}\t";
        $loc =~ s/\t/_/g;
        print FO "$loc\n";
     }
  }
  close(FO);

  my $r=0;
  if(scalar(keys %clusters) > 0){
     $r = Utilities::round($cnt*100/scalar(keys %clusters),2);
  } 
  print " Total $cnt putative $te insertions out of total ",scalar(keys %clusters)," insertion candidates ($r%)\n";
  if($cnt==0){
     system("rm $fo");
     return(\%clusters);
  }

  my $cnt1=0;
  open(IN,"bedtools closest -t first -d -a $fo -b $file|") or die $!;
  while(<IN>){
       next if(/^(\#)/);
       next if(/^(\@)/); 
       s/\n//; 
       s/\r//;  
       my @temp = split(/\t+/);
       my $loc = $temp[6];
       $loc=~ s/_/\t/g;
       if($clusters{$loc}){
          $cnt1++; 
          $clusters{$loc}{final}{bgTE} = $temp[10];
          $clusters{$loc}{final}{bgTEDist} = $temp[13];
        }else{
          print " ERROR in parsing the clusters: sub backgroundTE_overlaps \n";
          exit 1;
        }
  }
  if($cnt != $cnt1){
       print " WARNING: Numbers in background cluster proximity test do not add up for $te.". ($cnt-$cnt1)." short. Will be labelled as \"-\"\n";
  }
  system("rm $fo");
  return(\%clusters);
} 

sub perform_bg_enrichment{
  my ($clust) = shift;
  my %clusters = %{$clust};
  my $clustnum=0;
  my $file = $outprefix.".temp.insertions.bed";
  open FO ,"> $file" or die $!;
  foreach my $loc(sort keys %clusters){
     $clustnum++;
     my $l = $loc;
     $l =~ s/\t/_/g;
     print FO "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
     print FO "$clustnum\t$clusters{$loc}{final}{TE_mapReads}\t";
     print FO "$clusters{$loc}{final}{strand}\t";
     print FO "$clusters{$loc}{final}{TE}\t$l\t";
     my @w = split("\t",$loc);
     $w[1] = int($w[1]/1000);
     my $valnum=0;
     my $val=0;
     for(my $e=-2;$e<=2;$e++){
        if($cov{$w[0]}{$w[1]+$e}){
           $val += $cov{$w[0]}{$w[1]+$e};
           $valnum++;
        }
     }
     $val = Utilities::round($val/(10*$valnum),2) if($valnum>0);   ## readlength 100 vs bin size 1000, hence division by 10
     print FO "$val\n"; 
  }
  close(FO);
  
  my $testring = join(",",keys %transposons);
  $testring =~ s/PolyA//;
  $testring=~ s/,,/,/;
  $testring=~ s/^,//;
  $testring=~ s/,$//;
  #$testring=~ s/\/\S+?,/,/g;
  #$testring=~ s/\/\S+$//g;
  
  system("Rscript countmodel.R $baseoutprefix $wd $pval_cutoff $testring $refMapqQ $TEMapScore")==0 or die qq[ ERROR while running background enrichment model. Exiting !! \n];     
  
  open(I,"<",$wd."/".$baseoutprefix.".temp.bgModeledInsertions.txt") or die "Can't open bgModeled cluster file\n";
  while(<I>){
     chomp;
     next if(/^(\#)/); 
     next if(/^(\@)/); 
     s/\r//;  
     my @temp = split(/\t/);
     #chr8  113780680 113780692 13  * 1338  9 Alu (8)chr8_113780680_a  Alu_1338  BG=HERVK,0.00e+00,29;Rest_BG=SVA,1.39e-248,41;
     my $loc = $temp[8];
     $loc=~ s/_/\t/g;

     if($clusters{$loc}){
       if($temp[11] =~ "BG=;Rest_BG"){
         $clusters{$loc}{final}{bgModel} = "BG=;Rest_BG=;";
         $clusters{$loc}{final}{RAMCount} = 0;
         $clusters{$loc}{final}{pvalue} = 1;
       }else{
         $clusters{$loc}{final}{bgModel} = $temp[11];
         (my $l)  = $temp[11] =~ m/BG=(.*?);Rest_BG/;
         my @l = split(",",$l);
         $clusters{$loc}{final}{RAMCount} = $l[2];
         $clusters{$loc}{final}{pvalue} = $l[1];
         if($clusters{$loc}{final}{TE} eq "*"){
           $clusters{$loc}{final}{TE} = $l[0];
         }
       }
     }else{
       print " ERROR: Location $temp[8] was not found in the clusters object \n";
       exit 1;
     }
  }
  close(I);
  return(\%clusters); 
}



