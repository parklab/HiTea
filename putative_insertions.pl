#!/usr/bin/perl
# HiTEA

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

my %redb;
$redb{"MboI"}{motif} = "GATC|CTAG";
$redb{"MboI"}{ligmotif} = "GATCGATC|CTAGCTAG";
$redb{"MboI"}{mlen} = "4";
$redb{"MboI"}{ligmotifln} = "8";
$redb{"DpnII"}{motif} = "GATC|CTAG";
$redb{"DpnII"}{ligmotif} = "GATCGATC|CTAGCTAG";
$redb{"DpnII"}{mlen} = "4";
$redb{"DpnII"}{ligmotifln} = "8";
$redb{"HindIII"}{motif} = "[A]{0,1}AGCTT|AAGCT[T]{0,1}|[T]{0,1}TCGAA|TTCGA[A]{0,1}";
$redb{"HindIII"}{ligmotif} = "AAGCTAGCTT|TTCGATCGAA";
$redb{"HindIII"}{mlen} = "5";
$redb{"HindIII"}{ligmotifln} = "10";
$redb{"Arima"}{motif} = "GATC|CTAG|GA[ATGC]{1}T[C]{0,1}|[G]{0,1}A[ATGC]{1}TC|[C]{0,1}T[ATGC]{1}AG|CT[ATGC]{1}A[G]{0,1}";
$redb{"Arima"}{ligmotif} = "GATCGATC|CTAGCTAG|GA[ATGC]{1}TA[ATGC]{1}TC|GATCA[ATGC]{1}TC|GA[ATGC]{1}TGATC|CT[GATC]{1}AT[AGTC]{1}AG|CT[ATGC]{1}ACTAG|CTAGT[ATGC]{1}AG";
$redb{"Arima"}{mlen} = "4";
$redb{"Arima"}{ligmotifln} = "8";
$redb{"NcoI"}{motif} = "[C]{0,1}CATGG|CCATG[G]{0,1}|[G]{0,1}GTACC|GGTAC[C]{0,1}";
$redb{"NcoI"}{ligmotif} = "CCATGCATGG|GGTACGTACC";
$redb{"NcoI"}{mlen} = "5";
$redb{"NcoI"}{ligmotifln} = "10";

my $in="";
my $index ="";
my $polym ="";
my $repbase ="";
my $bgAnnotations ="";
my $wd = "";
my $outprefix="project";
my $sites ="";

my $enzyme ="MboI";
my $ncores=4;
my $refMapqQ = 28;
my $TEMapScore = 30;
my $clip = 20;
my $help=0;
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'index=s'          => \$index,
  'polym=s'          => \$polym,
  'repbase:s'        => \$repbase,
  'bgAnnotations:s'  => \$bgAnnotations,
  'outprefix=s'      => \$outprefix,
  'enzyme=s'          => \$enzyme,
  'q=s'        => \$refMapqQ,
  'clip=s'           => \$clip,
  'sites=s'          => \$sites,
  'wd:s'             => \$wd,
  'ncores:s'         => \$ncores,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl putative_insertions.pl -in [FILE_PATH] -outprefix [STRING] -polym [PATH] -q [refMapqQ] -index [PATH] -repbase [PATH] -bgAnnotations [DIRPATH] -clip [INT] -sites [PATH] -enzyme [STRING] -ncores [INT] -wd [DIR_PATH]\n\n";
  print "This program annotates finalized breaks using the bam file\n\n";
  print "Options:\n\n";
  print "***required:\n";
  print "  -in                    Clusters file in bed format \n";
  print "  -index                 TE consensus index (generated using BWA)\n";
  print "  -repbase               RepBase index (generated using BWA)\n";
  print "  -bgAnnotations         Reference TE annotations\n";
  print "  -q                     reference mapq score threshold \n";
  print "  -sites                 RE enzyme genomic annotations (Juicer format)\n";
  print "  -enzyme                RE enzyme (default:MboI)\n";
  print "***optional:\n";
  print "  -clip                  Clip length cutoff [default: 20] \n";
  print "  -outprefix             Outputprefix for generating 2 output files [default: project] \n";
  print "  -polym                 TE polymorphic example cases index (generated using BWA)\n";
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
print " ERROR: Restriction enzyme index is not provided\n" unless -e $sites;
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


my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "[putative_insertions] START:\t $run_time seconds\n";
print " Command: perl putative_insertions.pl -in $in -sites $sites -index $index -repbase $repbase -bgAnnotations $bgAnnotations -enzyme $enzyme -ncores $ncores -outprefix $outprefix -wd $wd \n";


## get transposons
%transposons = Utilities::get_fasta($index);
print " transposons being considered in the analyses: " , join(",",keys %transposons),"; sizes in bp (",join (",",values %transposons),")\n";
if($polym ne ""){
  %polymorphs = Utilities::get_fasta($polym);
  foreach my $i (keys %polymorphs){
   my @temp = split("~",$i); ## Family~Subfamily
   print " $i\t$temp[1]\t$polymorphs{$i}","bp\n"; 
   $temp[1] =~ s/\n//;
   $temp[1] =~ s/\r//;
   $crossref{$temp[1]}=$temp[0];
   $polymorphs{$temp[1]}=$polymorphs{$i};
   delete $polymorphs{$i};
  }
}
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
my @tes = grep { $_ ne 'PolyA' } (keys %transposons);
foreach my $te (@tes){
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

#my $oufilex=$outprefix.'.completed.clusters'; #Save 
#open FOO,">$oufilex" or die $!;
#print FOO Dumper %clusters;
#close(FOO);
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
    #define locus
    my @temp = split("\t",$loc);
    $clusters{$loc}{pos1}{chr} = $temp[0] ;
    $clusters{$loc}{pos1}{cliploc} = $temp[1];
    $clusters{$loc}{pos1}{strand} = "*";
    $clusters{$loc}{pos1}{clip_sidedness} = $temp[2];
    if(exists $clusters{$loc}{smate} and $temp[2] eq "b"){
      $clusters{$loc}{smate}{clip_sidedness} = "a";
    }elsif(exists $clusters{$loc}{smate} and $temp[2] eq "a"){
      $clusters{$loc}{smate}{clip_sidedness} = "b";
    }


    if(!defined($clusters{$loc}{pos1}{maxscore}) or !exists($clusters{$loc}{pos1}{dist2re}) or !exists($clusters{$loc}{pos1}{seqstring})){
      print qq[ ERROR: Can't find sequence maxscore/dist2re/seqstring for: $loc\n];
      print Dumper $clusters{$loc};
      exit 1;
    }
    ## counts
    my @features = ("fuzzyclipped");  
    foreach my $feature (@features){
      if(exists($clusters{$loc}{pos1}{$feature})){
         $clusters{$loc}{pos1}{$feature} = Utilities::unique($clusters{$loc}{pos1}{$feature}, ",","true"); 
       }
    }
    undef(@features);
    @features = ("clipseqs"); 
    foreach my $feature (@features){
      if(exists($clusters{$loc}{pos1}{$feature})){
         $clusters{$loc}{pos1}{$feature} = Utilities::unique($clusters{$loc}{pos1}{$feature},",","false");
      }
      if(exists $clusters{$loc}{smate} and defined $clusters{$loc}{smate}{$feature}){
         $clusters{$loc}{smate}{$feature} = Utilities::unique($clusters{$loc}{smate}{$feature},",","false");
      }
    }
    undef(@features);
    ##Unmapped cluster mate
    if(exists $clusters{$loc}{smate} and defined $clusters{$loc}{smate}{smate}){
      $clusters{$loc}{smate}{smate}=~ s/,$//;
      my @l = split(",",Utilities::unique($clusters{$loc}{smate}{smate},",","false") );
      my %fq;
      my %uq;
      foreach my $i (@l) {
        my $j = $i;
        $i=~ s/\|.*$//;
        $fq{$i}++;
        $j=~ s/^.*\|//;
        $uq{$j}++;
      }
      my @k = keys %fq;
      @k = sort{$fq{$b} <=> $fq{$a}} keys %fq if(scalar @k>1);
      $clusters{$loc}{smate}{smate} = $k[0]."|".$fq{$k[0]}; ## most represented location, total clip reads
      $clusters{$loc}{smate}{fuzzyclipped} = scalar (keys %uq) - $fq{$k[0]};
      ($clusters{$loc}{smate}{cliploc}) = $k[0]=~ m/\S*_(\d*)_\S*_\S*/; #chr10_68814835_a_*
    }
    
    ## reciprocal clusters
    @features = keys %transposons; 
    if(scalar keys %polymorphs > 0 ){
      push(@features,keys %polymorphs);
    }
    foreach my $feature (@features){
      if(defined($clusters{$loc}{pos1}{$feature} )){
        $clusters{$loc}{pos1}{$feature}{TESc} =~ s/,$//;
        $clusters{$loc}{pos1}{$feature}{TESc} = Utilities::max( split(",",$clusters{$loc}{pos1}{$feature}{TESc}) );
        my @array = @{Utilities::get_reciprocalClust($clusters{$loc}{pos1}{$feature}{TEClust}, $feature)}; ## position,freq
        $clusters{$loc}{pos1}{$feature}{TEClustPos} = $array[0];
        $clusters{$loc}{pos1}{$feature}{TEClustReads} = $array[1];
        $clusters{$loc}{pos1}{$feature}{TEClustPerc} = $array[2];
        $clusters{$loc}{pos1}{$feature}{TEClust} = $array[0].":".$array[1];
        $clusters{$loc}{pos1}{$feature}{mapq} = $clusters{$loc}{pos1}{mapq};
      }
      if(exists $clusters{$loc}{smate} and defined($clusters{$loc}{smate}{$feature} )){
        $clusters{$loc}{smate}{$feature}{TESc} =~ s/,$//;
        $clusters{$loc}{smate}{$feature}{TESc} = Utilities::max( split(",",$clusters{$loc}{smate}{$feature}{TESc}) );
        my @array = @{Utilities::get_reciprocalClust($clusters{$loc}{smate}{$feature}{TEClust}, $feature)}; ## position,freq
        $clusters{$loc}{smate}{$feature}{TEClustPos} = $array[0];
        $clusters{$loc}{smate}{$feature}{TEClustReads} = $array[1];
        $clusters{$loc}{smate}{$feature}{TEClustPerc} = $array[2];
        $clusters{$loc}{smate}{$feature}{TEClust} = $array[0].":".$array[1];
        $clusters{$loc}{smate}{$feature}{mapq} = $clusters{$loc}{smate}{mapq};
      }
    }
    undef(@features);

    ## if polymoph mapping is allowed, 
    if(scalar keys %polymorphs > 0 ){
     @features = keys %crossref;  
     my %chq;
     foreach my $feature (@features){
       if(defined($clusters{$loc}{pos1}{$feature} )){
         if(!exists $chq{pos} or (exists $chq{pos} and $chq{pos}{clusters}<$clusters{$loc}{pos1}{$feature}{TEClustReads})){
          $chq{pos}{te} = $feature;
          $chq{pos}{clusters} = $clusters{$loc}{pos1}{$feature}{TEClustReads};
         }
       }
       if(exists $clusters{$loc}{smate} and defined($clusters{$loc}{smate}{$feature} )){
         if(!exists $chq{smate} or (exists $chq{smate} and $chq{smate}{clusters}<$clusters{$loc}{smate}{$feature}{TEClustReads})){
          $chq{smate}{te} = $feature;
          $chq{smate}{clusters} = $clusters{$loc}{smate}{$feature}{TEClustReads};
         }
       }
     }
     foreach my $feature (@features){
       if(defined($clusters{$loc}{pos1}{$feature} )){
        if($feature ne $chq{pos}{te}){
          delete $clusters{$loc}{pos1}{$feature};
        }
       }
       if(exists $clusters{$loc}{smate} and defined($clusters{$loc}{smate}{$feature} )){
         if($feature ne $chq{smate}{te}){
          delete $clusters{$loc}{smate}{$feature};
         }
       }
     }
    }
    
    ## only keep those clusters where TE-consensus mapping are not present 
    if(scalar keys %polymorphs > 0){
      my @features = grep { $_ ne 'PolyA' } (keys %transposons);
      my $cnt=0;
      foreach my $feature (@features){
        $cnt++ if(defined($clusters{$loc}{pos1}{$feature} ));
        $cnt++ if(exists $clusters{$loc}{smate} and defined($clusters{$loc}{smate}{$feature} ));   
      }
      if($cnt>0){
        my @features = keys %polymorphs;  
        foreach my $feature (@features){
          if(defined($clusters{$loc}{pos1}{$feature} )){
            #delete $clusters{$loc}{pos1}{$feature};
          }
          if(exists $clusters{$loc}{smate} and defined($clusters{$loc}{smate}{$feature} )){
            #delete $clusters{$loc}{smate}{$feature};
          }
        }
      }  
    }

    ## print unmapped clipeseqs
    if(exists $clusters{$loc}{smate} and defined $clusters{$loc}{smate}{clipseqs}){
      my $l = $loc;
      $l =~ s/\t/_/g;
      my $cnt=0;
      $clusters{$loc}{smate}{clipseqs} =~ s/,$//;
      my @t = split(",",$clusters{$loc}{smate}{clipseqs});
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
      if(exists $clusters{$j}{smate}{AllClip}){
        %{$new{$j}{pos2}} = %{$clusters{$j}{smate}};
        $new{$j}{pos2}{maxscore} =0;
        $new{$j}{pos2}{strand} ="*";
      } 
      delete $clusters{$j};
      $new{$j}{pos1}{bgclusters} = $size;
    }elsif($size == 2){
      %{$new{$temp[0]}{pos1}} = %{$clusters{$temp[0]}{pos1}};
      %{$new{$temp[0]}{pos2}} = %{$clusters{$temp[1]}{pos1}};
      $new{$temp[0]}{pos1}{bgclusters} = $size;
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
     ## (1) locus,mapq
     $clusters{$loc}{final}{chr} = $clusters{$loc}{pos1}{chr};
     if(exists $clusters{$loc}{pos2}){
       $clusters{$loc}{final}{start} = Utilities::min(($clusters{$loc}{pos1}{cliploc},$clusters{$loc}{pos2}{cliploc})) ; 
       $clusters{$loc}{final}{end}   = Utilities::max(($clusters{$loc}{pos1}{cliploc},$clusters{$loc}{pos2}{cliploc}));     
       $clusters{$loc}{final}{dist2re} = Utilities::max(($clusters{$loc}{pos1}{dist2re},$clusters{$loc}{pos2}{dist2re}));
       $clusters{$loc}{final}{maxscore} = Utilities::max(($clusters{$loc}{pos1}{maxscore},$clusters{$loc}{pos2}{maxscore}));       
     }else{
       $clusters{$loc}{final}{start} = $clusters{$loc}{pos1}{cliploc}; 
       $clusters{$loc}{final}{end}   = $clusters{$loc}{pos1}{cliploc};     
       $clusters{$loc}{final}{dist2re} = $clusters{$loc}{pos1}{dist2re};
       $clusters{$loc}{final}{maxscore} = $clusters{$loc}{pos1}{maxscore};
     }
     
     ## (2) PolyA tails?
     $clusters{$loc}{final}{isPolyA}="F";   
     $clusters{$loc}{pos1}{isPolyA}="F";
     my $a = Utilities::getstrand($clusters{$loc}{pos1}{seqstring});
     $clusters{$loc}{pos1}{isPolyA}="T" if($a ne "NA" or exists $clusters{$loc}{pos1}{PolyA});
     $clusters{$loc}{pos1}{strand} = "+" if($a eq "PolyA");
     $clusters{$loc}{pos1}{strand} = "-" if($a eq "PolyT");
     $clusters{$loc}{final}{tailinfo} = $a; 
     $clusters{$loc}{final}{isPolyA}="T" if($a ne "NA" or exists $clusters{$loc}{pos1}{PolyA});
     if(exists $clusters{$loc}{pos2}){
       $clusters{$loc}{pos2}{isPolyA}="F";
       my $b = Utilities::getstrand($clusters{$loc}{pos2}{seqstring});
       $clusters{$loc}{pos2}{isPolyA}="T" if($b ne "NA" or exists $clusters{$loc}{pos2}{PolyA});
       $clusters{$loc}{pos2}{strand} = "+" if($b eq "PolyA");
       $clusters{$loc}{pos2}{strand} = "-" if($b eq "PolyT"); 
       $clusters{$loc}{final}{isPolyA}="T" if($b ne "NA" or exists $clusters{$loc}{pos2}{PolyA});
       $clusters{$loc}{final}{tailinfo} .= $b.",";
     }
     $clusters{$loc}{final}{tailinfo}=~ s/NA//;
     $clusters{$loc}{final}{strand} = "*"; 
     if($clusters{$loc}{final}{tailinfo} eq "PolyA"){
       $clusters{$loc}{final}{strand} = "+"; 
      }elsif($clusters{$loc}{final}{tailinfo} eq "PolyT"){
       $clusters{$loc}{final}{strand} = "-"; 
     }
     
     ## (2.1) if smate is polyA without any TE assignment, assign it to PolyA
     my @features = keys %transposons;
     if(scalar keys %polymorphs > 0){
       push(@features,keys %polymorphs);
     }
     $a=0;    
     foreach my $feature (@features){
       $a++ if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature});
     }
     if($a==0 and exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{smate} and $clusters{$loc}{pos2}{isPolyA} eq "T"){
       my @xx = ($clusters{$loc}{pos2}{smate} =~ m/\S*_\S*_\S*_(\S*)\|(\d*)/);
       if($xx[0] eq "*"){
         $clusters{$loc}{pos2}{PolyA}{TEClustPos} = "1";
         $clusters{$loc}{pos2}{PolyA}{TEClustReads} = $xx[1];
         $clusters{$loc}{pos2}{PolyA}{TEClustPerc} = "1.00";
         $clusters{$loc}{pos2}{PolyA}{TEClust} = "1:".$xx[1];
         $clusters{$loc}{pos2}{PolyA}{TESc} = "30";
         $clusters{$loc}{pos2}{PolyA}{num} = $xx[1];
         $clusters{$loc}{pos2}{PolyA}{mapq} = $clusters{$loc}{pos2}{mapq};
       }    
     }
     undef(@features);
     undef $a;
     
     ## (2.2) if clip read maps to PolyA or has an inferred PolyA mapping, assign the Mapscore to minimum 30 
      ## PolyA mapping is inferential and hence should not decide cluster ommission 
     $clusters{$loc}{pos1}{PolyA}{TESc} = 30 if(exists $clusters{$loc}{pos1}{PolyA} and $clusters{$loc}{pos1}{PolyA}{TESc}<30);
     if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{PolyA} and $clusters{$loc}{pos2}{PolyA}{TESc}<30){
      $clusters{$loc}{pos2}{PolyA}{TESc} = 30;
     }

     ## (3) TSD
     my @l = split(",",($clusters{$loc}{pos1}{seqstring}.",".$clusters{$loc}{final}{start}.",".$clusters{$loc}{final}{end}));
     my $tsd = Utilities::tsd(\@l);
     if($tsd eq "" or $tsd eq "ERROR-TSD" and exists $clusters{$loc}{pos2}){
         @l = split(",",($clusters{$loc}{pos2}{seqstring}.",".$clusters{$loc}{final}{start}.",".$clusters{$loc}{final}{end}));
         $tsd = Utilities::tsd(\@l);
     }
     if($tsd eq "" or $tsd eq "ERROR-TSD"){
         $clusters{$loc}{final}{TSD} = "*"; 
      }else{
        $clusters{$loc}{final}{TSD} = $tsd;
     }

     ## (4) distance to annotated RE motif
     my $dist2remotif="NA";
     if(exists $clusters{$loc}{pos1}){
       my $x = 0;
       my $a = &binarysearch($clusters{$loc}{final}{start}, $chromosomes{$clusters{$loc}{final}{chr}});
       if($a <=$clusters{$loc}{final}{start} and ($a+($redb{$enzyme}{mlen})) >= $clusters{$loc}{final}{start}){
         $x=0;
       }elsif($a < $clusters{$loc}{final}{start}){
          $x = abs($a-$clusters{$loc}{final}{start}) - ($redb{$enzyme}{mlen});
       }else{
          $x = abs($a-$clusters{$loc}{final}{start});
       }
       my $y = 0;
       my $b = &binarysearch($clusters{$loc}{final}{end}, $chromosomes{$clusters{$loc}{final}{chr}});
       if($b <=$clusters{$loc}{final}{end} and ($b+($redb{$enzyme}{mlen})) >= $clusters{$loc}{final}{end}){
         $y=0;
       }elsif($a < $clusters{$loc}{final}{end}){
          $y = abs($a-$clusters{$loc}{final}{end}) - $redb{$enzyme}{mlen} ;
       }else{
          $y = abs($a-$clusters{$loc}{final}{end});
       }
       $dist2remotif= Utilities::min( ($x,$y) );
     }
     #$clusters{$loc}{final}{dist2remotif} = $dist2remotif;
     $clusters{$loc}{final}{dist2remotif} = 100;

     ## (5) All clip read counts
     if(exists $clusters{$loc}{pos2}){
      $clusters{$loc}{final}{commonClips} = Utilities::common($clusters{$loc}{pos1}{AllClip},$clusters{$loc}{pos2}{AllClip},";","true");
     }else{
       $clusters{$loc}{final}{commonClips} =0;
     }     

     ## (6) All clip read counts
     undef $a;
     $a = $clusters{$loc}{pos1}{AllClip};
     $a .= $clusters{$loc}{pos2}{AllClip} if(exists $clusters{$loc}{pos2});
     $a =~ s/;$//;
     $clusters{$loc}{final}{AllClip} = Utilities::unique($a, ";","true");
     $clusters{$loc}{pos1}{AllClip} =~ s/;$//; 
     $clusters{$loc}{pos1}{AllClip} = Utilities::unique($clusters{$loc}{pos1}{AllClip}, ";","true");
     if(exists $clusters{$loc}{pos2}){
       $clusters{$loc}{pos2}{AllClip} =~ s/;$//; 
       $clusters{$loc}{pos2}{AllClip} = Utilities::unique($clusters{$loc}{pos2}{AllClip}, ";","true");
     }

     ## (7) Finalize TE 
     @features = grep { $_ ne 'PolyA' } (keys %transposons);
     if(scalar keys %polymorphs > 0){
       push(@features,keys %polymorphs);
     }
     foreach my $feature (@features){
       my %chq;
       $chq{lte}=""; $chq{rte}="";
       $chq{lpoly}=""; $chq{rpoly}="";
       $chq{TESc}=0; $chq{TEClustReads}=0;
       $chq{lpos}="-"; $chq{rpos}="-";
       $chq{TEClustPerc}=0;
       $chq{tmapq}=0;$chq{pmapq}=0;
       if(exists $clusters{$loc}{pos1}{$feature}){
         $chq{lte} .= $clusters{$loc}{pos1}{$feature}{num}.",";
         $chq{TESc}=$clusters{$loc}{pos1}{$feature}{TESc};
         $chq{lpos}=$clusters{$loc}{pos1}{$feature}{TEClustPos};
         $chq{TEClustReads}=$clusters{$loc}{pos1}{$feature}{TEClustReads}; 
         $chq{tmapq}=$clusters{$loc}{pos1}{$feature}{mapq};
         #$chq{TEClustReads}+=$clusters{$loc}{pos1}{PolyA}{TEClustReads} if(exists$clusters{$loc}{pos1}{PolyA});
       }
       if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature}){
         $chq{rte} .= $clusters{$loc}{pos2}{$feature}{num}.",";
         $chq{TESc}=$clusters{$loc}{pos2}{$feature}{TESc} if($clusters{$loc}{pos2}{$feature}{TESc}> $chq{TESc});
         $chq{rpos}=$clusters{$loc}{pos2}{$feature}{TEClustPos};
         $chq{TEClustReads}+=$clusters{$loc}{pos2}{$feature}{TEClustReads}; 
         $chq{tmapq}=$clusters{$loc}{pos2}{$feature}{mapq}  if($clusters{$loc}{pos2}{$feature}{mapq}> $chq{tmapq});
        #$chq{TEClustReads}+=$clusters{$loc}{pos2}{PolyA}{TEClustReads} if(exists$clusters{$loc}{pos2}{PolyA});
       }
       if(exists$clusters{$loc}{pos1}{PolyA}){
         $chq{lpoly} .= $clusters{$loc}{pos1}{PolyA}{num}."," ;
         $chq{pmapq}=$clusters{$loc}{pos1}{PolyA}{mapq};
       }
       if(exists$clusters{$loc}{pos2} and exists$clusters{$loc}{pos2}{PolyA}){
        $chq{rpoly} .= $clusters{$loc}{pos2}{PolyA}{num}.",";
        $chq{pmapq}=$clusters{$loc}{pos2}{PolyA}{mapq} if($clusters{$loc}{pos2}{PolyA}{mapq} > $chq{pmapq});
       }
       my $a = $chq{lte}.$chq{rte}.$chq{lpoly}.$chq{rpoly};
       $a=~s/,$//;
       $a = scalar(split(",",$a));
       $chq{tnum} = Utilities::unique($chq{lte}.$chq{rte}, ",","true");
       $chq{anum} = Utilities::unique($chq{lte}.$chq{rte}.$chq{lpoly}.$chq{rpoly}, ",","true");
       $chq{pnum} = $chq{anum} -$chq{tnum};
       $chq{umap} = $clusters{$loc}{final}{AllClip} - $chq{anum};
       if(defined $chq{TEClustReads} and defined $chq{pnum}){
         $chq{TEClustReads}+=$chq{pnum};
       }
       $chq{TEClustPerc}=Utilities::round($chq{TEClustReads}/$chq{anum},2) if($chq{anum}>0);
       if($chq{TEClustPerc}>0 and $chq{tnum}>0 and ($chq{tmapq}>$refMapqQ or $chq{pmapq}>$refMapqQ )){
        %{$clusters{$loc}{final}{$feature}} = %chq;
       }
     }
     undef @features;

     ## (7A) If polymophic mapping is allowed, rename them to Consensus
     if(scalar keys %polymorphs > 0){
      my %inv;
      foreach my $ele (keys %crossref){ $inv{$crossref{$ele}} .= $ele.","; }
      foreach my $ele (keys %inv){
        $inv{$ele} =~ s/,$//;
        my @features = split(",",$inv{$ele});
        my %chq;
        $chq{lte}=""; $chq{rte}="";
        $chq{lpoly}=""; $chq{rpoly}="";
        $chq{TESc}=0; $chq{TEClustReads}=0;
        $chq{lpos}="-"; $chq{rpos}="-";
        $chq{TEClustPerc}=0;
        $chq{tmapq}=0;$chq{pmapq}=0;
        foreach my $feature(@features){
          if(exists$clusters{$loc}{final}{$feature} and !defined$clusters{$loc}{final}{$ele}){
            $chq{lte} .= $clusters{$loc}{final}{$feature}{lte}."," if($clusters{$loc}{final}{$feature}{lte} ne "");
            $chq{rte} .= $clusters{$loc}{final}{$feature}{rte}."," if($clusters{$loc}{final}{$feature}{rte} ne "");
            
            if($clusters{$loc}{final}{$feature}{rpoly} eq "" and exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{PolyA}){
              $chq{rpoly} .= $clusters{$loc}{pos2}{PolyA}{num}.",";
              $chq{pmapq} = $clusters{$loc}{pos2}{PolyA}{mapq} if($clusters{$loc}{pos2}{PolyA}{mapq}>$chq{pmapq});
            }elsif($clusters{$loc}{final}{$feature}{rpoly} ne ""){
              $chq{rpoly} .= $clusters{$loc}{final}{$feature}{rpoly}.",";
            }
            if($clusters{$loc}{final}{$feature}{lpoly} eq "" and defined $clusters{$loc}{pos1}{PolyA}){
              $chq{lpoly} .= $clusters{$loc}{pos1}{PolyA}{num}.",";
              $chq{pmapq} = $clusters{$loc}{pos1}{PolyA}{mapq} if($clusters{$loc}{pos1}{PolyA}{mapq}>$chq{pmapq});
            }elsif($clusters{$loc}{final}{$feature}{lpoly} ne "" ){
              $chq{lpoly} .= $clusters{$loc}{final}{$feature}{lpoly}.",";
            }
            $chq{TESc} = $clusters{$loc}{final}{$feature}{TESc} if($clusters{$loc}{final}{$feature}{TESc}>$chq{TESc});
            $chq{lpos} = $clusters{$loc}{final}{$feature}{lpos} if($clusters{$loc}{final}{$feature}{lpos} ne "-");
            $chq{rpos} = $clusters{$loc}{final}{$feature}{rpos} if($clusters{$loc}{final}{$feature}{rpos} ne "-");
            $chq{TEClustReads} += $clusters{$loc}{final}{$feature}{TEClustReads};            
            
            if(exists$clusters{$loc}{pos1}{$feature} and $clusters{$loc}{pos1}{$feature}{mapq}> $chq{tmapq}){
              $chq{tmapq}=$clusters{$loc}{pos1}{$feature}{mapq};
            }
            if(exists$clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature} and $clusters{$loc}{pos2}{$feature}{mapq}> $chq{tmapq}){
              $chq{tmapq}=$clusters{$loc}{pos2}{$feature}{mapq};       
            }
            delete $clusters{$loc}{final}{$feature};
          }
        }
        my $a = $chq{lte}.$chq{rte}.$chq{lpoly}.$chq{rpoly};
        $a=~s/,$//;
        $a = scalar(split(",",$a));
        $chq{tnum} = Utilities::unique($chq{lte}.$chq{rte}, ",","true");
        $chq{anum} = Utilities::unique($chq{lte}.$chq{rte}.$chq{lpoly}.$chq{rpoly}, ",","true");
        $chq{pnum} = $chq{anum} -$chq{tnum};
        $chq{umap} = $clusters{$loc}{final}{AllClip} - $chq{anum};
        $chq{TEClustPerc}=Utilities::round($chq{TEClustReads}/$chq{anum},2) if($chq{anum}>0);
        if($chq{TEClustPerc}>0 and $chq{tnum}>0 and ($chq{tmapq}>$refMapqQ or $chq{pmapq}>$refMapqQ)) {
         %{$clusters{$loc}{final}{$ele}} = %chq;
        }
      }
     }

     ### (8) set the number counts
     @features = keys %transposons;
     if(scalar keys %polymorphs > 0){
       push(@features,keys %polymorphs);
     }
     foreach my $feature (@features){
      if(exists $clusters{$loc}{pos1}{$feature}){
        $clusters{$loc}{pos1}{$feature}{num} =~ s/,$//;
        $clusters{$loc}{pos1}{$feature}{num} = Utilities::unique($clusters{$loc}{pos1}{$feature}{num}, ",","true");           
      }
      if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature}){
        $clusters{$loc}{pos2}{$feature}{num} =~ s/,$//;
        $clusters{$loc}{pos2}{$feature}{num} = Utilities::unique($clusters{$loc}{pos2}{$feature}{num}, ",","true");           
      }
      if(exists $clusters{$loc}{final}{$feature}){
        $clusters{$loc}{final}{$feature}{lte} = Utilities::unique($clusters{$loc}{final}{$feature}{lte}, ",","true");          
        $clusters{$loc}{final}{$feature}{rte} = Utilities::unique($clusters{$loc}{final}{$feature}{rte}, ",","true");          
        $clusters{$loc}{final}{$feature}{lpoly} = Utilities::unique($clusters{$loc}{final}{$feature}{lpoly}, ",","true");          
        $clusters{$loc}{final}{$feature}{rpoly} = Utilities::unique($clusters{$loc}{final}{$feature}{rpoly}, ",","true");          
      }
     }
     
     ## repurpose the reads with PolyA tail but without TE consensus mapping to aid the discovery process
     $a=0;
     foreach my $feature (@features){
      $a++ if(exists $clusters{$loc}{pos2} and $clusters{$loc}{pos2}{$feature});
     }
     undef @features;
     @features = grep { $_ ne 'PolyA' } (keys %transposons);
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature} and exists $clusters{$loc}{pos2} and $a==0 and $clusters{$loc}{pos2}{isPolyA} eq "T"){
       $clusters{$loc}{final}{$feature}{rpoly} += $clusters{$loc}{pos2}{AllClip};
       $clusters{$loc}{final}{$feature}{anum} += $clusters{$loc}{pos2}{AllClip};
       $clusters{$loc}{final}{$feature}{TEClustReads} += $clusters{$loc}{pos2}{AllClip};
       $clusters{$loc}{final}{$feature}{TEClustPerc} = Utilities::round($clusters{$loc}{final}{$feature}{TEClustReads}/$clusters{$loc}{final}{$feature}{anum},2);       
      }
     }
     undef @features;
     
     ## (9) Purity ratio
     @features = grep { $_ ne 'PolyA' } (keys %transposons);
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature}){
       $clusters{$loc}{final}{$feature}{Purity} = Utilities::round($clusters{$loc}{final}{$feature}{anum}/$clusters{$loc}{final}{AllClip},2);
      }
     }
     
     ## (10) insert size
     @features = grep { $_ ne 'PolyA' } (keys %transposons);
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature}){
        $clusters{$loc}{final}{$feature}{insert_size} = "NA";
        my $min = Utilities::min(($clusters{$loc}{final}{$feature}{lpos},$clusters{$loc}{final}{$feature}{rpos}));
        my $max = Utilities::min(($clusters{$loc}{final}{$feature}{lpos},$clusters{$loc}{final}{$feature}{rpos}));
        if($clusters{$loc}{final}{isPolyA} eq "T" and defined $min and $min ne "" and exists $clusters{$loc}{pos2}){
          $clusters{$loc}{final}{$feature}{insert_size} = ($transposons{$feature} - $min);
        }elsif($clusters{$loc}{final}{isPolyA} eq "F" and defined $min and $min ne "" and defined $max and $max ne ""){
          $clusters{$loc}{final}{$feature}{insert_size} = ($max - $min);   
        }
        if($clusters{$loc}{final}{$feature}{insert_size} ne "NA" and $clusters{$loc}{final}{$feature}{insert_size} <1){
         $clusters{$loc}{final}{$feature}{insert_size} = 'NA';
        }
      }
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
     my $l = $loc;
     $l=~ s/\t/_/g;
     $clusters{$loc}{final}{LHS_cons} = Utilities::get_longestseq($clusters{$loc}{pos1}{clipseqs});
     print FO1 ">$l","_RC\n$clusters{$loc}{final}{LHS_cons}\n";
     if(exists $clusters{$loc}{pos2}){
      $clusters{$loc}{final}{RHS_cons} = Utilities::get_longestseq($clusters{$loc}{pos2}{clipseqs});  
      print FO1 ">$l","_LC\n$clusters{$loc}{final}{RHS_cons}\n";   
     }  
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
  
       foreach my $feature (keys %transposons){
         if($feature ne "PolyA" and exists $clusters{$loc}{final}{$feature}){
           if(defined $clusters{$loc}{final}{$feature}{subfamily}){
             if($rat >0.5 and $rat > $qw and $temp[2]=~ m/$feature/){
                $clusters{$loc}{final}{$feature}{subfamily} = $temp[2];
             }
           }elsif($rat >0.5 and $temp[2]=~ m/$feature/){
                $clusters{$loc}{final}{$feature}{subfamily} = $temp[2];
                $qw = $rat;
            }else{
                $clusters{$loc}{final}{$feature}{subfamily} = "*";
                $qw = 0.4;
            }
         }
       }
     }
     close(IN);
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
    if($te ne "PolyA" and exists $clusters{$loc}{final}{$te}){
      $cnt++;
      print FO "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
      print FO ".\t.\t$clusters{$loc}{final}{strand}\t";
      $loc =~ s/\t/_/g;
      print FO "$loc\n";
    }
  }
  close(FO);

  my $r=0;
  $r = Utilities::round($cnt*100/scalar(keys %clusters),2) if(scalar(keys %clusters) > 0);   
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
      $clusters{$loc}{final}{$te}{bgTE} = $temp[10];
      $clusters{$loc}{final}{$te}{bgTEDist} = $temp[13];
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
  my @features = grep { $_ ne 'PolyA' } (keys %transposons);
  foreach my $loc(sort keys %clusters){
    foreach my $feature(@features){
      if(exists $clusters{$loc}{final}{$feature}){
       $clustnum++;
       my $l = $loc;
       $l =~ s/\t/_/g;
       print FO "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
       print FO "$clustnum\t$clusters{$loc}{final}{$feature}{anum}\t";
       print FO "$clusters{$loc}{final}{strand}\t";
       print FO "$feature\t$l\t";
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
    }
  }
  close(FO);
  
  my $testring = join(",",@features);
  system("Rscript countmodel.R $baseoutprefix $wd $pval_cutoff $testring $refMapqQ $TEMapScore")==0 or die qq[ ERROR while running background enrichment model. Exiting !! \n];     
  
  open(I,"<",$wd."/".$baseoutprefix.".temp.bgModeledInsertions.txt") or die "Can't open bgModeled cluster file\n";
  while(<I>){
     chomp;
     next if(/^(\#)/); 
     next if(/^(\@)/); 
     s/\r//;  
     my @temp = split(/\t/);
     #(0)chr8  113780680 113780692 13  * 1338  9 Alu (8)chr8_113780680_a  (9)28.5 (10)Alu_1338 (11)fzcount (12)BG=HERVK,0.00e+00,29;Rest_BG=SVA,1.39e-248,41;
     my $loc = $temp[8];
     $loc=~ s/_/\t/g;

     if($clusters{$loc}){
       if($temp[12] =~ "BG=;Rest_BG"){
         $clusters{$loc}{final}{$temp[7]}{bgModel} = "BG=;Rest_BG=;";
         $clusters{$loc}{final}{$temp[7]}{RAMCount} = 0;
         $clusters{$loc}{final}{$temp[7]}{pvalue} = 1;
         $clusters{$loc}{final}{$temp[7]}{AllFuzzy} = $temp[11];
       }else{
         $clusters{$loc}{final}{$temp[7]}{bgModel} = $temp[12];
         (my $l)  = $temp[12] =~ m/BG=(.*?);Rest_BG/;
         my @l = split(",",$l);
         $clusters{$loc}{final}{$temp[7]}{RAMCount} = $l[2];
         $clusters{$loc}{final}{$temp[7]}{pvalue} = $l[1];
         $clusters{$loc}{final}{$temp[7]}{AllFuzzy} = $temp[11];
       }
     }else{
       print " ERROR: Location $temp[8] was not found in the clusters object \n";
       exit 1;
     }
  }
  close(I);
  return(\%clusters); 
}
