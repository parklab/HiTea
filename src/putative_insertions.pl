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
use src::Utilities;
use open qw(:std :utf8);
require "src/vars.pl";
our (%redb);
our (%chrs);
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
BEGIN { our $start_run = time(); }

my $in="";
my $index ="";
my $repbase ="";
my $bgAnnotations ="";
my $enzyme ="";
my $refMapqQ = 28;
my $polym ="";
my $wd = "";
my $outprefix="project";
my $ncores=4;
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
  'enzyme=s'         => \$enzyme,
  'q=s'              => \$refMapqQ,
  'clip=s'           => \$clip,
  'wd:s'             => \$wd,
  'ncores:s'         => \$ncores,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
sub help{
  my $j = shift;
  if($j){
   print "\nUsage: perl putative_insertions.pl -in [FILE_PATH] -outprefix [STRING] -polym [PATH] -q [refMapqQ] -index [PATH] -repbase [PATH] -bgAnnotations [DIRPATH] -clip [INT] -enzyme [STRING] -ncores [INT] -wd [DIR_PATH]\n\n";
   print "This program annotates finalized breaks using the bam file\n\n";
   print "Options:\n\n";
   print "***required:\n";
   print "  -in                    Clusters file in bed format \n";
   print "  -index                 TE consensus index (generated using BWA)\n";
   print "  -repbase               RepBase index (generated using BWA)\n";
   print "  -bgAnnotations         Reference TE annotations\n";
   print "  -q                     reference mapq score threshold \n";
   print "  -enzyme                RE enzyme \n";
   print "***optional:\n";
   print "  -clip                  Clip length cutoff (default: 20) \n";
   print "  -outprefix             Outputprefix for generating 2 output files (default: project) \n";
   print "  -polym                 TE polymorphic example cases index (generated using BWA)\n";
   print "  -wd                    Working directory (default: ~)\n";
   print "  -ncores                Number of threads while reading bam input (default: 4)\n";
   print "  -help|-h               Display usage information.\n\n\n";
   exit 1;
  }
}

$wd =~ s/\/$//;
my $baseoutprefix =$outprefix;
$outprefix = $wd."/".$outprefix if($wd ne "");
if($help or $repbase eq "" or $index eq "" or $bgAnnotations eq "" or $enzyme eq ""){
  print " One or more inputs are not recognized***\n";
  help(1);
}
print " ERROR: countmodel.R does not exist in the same directory!\n" unless -e "countmodel.R";
print " ERROR: pairLocations.R does not exist in the same directory!\n" unless -e "pairLocations.R";


#--------------------------------------------------------------------------------------------------------------
# I/O
#--------------------------------------------------------------------------------------------------------------
my %cov;
my $pval_cutoff = 0.05;
my %transposons;
my %polymorphs;
my %crossref;
my %clusters;
my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "[putative_insertions] START:\t $run_time seconds\n";
print " Command: perl putative_insertions.pl -in $in -index $index -polym $polym -repbase $repbase -bgAnnotations $bgAnnotations -enzyme $enzyme -ncores $ncores -outprefix $baseoutprefix -wd $wd \n";


## get transposons
%transposons = src::Utilities::get_fasta($index);
print " transposons being considered in the analyses: " , join(",",keys %transposons),"; sizes in bp (",join (",",values %transposons),")\n";
if($polym ne ""){
  %polymorphs = src::Utilities::get_fasta($polym);
  foreach my $i (keys %polymorphs){
   my @temp = split("~",$i); ## Family~Subfamily
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

## retrieve coverage information
%cov = %{retrieve($outprefix.".coverage.ph")};

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

exit 0;

#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub complete_clust_object {
  my ($clust) = shift;
  my %clusters = %{$clust};
   
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
    if(exists $clusters{$loc}{smate} and !exists $clusters{$loc}{smate}{identClip}){
      $clusters{$loc}{smate}{identClip} = "NA";
    }
    if(exists $clusters{$loc}{smate} and !exists $clusters{$loc}{smate}{OriRatio}){
      $clusters{$loc}{smate}{OriRatio} = "NA";
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
         $clusters{$loc}{pos1}{$feature} = src::Utilities::unique($clusters{$loc}{pos1}{$feature}, ",","true"); 
       }
    }
    undef(@features);
    @features = ("clipseqs"); 
    foreach my $feature (@features){
      if(exists($clusters{$loc}{pos1}{$feature})){
         $clusters{$loc}{pos1}{$feature} = src::Utilities::unique($clusters{$loc}{pos1}{$feature},",","false");
      }
      if(exists $clusters{$loc}{smate} and defined $clusters{$loc}{smate}{$feature}){
         $clusters{$loc}{smate}{$feature} = src::Utilities::unique($clusters{$loc}{smate}{$feature},",","false");
      }
    }
    undef(@features);
    ##Unmapped cluster mate
    if(exists $clusters{$loc}{smate} and defined $clusters{$loc}{smate}{smate}){
      $clusters{$loc}{smate}{smate}=~ s/,$//;
      my @l = split(",",src::Utilities::unique($clusters{$loc}{smate}{smate},",","false") );
      my %fq;
      my @knownte;
      foreach my $i (@l) {
        my $j = $i;
        $i=~ s/\|.*$//; ## locus
        my $k = $i;
        $k =~ s/.*_//; ## TE
        ($i) = $i=~ m/(\S*_\d*_\S*)_\S*/; #chr10_68814835_a_*
        $j=~ s/^.*\|//; ## readid
        $fq{$i}{num}++;
        if($k ne "*"){
          $fq{$i}{te}++;
          $fq{$i}{teO} = $k;
        }
        if($k eq "*" and !defined $fq{$i}{te}){
          $fq{$i}{te} =0;
          $fq{$i}{teO} = "*";
        }       
        $fq{$i}{teel}=$k;
        $fq{$i}{id}.= $j.";";
        if(exists $clusters{$loc}{smate}{$k}){
          push @knownte,$i;
        }
        ($fq{$i}{side}) = $i=~ m/\S*_\d*_(\S*)/; #chr10_68814835_a_*
        ($fq{$i}{loc}) = $i=~ m/\S*_(\d*)_\S*/; #chr10_68814835_a_*
      }
      #print Dumper %fq;
      my $cll= $clusters{$loc}{pos1}{cliploc};
      my @k = keys %fq;
      @k = sort{$fq{$b}{num} <=> $fq{$a}{num} or $fq{$b}{te} <=> $fq{$a}{te} or abs($fq{$a}{loc}-$cll) <=> abs($fq{$b}{loc}-$cll)} keys %fq if(scalar @k>1);
      if(scalar @knownte >0){
        foreach my $iid(@knownte){
          if(!$clusters{$loc}{smate}{cliploc} or $clusters{$loc}{smate}{tmp}< $fq{$iid}{num}){
           $clusters{$loc}{smate}{smate} = $iid."_".$fq{$iid}{teO}."|".$fq{$iid}{te}; ## most represented location, total clip reads
           $clusters{$loc}{smate}{tmp} = $fq{$iid}{num}; ## tmp
           $clusters{$loc}{smate}{fuzzyclipped} = scalar(@l) - $fq{$iid}{num};
           $clusters{$loc}{smate}{cliploc} = $fq{$iid}{loc}; 
           $clusters{$loc}{smate}{AllClip} = $fq{$iid}{id};
           if(defined $clusters{$loc}{smate}{clipseqs}){
             $clusters{$loc}{smate}{clipseqs} =~ s/,$//;
             my @tmm = grep {$_ =~ m/\|$fq{$iid}{loc}/; s/\|$fq{$iid}{loc}//; } split(",",$clusters{$loc}{smate}{clipseqs});
             $clusters{$loc}{smate}{clipseqs} = src::Utilities::getCorrectClipSeq(join(",",@tmm),$fq{$iid}{side});
           }
          }
        }
      }else{
          $clusters{$loc}{smate}{smate} = $k[0]."_".$fq{$k[0]}{teel}."|".$fq{$k[0]}{num}; ## most represented location, total clip reads
          $clusters{$loc}{smate}{tmp} = $fq{$k[0]}{num}; ## tmp
          $clusters{$loc}{smate}{fuzzyclipped} = scalar(@l) - $fq{$k[0]}{num};
          $clusters{$loc}{smate}{cliploc} = $fq{$k[0]}{loc}; 
          $clusters{$loc}{smate}{AllClip} = $fq{$k[0]}{id};
          if(defined $clusters{$loc}{smate}{clipseqs}){
            $clusters{$loc}{smate}{clipseqs} =~ s/,$//;
            my @tmm = grep {$_ =~ m/\|$fq{$k[0]}{loc}/; s/\|$fq{$k[0]}{loc}//; } split(",",$clusters{$loc}{smate}{clipseqs});
            $clusters{$loc}{smate}{clipseqs} = src::Utilities::getCorrectClipSeq(join(",",@tmm),$fq{$k[0]}{side});
          }          
      }
      my $ty = $clusters{$loc}{smate}{AllClip};
      $ty =~ s/=/chr/g; $ty =~ s/\*/chr/g; $ty =~ s/chr/ /g;
      my @ty = split(" ",$ty);
      @ty = grep { $_ !~ ';' } @ty;
      @ty = grep { $_ ne ' ' } @ty;
      my %ty; $ty{$_}++ for(@ty);
      my @keys = sort { $ty{$b} <=> $ty{$a} } keys %ty;
      $clusters{$loc}{smate}{identClip}=src::Utilities::getClonalClipPercent($clusters{$loc}{smate}{AllClip},";");
    }
    
    ## reciprocal clusters
    @features = keys %transposons; 
    if(scalar keys %polymorphs > 0 ){
      push(@features,keys %polymorphs);
    }
    foreach my $feature (@features){
      if(defined($clusters{$loc}{pos1}{$feature} )){
        $clusters{$loc}{pos1}{$feature}{TESc} =~ s/,$//;
        $clusters{$loc}{pos1}{$feature}{TESc} = src::Utilities::max( split(",",$clusters{$loc}{pos1}{$feature}{TESc}) );
        my @array = @{src::Utilities::get_reciprocalClust($clusters{$loc}{pos1}{$feature}{TEClust}, $feature,$clusters{$loc}{pos1}{$feature}{num})}; ## position,freq
        $clusters{$loc}{pos1}{$feature}{TEClustPos} = $array[0];
        $clusters{$loc}{pos1}{$feature}{TEClustReads} = $array[1];
        $clusters{$loc}{pos1}{$feature}{TEClustPerc} = $array[2];
        $clusters{$loc}{pos1}{$feature}{TEClust} = $array[0].":".$array[1];
        $clusters{$loc}{pos1}{$feature}{ClustIds} = $array[3];
        $clusters{$loc}{pos1}{$feature}{mapq} = $clusters{$loc}{pos1}{mapq};
      }
      if(exists $clusters{$loc}{smate} and defined($clusters{$loc}{smate}{$feature} )){
        $clusters{$loc}{smate}{$feature}{TESc} =~ s/,$//;
        $clusters{$loc}{smate}{$feature}{TESc} = src::Utilities::max( split(",",$clusters{$loc}{smate}{$feature}{TESc}) );
        my @array = @{src::Utilities::get_reciprocalClust($clusters{$loc}{smate}{$feature}{TEClust}, $feature, $clusters{$loc}{smate}{$feature}{num} )}; ## position,freq
        $clusters{$loc}{smate}{$feature}{TEClustPos} = $array[0];
        $clusters{$loc}{smate}{$feature}{TEClustReads} = $array[1];
        $clusters{$loc}{smate}{$feature}{TEClustPerc} = $array[2];
        $clusters{$loc}{smate}{$feature}{TEClust} = $array[0].":".$array[1];
        $clusters{$loc}{smate}{$feature}{ClustIds} = $array[3];
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
  }
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
    $j =src::Utilities::unique($j,",","false");
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
       $clusters{$loc}{final}{start} = src::Utilities::min(($clusters{$loc}{pos1}{cliploc},$clusters{$loc}{pos2}{cliploc})) ; 
       $clusters{$loc}{final}{end}   = src::Utilities::max(($clusters{$loc}{pos1}{cliploc},$clusters{$loc}{pos2}{cliploc}));     
       $clusters{$loc}{final}{dist2re} = src::Utilities::min(($clusters{$loc}{pos1}{dist2re},$clusters{$loc}{pos2}{dist2re}));
       $clusters{$loc}{final}{maxscore} = src::Utilities::max(($clusters{$loc}{pos1}{maxscore},$clusters{$loc}{pos2}{maxscore}));       
       $clusters{$loc}{final}{identClip} = src::Utilities::max(($clusters{$loc}{pos1}{identClip},$clusters{$loc}{pos2}{identClip}));       
       $clusters{$loc}{final}{OriRatio} = src::Utilities::min(($clusters{$loc}{pos1}{OriRatio},$clusters{$loc}{pos2}{OriRatio}));       
     }else{
       $clusters{$loc}{final}{start} = $clusters{$loc}{pos1}{cliploc}; 
       $clusters{$loc}{final}{end}   = $clusters{$loc}{pos1}{cliploc};     
       $clusters{$loc}{final}{dist2re} = $clusters{$loc}{pos1}{dist2re};
       $clusters{$loc}{final}{maxscore} = $clusters{$loc}{pos1}{maxscore};
       $clusters{$loc}{final}{identClip} = $clusters{$loc}{pos1}{identClip};
       $clusters{$loc}{final}{OriRatio} = $clusters{$loc}{pos1}{OriRatio};
     }
     
     ## (2) PolyA tails?
     $clusters{$loc}{pos1}{isPolyA}="F";
     $clusters{$loc}{final}{isPolyA}="F";   
     my $aa = src::Utilities::getstrand($clusters{$loc}{pos1}{seqstring});
     $clusters{$loc}{pos1}{isPolyA}="T" if($aa ne "NA" or exists $clusters{$loc}{pos1}{PolyA});
     $clusters{$loc}{final}{isPolyA}="T" if($aa ne "NA" or exists $clusters{$loc}{pos1}{PolyA});
     $clusters{$loc}{pos1}{strand} = "+" if($aa eq "PolyA");
     $clusters{$loc}{pos1}{strand} = "-" if($aa eq "PolyT");
     $clusters{$loc}{final}{tailinfo} = $aa; 
     if(exists $clusters{$loc}{pos2}){
       $clusters{$loc}{pos2}{isPolyA}="F";
       my $b = src::Utilities::getstrand($clusters{$loc}{pos2}{seqstring});
       if($b ne "NA"){
         $clusters{$loc}{pos2}{strand} = "+" if($b eq "PolyA");
         $clusters{$loc}{pos2}{strand} = "-" if($b eq "PolyT"); 
       }elsif(defined $clusters{$loc}{pos2}{clipseqs}){
         my $sq = $clusters{$loc}{pos2}{clipseqs};
         ## Seqstring instances do not necessarily support smate based strand identification
         if($clusters{$loc}{pos2}{clip_sidedness} eq "a"){
            my $cig = (length($sq)-10)."M10S";
            $b = src::Utilities::getstrand($sq.",".$cig.",-,1"); 
         }elsif($clusters{$loc}{pos2}{clip_sidedness} eq "b"){
            my $cig = "10S".(length($sq)-10)."M";
            $b = src::Utilities::getstrand($sq.",".$cig.",-,11"); 
         }
         $clusters{$loc}{pos2}{strand} = "+" if($b eq "PolyA");
         $clusters{$loc}{pos2}{strand} = "-" if($b eq "PolyT"); 
       }
       $clusters{$loc}{pos2}{isPolyA}="T" if($b ne "NA" or exists $clusters{$loc}{pos2}{PolyA});
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
     $aa=0;    
     foreach my $feature (@features){
       $aa++ if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature});
     }
     if($aa==0 and exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{smate} and $clusters{$loc}{pos2}{isPolyA} eq "T"){
       my @xx = ($clusters{$loc}{pos2}{smate} =~ m/\S*_\S*_\S*_(\S*)\|(\d*)/);
       if($xx[0] eq "*"){
         $clusters{$loc}{pos2}{PolyA}{TEClustPos} = "1";
         $clusters{$loc}{pos2}{PolyA}{TEClustReads} = $xx[1];
         $clusters{$loc}{pos2}{PolyA}{TEClustPerc} = "1.00";
         $clusters{$loc}{pos2}{PolyA}{TEClust} = "1:".$xx[1];
         $clusters{$loc}{pos2}{PolyA}{ClustIds} = $xx[1];
         $clusters{$loc}{pos2}{PolyA}{TESc} = "30";
         $clusters{$loc}{pos2}{PolyA}{num} = $xx[1];
         $clusters{$loc}{pos2}{PolyA}{mapq} = $clusters{$loc}{pos2}{mapq};
       }    
     }
     undef(@features);
     undef $aa;
     
     ## (2.2) if clip read maps to PolyA or has an inferred PolyA mapping, assign the Mapscore to minimum 30 
      ## PolyA mapping is inferential and hence should not decide cluster ommission 
     $clusters{$loc}{pos1}{PolyA}{TESc} = 30 if(exists $clusters{$loc}{pos1}{PolyA} and $clusters{$loc}{pos1}{PolyA}{TESc}<30);
     if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{PolyA} and $clusters{$loc}{pos2}{PolyA}{TESc}<30){
      $clusters{$loc}{pos2}{PolyA}{TESc} = 30;
     }

     ## (3) TSD
     my @l = split(",",($clusters{$loc}{pos1}{seqstring}.",".$clusters{$loc}{final}{start}.",".$clusters{$loc}{final}{end}));
     my $tsd = src::Utilities::tsd(\@l);
     if($tsd eq "" or $tsd eq "ERROR-TSD" and exists $clusters{$loc}{pos2}){
         @l = split(",",($clusters{$loc}{pos2}{seqstring}.",".$clusters{$loc}{final}{start}.",".$clusters{$loc}{final}{end}));
         $tsd = src::Utilities::tsd(\@l);
     }
     if($tsd eq "" or $tsd eq "ERROR-TSD"){
         $clusters{$loc}{final}{TSD} = "*"; 
      }else{
        $clusters{$loc}{final}{TSD} = $tsd;
     }

     ## (4) distance to annotated RE motif
     $clusters{$loc}{final}{dist2remotif} = $clusters{$loc}{final}{dist2re};
     #$clusters{$loc}{final}{dist2remotif} = ">100" if($clusters{$loc}{final}{dist2remotif}==500);

     ## (5) All clip read counts
     if(exists $clusters{$loc}{pos2}){
      $clusters{$loc}{final}{commonClips} = src::Utilities::common($clusters{$loc}{pos1}{AllClip},$clusters{$loc}{pos2}{AllClip},";","true");
     }else{
       $clusters{$loc}{final}{commonClips} =0;
     }     

     ## (6) All clip read counts
     undef $aa;
     $aa = $clusters{$loc}{pos1}{AllClip};
     $aa .= $clusters{$loc}{pos2}{AllClip} if(exists $clusters{$loc}{pos2});
     $aa =~ s/;$//;
     $clusters{$loc}{final}{AllClip} = src::Utilities::unique($aa, ";","true");
     $clusters{$loc}{pos1}{AllClip} =~ s/;$//; 
     $clusters{$loc}{pos1}{AllClip} = src::Utilities::unique($clusters{$loc}{pos1}{AllClip}, ";","true");
     if(exists $clusters{$loc}{pos2}){
       $clusters{$loc}{pos2}{AllClip} =~ s/;$//; 
       $clusters{$loc}{pos2}{AllClip} = src::Utilities::unique($clusters{$loc}{pos2}{AllClip}, ";","true");
     }

     ## (7) Finalize TE 
     my %inv;
     if($polym ne ""){
        foreach my $ele (keys %crossref){ $inv{$crossref{$ele}} .= $ele.",".$crossref{$ele}.","; } 
        foreach my $ele (keys %inv){ $inv{$ele} = src::Utilities::unique($inv{$ele},",","false");}
        delete $inv{'PolyA'} if(exists $inv{'PolyA'});
     }else{
        foreach my $ele (keys %transposons){ $inv{$ele} .= $ele.","; } 
        delete $inv{'PolyA'} if(exists $inv{'PolyA'});
     }
     if(scalar keys %inv > 0){
       foreach my $ele (keys %inv){
        $inv{$ele} =~ s/,$//;
        my @features = split(",",$inv{$ele});
        @features = grep {$_ ne 'PolyA'} @features;
        my %chq; my %cc;
        $chq{lpos}="-"; $chq{lte}=""; $chq{lpoly}=""; $chq{LTEClustReads}=0;
        $chq{rte}=""; $chq{rpoly}=""; $chq{rpos}="-"; $chq{RTEClustReads}=0;
        $chq{TESc}=0; $chq{TEClustPerc}=0;
        $chq{tmapq}=0;$chq{pmapq}=0;
        $chq{rTESc} =0; $chq{rtmapq}=0;
        $chq{TEsize}=0; $chq{OriginalTE}="";
        foreach my $feature(@features){
          if(exists $clusters{$loc}{pos1}{$feature}){
            $cc{$feature}{clust} += $clusters{$loc}{pos1}{$feature}{TEClustReads} ;
            $cc{$feature}{TESc} .= $clusters{$loc}{pos1}{$feature}{TESc}.",";
            $cc{$feature}{mapq} .= $clusters{$loc}{pos1}{$feature}{mapq}.",";
            $chq{lte} .= $clusters{$loc}{pos1}{$feature}{num}.",";
            if($chq{LTEClustReads} ==0 or $clusters{$loc}{pos1}{$feature}{TEClustReads} > $chq{LTEClustReads}){
              $chq{lpos}=$clusters{$loc}{pos1}{$feature}{TEClustPos};
              $chq{LTEClustReads}=$clusters{$loc}{pos1}{$feature}{TEClustReads};  
            }elsif($chq{LTEClustReads} >0 and $clusters{$loc}{pos1}{$feature}{TEClustReads} == $chq{LTEClustReads}){
              if($chq{TESc} < $clusters{$loc}{pos1}{$feature}{TESc}){
                $chq{lpos}=$clusters{$loc}{pos1}{$feature}{TEClustPos};
                $chq{LTEClustReads}=$clusters{$loc}{pos1}{$feature}{TEClustReads};  
              }elsif($chq{tmapq} < $clusters{$loc}{pos1}{$feature}{mapq}){
                $chq{lpos}=$clusters{$loc}{pos1}{$feature}{TEClustPos};
                $chq{LTEClustReads}=$clusters{$loc}{pos1}{$feature}{TEClustReads};  
              }     
            }
            $chq{TESc}=$clusters{$loc}{pos1}{$feature}{TESc} if($chq{TESc}==0 or $clusters{$loc}{pos1}{$feature}{TESc} > $chq{TESc});
            $chq{tmapq}=$clusters{$loc}{pos1}{$feature}{mapq} if($chq{tmapq}==0 or $clusters{$loc}{pos1}{$feature}{mapq} > $chq{tmapq});
          }
          if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature}){
            $cc{$feature}{clust} += $clusters{$loc}{pos2}{$feature}{TEClustReads} ;
            $cc{$feature}{TESc} .= $clusters{$loc}{pos2}{$feature}{TESc}.",";
            $cc{$feature}{mapq} .= $clusters{$loc}{pos2}{$feature}{mapq}.",";
            $chq{rte} .= $clusters{$loc}{pos2}{$feature}{num}.",";
            if($chq{RTEClustReads} ==0 or $clusters{$loc}{pos2}{$feature}{TEClustReads} > $chq{RTEClustReads}){
              $chq{rpos}=$clusters{$loc}{pos2}{$feature}{TEClustPos};
              $chq{RTEClustReads}=$clusters{$loc}{pos2}{$feature}{TEClustReads}; 
            }elsif($chq{RTEClustReads} >0 and $clusters{$loc}{pos2}{$feature}{TEClustReads} == $chq{RTEClustReads}){
              if($chq{rTESc} < $clusters{$loc}{pos2}{$feature}{TESc}){
                $chq{rpos}=$clusters{$loc}{pos2}{$feature}{TEClustPos};
                $chq{RLTEClustReads}=$clusters{$loc}{pos2}{$feature}{TEClustReads};  
              }elsif($chq{rtmapq} < $clusters{$loc}{pos2}{$feature}{mapq}){
                $chq{rpos}=$clusters{$loc}{pos2}{$feature}{TEClustPos};
                $chq{RTEClustReads}=$clusters{$loc}{pos2}{$feature}{TEClustReads};  
              }     
            }
            $chq{rTESc}=$clusters{$loc}{pos2}{$feature}{TESc} if($chq{rTESc} ==0 or $clusters{$loc}{pos2}{$feature}{TESc}> $chq{rTESc});
            $chq{rtmapq}=$clusters{$loc}{pos2}{$feature}{mapq}  if($chq{rtmapq}==0 or $clusters{$loc}{pos2}{$feature}{mapq}> $chq{rtmapq});         
          }
        }
        if(scalar (keys %cc) >0){
          foreach my $i (keys %cc){
            $cc{$i}{TESc} = src::Utilities::array_sort($cc{$i}{TESc});
            $cc{$i}{mapq} = src::Utilities::array_sort($cc{$i}{mapq});
          }
          my @keyscc = keys %cc;
          @keyscc = sort { no warnings;
                           my $aC = $cc{$a}{clust};
                           my $bC = $cc{$b}{clust};
                           my $aT = $cc{$a}{TESc};
                           my $bT = $cc{$b}{TESc};
                           my $aM = $cc{$a}{mapq};
                           my $bM = $cc{$b}{mapq}; 
                           $bC <=> $aC or $bT <=> $aT or $bM <=> $aM or $b cmp $a; } @keyscc if(scalar @keyscc >1);
          $chq{OriginalTE} = $keyscc[0];
          $chq{TEsize} = $transposons{$keyscc[0]} if($transposons{$keyscc[0]});
          $chq{TEsize} = $polymorphs{$keyscc[0]} if($polymorphs{$keyscc[0]});           
        }
        
        if(exists$clusters{$loc}{pos1}{PolyA}){
          $chq{lpoly} .= $clusters{$loc}{pos1}{PolyA}{num}."," ;
          $chq{pmapq}=$clusters{$loc}{pos1}{PolyA}{mapq};
        }
        if(exists$clusters{$loc}{pos2} and exists$clusters{$loc}{pos2}{PolyA}){
          $chq{rpoly} .= $clusters{$loc}{pos2}{PolyA}{num}.",";
          $chq{pmapq}=$clusters{$loc}{pos2}{PolyA}{mapq} if($clusters{$loc}{pos2}{PolyA}{mapq} > $chq{pmapq});
        }
        $chq{tnum} = src::Utilities::unique($chq{lte}.$chq{rte}, ",","true");
        $chq{anum} = src::Utilities::unique($chq{lte}.$chq{rte}.$chq{lpoly}.$chq{rpoly}, ",","true");
        $chq{pnum} = $chq{anum} -$chq{tnum};
        $chq{umap} = $clusters{$loc}{final}{AllClip} - $chq{anum};
        $chq{TEClustReads} = $chq{RTEClustReads}+$chq{LTEClustReads};
        $chq{TESc} = src::Utilities::max(($chq{TESc},$chq{rTESc}));
        $chq{tmapq} = src::Utilities::max(($chq{tmapq},$chq{rtmapq}));
        $chq{lcloneclip} = src::Utilities::getClonalClipPercent($chq{lte}.$chq{lpoly},",");
        $chq{rcloneclip} = src::Utilities::getClonalClipPercent($chq{rte}.$chq{rpoly},",");        
        delete($chq{rTESc});
        delete($chq{rtmapq});
        if(defined $chq{TEClustReads} and defined $chq{pnum}){
         $chq{TEClustReads}+=$chq{pnum};
        }
        $chq{TEClustPerc}=src::Utilities::round($chq{TEClustReads}/$chq{anum},2) if($chq{anum}>0);
        if($chq{rpos} eq "-" and $chq{lpos} eq "-"){
          $chq{TEClustPerc} = 0;
          $chq{TEClustReads} = 0;
        }
        if($chq{TEClustPerc}>0 and $chq{tnum}>0 and ($chq{tmapq}>$refMapqQ or $chq{pmapq}>$refMapqQ )){
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
        $clusters{$loc}{pos1}{$feature}{num} = src::Utilities::unique($clusters{$loc}{pos1}{$feature}{num}, ",","true");           
      }
      if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$feature}){
        $clusters{$loc}{pos2}{$feature}{num} = src::Utilities::unique($clusters{$loc}{pos2}{$feature}{num}, ",","true");           
      }
      if(exists $clusters{$loc}{final}{$feature}){
        $clusters{$loc}{final}{$feature}{lte} = src::Utilities::unique($clusters{$loc}{final}{$feature}{lte}, ",","true");          
        $clusters{$loc}{final}{$feature}{rte} = src::Utilities::unique($clusters{$loc}{final}{$feature}{rte}, ",","true");          
        $clusters{$loc}{final}{$feature}{lpoly} = src::Utilities::unique($clusters{$loc}{final}{$feature}{lpoly}, ",","true");          
        $clusters{$loc}{final}{$feature}{rpoly} = src::Utilities::unique($clusters{$loc}{final}{$feature}{rpoly}, ",","true");          
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
       $clusters{$loc}{final}{$feature}{TEClustPerc} = src::Utilities::round($clusters{$loc}{final}{$feature}{TEClustReads}/$clusters{$loc}{final}{$feature}{anum},2);       
      }
     }
     undef @features;
     
     ## (9) Purity ratio
     @features = grep { $_ ne 'PolyA' } (keys %transposons);
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature}){
       $clusters{$loc}{final}{$feature}{Purity} = src::Utilities::round($clusters{$loc}{final}{$feature}{anum}/$clusters{$loc}{final}{AllClip},2);
      }
     }
     
     ## (10) insert size
     @features = grep { $_ ne 'PolyA' } (keys %transposons);
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature}){
        $clusters{$loc}{final}{$feature}{insert_size} = "NA";
        my $min = src::Utilities::min(($clusters{$loc}{final}{$feature}{lpos},$clusters{$loc}{final}{$feature}{rpos}));
        my $max = src::Utilities::max(($clusters{$loc}{final}{$feature}{lpos},$clusters{$loc}{final}{$feature}{rpos}));
        if(defined $max and defined $min and $max != $min){
          $clusters{$loc}{final}{$feature}{insert_size} = ($max - $min);   
        }elsif(defined $max and defined $min and $max == $min){
           if($clusters{$loc}{final}{$feature}{rpoly}>0 or $clusters{$loc}{final}{$feature}{lpoly}>0){
             $clusters{$loc}{final}{$feature}{insert_size} = ($clusters{$loc}{final}{$feature}{TEsize} - $min);  
           }else{
              my $yy = $clusters{$loc}{final}{$feature}{OriginalTE};
              if(exists $clusters{$loc}{pos1}{$yy} and exists $clusters{$loc}{pos2} and $clusters{$loc}{pos2}{isPolyA} eq "T" and $clusters{$loc}{pos1}{isPolyA} eq "F"){
                $clusters{$loc}{final}{$feature}{insert_size} = ($clusters{$loc}{final}{$feature}{TEsize} - $min);
              }elsif(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$yy} and $clusters{$loc}{pos1}{isPolyA} eq "T" and $clusters{$loc}{pos2}{isPolyA} eq "F"){
                $clusters{$loc}{final}{$feature}{insert_size} = ($clusters{$loc}{final}{$feature}{TEsize} - $min);
              }
           }
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
     $clusters{$loc}{final}{LHS_cons} = src::Utilities::get_longestseq($clusters{$loc}{pos1}{clipseqs});
     print FO1 ">$l","_RC\n$clusters{$loc}{final}{LHS_cons}\n";
     if(exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{clipseqs}){
      $clusters{$loc}{final}{RHS_cons} = src::Utilities::get_longestseq($clusters{$loc}{pos2}{clipseqs});  
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
       my $rat = src::Utilities::round($score/length($temp[9]),2);
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
  $r = src::Utilities::round($cnt*100/scalar(keys %clusters),2) if(scalar(keys %clusters) > 0);   
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
      $clusters{$loc}{final}{$te}{bgTEDist} = $temp[14];
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
       $w[1] = int(($w[1]/1000)+0.5);
       my $valnum=0;
       my $val=0;
       for(my $e=-1;$e<=1;$e++){
         if($cov{$w[0]}{$w[1]+$e}){
           $val += $cov{$w[0]}{$w[1]+$e};
           $valnum++;
         }
       }
       $val = src::Utilities::round($val/(10*$valnum),2) if($valnum>0);   ## readlength 100 vs bin size 1000, hence division by 10
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
