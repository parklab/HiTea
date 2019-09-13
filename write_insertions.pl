#!/usr/bin/perl
# HiTEA

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
my $index="";
my $outprefix="project";
my $gap = 2;
my $readlen=100;
my $min_reads=3;
my $refMapqQ="";
my $help=0;
my $wd = "";
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'index=s'          => \$index,
  'gap=s'            => \$gap,
  'min_reads=s'      => \$min_reads,
  'readlen=s'        => \$readlen,
  'q=s'        => \$refMapqQ,
  'outprefix=s'      => \$outprefix,
  'wd:s'             => \$wd,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl write_insertions.pl -in [FILE_PATH] -outprefix [STRING] -gap [INT] -readlen [INT] -q [referenceMAPQ] -min_reads [INT] -index [TE_MergedIndex] -wd [DIR_PATH]\n\n";
  print "This program annotates finalized breaks using the bam file\n\n";
  print "Options:\n\n";
  print "***required:\n";
  print "  -in                    Clusters file in bed format \n";
  print "  -index                 transposable element reference assembly (fasta path) \n";
  print "  -q                     reference mapq score threshold \n";
  print "***optional:\n";
  print "  -min_reads             minimum reads that need to map at an insertion site [default: 3] \n";
  print "  -gap                   gap to merge the breakpoints [default: 2] \n";
  print "  -readlen               sequencing read length [default: 100] \n";
  print "  -outprefix             Outputprefix for generating 2 output files [default: project] \n";
  print "  -wd                    Working directory [default: ~]\n";
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
my $polyclust=0;
my $pval_cutoff = 0.05;
my $unassigned = 0;
my $filtered = 0;
my $RAMCountCUTOFF = 3;
my $PURITY_CUTOFF = 0.25;
my %clusters;
my %transposons;
my %cov;
my %cntcalls;
print " ERROR: makeMat.R does not exist in the same directory!\n" unless -e "makeMat.R";

#--------------------------------------------------------------------------------------------------------------
# I/O
#--------------------------------------------------------------------------------------------------------------
my $olog = $outprefix.".skippedclusters.logs.gz";
open(LOGS,"| gzip -c - > $olog") or die $!;

my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "[write_insertions] START:\t $run_time seconds\n";
print " Command: perl write_insertions.pl -in $in -gap $gap -min_reads $min_reads -readlen $readlen -outprefix $baseoutprefix -wd $wd \n";

%cov = %{retrieve($outprefix.".coverage.ph")};

## get transposons
%transposons = Utilities::get_fasta($index);
print " transposons being considered in the analyses: " , join(",",keys %transposons),"; sizes in bp (",join (",",values %transposons),")\n";

print " reading putative insertions and finalizing the calls\n";
%clusters = %{retrieve($in)};
print "  -> total putative insertions in completed object: ", scalar(keys%clusters),"\n";

## fitler the insertion breakpoints
my $clust = filter_clusters(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " filtered breakpoint object:\t $run_time seconds\n";
print "  -> total insertions after filtering: ", scalar(keys%clusters),"\n";


## write report
$clust = write_clusters(\%clusters);
%clusters = %{$clust};
$watch_run = time();
$run_time = $watch_run -  $start_run;
print " wrote insertion breakpoint report:\t $run_time seconds\n";

my $oufile=$outprefix.'.finalFilt'; #Save 
open FOO,">$oufile" or die $!;
print FOO Dumper %clusters;
close(FOO);
exit 1;

## generate matrix
my $testring = join(",",keys %transposons);
$testring =~ s/PolyA//;
$testring=~ s/,,/,/;
$testring=~ s/^,//;
$testring=~ s/,$//;
system( qq[Rscript makeMat.R $baseoutprefix $wd $testring] ) == 0 or die qq[ Cound not generate count matrices for coverage plots\n];
$watch_run = time();
$run_time = $watch_run - $start_run;
print " generated coverage matrices of RAM and split RAM around the insertion breakpoints:\t $run_time seconds\n";


$watch_run = time();
$run_time = $watch_run -  $start_run;
print "[write_insertions] END:\t $run_time seconds\n\n";
close(LOGS);

#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub filter_clusters{
  my ($clust) = shift;
  my %clusters = %{$clust};
  my @features = grep { $_ ne 'PolyA' } (keys %transposons);
  my $TESCORE=20;   
  
  foreach my $loc(sort keys %clusters){
   if(!defined($clusters{$loc}{final}{chr}) or !defined $clusters{$loc}{final}{start} ) {
     print  " ERROR: Something is not defined in the insertion calls\n";
     exit 1;
   }

   ## Dump clusters due to low coverage, STR or proximity to the RE motif
   my $check=0;
   foreach my $feature (@features){
    $check++ if(exists $clusters{$loc}{final}{$feature});
   }
   if($check ==0){
     my $a=0;
     $a += $clusters{$loc}{pos1}{PolyA}{num} if(exists$clusters{$loc}{pos1}{PolyA});
     my $b=0;
     $b += $clusters{$loc}{pos2}{PolyA}{num} if(exists$clusters{$loc}{pos2} and exists$clusters{$loc}{pos2}{PolyA});
     print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
     print LOGS "no assignment (write step)\tPolyAreads:$a;$b\tAllClip:$clusters{$loc}{final}{AllClip}\n";
     $filtered++;
     $unassigned++;       
     delete $clusters{$loc};
     next;
   }
   if($clusters{$loc}{final}{dist2re}<2){
    print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
    print LOGS "<proximity to RE (write step)\tAllreads:$clusters{$loc}{final}{AllClip}\tD2RE:$clusters{$loc}{final}{dist2re}\n";
    $filtered++;       
    delete $clusters{$loc};
    next;
   }
   foreach my $feature (@features){
    if(exists $clusters{$loc}{final}{$feature} and $clusters{$loc}{final}{$feature}{anum}<2 ){
       print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
       print LOGS "<2TE-mapping reads (write step)\tMapReads:$clusters{$loc}{final}{$feature}{anum};$clusters{$loc}{final}{AllClip}\t$feature\n";
      $clusters{$loc}{final}{$feature}{status} = 0;
      $clusters{$loc}{final}{$feature}{remark} = "<2anum;";
       $filtered++;       
       delete $clusters{$loc}{final}{$feature};
       next;
    }
   }
   ## STR?
   $clusters{$loc}{final}{STR} = "F";
   if($clusters{$loc}{final}{tailinfo} eq "PolyAPolyA" or $clusters{$loc}{final}{tailinfo} eq "PolyAPolyT" or $clusters{$loc}{final}{tailinfo} eq "PolyTPolyA" or $clusters{$loc}{final}{tailinfo} eq "PolyTPolyT"){
     $clusters{$loc}{final}{STR} = "T";
   }
 
   # Filtering the clusters ==Assign confidence to the clsuter
   ## Score 3: high conf calls,  lhs, rhs clip reads present and >10% clip reads map to TE
   ## Score 2: optimal,          either lhs or rhs read support and >10% clip reads map to TE
   ## Score 1: Overlaps with bg, The copy overlaps with background TE copy 
   ## Score 0: poor evidence,    Based on several filtering criteria 
   my $cstart = int ($clusters{$loc}{final}{start}/1000);
   my $cchr = $clusters{$loc}{final}{chr};
   if(defined $cstart){
     my $bgcov=0;
     for(my $i=-1;$i<=1;$i++){
       $cov{$cchr}{($cstart+$i)} = 0 if(!$cov{$cchr}{($cstart+$i)});
       $bgcov +=$cov{$cchr}{($cstart+$i)};
     }
     $bgcov = $bgcov/3; 
     $bgcov = Utilities::round( ($bgcov*$readlen)/1000,2); 
     $clusters{$loc}{final}{bgcoverage} = $bgcov;                         
   }  
   $clusters{$loc}{final}{sort}=0;   
   
   foreach my $feature (@features){
    if(exists $clusters{$loc}{final}{$feature}){     
      my $purity_ratio = $clusters{$loc}{final}{$feature}{Purity};
      my $tsp =  $feature;
      my $pval = $clusters{$loc}{final}{$feature}{pvalue};
      $clusters{$loc}{final}{$feature}{status} = 0;
      $clusters{$loc}{final}{$feature}{remark} = "-";
      if($purity_ratio >=$PURITY_CUTOFF and $pval<$pval_cutoff){
        if(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{lte}>=1 and $clusters{$loc}{final}{$tsp}{rte}>=1){
          $clusters{$loc}{final}{$tsp}{status} = 3; 
        }elsif(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{lte}>=1 and $clusters{$loc}{final}{$tsp}{rpoly}>=1){
          $clusters{$loc}{final}{$tsp}{status} = 3; 
        }elsif(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{rte}>=1 and $clusters{$loc}{final}{$tsp}{lpoly}>=1){
          $clusters{$loc}{final}{$tsp}{status} = 3; 
        }elsif($clusters{$loc}{final}{$tsp}{tnum}>=2 and $clusters{$loc}{final}{$tsp}{TESc}>=$TESCORE){
          $clusters{$loc}{final}{$tsp}{status} = 2; 
        }
      }
      
      ## STR 
      if($clusters{$loc}{final}{STR} eq "T"){
        $clusters{$loc}{final}{$feature}{status} = 0;
        $clusters{$loc}{final}{$feature}{remark} = "STR;";
      }
      ## TE map score -if it is <30 (arbitrary cutoff, check if PolyA has >30 mapping score, else assign the cluster as low conf)
      my $ATcnt = 0.0;
      if(exists $clusters{$loc}{final}{TSD} and $clusters{$loc}{final}{TSD} ne "-" and $clusters{$loc}{final}{TSD} ne "*"){
       my $sq1= $clusters{$loc}{final}{TSD};
       my $sq2= $clusters{$loc}{final}{TSD};
       $sq1 =~ s/[CGTN]/Q/g;
       $sq1 =~ s/Q$//;
       $sq2 =~ s/[CGAN]/Q/g;
       $sq2 =~ s/Q$//;
       $sq1 = $sq1."Q".$sq2;
       $sq1 =~ s/[Q]+/,/g;
       my %hash = map { $_ => length($_) } split(",",$sq1);
       my @keys = sort { $hash{$b} <=> $hash{$a} } keys %hash;
       $ATcnt = Utilities::round(length($keys[0])/length($clusters{$loc}{final}{TSD}),5) if(defined $keys[0]);
      }
      $clusters{$loc}{final}{ATpercent} = $ATcnt;
      if($clusters{$loc}{final}{$tsp}{TESc} <30 or $clusters{$loc}{final}{$tsp}{tmapq} < $refMapqQ){
        my $cq = 0;
        if($clusters{$loc}{final}{start}==$clusters{$loc}{final}{end}){
         $cq =1;
        }elsif(exists $clusters{$loc}{pos1}{PolyA} and $clusters{$loc}{pos1}{PolyA}{TESc}>=30 and $ATcnt < 0.7){
        }elsif(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{PolyA} and $clusters{$loc}{pos2}{PolyA}{TESc}>=30 and $ATcnt < 0.7){
        }elsif(exists $clusters{$loc}{pos1}{isPolyA} and $clusters{$loc}{pos1}{isPolyA} eq "T" and $ATcnt < 0.7){
        }elsif(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{isPolyA} and $clusters{$loc}{pos2}{isPolyA} eq "T" and $ATcnt < 0.7){
        }else{
          $cq =1;
        }
        if($cq >0){
         $clusters{$loc}{final}{$tsp}{status} = 0; 
         $clusters{$loc}{final}{$tsp}{remark} .= "PoorTEMapScore;";               
        }
      }
      ## double clipped reads
        ## check if the read defining insertion are double clipped and both ends are used for TSD 
      if($clusters{$loc}{final}{commonClips} >0){
        $clusters{$loc}{final}{$feature}{status} = 0;
        $clusters{$loc}{final}{$feature}{remark} = "doubleClip;";
      }elsif(exists $clusters{$loc}{pos2}){
        if($clusters{$loc}{pos1}{both} =~ m/$clusters{$loc}{pos2}{cliploc}/ or $clusters{$loc}{pos2}{both} =~ m/$clusters{$loc}{pos1}{cliploc}/ ){
         $clusters{$loc}{final}{$tsp}{status} = 0; 
         $clusters{$loc}{final}{$tsp}{remark} .= "doubleClip;";                         
        }
      }
      ## Pval
      if($pval >=$pval_cutoff){
       $clusters{$loc}{final}{$tsp}{status} = 0; 
       $clusters{$loc}{final}{$tsp}{remark} .= "NoEnrichment;";         
      }
      ## Fuzzy
      my $a = 0;
      $clusters{$loc}{final}{$tsp}{FbyC} = "0.0";
      if($clusters{$loc}{final}{$tsp}{AllFuzzy}>$clusters{$loc}{final}{AllClip}){
        my $a = $clusters{$loc}{final}{$tsp}{AllFuzzy};
        my $b = $clusters{$loc}{final}{AllClip};
        $b = $clusters{$loc}{final}{$tsp}{anum} if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{smate}); 
          ## as smate feature does nothave possibility of finding fuzzy reads
        $a = $a-$b;
        $clusters{$loc}{final}{$tsp}{AllFuzzy} = $a;
        $clusters{$loc}{final}{$tsp}{FbyC} = Utilities::round($a/($b+$a),2);
      }
      if($clusters{$loc}{final}{$tsp}{FbyC} >= 0.5){
        $clusters{$loc}{final}{$tsp}{status} = 0 ;  #if there are >50% fuzzy clipped reads in the vicinity, flag a cluster as Category 0 
        $clusters{$loc}{final}{$tsp}{remark} .= "fuzzy;";         
      }
      #if($clusters{$loc}{final}{$tsp}{FbyC} >= 0.2 and $a>=3){
      #  $clusters{$loc}{final}{$tsp}{status} = 0 ;  #if there are >20% fuzzy clipped reads in the vicinity, flag a cluster as Category 0 
      #  $clusters{$loc}{final}{$tsp}{remark} .= "fuzzy;";         
      #}elsif($clusters{$loc}{pos1}{bgclusters} > 2){ ## if the region has more than 2 clusters
      # $clusters{$loc}{final}{$tsp}{status} = 0;
      # $clusters{$loc}{final}{$tsp}{strand} ="*";
      # $clusters{$loc}{final}{$tsp}{remark} .= "fuzzy;";         
      #}
      ## Purity cutoff
      if($purity_ratio <$PURITY_CUTOFF){
        $clusters{$loc}{final}{$tsp}{status} = 0; 
        $clusters{$loc}{final}{$tsp}{remark} .= "InsuffClipCov;";               
      }
      ## distance2RE motif
      if($clusters{$loc}{final}{dist2remotif}<$gap){
         $clusters{$loc}{final}{$tsp}{status} = 0;
         $clusters{$loc}{final}{$tsp}{remark} .= "RE;";                             
      }
      ## <50% reciprocal coverage 
      if($clusters{$loc}{final}{$tsp}{TEClustPerc} < 0.5){
         $clusters{$loc}{final}{$tsp}{status} = 0; 
         $clusters{$loc}{final}{$tsp}{remark} .= "poorRecClust;";
      }
      ## RAM <1
      if( ($clusters{$loc}{final}{$tsp}{RAMCount}) ==0){
        $clusters{$loc}{final}{$tsp}{status} = 0;
        $clusters{$loc}{final}{$tsp}{remark} .= "RAM<1;";         
      }
      ## 10% of the bgCov reads mapping to TE
      if($clusters{$loc}{final}{$tsp}{anum} < Utilities::round($clusters{$loc}{final}{bgcoverage}/10,2) ){
        $clusters{$loc}{final}{$tsp}{status} = 0;
        $clusters{$loc}{final}{$tsp}{remark} .= "<10xClipCov;";   
      } 
      if($clusters{$loc}{final}{$tsp}{tnum} < Utilities::round($clusters{$loc}{final}{bgcoverage}/20,2) ){
        $clusters{$loc}{final}{$tsp}{status} = 0;
        $clusters{$loc}{final}{$tsp}{remark} .= "<10xClipTECov;";   
      } 
      if($clusters{$loc}{final}{$tsp}{RAMCount} < Utilities::round($clusters{$loc}{final}{bgcoverage}/20,2) ){ ## 20th of bgcov for RAM
        $clusters{$loc}{final}{$tsp}{status} = 0;
        $clusters{$loc}{final}{$tsp}{remark} .= "<5xRAMCov;";         
      }
      ## Overlap with reference copy
      if($clusters{$loc}{final}{$tsp}{bgTEDist}<=2){
        $clusters{$loc}{final}{$tsp}{remark} .= "bgOlap;";
       if($clusters{$loc}{final}{$tsp}{status}==3){  ## make status 1 of calls if original status is 3 and it overlaps with bgTE
         $clusters{$loc}{final}{$tsp}{status} = 1;  
       }elsif($clusters{$loc}{final}{$tsp}{status}==2){  ## make status 1 of calls if original status is 3 and it overlaps with bgTE
         $clusters{$loc}{final}{$tsp}{status} = 0;  
       }else{
         $clusters{$loc}{final}{$tsp}{status} = 0;  
       }
      }

      ## LHS/RHS map within 20bp on TE consensus
      if($clusters{$loc}{final}{start}!=$clusters{$loc}{final}{end}){
       if($clusters{$loc}{final}{$tsp}{rpos} ne "-" and $clusters{$loc}{final}{$tsp}{lpos} ne "-" and $clusters{$loc}{final}{isPolyA} eq "F"){
         if(abs($clusters{$loc}{final}{$tsp}{rpos}-$clusters{$loc}{final}{$tsp}{lpos})<20){
          $clusters{$loc}{final}{$tsp}{status} = 0;  
          $clusters{$loc}{final}{$tsp}{remark} .= "badPair;";
         }      
       }
      }

     ## if no size but PolyA and TSD present
     if($clusters{$loc}{final}{start}!=$clusters{$loc}{final}{end}){
      if($clusters{$loc}{final}{tailinfo} =~ /PolyAPolyA/ or $clusters{$loc}{final}{tailinfo} =~ /PolyAPolyT/ or $clusters{$loc}{final}{tailinfo} =~ /PolyTPolyA/ or $clusters{$loc}{final}{tailinfo} =~ /PolyTPolyT/){
        if($clusters{$loc}{final}{$tsp}{insert_size}=~ /NA/ or $clusters{$loc}{final}{$tsp}{insert_size} < 20){
         $clusters{$loc}{final}{$tsp}{status} = 0;  
         $clusters{$loc}{final}{$tsp}{remark} .= "badPair;";     
        }
      }elsif($clusters{$loc}{pos1}{isPolyA} eq "T" and exists $clusters{$loc}{pos2} and $clusters{$loc}{pos2}{isPolyA} eq "T" ){
        if($clusters{$loc}{final}{$tsp}{insert_size}=~ /NA/ or $clusters{$loc}{final}{$tsp}{insert_size} < 20){
         $clusters{$loc}{final}{$tsp}{status} = 0;  
         $clusters{$loc}{final}{$tsp}{remark} .= "badPair;";     
        }
      }
     }
     ## remove false copies of L1
     if($tsp eq "L1"){
      ## if no polyA tail but TSD present
      if($clusters{$loc}{final}{start}!=$clusters{$loc}{final}{end}){
       if($clusters{$loc}{final}{$tsp}{rpos} ne "-" and $clusters{$loc}{final}{$tsp}{lpos} ne "-" and $clusters{$loc}{final}{isPolyA} eq "F"){
        my $yx= Utilities::max(($clusters{$loc}{final}{$tsp}{rpos},$clusters{$loc}{final}{$tsp}{lpos}));
        if(defined $yx and abs($yx-$transposons{$tsp})>100){
         $clusters{$loc}{final}{$tsp}{status} = 0; 
         $clusters{$loc}{final}{$tsp}{remark} .= "badL1;";            
        }
       }
      }
      ## If both ends map in PolyA/T sequences
      if($clusters{$loc}{final}{$tsp}{insert_size}!~ /NA/ and $clusters{$loc}{final}{$tsp}{insert_size} < 50){
         #$clusters{$loc}{final}{$tsp}{status} = 0; 
         #$clusters{$loc}{final}{$tsp}{remark} .= "badL1;";            
      }
     }




    }
   }

   foreach my $feature (@features){
    if(exists $clusters{$loc}{final}{$feature}){
      $cntcalls{$clusters{$loc}{final}{$feature}{status}}++;
      if($clusters{$loc}{final}{$feature}{status}<$clusters{$loc}{final}{sort}){
        $clusters{$loc}{final}{sort}= $clusters{$loc}{final}{$feature}{status};
      }
    }
   }

  }
 if(!defined $cntcalls{3}){$cntcalls{3}=0;}
 if(!defined $cntcalls{2}){$cntcalls{2}=0;}
 if(!defined $cntcalls{1}){$cntcalls{1}=0;}
 if(!defined $cntcalls{0}){$cntcalls{0}=0;}
   
 return(\%clusters);
}

sub write_clusters{
   my ($clust) = shift;
   my %clusters = %{$clust};
   my @features = grep { $_ ne 'PolyA' } (keys %transposons);
  
   my @keys = sort { no warnings;
                     my ($side1) = $clusters{$a}{final}{sort};
                     my ($side2) = $clusters{$b}{final}{sort};
                     $side2 <=> $side1; } keys %clusters;

   ## (6) : Write insertion output to a file
   my $file = $outprefix.".candidate.insertions.bed";
   open FO ,"> $file" or die $!;
   my $clustnum=0;
   print FO "########################################################################################\n";
   print FO "##HiTEA report: V1.0\n";
   print FO "## Input Clusters files: ".$in."\n";
   print FO "## Description:\n";
   print FO "## Status 3: High confidence calls with clip reads on both ends ($cntcalls{3} calls) \n";
   print FO "## Status 2: Calls with either right or left hand side clip support ($cntcalls{2} calls) \n"; 
   print FO "## Status 1: Calls overlapping with known genomic copies of transposable element ($cntcalls{1} calls)\n";
   print FO "## Status 0: Poor quality/ low confidence calls ($cntcalls{0} calls) \n";
   print FO "########################################################################################\n";
   print FO "#chr\tstart\tend\tid\tscore\tstrand\tTE\tstatus\tdescription\tremark\n";
   
   foreach my $loc (@keys){
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature}){
        $clustnum++;
        ## bed entry
        print FO "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
        print FO "$clustnum\t$clusters{$loc}{final}{$feature}{anum}\t$clusters{$loc}{final}{strand}\t";
        ## type of insertion
        print FO "$feature\t";
        ## status
        print FO "$clusters{$loc}{final}{$feature}{status}\t";
        ## info related to cluster
        print FO "TSD=$clusters{$loc}{final}{TSD};SUBFAMILY=$clusters{$loc}{final}{$feature}{subfamily};";

        print FO "TE=$clusters{$loc}{final}{$feature}{lte},$clusters{$loc}{final}{$feature}{rte};";
        print FO "POLYA=$clusters{$loc}{final}{$feature}{lpoly},$clusters{$loc}{final}{$feature}{rpoly};";
        print FO "TAIL=$clusters{$loc}{final}{tailinfo};isPolyA=$clusters{$loc}{final}{isPolyA};";
        print FO "SIZE=$clusters{$loc}{final}{$feature}{insert_size};";
        ## info related to insertion
        $a=0;
        $a += $clusters{$loc}{pos1}{fuzzyclipped} if(defined $clusters{$loc}{pos1}{fuzzyclipped});
        $a += $clusters{$loc}{pos2}{fuzzyclipped} if(exists$clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{fuzzyclipped});
        print FO "UN=$clusters{$loc}{final}{$feature}{umap};FUZZY=$a;TOTALCLIP=$clusters{$loc}{final}{AllClip};";
        print FO "PURITY=$clusters{$loc}{final}{$feature}{Purity};";
        print FO "RAM=$clusters{$loc}{final}{$feature}{RAMCount};Pval=$clusters{$loc}{final}{$feature}{pvalue};";
        print FO "COVERAGE=$clusters{$loc}{final}{bgcoverage};";
        print FO "gTE=$clusters{$loc}{final}{$feature}{bgTE};D2gTE=$clusters{$loc}{final}{$feature}{bgTEDist};";
        #print FO "TECLUSTER1=$clusters{$loc}{final}{$feature}{TEClust};";
        print FO "MAXMAP=$clusters{$loc}{final}{$feature}{TESc};";
        print FO "D2RE=$clusters{$loc}{final}{dist2remotif};%FUZZY=$clusters{$loc}{final}{$feature}{FbyC};";
        print FO "N_CLUSTERS=$clusters{$loc}{pos1}{bgclusters};RECI%=$clusters{$loc}{final}{$feature}{TEClustPerc};";
        print FO "$clusters{$loc}{final}{$feature}{bgModel}\t";
        print FO "$clusters{$loc}{final}{$feature}{remark}\n";
      }
     }
   }  
   close(FO);


   print "\n";
   print "   Total clusters: ".$clustnum."\n";
   my $polyr = 0;
   if($clustnum >0){
     $polyr = Utilities::round($polyclust*100/($clustnum+$filtered),2);
   }
   print "  Total $polyclust putative polyA/T expansions out of total ",($clustnum+$filtered)," insertion candidates ($polyr%)\n";
   $polyr = 0;
   if($clustnum >0){
     $polyr = Utilities::round($unassigned*100/($clustnum+$filtered),2);
   }
   print "  Total unassigned $unassigned clusters out of total ",($clustnum+$filtered)," insertion candidates ($polyr%)\n";
   print "  Total filtered clusters by HiTEA (e.g. status 0): $filtered\n";
   print "  Clusters written in the reporting bed file: $clustnum\n";
   print "  Number of HiTEA insertions:",($cntcalls{3}+$cntcalls{2}+$cntcalls{1}),"\n";
   
   return(\%clusters);
}

