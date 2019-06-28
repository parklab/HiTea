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
use open qw(:std :utf8);
BEGIN { our $start_run = time(); }

my $in="";
my $index="";
my $outprefix="project";
my $gap = 2;
my $readlen=100;
my $min_reads=3;
my $help=0;
my $wd = "";
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'index=s'          => \$index,
  'gap=s'            => \$gap,
  'min_reads=s'      => \$min_reads,
  'readlen=s'        => \$readlen,
  'outprefix=s'      => \$outprefix,
  'wd:s'             => \$wd,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl write_insertions.pl -in [FILE_PATH] -outprefix [STRING] -gap [INT] -readlen [INT] -min_reads [INT] -index [TE_MergedIndex] -wd [DIR_PATH]\n\n";
  print "This program annotates finalized breaks using the bam file\n\n";
  print "Options:\n\n";
  print "***required:\n";
  print "  -in                    Clusters file in bed format \n";
  print "  -index                 transposable element reference assembly (fasta path) \n";
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
print " ERROR: makeMat.R does not exist in the same directory!\n" unless -e "makeMat.R";

#--------------------------------------------------------------------------------------------------------------
# I/O
#--------------------------------------------------------------------------------------------------------------

my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "[write_insertions] START:\t $run_time seconds\n";
print " Command: perl write_insertions.pl -in $in -gap $gap -min_reads $min_reads -readlen $readlen -outprefix $baseoutprefix -wd $wd \n";

%cov = %{retrieve($outprefix.".coverage.ph")};

## get transposons
%transposons = get_fasta($index);
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


#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub filter_clusters{
  my ($clust) = shift;
  my %clusters = %{$clust};
  
  foreach my $loc(sort keys %clusters){

     if(!defined($clusters{$loc}{pos2}{chr}) or !defined $clusters{$loc}{final}{TE_mapReads} ) {
        print  " ERROR: Something is not defined in the insertion calls\n";
        exit 1;
     }

     ## Clusters without TE assignment and STRs
     if(!defined$clusters{$loc}{final}{bgTE} or !defined$clusters{$loc}{final}{bgTEDist}){
        $clusters{$loc}{final}{bgTE}="-";
        $clusters{$loc}{final}{bgTEDist}=0;       
        if($clusters{$loc}{final}{polyAT} eq "PolyA,PolyA" or $clusters{$loc}{final}{polyAT} eq "PolyT,PolyT"){
          #print " PolyA cluster: $loc\t$clusters{$loc}{final}{PolyA}\n";  
          if($clusters{$loc}{final}{PolyA}>0){
            my $check=0;
            foreach my $te(keys %transposons){
              if($te ne "PolyA" and $te ne "*" and $clusters{$loc}{final}{$te}>0){
               $check++;
              }
            }
            if($check ==0){
              $clusters{$loc}{final}{TE} = "PolyAT-Expansion";
              $polyclust++;
            }
          }                 
        }else{
          $unassigned++;
        }
     }     
     
     ## Dump cluster due to low coverage 
     if($clusters{$loc}{final}{TE_mapReads} < 1){
       print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\tNo assignment (write step)\tMapReads:$clusters{$loc}{final}{subfamily};$clusters{$loc}{final}{AllClip}\t$clusters{$loc}{final}{subfamily}\n";
       $filtered++;
       delete($clusters{$loc});
       next;
     }   
     ## Dump cluster due to proximity to RE 
     if($clusters{$loc}{final}{dist2re} < 2){
        print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\tTE-proximity (write step)\tdist2RE:$clusters{$loc}{final}{dist2re}\t$clusters{$loc}{final}{subfamily}\n";
        delete($clusters{$loc});
        $filtered++;
        next;
     }
 
     # Filtering the clusters ==Assign confidence to the clsuter
     #h <- h[h$V21>10 & h$V13>0.1 & h$V11>4 & h$fuzzy < 0.5*h$score & h$V22 <0.5 & ((h$rhs>1 & h$lhs>1)| (h$V12)>=35)]
     ## Score 3: high conf calls,  lhs, rhs clip reads present and >10% clip reads map to TE
     ## Score 2: optimal,          either lhs or rhs read support and >10% clip reads map to TE
     ## Score 1: Overlaps with bg, The copy overlaps with background TE copy 
     ## Score 0: poor evidence,    Based on several filtering criteria 
     my $purity_ratio = $clusters{$loc}{final}{purity};
     my $tsp =  $clusters{$loc}{final}{TE};
     my $pval = $clusters{$loc}{final}{pvalue};
     $clusters{$loc}{final}{status} = 0;
     $clusters{$loc}{final}{remark} = "-";
     if($clusters{$loc}{final}{TE} eq "PolyAT-Expansion"){
        $clusters{$loc}{final}{status} = 0 ;  #
        $clusters{$loc}{final}{remark} .= "STR;";         
        next;
     }
     
     if($purity_ratio>=$PURITY_CUTOFF and $pval<$pval_cutoff){
       if($clusters{$loc}{pos1}{$tsp}>=1 and $clusters{$loc}{pos2}{$tsp}>=1){
         $clusters{$loc}{final}{status} = 3; 
       }elsif((($clusters{$loc}{pos1}{$tsp}>=1 and $clusters{$loc}{pos2}{PolyA}>=1) or ($clusters{$loc}{pos2}{$tsp}>=1 and $clusters{$loc}{pos1}{PolyA}>=1))){
         $clusters{$loc}{final}{status} = 3; 
       }elsif($clusters{$loc}{final}{$tsp}>=2 and $clusters{$loc}{final}{maxscore}>=30){
         $clusters{$loc}{final}{status} = 2; 
       }elsif($clusters{$loc}{final}{$tsp}>=1 and $clusters{$loc}{final}{maxscore}>=30 and $clusters{$loc}{final}{PolyA}>=2){
         $clusters{$loc}{final}{status} = 2; 
       }
     }
     
     if($clusters{$loc}{final}{maxscore} <30){
       $clusters{$loc}{final}{status} = 0; 
       $clusters{$loc}{final}{remark} .= "PoorTEMapScore;";               
     }

     if($pval >=$pval_cutoff){
       $clusters{$loc}{final}{status} = 0; 
       $clusters{$loc}{final}{remark} .= "NoEnrichment;";         
     }

     $clusters{$loc}{final}{FbyC} = Utilities::round($clusters{$loc}{final}{fuzzyclipped}/$clusters{$loc}{final}{AllClip},2);
     if($clusters{$loc}{final}{FbyC} >= 0.25 and $clusters{$loc}{final}{fuzzyclipped}>5){
        $clusters{$loc}{final}{status} = 0 ;  #if there are >50% fuzzy clipped reads in the vicinity, flag a cluster as Category 0 
        $clusters{$loc}{final}{remark} .= "fuzzy;";         
     }elsif($clusters{$loc}{pos1}{bgclusters} > 5){ ## some random threshold suggesting that a cluster is bad if there are atleast 5 other clusters in the vicinity. In otehr words, this is a bad locus
       $clusters{$loc}{final}{status} = 0;
       $clusters{$loc}{final}{strand} ="*";
       $clusters{$loc}{final}{remark} .= "fuzzy;";         
     }

     if($tsp eq "*"){
       $clusters{$loc}{final}{status} = 0; 
       $clusters{$loc}{final}{remark} .= "NoAssignment;";               
     }
     if($purity_ratio <$PURITY_CUTOFF){
       $clusters{$loc}{final}{status} = 0; 
       $clusters{$loc}{final}{remark} .= "InsuffClipCov;";               
     }
     if($clusters{$loc}{final}{dist2remotif}<$gap){
         $clusters{$loc}{final}{status} = 0;
         $clusters{$loc}{final}{remark} .= "RE;";                             
     }

     ## TEcluster should have at least 50% of the reads forming a cluster (2/3 or 6/10)
     if($clusters{$loc}{final}{RecrpClustPercent} < 50){
        $clusters{$loc}{final}{status} = 0; 
        $clusters{$loc}{final}{remark} .= "poorTECluster;";
     }
     ## if the coverage of clip reads at hte cluster is <5% coverage in a 5kb window, assign the insertion as class 0
     my $cstart = int ($clusters{$loc}{final}{start}/1000);
     my $cchr = $clusters{$loc}{final}{chr};
     if(defined $cstart){
        my $bgcov=0;
        for(my $i=-2;$i<=2;$i++){
           $cov{$cchr}{($cstart+$i)} = 0 if(!$cov{$cchr}{($cstart+$i)});
           $bgcov +=$cov{$cchr}{($cstart+$i)};
        }
        $bgcov = $bgcov/5; 
        $bgcov = Utilities::round( ($bgcov*$readlen)/10000,2); # 10th of coverage 
        $clusters{$loc}{final}{bgcoverage} = $bgcov;                         
        if($clusters{$loc}{final}{TE_mapReads} < ($bgcov) ){
             $clusters{$loc}{final}{status} = 0;
             $clusters{$loc}{final}{remark} .= "<5xcov;";   
        } 
        if( ($clusters{$loc}{final}{RAMCount}) < int($bgcov/2)){ ## 20th of bgcov for RAM
          $clusters{$loc}{final}{status} = 0;
          $clusters{$loc}{final}{remark} .= "RAM<".$RAMCountCUTOFF.";";         
        }
     }
     if( ($clusters{$loc}{final}{RAMCount}) ==0){
       $clusters{$loc}{final}{status} = 0;
       $clusters{$loc}{final}{remark} .= "RAM<1;";         
     }


     
     if($clusters{$loc}{final}{bgTEDist}<=2){
       $clusters{$loc}{final}{remark} .= "bgOlap;";
       if($clusters{$loc}{final}{status}==3){  ## make status 1 of calls if original status is 3 and it overlaps with bgTE
         $clusters{$loc}{final}{status} = 1;  
       }elsif($clusters{$loc}{final}{status}==2){  ## make status 1 of calls if original status is 3 and it overlaps with bgTE
         $clusters{$loc}{final}{status} = 0;  
       }else{
         $clusters{$loc}{final}{status} = 0;  
       }
     }

     ## if LHS and RHS map within 10bp of each other on TE consensus
     if($clusters{$loc}{pos1}{cliploc}!=$clusters{$loc}{pos2}{cliploc} and $clusters{$loc}{pos1}{clip_sidedness} ne $clusters{$loc}{pos2}{clip_sidedness}){
        if($clusters{$loc}{pos1}{RecrpPolyA} ne "TRUE" and $clusters{$loc}{pos2}{RecrpPolyA} ne "TRUE"){
          if($tsp eq $clusters{$loc}{pos1}{RecrpTE}){
            if($clusters{$loc}{pos1}{RecrpClustLoc} ne "-" and $clusters{$loc}{pos2}{RecrpClustLoc} ne "-"){
              if(abs($clusters{$loc}{pos1}{RecrpClustLoc}-$clusters{$loc}{pos2}{RecrpClustLoc})<10){
                $clusters{$loc}{final}{status} = 0; 
                $clusters{$loc}{final}{remark} .= "singularMap;";
              }
            }
          }  
        }
     }
     
     ## L1: If LHS n RHS present and both mapaway from 3' end, filter
     if($tsp eq "L1"){
      if($clusters{$loc}{pos1}{cliploc}!=$clusters{$loc}{pos2}{cliploc} and $clusters{$loc}{pos1}{clip_sidedness} ne $clusters{$loc}{pos2}{clip_sidedness}){
       if($clusters{$loc}{pos1}{RecrpPolyA} ne "TRUE" and $clusters{$loc}{pos2}{RecrpPolyA} ne "TRUE"  and $clusters{$loc}{final}{polyAT} eq "NA"){
         if($tsp eq $clusters{$loc}{pos1}{RecrpTE} and $tsp eq $clusters{$loc}{pos2}{RecrpTE}){
           my $yx= Utilities::max(($clusters{$loc}{pos1}{RecrpClustLoc},$clusters{$loc}{pos2}{RecrpClustLoc}));
           if(abs($yx-$transposons{$tsp})>100){
               $clusters{$loc}{final}{status} = 0; 
               $clusters{$loc}{final}{remark} .= "badL1;";            
           }
         }
        }
       }
     }     
  }
  return(\%clusters);
}

sub write_clusters{
   my ($clust) = shift;
   my %clusters = %{$clust};
  
   my %cntcalls;
   foreach my $loc(sort keys %clusters){ $cntcalls{$clusters{$loc}{final}{status}}++; }
   if(!defined $cntcalls{3}){$cntcalls{3}=0;}
   if(!defined $cntcalls{2}){$cntcalls{2}=0;}
   if(!defined $cntcalls{1}){$cntcalls{1}=0;}
   if(!defined $cntcalls{0}){$cntcalls{0}=0;}
   my @keys = sort { no warnings;
                     my ($side1) = $clusters{$a}{final}{status};
                     my ($side2) = $clusters{$b}{final}{status};
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
      
   foreach my $loc(@keys){
     my $mei = $clusters{$loc}{final}{TE};
     if($mei eq "*" or $mei eq "PolyAT-Expansion"){
      $mei = "PolyA";
     }
     $clustnum++;
     ## bed entry
     print FO "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
     print FO "$clustnum\t$clusters{$loc}{final}{TE_mapReads}\t$clusters{$loc}{final}{strand}\t";
     ## type of insertion
     print FO "$clusters{$loc}{final}{TE}\t";
     ## status
     print FO "$clusters{$loc}{final}{status}\t";
     ## info related to cluster
     print FO "TSD=$clusters{$loc}{final}{TSD};SUBFAMILY=$clusters{$loc}{final}{subfamily};";
     if($mei eq "PolyA"){
       print FO "TE=0,0;";    
     }else{
       print FO "TE=$clusters{$loc}{pos1}{$mei},$clusters{$loc}{pos2}{$mei};";
     }
     print FO "POLYA=$clusters{$loc}{pos1}{PolyA},$clusters{$loc}{pos2}{PolyA};";
     print FO "TAIL=$clusters{$loc}{final}{tailinfo};";
     print FO "INSERT_SIZE=$clusters{$loc}{final}{insert_size};";
     print FO "UN=$clusters{$loc}{final}{unmapped};FUZZY=$clusters{$loc}{final}{fuzzyclipped};TOTALCLIP=$clusters{$loc}{final}{AllClip};";
     print FO "PURITY=$clusters{$loc}{final}{purity};";
     print FO "RAM=$clusters{$loc}{final}{RAMCount};Pval=$clusters{$loc}{final}{pvalue};";
     $clusters{$loc}{final}{bgcoverage} = $clusters{$loc}{final}{bgcoverage}*10;
     print FO "COVERAGE=$clusters{$loc}{final}{bgcoverage};";
     ## info related to insertion
     print FO "gTE=$clusters{$loc}{final}{bgTE};D2gTE=$clusters{$loc}{final}{bgTEDist};";
     print FO "TECLUSTER1=$clusters{$loc}{pos1}{TEClust};TECLUSTER2=$clusters{$loc}{pos2}{TEClust};MAXMAP=$clusters{$loc}{final}{maxscore};"; 
     print FO "D2RE=$clusters{$loc}{final}{dist2remotif};%FUZZY=$clusters{$loc}{final}{FbyC};";
     print FO "N_CLUSTERS=$clusters{$loc}{pos1}{bgclusters};RECI%=$clusters{$loc}{final}{RecrpClustPercent};";
     print FO "CLUSTReads=$clusters{$loc}{pos1}{RecrpClustReads},$clusters{$loc}{pos2}{RecrpClustReads},$clusters{$loc}{pos1}{RecrpPolyA},$clusters{$loc}{pos2}{RecrpPolyA};";
     print FO "$clusters{$loc}{final}{bgModel}\t";
     print FO "$clusters{$loc}{final}{remark}\n";
   }
   close(FO);


   print "\n";
   #print "  ###----------------------------------------------\n";
   #print "  #Summary\n";
   #print "  ###----------------------------------------------\n";
   print "   Total clusters: ".scalar(keys%clusters)."\n";
   my $polyr = 0;
   if(scalar(keys %clusters) >0){
     $polyr = Utilities::round($polyclust*100/scalar(keys %clusters),2);
   }
   #print "  Total $polyclust putative polyA/T expansions out of total ",scalar(keys %clusters)," insertion candidates ($polyr%)\n";
   $polyr = 0;
   if(scalar(keys %clusters) >0){
     $polyr = Utilities::round($unassigned*100/scalar(keys %clusters),2);
   }
   print "   Total $unassigned non-assigned clusters out of total ",scalar(keys %clusters)," insertion candidates ($polyr%)\n";
   print "   Filtered based on coverage/RE proximity etc: $filtered\n";
   print "   Final reported clusters: $clustnum\n";

   ## (7): cleanup and Return final clusters object and save
   return(\%clusters);
}

