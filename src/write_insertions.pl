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
use src::Utilities;
require "src/vars.pl";
our (%redb);
our (%chrs);
our $PVAL_CUTOFF;
our $RAMCountCUTOFF;
our $PURITY_CUTOFF;
my $INDEL_CUTOFF = 50;
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
use open qw(:std :utf8);
BEGIN { our $start_run = time(); }

my $in="";
my $index="";
my $polym="";
my $refMapqQ="";
my $outprefix="project";
my $gap = 2;
my $enzyme = "";
my $readlen=100;
my $min_reads=3;
my $help=0;
my $wd = "";
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'index=s'          => \$index,
  'polym=s'          => \$polym,
  'gap=s'            => \$gap,
  'min_reads=s'      => \$min_reads,
  'readlen=s'        => \$readlen,
  'e=s'              => \$enzyme,
  'q=s'              => \$refMapqQ,
  'outprefix=s'      => \$outprefix,
  'wd:s'             => \$wd,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
sub help{
  my $j = shift;
  if($j){
   print "\nUsage: perl write_insertions.pl -in [FILE_PATH] -outprefix [STRING] -gap [INT] -readlen [INT] -e [STRING] -q [referenceMAPQ] -min_reads [INT] -index [TE_MergedIndex] -polym [TE_MergedIndex] -wd [DIR_PATH]\n\n";
   print "This program annotates finalized breaks using the bam file\n\n";
   print "Options:\n\n";
   print "***required:\n";
   print "  -in                    Clusters file in bed format \n";
   print "  -index                 transposable element reference assembly (fasta path) \n";
   print "  -q                     reference mapq score threshold \n";
   print "  -e                     RE enzyme \n";
   print "***optional:\n";
   print "  -polym                 transposable element reference assembly (fasta path) \n";
   print "  -min_reads             minimum reads that need to map at an insertion site [default: 3] \n";
   print "  -gap                   gap to merge the breakpoints [default: 2] \n";
   print "  -readlen               sequencing read length [default: 100] \n";
   print "  -outprefix             Outputprefix for generating 2 output files [default: project] \n";
   print "  -wd                    Working directory [default: ~]\n";
   print "  -help|-h               Display usage information.\n\n\n";
   exit 1;
  }
}

$wd =~ s/\/$//;
my $baseoutprefix =$outprefix;
$outprefix = $wd."/".$outprefix if($wd ne "");
if($help or $in eq "" or $index eq "" or $refMapqQ eq "" or $enzyme eq ""){
  print " One or more required inputs are not recognized***\n";
  help(1);
}
print " ERROR: makeMat.R does not exist in the same directory!\n" unless -e "src/makeMat.R";

#--------------------------------------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------------------------------------
my $polyclust=0;
my $unassigned = 0;
my $filtered = 0;
my %clusters;
my %transposons;
my %cov;
my %cntcalls;
my %TELOCS;

my %tmp;
%tmp = src::Utilities::get_fasta_seqs($index);
foreach my $te (keys %tmp){
  my @dist;
  while ($tmp{$te} =~ /$redb{$enzyme}{motif}/g){
    push @dist, $-[0];
    push @dist, $+[0];
  }
  @dist = sort {$a <=> $b} @dist;
  $TELOCS{$te} = \@dist if(scalar @dist>1);
}
if($polym ne ""){
  %tmp = src::Utilities::get_fasta_seqs($polym);
  foreach my $te (keys %tmp){
    my @dist;
    while ($tmp{$te} =~ /$redb{$enzyme}{motif}/g){
      push @dist, $-[0];
      push @dist, $+[0];
    }
    @dist = sort {$a <=> $b} @dist;
    $te=~ s/^\S*~//;
    $TELOCS{$te} = \@dist if(scalar @dist>1);
  }
}
#--------------------------------------------------------------------------------------------------------------
# I/O
#--------------------------------------------------------------------------------------------------------------
my $olog = $outprefix.".skippedclusters.logs.gz";
open(LOGS,"| gzip -c - > $olog") or die $!;

my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print "\n\n[write_insertions] START:\t $run_time seconds\n";
print " Command: perl write_insertions.pl -in $in -index $index -polym $polym -q $refMapqQ -gap $gap -e $enzyme -min_reads $min_reads -readlen $readlen -outprefix $baseoutprefix -wd $wd \n\n";


%cov = %{retrieve($outprefix.".coverage.ph")};
print " RE index generated on following TE sequences: ", join(",",keys %TELOCS),"\n";

## get transposons
%transposons = src::Utilities::get_fasta($index);
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

my $oufile=$outprefix.'.finalFilt'; #Save 
open FOO,">$oufile" or die $!;
print FOO Dumper %clusters;
close(FOO);
print " wrote insertion breakpoint report:\t $run_time seconds\n";


## generate matrix
my $testring = join(",",keys %transposons);
$testring =~ s/PolyA//;
$testring=~ s/,,/,/;
$testring=~ s/^,//;
$testring=~ s/,$//;
system( qq[Rscript src/makeMat.R $baseoutprefix $wd $testring] ) == 0 or die qq[ Cound not generate count matrices for coverage plots\n];
$watch_run = time();
$run_time = $watch_run - $start_run;
print " generated coverage matrices of RAM and split RAM around the insertion breakpoints:\t $run_time seconds\n";


## genrate HTML report
if(system( qq[Rscript src/createReport.R $baseoutprefix $wd] ) == 0){
  print " generated HTML report successfully\n";
}else{
  print " Cound not generate HTML report. One or more R-packages are missing\n";
}

$watch_run = time();
$run_time = $watch_run -  $start_run;
print "[write_insertions] END:\t $run_time seconds\n\n";
close(LOGS);

exit 0;

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

    ## Dump clusters due polyA mapping alone
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
    ## Dump clusters with <2 AllClip reads
    if($clusters{$loc}{final}{AllClip}<2){
      print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
      print LOGS "<2Clips\tAllreads:$clusters{$loc}{final}{AllClip}\tD2RE:$clusters{$loc}{final}{AllClip}\n";
      $filtered++;       
      delete $clusters{$loc};
      next;      
    }
    ## Dump clusters due proximity to RE
    if($clusters{$loc}{final}{dist2re}<2){
      print LOGS "$clusters{$loc}{final}{chr}\t$clusters{$loc}{final}{start}\t$clusters{$loc}{final}{end}\t";
      print LOGS "<proximity to RE (write step)\tAllreads:$clusters{$loc}{final}{AllClip}\tD2RE:$clusters{$loc}{final}{dist2re}\n";
      $filtered++;       
      delete $clusters{$loc};
      next;
    }
    ## Dump clusters if TE-mapping clip reads are <2
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

    ## background coverage
    my $cstart = int (($clusters{$loc}{final}{start}/1000)+0.5);
    my $cchr = $clusters{$loc}{final}{chr};
    if(defined $cstart){
      my $bgcov=0;
      for(my $i=-1;$i<=1;$i++){
        $cov{$cchr}{($cstart+$i)} = 0 if(!$cov{$cchr}{($cstart+$i)});
        $bgcov +=$cov{$cchr}{($cstart+$i)};
      }
      $bgcov = $bgcov/3; 
      $bgcov = src::Utilities::round( ($bgcov*$readlen)/1000,2); 
      $clusters{$loc}{final}{bgcoverage} = $bgcov;                         
    }  
    $clusters{$loc}{final}{sort}=0;   
   
    ## STR?
    $clusters{$loc}{final}{STR} = "F";
    if($clusters{$loc}{final}{tailinfo} eq "PolyAPolyA" or $clusters{$loc}{final}{tailinfo} eq "PolyAPolyT" or $clusters{$loc}{final}{tailinfo} eq "PolyTPolyA" or $clusters{$loc}{final}{tailinfo} eq "PolyTPolyT"){
       $clusters{$loc}{final}{STR} = "T";
    }
 
    ## only report the clsuter with maximum number of read support at a given locus
    my %tmp;
    foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature}){
        $tmp{$feature}{RAM} = $clusters{$loc}{final}{$feature}{RAMCount}+$clusters{$loc}{final}{$feature}{tnum};
        $tmp{$feature}{tnum} = $clusters{$loc}{final}{$feature}{tnum};
        $tmp{$feature}{pval} = $clusters{$loc}{final}{$feature}{pvalue};
        $tmp{$feature}{Purity} = $clusters{$loc}{final}{$feature}{Purity};
      }
    }
    my @tmpkeys = sort { $tmp{$b}{RAM} <=> $tmp{$a}{RAM} or $tmp{$b}{tnum} <=> $tmp{$a}{tnum} or $tmp{$a}{pval} <=> $tmp{$b}{pval} or $tmp{$b}{Purity} <=> $tmp{$a}{Purity}} keys %tmp;
    if(scalar @tmpkeys >0){
      $clusters{$loc}{final}{ReportedTE} = $tmpkeys[0];
      $clusters{$loc}{final}{ReportedTEScore} = $tmp{$tmpkeys[0]};
    }

    ### Providing status calls the clusters
    # Filtering the clusters ==Assign confidence to the clsuter
    ## Score 3: high conf calls,  lhs, rhs clip reads present and >10% clip reads map to TE
    ## Score 2: optimal,          either lhs or rhs read support and >10% clip reads map to TE
    ## Score 1: Overlaps with bg, The copy overlaps with background TE copy 
    ## Score 0: poor evidence,    Based on several filtering criteria 
    foreach my $feature (@features){
      
      if(exists $clusters{$loc}{final}{$feature}){     
        my $purity_ratio = $clusters{$loc}{final}{$feature}{Purity};
        my $tsp =  $feature;
        my $pval = $clusters{$loc}{final}{$feature}{pvalue};
        $clusters{$loc}{final}{$feature}{status} = 0;
        $clusters{$loc}{final}{$feature}{remark} = "-";
      
       ## If LHS and RHS map to two distinct TEs (without PolyA), change the start/end of the cluster
       if($clusters{$loc}{final}{start}!=$clusters{$loc}{final}{end}){
          my %s;
          my %e;
          foreach my $f ( (@features,"PolyA") ){
            $s{$f} = $clusters{$loc}{pos1}{cliploc} if(exists $clusters{$loc}{pos1}{$f});
            $e{$f} = $clusters{$loc}{pos2}{cliploc} if(exists $clusters{$loc}{pos2}{$f});
          }
          if(!$s{"PolyA"} and !$e{"PolyA"}){
            if($s{$feature} and !$e{$feature} and scalar(keys %e)>0){
              $clusters{$loc}{final}{start}=$s{$feature};
              $clusters{$loc}{final}{end}= $s{$feature};     
              $clusters{$loc}{final}{identClip}=$clusters{$loc}{pos1}{identClip};
            }elsif($e{$feature} and !$s{$feature} and scalar(keys %s)>0){
              $clusters{$loc}{final}{start}=$e{$feature};
              $clusters{$loc}{final}{end}= $e{$feature};     
              $clusters{$loc}{final}{identClip}=$clusters{$loc}{pos2}{identClip};
            }
          }
       } 
       
       ## status assignment
       if($purity_ratio >=$PURITY_CUTOFF and $pval<$PVAL_CUTOFF){
          if(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{lte}>=1 and $clusters{$loc}{final}{$tsp}{rte}>=1){
            $clusters{$loc}{final}{$tsp}{status} = 3; 
          }elsif(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{lte}>=1 and $clusters{$loc}{final}{$tsp}{rpoly}>=1){
            $clusters{$loc}{final}{$tsp}{status} = 3; 
          }elsif(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{rte}>=1 and $clusters{$loc}{final}{$tsp}{lpoly}>=1){
            $clusters{$loc}{final}{$tsp}{status} = 3; 
          }elsif(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{lte}>=1 and exists $clusters{$loc}{pos2} and $clusters{$loc}{pos2}{isPolyA} eq "T"){
            $clusters{$loc}{final}{$tsp}{status} = 3; 
          }elsif(exists $clusters{$loc}{final}{$tsp} and $clusters{$loc}{final}{$tsp}{rte}>=1 and $clusters{$loc}{pos1}{isPolyA} eq "T"){
            $clusters{$loc}{final}{$tsp}{status} = 3; 
          }elsif($clusters{$loc}{final}{$tsp}{tnum}>=2 and $clusters{$loc}{final}{$tsp}{TESc}>=$TESCORE){
            $clusters{$loc}{final}{$tsp}{status} = 1; 
          }
       }      
       ## Overlap with reference copy
       if($clusters{$loc}{final}{$tsp}{bgTEDist}<=2){
          $clusters{$loc}{final}{$tsp}{remark} .= "bgOlap;";
          if($clusters{$loc}{final}{$tsp}{status}==3){  ## make status 1 of calls if original status is 3 and it overlaps with bgTE
            $clusters{$loc}{final}{$tsp}{status} = 2;  
          }elsif($clusters{$loc}{final}{$tsp}{status}==1){  ## make status 1 of calls if original status is 3 and it overlaps with bgTE
            $clusters{$loc}{final}{$tsp}{status} = 0;  
          }else{
            $clusters{$loc}{final}{$tsp}{status} = 0;  
          }
       }
       ## TE map score and refmapscore of clipped reads
       if($clusters{$loc}{final}{AllClip}>=2){
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
            $ATcnt = src::Utilities::round(length($keys[0])/length($clusters{$loc}{final}{TSD}),5) if(defined $keys[0]);
          }
          $clusters{$loc}{final}{ATpercent} = $ATcnt;
          my $x = $clusters{$loc}{final}{$tsp}{OriginalTE};
          if(exists $clusters{$loc}{pos1}{$x} and $clusters{$loc}{pos1}{$x}{TESc}>=30){
          }elsif(exists $clusters{$loc}{pos2} and exists$clusters{$loc}{pos2}{$x} and $clusters{$loc}{pos2}{$x}{TESc}>=30 ){ 
          }elsif(exists $clusters{$loc}{pos1}{PolyA} and $clusters{$loc}{pos1}{PolyA}{TESc}>=30 and $ATcnt < 0.7){
          }elsif(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{PolyA} and $clusters{$loc}{pos2}{PolyA}{TESc}>=30 and $ATcnt < 0.7){
          }else{
            $clusters{$loc}{final}{$tsp}{status} = 0; 
            $clusters{$loc}{final}{$tsp}{remark} .= "<30MapSc;";               
          }

          if(exists $clusters{$loc}{pos1}{$x} and $clusters{$loc}{pos1}{$x}{mapq} >= 5){
          }elsif(exists $clusters{$loc}{pos2} and exists$clusters{$loc}{pos2}{$x} and $clusters{$loc}{pos2}{$x}{mapq} >= 5 ){ 
          }else{ 
            $clusters{$loc}{final}{$tsp}{status} = 0; 
            $clusters{$loc}{final}{$tsp}{remark} .= "<5-refMapq;";               
          }        
       }
       ## Pval
       if($pval >=$PVAL_CUTOFF){
         $clusters{$loc}{final}{$tsp}{status} = 0; 
         $clusters{$loc}{final}{$tsp}{remark} .= "NoEnrichment;";         
       }
       ## Purity cutoff
       if($purity_ratio <$PURITY_CUTOFF){
         $clusters{$loc}{final}{$tsp}{status} = 0; 
         $clusters{$loc}{final}{$tsp}{remark} .= "InsuffClipCov;";               
       }
       ## 1 or less TE-consensus mapping clusters
       if($clusters{$loc}{final}{$tsp}{tnum}>=1){
         my $tot=0;
         $tot += $clusters{$loc}{final}{$tsp}{lte} if($clusters{$loc}{final}{$tsp}{lpoly}==0);
         $tot += $clusters{$loc}{final}{$tsp}{rte} if($clusters{$loc}{final}{$tsp}{rpoly}==0);
         if($tot<1 or ($tot<2 and $clusters{$loc}{final}{$tsp}{anum}>20)){
            $clusters{$loc}{final}{$tsp}{status} = 0;
            $clusters{$loc}{final}{$tsp}{remark} .= "<2TE-mapping;";   
         }
       }       
       ############################################
       ## check proximity to the RE motif (not ligation motif, but genomic motif)
         ### (0) based on counted distance to ligation motif
       if($clusters{$loc}{final}{dist2remotif}<$gap){
         $clusters{$loc}{final}{$tsp}{status} = 0;
         $clusters{$loc}{final}{$tsp}{remark} .= "+LigMotif;";                             
       }    
         ### (1) check based on the provided read-sequence
       if($clusters{$loc}{final}{AllClip}>=2){
         my $dist =500;
         my $hpoly = "*";
         my $HPOLYCUTOFF = 95;
           
         if(defined $clusters{$loc}{pos1}{seqstring}){
           my @rr = split(",",$clusters{$loc}{pos1}{seqstring});
           my ($a) = $rr[1]=~ /^(\d+)S\S+/;
           my ($b) = $rr[1]=~ /\D(\d+)S$/;
           $a=0 if(!defined($a) || $a eq"");  ## left hand side
           $b=0 if(!defined($b) || $b eq"");  ## right hand side   
           my $xstart = $clusters{$loc}{final}{start};
           my $xend = $clusters{$loc}{final}{end}-length($rr[0])+$b;
           if( ($rr[3]-$gap) <= $xstart and ($rr[3]+$gap) >= $xstart   and $clusters{$loc}{pos1}{clip_sidedness} eq "a"){
             $dist = src::Utilities::checkREatClip($rr[0],$redb{$enzyme}{motif},$a);
             $hpoly = src::Utilities::hompolymer_check($rr[0],$a,$HPOLYCUTOFF);
           }elsif( ($rr[3]-$gap) <= $xend and ($rr[3]+$gap) >= $xend and $clusters{$loc}{pos1}{clip_sidedness} eq "b"){
             $dist = src::Utilities::checkREatClip($rr[0],$redb{$enzyme}{motif},(length($rr[0])-$b)); 
             $hpoly = src::Utilities::hompolymer_check($rr[0],(length($rr[0])-$b),$HPOLYCUTOFF); 
           }
         }
         if(exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{seqstring}){
           my $diff = 500;
           my $dhpoly="*";
           my @rr = split(",",$clusters{$loc}{pos2}{seqstring});
           my ($a) = $rr[1]=~ /^(\d+)S\S+/;
           my ($b) = $rr[1]=~ /\D(\d+)S$/;
           $a=0 if(!defined($a) || $a eq"");  ## left hand side
           $b=0 if(!defined($b) || $b eq"");  ## right hand side   
           my $xstart = $clusters{$loc}{final}{start};
           my $xend = $clusters{$loc}{final}{end}-length($rr[0])+$b;
           if( ($rr[3]-$gap) <= $xend and ($rr[3]+$gap) >= $xend and $clusters{$loc}{pos2}{clip_sidedness} eq "b"){
             $diff = src::Utilities::checkREatClip($rr[0],$redb{$enzyme}{motif},(length($rr[0])-$b));
             $dhpoly = src::Utilities::hompolymer_check($rr[0],(length($rr[0])-$b),$HPOLYCUTOFF);
           }elsif( ($rr[3]-$gap) <= $xstart and ($rr[3]+$gap) >= $xstart and $clusters{$loc}{pos2}{clip_sidedness} eq "a"){
             $diff = src::Utilities::checkREatClip($rr[0],$redb{$enzyme}{motif},$a);
             $dhpoly = src::Utilities::hompolymer_check($rr[0],$a,$HPOLYCUTOFF);
           }
           $dist = $diff if($diff > $dist);
           $hpoly = $dhpoly if($hpoly eq "*");
         }
         $clusters{$loc}{final}{TEproximity} = $dist;
         $clusters{$loc}{final}{hompolymer} = $hpoly;
         if($clusters{$loc}{final}{TEproximity} <= 3){
           $clusters{$loc}{final}{$tsp}{status} = 0; 
           $clusters{$loc}{final}{$tsp}{remark} .= "refRE+;";                   
         }
         if($clusters{$loc}{final}{hompolymer} ne "*"){
           $clusters{$loc}{final}{$tsp}{status} = 0; 
           $clusters{$loc}{final}{$tsp}{remark} .= $clusters{$loc}{final}{hompolymer}.";";                   
         }                  
       }   
         ### (2)check using the index generated from TE-consensus
       if(defined $clusters{$loc}{final}{$tsp}{OriginalTE}){
         my $x = $clusters{$loc}{final}{$tsp}{OriginalTE};
         $x=~ s/\r//; 
         $x=~ s/\s+//;
         my $dist=500;
         my $diff=500;
         if(exists $TELOCS{$x}){
           if(exists $clusters{$loc}{pos1}{$x} and $clusters{$loc}{pos1}{$x}{TEClustPos} ne "-"){
             my $lll = $clusters{$loc}{pos1}{$x}{TEClustPos};
             $dist = abs($lll - src::Utilities::binarysearch($lll,$x,$TELOCS{$x}) ) if(exists $TELOCS{$x}); 
           }
           if(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$x} and $clusters{$loc}{pos2}{$x}{TEClustPos} ne "-"){
             my $lll = $clusters{$loc}{pos2}{$x}{TEClustPos};
             $diff = abs($lll - src::Utilities::binarysearch($lll,$x,$TELOCS{$x}) ) if(exists $TELOCS{$x}); 
           }
           if( ($dist <= 3 or $diff <= 3) and $clusters{$loc}{final}{start}==$clusters{$loc}{final}{end}){
             $clusters{$loc}{final}{TEproximity} = src::Utilities::min(($dist,$diff));
             $clusters{$loc}{final}{$tsp}{status} = 0; 
             $clusters{$loc}{final}{$tsp}{remark} .= "consRE+;";                             
           }elsif($dist <= 3 and $diff <= 3 and $clusters{$loc}{final}{start}!=$clusters{$loc}{final}{end}){
             $clusters{$loc}{final}{TEproximity} = src::Utilities::min(($dist,$diff));
             $clusters{$loc}{final}{$tsp}{status} = 0; 
             $clusters{$loc}{final}{$tsp}{remark} .= "consRE+;";                                         
           }
         }
       }
       ############################################
       ## at least 1/3 reads mapping to TE/PolyA
       if($clusters{$loc}{final}{$tsp}{umap} >= src::Utilities::round($clusters{$loc}{final}{AllClip}*0.67,0) ){
          $clusters{$loc}{final}{$tsp}{status} = 0;
          $clusters{$loc}{final}{$tsp}{remark} .= ">2/3Unmap;";   
       }
       ## double clipped reads
       if(exists $clusters{$loc}{pos2}){
         my $chk=0;
         my $x = $clusters{$loc}{pos2}{cliploc};
         for(my $i=-1*$gap; $i<=$gap;$i++){
          my $y = $x+$i;
          $chk++ if($clusters{$loc}{pos1}{both} =~ m/$y/); 
         }
         $x = $clusters{$loc}{pos1}{cliploc};
         for(my $i=-1*$gap; $i<=$gap;$i++){
          my $y = $x+$i;
          $chk++ if($clusters{$loc}{pos2}{both} =~ m/$y/); 
         } 
         if($chk>0){ ## two distict TE,TE or TE,Poly mapping clusters
           if($clusters{$loc}{pos1}{bgclusters}>1 or $clusters{$loc}{final}{commonClips}>1){              
             $clusters{$loc}{final}{$tsp}{status} = 0; 
             $clusters{$loc}{final}{$tsp}{remark} .= "doubleClip;";                         
           }
         }
       }              
       ## <50% reciprocal coverage 
       if($clusters{$loc}{final}{AllClip}>=2){
         my $chk=0;
         if($clusters{$loc}{final}{$tsp}{TEClustPerc} <= 0.5){
           $chk++;
         }elsif($clusters{$loc}{final}{$tsp}{rpos} eq "-" and $clusters{$loc}{final}{$tsp}{lpos} eq "-"){
             $clusters{$loc}{final}{$tsp}{TEClustPerc} = 0; 
             $clusters{$loc}{final}{$tsp}{TEClustReads} = 0;
             $chk++;
         }else{
           my $x = $clusters{$loc}{final}{$tsp}{OriginalTE};
           if(exists $clusters{$loc}{pos1}{$x} and $clusters{$loc}{pos1}{$x}{TEClustPerc}> 0.5){
           }elsif(exists $clusters{$loc}{pos2} and exists $clusters{$loc}{pos2}{$x} and $clusters{$loc}{pos2}{$x}{TEClustPerc}> 0.5){            
           }else{
             $clusters{$loc}{final}{$tsp}{TEClustPerc} = 0; 
             $clusters{$loc}{final}{$tsp}{TEClustReads} = 0;
             $chk++;
           }
         }
         if($chk>0){
           $clusters{$loc}{final}{$tsp}{status} = 0; 
           $clusters{$loc}{final}{$tsp}{remark} .= "RecClust-;";
         }
       }
       ## Fuzzy/ Ambiguous breakpoints
       if($clusters{$loc}{final}{AllClip}>=2){
          $clusters{$loc}{final}{$tsp}{FbyC} = "0.0";
          if(defined $clusters{$loc}{final}{$tsp}{AllFuzzy}){
            $clusters{$loc}{final}{$tsp}{AllFuzzy} -= $clusters{$loc}{final}{AllClip};
          }else{
            $clusters{$loc}{final}{$tsp}{AllFuzzy}=0;
          }
          $clusters{$loc}{final}{$tsp}{AllFuzzy} += $clusters{$loc}{final}{$tsp}{umap} if(defined $clusters{$loc}{final}{$tsp}{umap});
          my $a = 0;
          my $b = $clusters{$loc}{final}{AllClip};
          $a = $clusters{$loc}{final}{$tsp}{AllFuzzy};
          $clusters{$loc}{final}{$tsp}{FbyC} = src::Utilities::round($a/($b+$a),2);
          
          if(defined $clusters{$loc}{pos1}{fuzzyclipped}){
            my $a = $clusters{$loc}{pos1}{fuzzyclipped};
            my $b = ($clusters{$loc}{pos1}{fuzzyclipped} + $clusters{$loc}{pos1}{AllClip});
            my $c = ($clusters{$loc}{pos1}{AllClip} - $clusters{$loc}{final}{$tsp}{lpoly} - $clusters{$loc}{final}{$tsp}{lte});
            $a += $c if($c>0);
            $clusters{$loc}{pos1}{FbyC} = src::Utilities::round( $a/$b,2);
          }  
          if(exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{fuzzyclipped}){
            my $a = $clusters{$loc}{pos2}{fuzzyclipped};
            my $b = ($clusters{$loc}{pos2}{fuzzyclipped}+$clusters{$loc}{pos2}{AllClip});
            my $c = ($clusters{$loc}{pos2}{AllClip} - $clusters{$loc}{final}{$tsp}{rpoly} - $clusters{$loc}{final}{$tsp}{rte});
            $a += $c if($c>0);
            $clusters{$loc}{pos2}{FbyC} = src::Utilities::round( $a/$b,2);
          } 
          #if there are >1/2rd of total fuzzy clipped reads in the vicinity, flag a cluster as Category 0           
          my $chk = 0;
          if($clusters{$loc}{final}{$tsp}{FbyC} >= 0.5){
            $chk++;
          }elsif(defined $clusters{$loc}{pos1}{FbyC} and $clusters{$loc}{pos1}{FbyC}>=0.8 and $clusters{$loc}{pos1}{fuzzyclipped}>=5){
            $chk++;
          }elsif(exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{FbyC} and$clusters{$loc}{pos2}{FbyC}>=0.8 and $clusters{$loc}{pos2}{fuzzyclipped}>=5){
            $chk++;
          }
          if($chk >0){
           $clusters{$loc}{final}{$tsp}{status} = 0 ;  
           $clusters{$loc}{final}{$tsp}{remark} .= "Amb;";                     
          }
       }      
       ## STR?
       if($clusters{$loc}{final}{AllClip}>=2){
         my $chk=0;
         if($clusters{$loc}{final}{STR} eq "T"){
           #$chk++;
         }elsif($clusters{$loc}{final}{tailinfo} =~ /PolyAPolyA/ or $clusters{$loc}{final}{tailinfo} =~ /PolyAPolyT/ or $clusters{$loc}{final}{tailinfo} =~ /PolyTPolyA/ or $clusters{$loc}{final}{tailinfo} =~ /PolyTPolyT/){
           #$chk++;
         }elsif($clusters{$loc}{final}{$tsp}{rpoly} >0 and $clusters{$loc}{final}{$tsp}{lpoly} >0 ){
           $chk++;
         }elsif(defined $clusters{$loc}{pos1}{PolyA} and exists $clusters{$loc}{pos2} and defined $clusters{$loc}{pos2}{PolyA}){
           $chk++;     
         }
         if($chk>0){
           #$clusters{$loc}{final}{$feature}{status} = 0;
           #$clusters{$loc}{final}{$feature}{remark} = "STR?;";
         } 
       }
       ## Indel?
       if($clusters{$loc}{final}{$tsp}{insert_size}!~ /NA/ and $clusters{$loc}{final}{$tsp}{insert_size} < $INDEL_CUTOFF){
          $clusters{$loc}{final}{$tsp}{status} = 0;  
          $clusters{$loc}{final}{$tsp}{remark} .= "indel?;";     
       }
       ## Clip-reads are clonal pileup?
       if($clusters{$loc}{final}{AllClip}>=2){
         $clusters{$loc}{pos2}{identClip} = "NA" if(exists $clusters{$loc}{pos2} and $clusters{$loc}{pos2}{AllClip}==1);
         $clusters{$loc}{pos1}{identClip} = "NA" if($clusters{$loc}{pos1}{AllClip}==1);
         if(exists $clusters{$loc}{pos2}){
           if($clusters{$loc}{pos1}{identClip} eq "NA" and $clusters{$loc}{pos2}{identClip} eq "NA"){
             $clusters{$loc}{final}{identClip} = "NA" 
           }else{
             $clusters{$loc}{final}{identClip} = src::Utilities::max(($clusters{$loc}{pos1}{identClip},$clusters{$loc}{pos2}{identClip}));       
           }        
         }elsif($clusters{$loc}{pos1}{identClip} eq "NA"){
           $clusters{$loc}{final}{identClip} = "NA"
         }elsif($clusters{$loc}{pos1}{identClip} ne "NA"){
           $clusters{$loc}{final}{identClip} = $clusters{$loc}{pos1}{identClip};
         }
         my $chk=0;  
         if($clusters{$loc}{final}{identClip} ne "NA" and $clusters{$loc}{final}{identClip}>=0.8){
           $chk++ if($clusters{$loc}{final}{start}==$clusters{$loc}{final}{end});
           if($clusters{$loc}{final}{$tsp}{lte}>0 and $clusters{$loc}{final}{$tsp}{rte}==0 ){
             $chk++ if($clusters{$loc}{final}{$tsp}{rpoly}==0 and exists $clusters{$loc}{pos2} and $clusters{$loc}{pos2}{isPolyA} eq 'F');
           }elsif($clusters{$loc}{final}{$tsp}{lte}==0 and $clusters{$loc}{final}{$tsp}{rte}>0 ){
             $chk++ if($clusters{$loc}{final}{$tsp}{lpoly}==0 and $clusters{$loc}{pos1}{isPolyA} eq 'F');
           }
         }
         if($clusters{$loc}{final}{$tsp}{lcloneclip}>=0.8 and $clusters{$loc}{final}{$tsp}{rcloneclip}>=0.8){
           if( ($clusters{$loc}{final}{$tsp}{lte}+$clusters{$loc}{final}{$tsp}{lpoly})>=2){
             if( ($clusters{$loc}{final}{$tsp}{rte}+$clusters{$loc}{final}{$tsp}{rpoly})>=2){
               $chk++;
             }
           }
         }elsif($clusters{$loc}{final}{$tsp}{lcloneclip}>=0.8 and $clusters{$loc}{final}{$tsp}{rcloneclip}==0){
              $chk++ if( ($clusters{$loc}{final}{$tsp}{lte}+$clusters{$loc}{final}{$tsp}{lpoly})>=2); 
         }elsif($clusters{$loc}{final}{$tsp}{lcloneclip}==0 and $clusters{$loc}{final}{$tsp}{rcloneclip}>=0.8){
              $chk++ if( ($clusters{$loc}{final}{$tsp}{rte}+$clusters{$loc}{final}{$tsp}{rpoly})>=2);  
         }
         if($chk>0){
           $clusters{$loc}{final}{$tsp}{status} = 0; 
           $clusters{$loc}{final}{$tsp}{remark} .= "Clonal?;";                   
         }
       }
       ## remove false copies of L1/SVA (if no polyA tail but TSD present)
       if($clusters{$loc}{final}{start}!=$clusters{$loc}{final}{end}){
         my $chk=0;
         if($clusters{$loc}{final}{$tsp}{rpos} ne "-" and $clusters{$loc}{final}{$tsp}{lpos} ne "-" and $clusters{$loc}{final}{isPolyA} eq "F"){
           my $yx= src::Utilities::max(($clusters{$loc}{final}{$tsp}{rpos},$clusters{$loc}{final}{$tsp}{lpos}));
           if(defined $yx and abs($yx-$clusters{$loc}{final}{$tsp}{TEsize})>$INDEL_CUTOFF){
             $chk++;
           }elsif(abs($clusters{$loc}{final}{$tsp}{rpos}-$clusters{$loc}{final}{$tsp}{lpos})<10){
             $chk++;
           }      
         }
         if($chk>0){
            $clusters{$loc}{final}{$tsp}{status} = 0; 
            $clusters{$loc}{final}{$tsp}{remark} .= "badCopy;";                      
         }
       }      
       ## Coverage based filters, reads >5% of the coverage or at least 4 clip and 5 RAM
       if($clusters{$loc}{final}{AllClip}>=2){
         my $XX=($clusters{$loc}{final}{bgcoverage}/20);
         if($clusters{$loc}{final}{AllClip} < $XX ){
            $clusters{$loc}{final}{$tsp}{status} = 0;
            $clusters{$loc}{final}{$tsp}{remark} .= "<5%ClipCov;";   
         }
         if($clusters{$loc}{final}{$tsp}{RAMCount} <2 ){ 
           $clusters{$loc}{final}{$tsp}{status} = 0;
           $clusters{$loc}{final}{$tsp}{remark} .= "RAM=0;";         
         }elsif($clusters{$loc}{final}{$tsp}{RAMCount} < $XX/2 ){
           $clusters{$loc}{final}{$tsp}{status} = 0;
           $clusters{$loc}{final}{$tsp}{remark} .= "<2.5%RAMCov;";         
         }
       }
       ## Clip at End of TE-consensys with polyA/T reads
       if(defined $clusters{$loc}{final}{$tsp}{TEsize} and $clusters{$loc}{final}{isPolyA} eq "T"){
         my $x = $clusters{$loc}{final}{$tsp}{TEsize}*0.99;
         my $chk=0;
         my $chy=0;
         if($clusters{$loc}{final}{$tsp}{lte}>0 and $clusters{$loc}{final}{$tsp}{rte}==0 and $clusters{$loc}{final}{$tsp}{lpos} ne "-"){
           if($clusters{$loc}{final}{$tsp}{lpos} >= $x or abs($clusters{$loc}{final}{$tsp}{lpos}-$clusters{$loc}{final}{$tsp}{TEsize})<=$INDEL_CUTOFF ){
              $chk++;
           }
         }
         if($clusters{$loc}{final}{$tsp}{rte}>0 and $clusters{$loc}{final}{$tsp}{lte}==0 and $clusters{$loc}{final}{$tsp}{rpos} ne "-"){
           if($clusters{$loc}{final}{$tsp}{rpos} >= $x or abs($clusters{$loc}{final}{$tsp}{rpos}-$clusters{$loc}{final}{$tsp}{TEsize})<=$INDEL_CUTOFF ){
             $chk++;
           }
         }
         if($chk>0){
           $clusters{$loc}{final}{$tsp}{status} = 0;
           $clusters{$loc}{final}{$tsp}{remark} .= "ClipAtEnd;";                   
         }
       }
       ## PolyA/diminished-PolyA tail at the locus?
       if($clusters{$loc}{final}{isPolyA} eq "F" and ($tsp eq "Alu" or $tsp eq "L1" or $tsp eq "SVA") ){
          ## last attempt to salvage the insertion calls
          ## check following
          ## 1. Scan the region of breakpoint for PolyA reads.
          ## 2. If unavailable, check whether the clip-reads at the locus carry PolyA/T features
          ##   PolyA/T features include
          ##     a. at least 2/3rd As or Ts
          ##     b. at least one homo-7 polymer of A/T
          ##     c. read with at least 7/10  As/Ts at the clip junction
          my $chk=0; 
          if($clusters{$loc}{final}{$tsp}{remark} eq "-"){
            my $ffile = $outprefix.".RAM.sort.bam";
            my $X = $clusters{$loc}{final}{$tsp}{OriginalTE};
            my $CLP = "";
            my $sidedness ="";
            my $floc = "";
            if(defined $clusters{$loc}{pos1}{$X}){
               $sidedness = $clusters{$loc}{pos1}{clip_sidedness};
               $CLP  = $clusters{$loc}{pos1}{cliploc};
               $floc = $clusters{$loc}{pos1}{chr}.":".($CLP-30)."-".($CLP+40) if($sidedness eq "a");
               $floc = $clusters{$loc}{pos1}{chr}.":".($CLP-40)."-".($CLP+30) if($sidedness eq "b");
            }
            if($floc ne ""){
              if($sidedness eq "a"){
                $sidedness = "b";
              }else{ $sidedness = "a";}
              my $result = src::Utilities::gridCheck_PolyA($CLP,$sidedness,$floc,$ffile,-30,40);
              if($result ne "n/a"){
                my ($coord) = $result =~ m/(\d*)\|/;
                $clusters{$loc}{final}{$feature}{remark} = $CLP-$coord;            
                print "   found potential degenerate-tail for $X within $floc .. precisely at $result\n";
                $chk++;
              }
            }
          }
          if($chk==0){ # and
            $clusters{$loc}{final}{$feature}{status} = 0;
            $clusters{$loc}{final}{$feature}{remark} = "noTail;";        
          }  
       }    

      }
    }

    ### provide sort order for reporting
    foreach my $feature (@features){
     if(exists $clusters{$loc}{final}{$feature} and $clusters{$loc}{final}{ReportedTE} eq $feature){
       $cntcalls{$clusters{$loc}{final}{$feature}{status}}++;
       if($clusters{$loc}{final}{$feature}{status}>$clusters{$loc}{final}{sort}){
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
   print FO "## Status 3: High confidence insertion with clip reads on both ends (total= $cntcalls{3}) \n";
   print FO "## Status 2: Insertions with either right or left hand side clip support (total= $cntcalls{2}) \n"; 
   print FO "## Status 1: Insertions overlapping with known genomic copies of transposable element (total= $cntcalls{1})\n";
   print FO "## Status 0: Poor quality/ low confidence insertions (total= $cntcalls{0}, should be ommitted) \n";
   print FO "########################################################################################\n";
   print FO "#chr\tstart\tend\tid\tscore\tstrand\tTE\tstatus\tdescription\tremark\n";
   
   foreach my $loc (@keys){
     foreach my $feature (@features){
      if(exists $clusters{$loc}{final}{$feature} and $feature eq $clusters{$loc}{final}{ReportedTE}){
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
        my $a=0;
        $a = $clusters{$loc}{final}{$feature}{AllFuzzy} if(defined $clusters{$loc}{final}{$feature}{AllFuzzy});
        print FO "UN=$clusters{$loc}{final}{$feature}{umap};FUZZY=$a;TOTALCLIP=$clusters{$loc}{final}{AllClip};";
        print FO "PURITY=$clusters{$loc}{final}{$feature}{Purity};";
        print FO "RAM=$clusters{$loc}{final}{$feature}{RAMCount};Pval=$clusters{$loc}{final}{$feature}{pvalue};";
        print FO "COVERAGE=$clusters{$loc}{final}{bgcoverage};";
        print FO "gTE=$clusters{$loc}{final}{$feature}{bgTE};D2gTE=$clusters{$loc}{final}{$feature}{bgTEDist};";
        #print FO "TECLUSTER1=$clusters{$loc}{final}{$feature}{TEClust};";
        if($clusters{$loc}{final}{identClip} ne "NA"){
          $clusters{$loc}{final}{identClip} = src::Utilities::round($clusters{$loc}{final}{identClip},2);
        }
        $clusters{$loc}{final}{OriRatio} = src::Utilities::round($clusters{$loc}{final}{OriRatio},2);
        print FO "MAXMAP=$clusters{$loc}{final}{$feature}{TESc};";
        print FO "D2RE=$clusters{$loc}{final}{dist2remotif};";
        print FO "%FUZZY=$clusters{$loc}{final}{$feature}{FbyC};";
        print FO "N_CLUSTERS=$clusters{$loc}{pos1}{bgclusters};RECI%=$clusters{$loc}{final}{$feature}{TEClustPerc};";
        print FO "MAXCLONECLIP=$clusters{$loc}{final}{identClip};ORIRATIO=$clusters{$loc}{final}{OriRatio};";
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
     $polyr = src::Utilities::round($polyclust*100/($clustnum+$filtered),2);
   }
   print "  Total $polyclust putative polyA/T expansions out of total ",($clustnum+$filtered)," insertion candidates ($polyr%)\n";
   $polyr = 0;
   if($clustnum >0){
     $polyr = src::Utilities::round($unassigned*100/($clustnum+$filtered),2);
   }
   print "  Total unassigned $unassigned clusters out of total ",($clustnum+$filtered)," insertion candidates ($polyr%)\n";
   print "  Total filtered clusters by HiTEA (e.g. status 0): $filtered\n";
   print "  Clusters written in the reporting bed file: $clustnum\n";
   print "  Number of HiTEA insertions:",($cntcalls{3}+$cntcalls{2}+$cntcalls{1}),"\n";
   
   return(\%clusters);
}

