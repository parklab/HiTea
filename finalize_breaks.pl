#!/usr/bin/perl
# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard

# This program finalizes breakpoint locations from sorted cluster.bed.gz file

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

my $in = "";
my $wd = "";
my $outprefix ="";
my $algnscore = 20;
my $help=0;
Getopt::Long::GetOptions(
  'in=s'             => \$in,
  'wd:s'             => \$wd,
  'outprefix:s'      => \$outprefix,
  'algnscore:s'      => \$algnscore,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl finalize_breaks.pl -in [FILE_PATH] -algnscore [INT] -outprefix [STRING] -wd [DIR_PATH]\n\n";
  print "This program finalizes breakpoints where the clustered breaks are scanned and most represented position is used as a breakpoint\n\n";
  print "Options:\n\n";
  print "***required:\n";
  print "  -in                    Clusters file in bed format \n";
  print "***optional:\n";
  print "  -algnscore             Minimum alignment score [default: 20]\n";
  print "  -outprefix             Outputprefix for generating 2 output files [default: project] \n";
  print "  -wd                    Working directory [default: ~]\n";
  print "  -help|-h               Display usage information.\n\n\n";
   exit 0;
}
$wd =~ s/\/$//;
my $baseoutprefix =$outprefix;
if($wd ne ""){
  $outprefix = $wd."/".$outprefix;
}
my $watch_run = time();
my $run_time_start = $watch_run - our $start_run;

my %clusters;
my %omtcl;
my $clust = generate_clusters_object($in,\%clusters);
%clusters =%{$clust};
store \%clusters, $outprefix.'.ClustObj.ph'; #Save
$watch_run = time();
my $run_time_end = $watch_run - $start_run;
if(scalar(keys%clusters)<1){
  print " Aborting as there are no clusters identified in the input file \n";
  exit 1;
}

print "[finalize_breaks] START:\t $run_time_start seconds\n";
print "  Command: perl finalize_breaks.pl -in $in -algnscore $algnscore -outprefix $baseoutprefix -wd $wd\n";
print "  Clusters ommitted due to poor alignment score/single clip read: ",$omtcl{"Poor Alignment/Single Clip"},"\n";
print "  Total cluster locations: ", scalar(keys%clusters)," \n";
if(exists $omtcl{"0 softclipped reads"}){
  print "  Warning: Clusters with 0 softclipped reads (check $outprefix.skippedclusters.logs.gz file for details): ",$omtcl{"0 softclipped reads"},"\n";
}
print "[finalize_breaks] END:\t $run_time_end seconds\n";

exit 0;
#--------------------------------------------------------------------------------------------------------------
# subroutines
#--------------------------------------------------------------------------------------------------------------
sub generate_clusters_object{
    my ($file,$clust) = @_;
    my %clusters = %{$clust};
    if($file =~ /.gz/){
       open IN,"zcat $file|" or next "Can't open file $file";  ## sorted clusters file, sort -k1,1 -k2,2n -k4,4 -u
    }else{
       open IN,"$file" or next "Can't open file $file";
    }
    
    my $outf=$outprefix."_HighConfidenceBreaks.txt.gz";
    open(FO1, "| gzip -c - > $outf") or die $!;
    print FO1 "#chr\tstart\tside\tmaxmapq\t#reads\tmapqs\n";
         
    while (<IN>) {
        next if(/^(\#)/);
        next if(/^(\@)/); 
        s/\n//; 
        s/\r//;  
        my @temp = split(/\t/);
        #chr1 247754064 247754068 - chr124775403134M61Schr468093925|54|c:247754065 chr124775403134M61Schr468093925|54|c:247754065 ,chr124775403334M61Schr5111307483|46|c:247754067 chr124775403334M61Schr5111307483|46|c:247754067 
        $temp[4] =~ s/ $//;
        $temp[4] =~ s/ ,/ /g;
        $temp[4] = Utilities::unique($temp[4]," ","false");
        my (@scores) = $temp[4]=~ m/\|(\d*)\|/g;
        @scores = sort{$b<=>$a} @scores;
        my (@clips) = $temp[4]=~ m/\|c:(\d*)/g;
        
        my %clcnt;
        foreach(@clips){$clcnt{$_}++;}
        @clips = sort { no warnings;
                        $clcnt{$b} <=> $clcnt{$a}; 
                      } keys %clcnt;

        my $side = "a";
        $side = "b" if($temp[3] eq "-");
        my $j = scalar(split(" ",$temp[4]));
        my $loc = $temp[0]."\t".$clips[0]."\t".$side; 

        if(scalar(@scores)==0){
            $omtcl{"0 softclipped reads"}++;
            print  " Warning: Cluster with 0 softclipped reads is reported: $loc\n";                   
        }elsif(($j >=2 and $scores[0]>=$algnscore)) { 
            $clusters{$loc}{pos1}{maxscore} = $scores[0];
            $clusters{$loc}{pos1}{clipreads} = $j;            
            print FO1 "$loc\t$clusters{$loc}{pos1}{maxscore}\t$clusters{$loc}{pos1}{clipreads}\t",join(",",@scores),"\n";
        }else{
            $omtcl{"Poor Alignment/Single Clip"}++;
        }
    }
    close(IN);
    close(FO1); 
    return(\%clusters);
}


