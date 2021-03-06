#!/usr/bin/perl
# HiTEA
# This program finalizes breakpoint locations from sorted cluster.bed.gz file

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use constant { true => 1, false => 0 };
use utf8;
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0);
use Storable;
use open qw(:std :utf8);
use src::Utilities;
BEGIN { our $start_run = time(); }
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

my $in = "";
my $wd = "";
my $outprefix ="project";
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

sub help{
  my $j = shift;
  if($j){
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
   exit 1;
  }
}
$wd =~ s/\/$//;
my $baseoutprefix =$outprefix;
$outprefix = $wd."/".$outprefix if($wd ne "");
if($help or $in eq ""){
  print " Inputs are not correctly specified\n";
  help(1);
}

#--------------------------------------------------------------------------------------------------------------
# I/O
#--------------------------------------------------------------------------------------------------------------
my %clusters;
my %omtcl;
my $watch_run = time();
my $run_time_start = $watch_run - our $start_run;

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
if(exists $omtcl{"Pile up cluster"}){
  print "  Clusters with >80% identical read-pile up at a given locus: ",$omtcl{"Pile up cluster"},"\n";
}
if(exists $omtcl{"0 softclipped reads"}){
  print "  Warning: Clusters with 0 softclipped reads: ",$omtcl{"0 softclipped reads"},"\n";
}
print "[finalize_breaks] END:\t $run_time_end seconds\n";

#my $oufile=$outprefix.'.Obj.clusters'; #Save 
#open FOO,">$oufile" or die $!;
#print FOO Dumper %clusters;
#close(FOO);
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
        #chr1    1657990 1657990 -       chr116579487S42M52Schr1038277337|21|c:1657990 chr116579487S42M52Schr1038277337|21|c:1657990
        $temp[4] =~ s/ $//;
        $temp[4] =~ s/ ,/ /g;
        $temp[4] = src::Utilities::unique($temp[4]," ","false");
        my (@scores) = $temp[4]=~ m/\|(\d*)\|/g;
        @scores = sort{$b<=>$a} @scores;
        my (@clips) = $temp[4]=~ m/\|c:(\d*)/g;
        
        ## remove clusters where identical mate piles up on a location
        my $ty = $temp[4];
        my $posstr = () = $ty=~ /\+/g;
        $posstr = 0 if(!$posstr);
        $ty =~ s/=/chr/g;
        $ty =~ s/\*/chr/g;
        $ty=~ s/,//g;
        $ty=~ s/  / /g;
        $ty =~ s/chr/ /g;
        my @ty = split(" ",$ty);
        @ty = grep { $_ !~ '\|' } @ty;
        @ty = grep { $_ ne ' ' } @ty;
        my %ty; $ty{$_}++ for(@ty);
        my @keys = sort { $ty{$b} <=> $ty{$a} }  keys %ty;
        my $tyout=0;
        $tyout = $ty{$keys[0]}/scalar(@ty) if(scalar @keys>0);
        my $strrat = src::Utilities::max(($posstr / scalar@ty), (scalar@ty - $posstr)/scalar@ty);


        if(scalar $tyout >= 0.8){
          $omtcl{"Pile up cluster"}++;
        }

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
            $clusters{$loc}{pos1}{identClip} = $tyout;            
            $clusters{$loc}{pos1}{OriRatio} = $strrat;
            print FO1 "$loc\t$clusters{$loc}{pos1}{maxscore}\t$clusters{$loc}{pos1}{clipreads}\t",join(",",@scores),"\n";
        }else{
            $omtcl{"Poor Alignment/Single Clip"}++;
        }
    }
    close(IN);
    close(FO1); 
    return(\%clusters);
}


