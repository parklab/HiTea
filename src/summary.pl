#!/usr/bin/perl
# HiTEA
use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use src::Utilities;
use open qw(:std :utf8);
BEGIN { our $start_run = time(); }

###-------------------------------------------------------------------------------------------------------------------------------
# inputs etc
###-------------------------------------------------------------------------------------------------------------------------------
my $outprefix = "";
my $wd = "";
my $help = 0;
Getopt::Long::GetOptions(
  'wd:s'             => \$wd,
  'outprefix:s'      => \$outprefix,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl summary.pl -outprefix [STRING] -wd [WORK_DIR]\n\n";
  print "This script reads the perl objects generated at workdirectory location and writes summary of input reads. \n";
  print "Options:\n";   
  print "   -outprefix        Output file PREFIX\n";
  print "   -wd               Working directory\n";
  print "   -help|-h          Display usage information.\n";
  print "Default outputs:\n";
  print "    Writes summary report in text format \n\n\n";
  exit 0;
}
$wd =~ s/\/$//;
if($wd ne ""){
  $outprefix = $wd."/".$outprefix;
}

my @files = glob($wd."/*.summary.log.ph");
if(scalar(@files) ==0){
	print " Summary is already displayed in a log file\n";
	exit 0;
}

my $outlog=$outprefix.".ReadSummary.logs";
open(LOGS, "> $outlog") or die "can't create $outlog";
print LOGS "##-----------------------------------------------------------------------------------------------------\n";

my %flags;
my %cov;
foreach my $f(@files){
   my $hashref = retrieve($f);
   my %r = %{$hashref};
   
   $f=~ s/.summary.log.ph//;
   #print LOGS "files: $f\n";   
   foreach my $dis (sort keys %{$r{"ori"}}){
      foreach my $ori (sort keys %{$r{"ori"}{$dis}}){
        $r{"ori"}{$dis}{$ori} = 0 if(!$r{"ori"}{$dis}{$ori});
        $flags{"ori"}{$dis}{$ori} += $r{"ori"}{$dis}{$ori};
      }
   }
   foreach my $cl (sort keys %{$r{"raw"}{"counts"}}){
     $r{"raw"}{"counts"}{$cl} = 0 if(!$r{"raw"}{"counts"}{$cl});
     $flags{"raw"}{"counts"}{$cl} += $r{"raw"}{"counts"}{$cl};	
   }
   foreach my $cl (sort keys %{$r{"raw"}{"class"}}){
     $flags{"raw"}{"class"}{$cl} += $r{"raw"}{"class"}{$cl};
     $r{"reported"}{"class"}{$cl} = 0 if(!$r{"reported"}{"class"}{$cl});	
     $flags{"reported"}{"class"}{$cl} += $r{"reported"}{"class"}{$cl};	
   }
   foreach my $cl (sort keys %{$r{"raw"}{"read1"}}){
     $flags{"raw"}{"read1"}{$cl} += $r{"raw"}{"read1"}{$cl};	
     $r{"reported"}{"read1"}{$cl} = 0 if(!$r{"reported"}{"read1"}{$cl});	
     $flags{"reported"}{"read1"}{$cl} += $r{"reported"}{"read1"}{$cl};	
   }
   foreach my $cl (sort keys %{$r{"raw"}{"read2"}}){
     $flags{"raw"}{"read2"}{$cl} += $r{"raw"}{"read2"}{$cl};	
     $r{"reported"}{"read2"}{$cl} = 0 if(!$r{"reported"}{"read2"}{$cl});	
     $flags{"reported"}{"read2"}{$cl} += $r{"reported"}{"read2"}{$cl};	
   }
   foreach my $chr(sort keys %{$r{"cov"}}){
     foreach my $pos(keys %{$r{"cov"}{$chr}}){
      $cov{$chr}{$pos} = 0 if(!$cov{$chr}{$pos});
      $cov{$chr}{$pos} += $r{"cov"}{$chr}{$pos};
     }
   }
}
## store the coverage hash
my $covObj=$outprefix.".coverage.ph";
store \%cov, $covObj;


$flags{"raw"}{"counts"}{"linear_rescued"} = 0 if(!$flags{"raw"}{"counts"}{"linear_rescued"});
$flags{"raw"}{"counts"}{"uninformative"} = 0 if(!$flags{"raw"}{"counts"}{"uninformative"});
$flags{"raw"}{"counts"}{"lowqual"} = 0 if(!$flags{"raw"}{"counts"}{"lowqual"});
$flags{"raw"}{"counts"}{"nn"} = 0 if(!$flags{"raw"}{"counts"}{"nn"});
$flags{"raw"}{"counts"}{"outofindex"} = 0 if(!$flags{"raw"}{"counts"}{"outofindex"});
$flags{"raw"}{"counts"}{"lonelymates"} = 0 if(!$flags{"raw"}{"counts"}{"lonelymates"});
$flags{"raw"}{"counts"}{"nonHicChimera"} = 0 if(!$flags{"raw"}{"counts"}{"nonHicChimera"});
$flags{"raw"}{"counts"}{"twojunctoins"} = 0 if(!$flags{"raw"}{"counts"}{"twojunctoins"});
$flags{"raw"}{"counts"}{"wgs_nonchimera"} = 0 if(!$flags{"raw"}{"counts"}{"wgs_nonchimera"});
$flags{"raw"}{"counts"}{"cis"} = 0 if(!$flags{"raw"}{"counts"}{"cis"});
$flags{"raw"}{"counts"}{"trans"} = 0 if(!$flags{"raw"}{"counts"}{"trans"});


print LOGS "Total read-pairs in the file\t".$flags{"raw"}{"counts"}{"total"}."\n";
my $perc = src::Utilities::round($flags{"raw"}{"counts"}{"cis"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Total \"cis\" read-pairs in the file\t".$flags{"raw"}{"counts"}{"cis"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"trans"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Total \"trans\" read-pairs in the file\t".$flags{"raw"}{"counts"}{"trans"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"linear_rescued"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Linearly Mapped and rescued read-pairs (UU/UR/RU)\t".$flags{"raw"}{"counts"}{"linear_rescued"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"uninformative"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Uninformative read-pairs\t".$flags{"raw"}{"counts"}{"uninformative"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"lowqual"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Low mapping quality read-pairs (MAPQ < 28)\t".$flags{"raw"}{"counts"}{"lowqual"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"nn"}*100/$flags{"raw"}{"counts"}{"total"},4);
print LOGS "Read-pairs with at least 2 Ns in both the read-sequences\t".$flags{"raw"}{"counts"}{"nn"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"outofindex"}*100/$flags{"raw"}{"counts"}{"total"},3);
print LOGS "Read-pairs mapping to non-selected chromosomes\t".$flags{"raw"}{"counts"}{"outofindex"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"lonelymates"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Reads with only single mate mapped\t".$flags{"raw"}{"counts"}{"lonelymates"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"wgs_nonchimera"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "WGS non-chimeric reads\t".$flags{"raw"}{"counts"}{"wgs_nonchimera"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"nonHicChimera"}*100/$flags{"raw"}{"counts"}{"total"},2);
print LOGS "Non-HiC chimeric read pairs\t".$flags{"raw"}{"counts"}{"nonHicChimera"}." (".$perc."%)\n";
$perc = src::Utilities::round($flags{"raw"}{"counts"}{"twojunctoins"}*100/$flags{"raw"}{"counts"}{"total"},3);
print LOGS "Pairs with 2 RE junctions\t",$flags{"raw"}{"counts"}{"twojunctoins"}." (".$perc."%)\n";


print LOGS "\n# Read/Flag-Class\t(Inp.)\t\t(Report)\t\n";
my $ptotal = 0;
foreach my $cl (sort keys %{$flags{"reported"}{"class"}}){
    my $percent = src::Utilities::round($flags{"reported"}{"class"}{$cl}*100/$flags{"raw"}{"class"}{$cl},2);
    $percent = "(".$percent."%)";
    print LOGS "# ",$cl,"\t",$flags{"raw"}{"class"}{$cl},"\t\t",$flags{"reported"}{"class"}{$cl}," $percent\t\n";
    $ptotal += $flags{"reported"}{"class"}{$cl};
}
print LOGS "\n# Flag-Class\t(R1,Inp.)\t(R2,Inp.)\t(R1,Report)\t(R2,Report)\n";
my $pt1 = 0;
my $pt2 = 0;
foreach my $cl (sort keys %{$flags{"raw"}{"read1"}}){
    my $p1;
    my $p2;
    $flags{"reported"}{"read1"}{$cl}=0 if(!$flags{"reported"}{"read1"}{$cl});
    $flags{"reported"}{"read2"}{$cl}=0 if(!$flags{"reported"}{"read2"}{$cl});
    $p1 = "(".src::Utilities::round($flags{"reported"}{"read1"}{$cl}*100/$flags{"raw"}{"read1"}{$cl},2)."%)";
    $p2 = "(".src::Utilities::round($flags{"reported"}{"read2"}{$cl}*100/$flags{"raw"}{"read2"}{$cl},2)."%)";
    print LOGS "# ",$cl,"\t",$flags{"raw"}{"read1"}{$cl},"\t",$flags{"raw"}{"read2"}{$cl},"\t";
    print LOGS $flags{"reported"}{"read1"}{$cl}," $p1\t",$flags{"reported"}{"read2"}{$cl}," $p2\n";
    $pt1 += $flags{"reported"}{"read1"}{$cl}+$flags{"reported"}{"read2"}{$cl};
    $pt2 += $flags{"raw"}{"read1"}{$cl}+$flags{"raw"}{"read2"}{$cl};
}
my $percent = src::Utilities::round($ptotal*100/$flags{"raw"}{"counts"}{"total"},2);
$percent = "(".$percent."%)";
print LOGS "Total pairs used in the analyses:\t $ptotal $percent\n";
$percent = "(".src::Utilities::round($pt1*100/$pt2,2)."%)";
print LOGS "Total %flags retained for reporting:\t$pt1\t$percent (out of $pt2)\n";
print LOGS "##-----------------------------------------------------------------------------------------------------\n";
close(LOGS);


my $outlog1=$outprefix.".ReadPairInsertSizeOriSummary.logs.gz";
open(LOGS, "| gzip -c - > $outlog1") or die "can't create $outlog1";
foreach my $dis (sort keys %{$flags{"ori"}}) {
  foreach my $ori (sort keys %{$flags{"ori"}{$dis}}){
    print LOGS $dis."\t".$ori."\t".$flags{"ori"}{$dis}{$ori}."\n";
  }
}
close(LOGS);

exit 0;
