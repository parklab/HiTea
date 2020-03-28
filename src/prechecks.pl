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
use open qw(:std :utf8);
require "src/vars.pl";
our (%redb);
our (%chrs);
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

my $index ="";
my $repbase ="";
my $bgAnnotations ="";
my $polym ="";
my $enzyme ="";
my $help=0;
Getopt::Long::GetOptions(
  'index=s'          => \$index,
  'polym=s'          => \$polym,
  'repbase:s'        => \$repbase,
  'bgAnnotations:s'  => \$bgAnnotations,
  'enzyme=s'         => \$enzyme,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
sub help{
  my $j = shift;
  if($j){
   print "\nUsage: perl prechecks.pl -polym [PATH] -index [PATH] -repbase [PATH] -bgAnnotations [DIRPATH] -enzyme [STRING]\n\n";
   print "This program annotates finalized breaks using the bam file\n\n";
   print "Options:\n\n";
   print "  -index                 TE consensus index (generated using BWA)\n";
   print "  -repbase               RepBase index (generated using BWA)\n";
   print "  -bgAnnotations         Reference TE annotations\n";
   print "  -polym                 TE polymorphic example cases index (generated using BWA)\n";
   print "  -enzyme                RE enzyme \n";
   print "  -help|-h               Display usage information.\n\n\n";
   exit 1;
  }
}

print "    (perl prechecks.pl -index $index -polym $polym -repbase $repbase -bgAnnotations bgAnnotations -enzyme $enzyme)\n";

if($help or $repbase eq "" or $index eq "" or $bgAnnotations eq "" or $enzyme eq ""){
  print " One or more inputs are not recognized***\n";
  help(1);
}

##--------------------------------------------------- I/O --------------------------------
my %polymorphs;
my %transposons;
my %files;
my @files;

if(!defined $redb{$enzyme}){
  print "       [prechck] Enzyme $enzyme is not supported at the moment\n";
  exit 1;
}

%transposons = src::Utilities::get_fasta($index);

@files = glob($bgAnnotations."/*.bed");
%files = map { $_ => 1 } @files;
foreach my $i (keys %transposons){
 next if($i eq "" or $i eq "PolyA" or $i eq "PolyT");	
 if(!$files{$bgAnnotations."/".$i.".bed"}){
	print "       Background Annotation file for $i not found at location: $bgAnnotations\n";
	exit 1;
 }
}


if($polym ne ""){
  %polymorphs = src::Utilities::get_fasta($polym);
  foreach my $i (keys %polymorphs){
   my @temp = split("~",$i); ## Family~Subfamily
   print "       $i\t$temp[1]\t$polymorphs{$i}","bp\n"; 
   $temp[1] =~ s/\n//;
   $temp[1] =~ s/\r//;
   if(!defined $transposons{$temp[0]}){
     print "       Error: The header of subfamily sequences should be in the format - Family~Subfaily\n";
     print "              The family $temp[0] could not be found\n";
     exit 1;
   }
  }
}




