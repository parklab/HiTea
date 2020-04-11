#!/usr/bin/perl
# HiTEA
# This script parses Hi-C bam file (duplicate-marked lossless bam) and to generate fastq file for mapping on to the TE assembly 

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0);
use src::Utilities;
require "src/vars.pl";
our (%redb);
our (%chrs);
use open qw(:std :utf8);
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
BEGIN { our $start_run = time(); }

###-------------------------------------------------------------------------------------------------------------------------------
# inputs etc
###-------------------------------------------------------------------------------------------------------------------------------
my $bam = "";
my $enzyme = ""; 
my $min_mapq = "";
my $clip = ""; 
my $wd = "";
my $outprefix = "project"; #default 
my $help = 0;
Getopt::Long::GetOptions(
  'bam=s'            => \$bam,
  'e=s'              => \$enzyme,
  'wd:s'             => \$wd,
  'outprefix:s'      => \$outprefix,
  'q=s'              => \$min_mapq,
  'clip=s'           => \$clip,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
sub help{
  my $j = shift;
  if($j){
   print "\nUsage: perl parse.pl -bam [FILE_PATH] -e [STRING] -q [INT] -clip [INT] -outprefix [STRING] -wd [WORK_DIR] \n\n";
   print "This script parses pairsam file to extract discordent reads for assessment of TE-insertions. \n\n";
   print "Options (required*):\n";
   print "   -bam              Input file (either lossless bam format or WGS bam)\n";
   print "   -e                RE used in the Hi-C experiment\n";
   print "   -q                Minimum mapping quality for the reference alignment. This value is used to determine repeat anchored mates in the genome (default:28) \n";
   print "   -clip             Minimum softclipped read length for mapping the reads to the TE assembly\n";
   print "                       (should be at least the length of the ligation motif. Too short lengths may lead to increased false positives)\n";
   print "                       (ideal= 20bp, suggested>=10bp)\n";
   print "Options (optional):\n";
   print "   -outprefix        Output file PREFIX (default: project)\n";
   print "   -wd               Working directory (default: ~)\n";
   print "   -help|-h          Display usage information.\n";
   print "Default outputs:\n";
   print "    Writes fastq file with discordent reads \n\n\n";
   exit 1;
  }
}
$wd =~ s/\/$//;
$outprefix = $wd."/".$outprefix if($wd ne "");
if($help or $min_mapq eq "" or $enzyme eq "" or $bam eq "" or $clip eq ""){
  print "One or more required inputs are not recognized***\n\n"; 
  help(1);
}
if($clip < $redb{$enzyme}{ligmotifln}){
  print " Clip length ($clip) can not be smaller than the ligation motif length ($redb{$enzyme}{ligmotifln})\n";
  help(1);
}

###-------------------------------------------------------------------------------------------------------------------------------
# I/O
###-------------------------------------------------------------------------------------------------------------------------------
my $watch_run=0;
my $run_time=0;
my $run_time_start=0;
my %flags;
my $lines =0;
my $outfq=$outprefix.".temp.fq.gz";  ## softclipped read sequences
open(O1, "| gzip -c - > $outfq") or die "can't create $outfq";
my $outfq2=$outprefix.".temp2.fq.gz";   ## non-clip reads
open(O2, "| gzip -c - > $outfq2") or die "can't create $outfq2";

$watch_run = time();
$run_time_start = $watch_run - our $start_run;
print " Input bam file: $bam\t  ";
$lines = bam_read($bam);
my $sumObj=$outprefix.".summary.log.ph"; ## store summary info
store \%flags, $sumObj;

$watch_run = time();
$run_time = $watch_run - $start_run;
print " $run_time_start .. $run_time seconds\t lines in the file: $lines\n";

close(O1);
close(O2);
exit 0;

###-------------------------------------------------------------------------------------------------------------------------------
# Subroutines
###-------------------------------------------------------------------------------------------------------------------------------
  ## What are the Hi-C non-conforming reads?
  # 1) Split-reads: clipped mate without ligation motif 
  # 2) Unsplit-reads: One of the mates with unique mapping(mapq>=$min_mapq), while other unmapped/multimapped
  # 3) Flags on the read:
  #   DE : >=$clip clipped bases
  #   IE : Indirect support towards insertion/no clipping of the read
  #   FP : Carrying ligation motif at clip 
  #   TP : No ligation motif at clip
  
sub bam_read{
  my ($file) = shift;
  if($file=~ m/.bam/){
     open IN,"samtools view -@ 8 $file|" or next "Can't open file $file";
  }else{
     print "input is pipe/sam\t";
     open IN,"$file" or next "Can't open file $file"; 
  }
  
  my (@read1,@read2);
  my $oid="";

  my $line = 0;
  while(<IN>) {
    next if(/^(\#)/); 
    next if(/^(\@)/); 
    chomp;
    s/\r//;  
    $line++;
    my ($flag) =  $_ =~ /^.*?\t(\d*)/;
    next if($flag>=256);
    my @sam = split(/\t/);
    if(scalar @sam <11){
      print " ERROR reading the input sam/bam file. Exiting!! \n";
      exit 1;
    }

    if($oid eq ""){ $oid = $sam[0];   }

    if($oid eq $sam[0]){
       if($sam[1] & 0x40){
         @read1 = @sam;
       }elsif($sam[1] & 0x80){
         @read2 = @sam;
       }
    }
    
    if($oid ne $sam[0] or eof){
       my $is_next=0;
       
       if($is_next==0){
         $flags{"raw"}{"counts"}{"total"}++;
         ### Filter out the reads
         if(!defined $read1[1] or !defined $read2[1]){
           $flags{"raw"}{"counts"}{"lonelymates"}++;
           $is_next=1;
         }elsif($read1[2] eq "*" and $read2[2] eq "*") {
           $flags{"raw"}{"counts"}{"uninformative"}++; 
           $is_next=1;
         }elsif(length($read1[9])<40 or length($read2[9])<40) {
           $flags{"raw"}{"counts"}{"uninformative"}++; 
           $is_next=1;
         }
         if($is_next==0){
           if($read1[2] eq $read2[2]){
              $flags{"raw"}{"counts"}{"cis"}++;
           }else{
              $flags{"raw"}{"counts"}{"trans"}++;
           }
           write_fq_from_bam(\@read1,\@read2);
         }
       }         
       undef (@read1);
       undef (@read2);
       $oid = $sam[0];
       if($sam[1] & 0x40){
         @read1 = @sam;
       }elsif($sam[1] & 0x80){
         @read2 = @sam;
       }
    }
    #if($flags{"raw"}{"counts"}{"total"}==10000 and $flags{"raw"}{"counts"}{"lonelymates"}>5000){
    #  print " Input bamfile is not sorted by names. Exiting!\n";
    #  exit 1;
    #}
  }
  close(IN);
  return($line);
}

sub write_fq_from_bam{
  my ($a,$b) = @_;
  my @read1 = @{$a};
  my @read2 = @{$b};
  
  ## get flag status
  my $r1_type = "WGS";
  if($read1[-1] =~ /Yt:Z:/){ 
     ($r1_type = $read1[-1])=~ s/Yt:Z://g;
  }
  if($r1_type eq "WGS"){
    push(@read1,"Yt:Z:WGS");
    push(@read2,"Yt:Z:WGS");
  }
  $flags{"raw"}{"class"}{$r1_type}++;
  
  ## filters 1
  if($r1_type eq "DD" or $r1_type eq "NN"){
    $flags{"raw"}{"counts"}{"uninformative"}++;
    return(1);    
  }elsif(!exists($chrs{$read1[2]}) and !exists($chrs{$read2[2]})){ ## both mates out of index
    $flags{"raw"}{"counts"}{"outofindex"}++;
    $flags{"raw"}{"read1"}{"IE,FP"}++;
    $flags{"raw"}{"read2"}{"IE,FP"}++;
    return(1);    
  }elsif($read1[1] ==4 and $read2[1] ==4){ ## both mates unmapped
    $flags{"raw"}{"counts"}{"uninformative"}++;
    $flags{"raw"}{"read1"}{"IE,FP"}++;
    $flags{"raw"}{"read2"}{"IE,FP"}++;
    return(1);      
  }

  ## find chimera
  my $is_chimera = "false";
  my $strand1 = "+";
  my $strand2 = "+";
  $strand1="-" if($read1[1] & 0x10);
  $strand2="-" if($read2[1] & 0x10);
  if($strand1 eq $strand2){
     $is_chimera = "true";
  }elsif($read1[2] ne $read2[2]){
     $is_chimera = "true";
  }elsif(( abs($read1[8]) > 500 and abs($read2[8])>500) or ($read1[8]==0 and $read2[8]==0)){
     $is_chimera = "true";
  }elsif($read1[1] & 0x4 or $read2[1] & 0x4){
     $is_chimera="true";
  }
  if($is_chimera eq "false"){
     $flags{"raw"}{"counts"}{"wgs_nonchimera"}++;
  }
  
  ## get coverage of all mapped reads
  $flags{"cov"}{$read1[2]}{ int (($read1[3]/1000)+0.5) }++ if($read1[3]!=0 and $read1[4] >= $min_mapq);
  $flags{"cov"}{$read2[2]}{ int (($read2[3]/1000)+0.5) }++ if($read2[3]!=0 and $read2[4] >= $min_mapq);
  
  ## filters 2
  if($r1_type eq "UR" or $r1_type eq "RU" or $r1_type eq "UU" or ($r1_type eq "WGS" and $is_chimera eq "false")){
    $flags{"raw"}{"counts"}{"linear_rescued"}++; 
    if($read1[2] eq $read2[2] and $read1[8] !=0){
      my $logdist = src::Utilities::round(log(abs($read1[8]))/log(10),1);
      if($read1[3]<$read2[3]){
         $flags{"ori"}{$logdist}{$strand1.$strand2}++;
      }else{
         $flags{"ori"}{$logdist}{$strand2.$strand1}++;
      }
    }
    return(1);
  } 

  ## filter out read-pairs
  if($read1[5]!~ /S/  and $read2[5]!~ /S/){  
    if($read1[4] >= $min_mapq and $read2[4] >= $min_mapq){ ## both mates without clip info and both have mapq > min_mapq
      $flags{"raw"}{"read1"}{"IE,FP"}++;
      $flags{"raw"}{"read2"}{"IE,FP"}++;
      $flags{"raw"}{"counts"}{"uninformative"}++;
      return(1);
    }elsif($read1[4] ==0 and $read2[4] == 0){  ## both mates without clip info and both are multimappers 
      $flags{"raw"}{"read1"}{"IE,TP"}++;
      $flags{"raw"}{"read2"}{"IE,TP"}++;
      $flags{"raw"}{"counts"}{"lowqual"}++;
      return(1);
    }elsif($is_chimera eq "false"){
      $flags{"raw"}{"read1"}{"IE,TP"}++;
      $flags{"raw"}{"read2"}{"IE,TP"}++;
      $flags{"raw"}{"counts"}{"uninformative"}++;
      return(1);    
    }
  }
  ## Filter out reads if they have 2Ns
  my $num_N1 = $read1[9] =~ tr/N//;
  my $num_N2 = $read2[9] =~ tr/N//;    
  if($num_N1 >=2 and $num_N2>=2){  ## both reads carrying at least N bases
      $flags{"raw"}{"counts"}{"nn"}++;
      $flags{"raw"}{"read1"}{"IE,FP"}++;
      $flags{"raw"}{"read2"}{"IE,FP"}++;
      return(1);
  }  
  ## count total reads used in the pipeline
  $flags{"raw"}{"counts"}{"total_read_pairs"}++; 
 
  
  ## gather evidence
  my @evi1 = @{check_REMotif_presence_and_gain($read1[9],$read1[10],$read1[5],$read1[2],$read1[3],$strand1,$read1[4])};
  my @evi2 = @{check_REMotif_presence_and_gain($read2[9],$read2[10],$read2[5],$read2[2],$read2[3],$strand2,$read2[4])};  
  #my ($seq,$qual,$cigar,$chr,$start,$strand, $mapq) = @_;
   
  if($evi1[0] eq "" and $evi1[1] eq "" and $evi2[0] eq "" and $evi2[1] eq ""){
     $flags{"raw"}{"counts"}{"uninformative"}++;   ## mate evidence is not defined
     return(1);
  }
  if($evi1[0] eq "DE" and $evi1[1] eq "FP" and $read1[4]<$min_mapq){
     @evi1 = ("IE","TP",0,"-","-",$read1[9],$read1[3],$read1[10]);
  }elsif($evi2[0] eq "DE" and  $evi2[1] eq "FP" and $read2[4]<$min_mapq ){
     @evi2 = ("IE","TP",0,"-","-",$read2[9],$read2[3],$read2[10]);     
  } 
  $flags{"raw"}{"read1"}{$evi1[0].",".$evi1[1]}++;
  $flags{"raw"}{"read2"}{$evi2[0].",".$evi2[1]}++;
      
  if($evi1[0] eq "IE" and $evi1[1] eq "FP" and $evi2[0] eq "IE" and $evi2[1] eq "FP"){
     $flags{"raw"}{"counts"}{"uninformative"}++;   ## both mates mapped uniquely 
     return(1);
  }elsif($evi1[0] eq "IE" and $evi1[1] eq "TP" and $evi2[0] eq "IE" and $evi2[1] eq "TP"){
     $flags{"raw"}{"counts"}{"uninformative"}++;   ## both mates with low qual
     return(1);
  }elsif($evi1[0] eq "DE" and $evi1[1] eq "FP" and $evi2[0] eq "DE" and $evi2[1] eq "FP"){
      ## each mate has >1 ligation junction (fragment with short insert, ligation artifact)
      $flags{"raw"}{"counts"}{"twojunctoins"}++;   
      $flags{"raw"}{"counts"}{"uninformative"}++;   
      return(1);
  }elsif($evi1[1] eq "FP" and $evi2[0] eq "IE" and $is_chimera eq "false"){
      $flags{"raw"}{"counts"}{"uninformative"}++;   ## WGS-like pair, not informative 
      return(1);        
  }elsif($evi1[0] eq "IE" and $evi2[1] eq "FP" and $is_chimera eq "false"){
      $flags{"raw"}{"counts"}{"uninformative"}++;   ## WGS-like pair, not informative 
      return(1);        
  }
      
  if(($evi1[1] eq "TP" or $evi2[1] eq "TP") and $is_chimera eq "true"){
     $flags{"raw"}{"counts"}{"nonHicChimera"}++;
  }
  
  ## Use bam entry in the read header 
  ## For DE,TP/ DE,FP/ IE/FP, use the bam entry for the mate
  ## For IE,TP, use the bam entry for the other mate
  my $addum1 = "-";
  my $addum2 = "-";
  if($evi1[0] eq "IE" and $evi1[1] eq "TP"){
      $addum1 = join("\x{019}",@read2);
      $addum1 .= "\x{019}OP:Z:evi=IE,FP,".join(",",@evi2[2..3]).";side=".$evi2[4].";clip=".$evi2[6].";both=".$evi2[8];
      $addum1 = "-" if(!$chrs{$read2[2]});
  }else{
      $addum1 = join("\x{019}",@read1);
      $addum1 .= "\x{019}OP:Z:evi=".join(",",@evi1[0..3]).";side=".$evi1[4].";clip=".$evi1[6].";both=".$evi1[8]; 
      $addum1 = "-" if(!$chrs{$read1[2]});    
  }

  if($evi2[0] eq "IE" and $evi2[1] eq "TP"){ 
      $addum2 = join("\x{019}",@read1);
      $addum2 .= "\x{019}OP:Z:evi=IE,FP,".join(",",@evi1[2..3]).";side=".$evi1[4].";clip=".$evi1[6].";both=".$evi1[8];
      $addum2 = "-" if(!$chrs{$read1[2]});
  }else{
      $addum2 = join("\x{019}",@read2);
      $addum2 .= "\x{019}OP:Z:evi=".join(",",@evi2[0..3]).";side=".$evi2[4].";clip=".$evi2[6].";both=".$evi2[8];
      $addum2 = "-" if(!$chrs{$read2[2]});
  }
  
  undef(@read1) if($addum1 eq "-");
  undef(@read2) if($addum2 eq "-");      
  undef(@read1) if(@read1 and $num_N1>1);
  undef(@read2) if(@read2 and $num_N2>1);
  ## also remove reads that has >map_qual map score , as re-mapping them to TE assembly is futile
  undef(@read1) if(@read1 and $evi1[0] eq "IE" and $evi1[1] eq "FP");
  undef(@read2) if(@read2 and $evi2[0] eq "IE" and $evi2[1] eq "FP");
  
  ## checks!! 
  if( (@read1 and !defined $addum1) or (@read2 and !defined $addum2) ){
    print "\nERROR while classifying the reads\nsubroutine: write_fq_from_raw_pairsam\n";
    print join("\t",@read1),"\n";
    print join("\t",@read2),"\nExiting !!\n";
    exit 1;      
  }
  if(@read1 or @read2){   $flags{"reported"}{"class"}{$r1_type}++;    }      
  
  if(@read1 and defined $addum1){
    #print join("\t",@read1),"\n",join("\t",@evi1),"\n",$addum1,"\n";
    $flags{"reported"}{"read1"}{$evi1[0].",".$evi1[1]}++;
    if($evi1[0] eq "DE" and $evi1[1] eq "TP"){
      print O1 "@".$addum1."\n".$evi1[5]."\n+\n".$evi1[7]."\n"; 
      if(scalar(@evi1)==18 and $evi1[9] eq "DE" and $evi1[10] eq "TP"){
        $addum1 = join("\x{019}",@read1);
        $addum1 .= "\x{019}OP:Z:evi=".join(",",@evi1[9..12]).";side=".$evi1[13].";clip=".$evi1[15].";both=".$evi1[17];
        print O1 "@".$addum1."\n".$evi1[14]."\n+\n".$evi1[16]."\n";            
      }
    }elsif($is_chimera eq "true"){
      print O2 "@".$addum1."\n".$evi1[5]."\n+\n".$evi1[7]."\n"; 
    }
  }
  if(@read2 and defined $addum2){
    #print join("\t",@read2),"\n",join("\t",@evi2),"\n",$addum2,"\n";  
    $flags{"reported"}{"read2"}{$evi2[0].",".$evi2[1]}++;
    if($evi2[0] eq "DE" and $evi2[1] eq "TP"){
      print O1 "@".$addum2."\n".$evi2[5]."\n+\n".$evi2[7]."\n";
      if(scalar(@evi2)==18 and $evi2[9] eq "DE" and $evi2[10] eq "TP"){
        $addum2 = join("\x{019}",@read2);
        $addum2 .= "\x{019}OP:Z:evi=".join(",",@evi2[9..12]).";side=".$evi2[13].";clip=".$evi2[15].";both=".$evi2[17];
        print O1 "@".$addum2."\n".$evi2[14]."\n+\n".$evi2[16]."\n";            
      }
    }elsif($is_chimera eq "true"){
      print O2 "@".$addum2."\n".$evi2[5]."\n+\n".$evi2[7]."\n";
    }
  }
  return(1);
}

sub check_REMotif_presence_and_gain {
  my ($seq,$qual,$cigar,$chr,$start,$strand,$mapq) = @_;
  my $single= "0";
  my $check="F";
  my @evi = ("","","","","","","","",$single);  
  ## (0)DE,(1)TP,(2)SClen,(3)dist2RE, (4)clip-sidedness, (5)seq, (6)cliploc, (7)qual, (8)double/single clip
  
  ## checks
  if(!defined $seq or !defined $qual){
    print " ISSUE: $cigar,$chr,$start,$strand,$mapq\n";
    return(\@evi);
  }elsif(length$seq != length $qual){
    print " ERROR: Sequence length and quality string length do not match\n";
    exit 1;
  }elsif($cigar eq "*" or $chr eq "*"){  ## unmapped read
     @evi = ("IE","TP",0,"-","-",$seq,$start,$qual,$single);
     return(\@evi);
  }
  
  my ($a) = $cigar=~ /^(\d+)S\S+/;
  my ($b) = $cigar=~ /\D(\d+)S$/;
  $a=0 if(!defined($a) || $a eq"");  ## left hand side
  $b=0 if(!defined($b) || $b eq"");  ## right hand side   
  if($a < $clip and $b < $clip){
   if($mapq < $min_mapq){   
     @evi = ("IE","TP",0,"-","-",$seq,$start,$qual,$single);  ## this will be kept
     return(\@evi);
   }else{
     @evi = ("IE","FP",0,"-","-",$seq,$start,$qual,$single);  ## this will be discarded ultimately
     return(\@evi); 
   }    
  } 
  
  my ($dista,$distb)="";
  my ($xa,$xb)=0;
  if($a>=10){
    $dista = src::Utilities::checkREatClip(substr($seq,$a-$redb{$enzyme}{mlen},2*$redb{$enzyme}{mlen}),$redb{$enzyme}{motif},-1);
    if($dista!=0){
      $dista = src::Utilities::checkREatClip($seq,$redb{$enzyme}{ligmotif},$a);
    }
    $xa = src::Utilities::get_clip_coordV1($start,$strand,$cigar,"lhs");
  }
  if($b>=10){
    $distb = src::Utilities::checkREatClip(substr($seq,(length($seq)-$b-$redb{$enzyme}{mlen}),2*$redb{$enzyme}{mlen}),$redb{$enzyme}{motif},-1);
    if($distb!=0){
      $distb = src::Utilities::checkREatClip($seq,$redb{$enzyme}{ligmotif},(length($seq)-$b));
    }
    $xb = src::Utilities::get_clip_coordV1($start,$strand,$cigar,"rhs");
  }
  $check="T" if($a>=10 and $b>=10 and $dista>=10 and $distb>=10);

  if($a>=$clip and $b>=$clip){ ## return both clip locations 
     my @aa = ("DE","FP",$a,$dista,"a",substr($seq,0,$a),$xa,substr($qual,0,$a),$single);
     $aa[1] = "TP" if($dista >= $clip);
     $aa[8] = $xb if($check eq "T" ); 
     
     my @ab = ("DE","FP",$b,$distb,"b",substr($seq,-$b),$xb,substr($qual,-$b),$single);
     $ab[1] = "TP" if($distb >= $clip);
     $ab[8] = $xa if($check eq "T" ); 
  
     if($aa[1] eq "TP" and $ab[1] eq "FP"){
        @evi = @aa;
     }elsif($aa[1] eq "FP" and $ab[1] eq "TP"){
        @evi = @ab;
     }elsif($dista>$distb){
        @evi = (@aa,@ab);
     }else{
       @evi = (@ab,@aa);
     }  
  }elsif($a>=$clip and $b<$clip){
    @evi = ("DE","FP",$a,$dista,"a",substr($seq,0,$a),$xa,substr($qual,0,$a),$single);
    $evi[1] = "TP" if($dista >= $clip);
    $evi[8] = $xb if($check eq 'T');

  }elsif($b>=$clip and $a<$clip){
    @evi = ("DE","FP",$b,$distb,"b",substr($seq,-$b),$xb,substr($qual,-$b),$single);
    $evi[1] = "TP" if($distb >= $clip);
    $evi[8] = $xa if($check eq 'T');
  }
  return(\@evi);
}
#####-------------------------------------------- END OF SCRIPT ------------------------------------------------------------------ 