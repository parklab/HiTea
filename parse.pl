#!/usr/bin/perl

# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard

# This script parses Hi-C bam file (duplicate-marked lossless bam) and to generate fastq file for mapping on to the TE assembly 

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use Utilities;

use open qw(:std :utf8);
BEGIN { our $start_run = time(); }

###-------------------------------------------------------------------------------------------------------------------------------
# inputs etc
###-------------------------------------------------------------------------------------------------------------------------------
my $psam = "";
my $bam = "";
my $wd = "";
my $outprefix = "project"; #default 
my $min_mapq = 28;  #default
my $clip = 20; #default 
my $dist2re =3;  #default
my $motif = "GATC"; #default 
my $help = 0;
my $subset_fraction=0.33;
Getopt::Long::GetOptions(
  'psam=s'           => \$psam,
  'bam=s'            => \$bam,
  'm=s'              => \$motif,
  'wd:s'             => \$wd,
  'outprefix:s'      => \$outprefix,
  'q=s'              => \$min_mapq,
  'clip=s'           => \$clip,
  'subset_fraction=s' => \$subset_fraction,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
if ($help) {
  print "\nUsage: perl parse.pl -psam [FILE_PATH] -bam [FILE_PATH] -m [STRING] -outprefix [STRING] -wd [WORK_DIR] -clip [INT] -subset_fraction [Fraction] -q [INT] \n\n";
  print "This script parses pairsam file to extract discordent reads for assessment of TE-insertions. \n\n";
  print "Options (required*):\n";
  print "   -psam             Input file (pairsam format)\n";
  print "   -bam              Input file (either lossless bam format or WGS bam)\n";
  print "Options (optional):\n";
  print "   -m                RE Motif sequence (default:GATC)\n";
  print "   -outprefix        Output file PREFIX (default: project)\n";
  print "   -wd               Working directory (default: ~)\n";
  print "   -q                Minimum mapping quality for the reference alignment. This value is used to determine repeat anchored mates in the genome (default:28) \n";
  print "   -clip             Minimum softclipped read length for mapping the reads to the TE assembly (default:20)\n";
  print "   -help|-h          Display usage information.\n";
  print "Default outputs:\n";
  print "    Writes fastq file with discordent reads \n\n\n";
  exit 0;
}
$wd =~ s/\/$//;
if($wd ne ""){
  $outprefix = $wd."/".$outprefix;
}
##################################################################################################################################
## Soft-Inputs
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
    "!" => '1',
);

###-------------------------------------------------------------------------------------------------------------------------------
# I/O
###-------------------------------------------------------------------------------------------------------------------------------
my $watch_run=0;
my $run_time=0;
my %flags;

############# (1) locate clipped and RAM readssssss, print them to fastq files
#############     print fragment ends and pair-orientations  
my $outfq=$outprefix.".temp.fq.gz";  ## softclipped read sequences
open(O1, "| gzip -c - > $outfq") or die "can't create $outfq";
my $outfq2=$outprefix.".temp2.fq.gz";   ## non-clip reads
open(O2, "| gzip -c - > $outfq2") or die "can't create $outfq2";
#my $out3=$outprefix.".temp3.txt.gz";   ## non-clip reads
#open(O3, "| gzip -c - > $out3") or die "can't create $out3";

$watch_run = time();
$run_time = $watch_run - our $start_run;
my $lines = 0;
if($bam ne "" and $psam eq ""){
  print " Input bam file: $bam\t  $run_time ...";
  $lines = bam_read($bam);
}else{
  print " ERROR: Incorrect input file specified. Exiting!! \n";
  exit 1;
}

my $sumObj=$outprefix.".summary.log.ph"; ## store summary info
store \%flags, $sumObj;

$watch_run = time();
$run_time = $watch_run - $start_run;
print " $run_time seconds\t lines in the file: $lines\n";

close(O1);
close(O2);
#close(O3);
exit 0;

###-------------------------------------------------------------------------------------------------------------------------------
# Subroutines
###-------------------------------------------------------------------------------------------------------------------------------
  ## What are the Hi-C discordent reads?
  # 1) Split-reads: 
  #     a. clipped mate without ligation motif 
  #     Consider these reads as single end read support
  # 2) Unsplit-reads: 
  #     a. One mate with unique mapping(mapq>=28), while other with ambiguous mapping (mapq<28)
  #     b. One mate with unique mapping(mapq>=28), while other mate is unmapped 
  #     c. Unique mapping mate can have Hi-C chimeric junction
  #  Flags on the read:
  #   DE : >=20bp (or clip-length) clipped bases
  #   IE : Indirect support towards insertion/no clipping of the read
  #   FP : Carrying ligation junction
  #   TP : No ligation junction
  
sub bam_read{
  my ($file) = shift;
  if($file=~ m/.bam/){
     open IN,"samtools view -@ 8 $file|" or next "Can't open file $file";
  }else{
     print "input is in sam format\n";
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
       $is_next=1 if(rand()>$subset_fraction);  ###################################### TEMPORARY
       
       if($is_next==0){
         $flags{"raw"}{"counts"}{"total"}++;
         ### Filter out the reads
         if(!defined $read1[1] or !defined $read2[1]){
           $flags{"raw"}{"counts"}{"lonelymates"}++;
           $is_next=1;
         }elsif($read1[2] eq "*" and $read2[2] eq "*") {
           $flags{"raw"}{"counts"}{"uninformative"}++; 
           $is_next=1;
         }elsif(length($read1[9])<50 or length($read2[9])<50) {
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
  $flags{"cov"}{$read1[2]}{ int ($read1[3]/1000) }++ if($read1[3]!=0 and $read1[4] >= $min_mapq);
  $flags{"cov"}{$read2[2]}{ int ($read2[3]/1000) }++ if($read2[3]!=0 and $read2[4] >= $min_mapq);
  
  ## filters 2
  if($r1_type eq "UR" or $r1_type eq "RU" or $r1_type eq "UU" or ($r1_type eq "WGS" and $is_chimera eq "false")){
    $flags{"raw"}{"counts"}{"linear_rescued"}++; 
    if($read1[2] eq $read2[2] and $read1[8] !=0){
      my $logdist = Utilities::round(log(abs($read1[8]))/log(10),1);
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
  if($num_N1 >1 and $num_N2>1){  ## both reads carrying more than 1 N bases
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
      $flags{"raw"}{"counts"}{"uninformative"}++;   ## WGS non-chinera pair, not informative 
      return(1);        
  }elsif($evi1[0] eq "IE" and $evi2[1] eq "FP" and $is_chimera eq "false"){
      $flags{"raw"}{"counts"}{"uninformative"}++;   ## WGS non-chinera pair, not informative 
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
      $addum1 .= "\x{019}OP:Z:evi=IE,FP,".join(",",@evi2[2..3]).";side=".$evi2[4].";clip=".$evi2[6];
      $addum1 = "-" if(!$chrs{$read2[2]});
  }else{
      $addum1 = join("\x{019}",@read1);
      $addum1 .= "\x{019}OP:Z:evi=".join(",",@evi1[0..3]).";side=".$evi1[4].";clip=".$evi1[6]; 
      $addum1 = "-" if(!$chrs{$read1[2]});    
  }

  if($evi2[0] eq "IE" and $evi2[1] eq "TP"){ 
      $addum2 = join("\x{019}",@read1);
      $addum2 .= "\x{019}OP:Z:evi=IE,FP,".join(",",@evi1[2..3]).";side=".$evi1[4].";clip=".$evi1[6];
      $addum2 = "-" if(!$chrs{$read1[2]});
  }else{
      $addum2 = join("\x{019}",@read2);
      $addum2 .= "\x{019}OP:Z:evi=".join(",",@evi2[0..3]).";side=".$evi2[4].";clip=".$evi2[6];
      $addum2 = "-" if(!$chrs{$read2[2]});
  }
  
  #print join("\t",@read1),"\n";
  #print join ("\t",@evi1),"\n";
  #print $addum1,"\n\n";
  #print join("\t",@read2),"\n";
  #print join ("\t",@evi2),"\n";
  #print $addum2,"\n\n";

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
           if(scalar(@evi1)==16 and $evi1[8] eq "DE" and $evi1[9] eq "TP"){
              $addum1 = join("\x{019}",@read1);
              $addum1 .= "\x{019}OP:Z:evi=".join(",",@evi1[8..11]).";side=".$evi1[4].";clip=".$evi1[6];
              print O1 "@".$addum1."\n".$evi1[13]."\n+\n".$evi1[15]."\n";            
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
          if(scalar(@evi2)==16 and $evi2[8] eq "DE" and $evi2[9] eq "TP"){
              $addum2 = join("\x{019}",@read2);
              $addum2 .= "\x{019}OP:Z:evi=".join(",",@evi2[8..11]).";side=".$evi2[4].";clip=".$evi2[6];
             print O1 "@".$addum2."\n".$evi2[13]."\n+\n".$evi2[15]."\n";            
          }
        }elsif($is_chimera eq "true"){
             print O2 "@".$addum2."\n".$evi2[5]."\n+\n".$evi2[7]."\n";
        }
  }
  return(1);
}

sub check_REMotif_presence_and_gain {
  my ($seq,$qual,$cigar,$chr,$start,$strand,$mapq) = @_;
  my @evi = ("","","","","","","","");  ## (0)DE,(1)TP,(2)SClen,(3)dist2RE, (4)clip-sidedness, (5)seq, (6)cliploc, (7)qual
  ## checks
  if(!defined $seq or !defined $qual){
    print " ISSUE: $cigar,$chr,$start,$strand,$mapq\n";
    return(\@evi);
  }elsif(length$seq != length $qual){
    print " ERROR: Sequence length and quality string length do not match\n";
    exit 1;
  }
  if($cigar eq "*" or $chr eq "*"){  ## unmapped read
     @evi = ("IE","TP",0,"-","-",$seq,$start,$qual);
     return(\@evi);
  }

  my ($a) = $cigar=~ /^(\d+)S\S+/;
  my ($b) = $cigar=~ /\D(\d+)S$/;
  $a=0 if(!defined($a) || $a eq"");  ## left hand side
  $b=0 if(!defined($b) || $b eq"");  ## right hand side   
  if($a < $clip and $b < $clip){
    if($mapq < $min_mapq){   
       @evi = ("IE","TP",0,"-","-",$seq,$start,$qual);  ## this will be kept
       return(\@evi);
    }else{
      @evi = ("IE","FP",0,"-","-",$seq,$start,$qual);  ## this will be discarded ultimately
      return(\@evi); 
    }    
  } 
  
  my $motif_len = length($motif);
  my $dista="";
  my $distb="";  
  my $xa=0;
  my $xb=0;
  if($a>=$clip){
    $dista = get_dist2RE_fromSeq($seq,$motif,$a,"a");
    $xa = Utilities::get_clip_coordV1($start,$strand,$cigar,"lhs");
  }
  if($b>=$clip){
    $distb = get_dist2RE_fromSeq($seq,$motif,(length($seq)-$b),"b");
    $xb = Utilities::get_clip_coordV1($start,$strand,$cigar,"rhs");
  }

  if($a>=$clip and $b>=$clip){ ## return both clip locations 
     my @aa;
     my @ab;
     
     if($a>=$clip){
        $aa[0] = "DE";
     }
     if($dista < $dist2re){
        $aa[1] = "FP";
      }else{
        $aa[1] = "TP";
     }
     $aa[2] = $a;
     $aa[3] = $dista;
     $aa[4] = "a";
     $aa[5] = substr($seq,0,$a);
     $aa[6] = $xa;
     $aa[7] = substr($qual,0,$a);

     if($b>=$clip){
        $ab[0] = "DE";
     }
     if($distb < $dist2re){
        $ab[1] = "FP";
      }else{
        $ab[1] = "TP";
     }
     $ab[2] = $b;
     $ab[3] = $distb;
     $ab[4] = "b";
     $ab[5] = substr($seq,-$b);;
     $ab[6] = $xb;
     $ab[7] = substr($qual,-$b);

     if($aa[1] eq "TP" and $ab[1] eq "FP"){
        @evi = @aa;
     }elsif($aa[1] eq "FP" and $ab[1] eq "TP"){
        @evi = @aa;
     }elsif($dista>$distb){
        @evi = (@aa,@ab);
     }else{
       @evi = (@ab,@aa);
     }  
  }elsif($a>=$clip and $b<$clip){
    if($a>=$clip){
        $evi[0] = "DE";
    }
    if($dista<$dist2re){
      $evi[1] = "FP";
     }else{
      $evi[1] = "TP";
    }
    $evi[2] = $a;
    $evi[3] = $dista;
    $evi[4] = "a";
    $evi[5] = substr($seq,0,$a);
    $evi[7] = substr($qual,0,$a);
    $evi[6] = $xa;
  }elsif($b>=$clip and $a<$clip){
    if($b>=$clip){
        $evi[0] = "DE";
    }
    if($distb<$dist2re){
      $evi[1] = "FP";
     }else{
      $evi[1] = "TP";
    }
    $evi[2] = $b;
    $evi[3] = $distb;
    $evi[4] = "b";
    $evi[5] = substr($seq,-$b);
    $evi[7] = substr($qual,-$b);
    $evi[6] = $xb;
  }
  return(\@evi);
}

sub get_dist2RE_fromSeq{
  my ($seq,$motif,$clip,$side) = @_;
  my $dist = 500;
  my $offset = 0;
  ## If clip seq contains RE motif, it should be atleast $clip bp away from the clipped position in the read
  ## This ensures removal of Hi-C chimera derived through unmapped clip reads
     
  if($motif eq "" or !defined $motif){
    print " Motif not defined in Dist2RE_fromSeq sub!! exiting! \n";
    exit 1;
  }elsif($motif eq "GATC"){
     my $cnt=() = $seq =~ /GATCGATC/g;
     return("0") if($cnt>1); 
     return($dist) if($cnt == 0);
     my $c = index($seq, "GATCGATC",0);
     return($dist) if($c == -1);
     $clip -=1 if($clip !=0); 
     if( ($c-$dist2re) <= $clip and ($c+8+$dist2re) >= $clip ){
        $dist = 0;     
     }elsif($c>$clip){
        $dist = $c - $clip;
        $offset = 1;
     }else{
        $dist = ($clip - $c - 8);
        $offset = 2;
     }
     my $mlen=length($motif);
     my $slen = length($seq);
     if( ($mlen+$clip) <= $slen and $clip>=$mlen){
        my $sq = substr($seq,($clip-$mlen),($clip+$mlen));
        if($sq =~ m/GATC/){
          $dist=0;
        }
     }
     if($side eq "a" and $offset == 2 and $dist<=$clip){
        $dist=0;
     }elsif($side eq "b" and $offset == 1 and $dist<=$clip){
        $dist=0;
     }
  }elsif($motif eq "AAGCTT"){
     my $cnt=() = $seq =~ /AAGCTAGCTT/g;
     return("0") if($cnt>1); 
     return($dist) if($cnt == 0);
     my $c = index($seq, "AAGCTAGCTT",0);
     return($dist) if($c == -1);
     $clip -=1 if($clip !=0); 
     if( ($c-$dist2re) <= $clip and ($c+10+$dist2re) >= $clip ){
        $dist = 0;     
     }elsif($c>$clip){
        $dist = $c - $clip;
        $offset = 1;
     }else{
        $dist = ($clip - $c - 10);
        $offset = 2;
     }     
     my $mlen=length($motif);
     my $slen = length($seq);
     if( ($mlen+$clip) <= $slen and $clip>=$mlen){
        my $sq = substr($seq,($clip-$mlen),($clip+$mlen));
        if($sq =~ m/AAGCTT/){
          $dist=0;
        }
     }
     if($side eq "a" and $offset == 2 and $dist<=$clip){
        $dist=0;
     }elsif($side eq "b" and $offset == 1 and $dist<=$clip){
        $dist=0;
     }
  }elsif($motif eq "GCGGCCGC"){
     my $cnt=() = $seq =~ /GCGGCCGGCCGC/g;
     return("0") if($cnt>1); 
     return($dist) if($cnt == 0);
     my $c = index($seq, "GCGGCCGGCCGC",0);
     return($dist) if($c == -1);
     $clip -=1 if($clip !=0); 
     if( ($c-$dist2re) <= $clip and ($c+12+$dist2re) >= $clip ){
        $dist = 0;     
     }elsif($c>$clip){
         $offset = 1;
         $dist = $c - $clip;
     }else{
         $offset = 2;
        $dist = ($clip - $c - 12);
     }     
     my $mlen=length($motif);
     my $slen = length($seq);
     if( ($mlen+$clip) <= $slen and $clip>=$mlen){
        my $sq = substr($seq,($clip-$mlen),($clip+$mlen));
        if($sq =~ m/GCGGCCGC/){
          $dist=0;
        }
     }
     if($side eq "a" and $offset == 2 and $dist<=$clip){
        $dist=0;
     }elsif($side eq "b" and $offset == 1 and $dist<=$clip){
        $dist=0;
     }
  }elsif($motif eq "CCATGG"){
     my $cnt=() = $seq =~ /CCATGCATGG/g;
     return("0") if($cnt>1); 
     return($dist) if($cnt == 0);
     my $c = index($seq, "CCATGCATGG",0);
     return($dist) if($c == -1);
     $clip -=1 if($clip !=0); 
     if( ($c-$dist2re) <= $clip and ($c+10+$dist2re) >= $clip ){
        $dist = 0;     
     }elsif($c>$clip){
        $dist = $c - $clip;
        $offset = 1;
     }else{
        $dist = ($clip - $c - 10);
        $offset = 2;
     }
     my $mlen=length($motif);
     my $slen = length($seq);
     if( ($mlen+$clip) <= $slen and $clip>=$mlen){
        my $sq = substr($seq,($clip-$mlen),($clip+$mlen));
        if($sq =~ m/CCATGG/){
          $dist=0;
        }
     }
     if($side eq "a" and $offset == 2 and $dist<=$clip){
        $dist=0;
     }elsif($side eq "b" and $offset == 1 and $dist<=$clip){
        $dist=0;
     }
  }
  return($dist);
}

sub parse_psam_reads {
  my ($side1, $side2) =@_;
  $side1=~ s/\n//; 
  $side1=~ s/\r//;  
  $side2=~ s/\n//; 
  $side2=~ s/\r//;  
  my @a = split("NEXT_SAM",$side1."NEXT_SAM".$side2);
  my @r1;
  my @r2;
  foreach my $l (@a){
    $l =~ s/^\x{019}*//;
    $l =~ s/\x{019}*$//;
    my ($flag) =  $l =~ /^.*?\x{019}(\d*)/;
    next if($flag>=256);
    @r1 = split("\x{019}",$l) if($flag & 0x40);
    @r2 = split("\x{019}",$l) if($flag & 0x80);
  }

  return(\@r1,\@r2); 
}
#####-------------------------------------------- END OF SCRIPT ------------------------------------------------------------------ 