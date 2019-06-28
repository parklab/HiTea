#!/usr/bin/perl
# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard
package Utilities;
use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use constant { true => 1, false => 0 };
use utf8;
use Storable;
use open qw(:std :utf8);

###-------------- subroutines
sub get_clip_coord{
  my ($start,$cigar) = @_;
  #print $cigar,"\n";

  if($cigar eq "*" or $cigar eq ""){
    return($start);
  }

  my @len1 = split (/\D+/,$cigar); # storing the length per operation
  my @ops1 = split (/\d+/,$cigar); # storing the operation
  shift @ops1; # remove the empty first element
  if(!defined $ops1[0]){ 
    print " Error: $start $cigar\n";
  }
  my $x=0;
  my $len1 =0;
  my $i=0;
  for($i=0; $i<scalar@ops1; $i++){
      if($ops1[$i] eq "M" or $ops1[$i] eq "D" ){
        $len1 += $len1[$i]; 
      }
  }
  
  if($ops1[0] eq "S" && $ops1[-1] ne "S"){
    $x=$start;
  }
  elsif($ops1[-1] eq "S" && $ops1[0] ne "S"){
    $x=$start+$len1;
  }
  elsif($ops1[0] eq "S" &&  $ops1[-1] eq "S"){  #$len1[0] > $len1[-1]
    $x=$start; 
  }
  elsif($ops1[0] eq "S" &&  $ops1[-1] eq "S"){ #$len1[-1] > $len1[0]
    $x=$start+$len1;
  }
  else{
    $x=$start;
  }
  return($x);
}

sub get_clip_coordV1{
  my ($start,$strand,$cigar,$side) = @_;
  #print $cigar,"\n";

  if(!($side eq "lhs" or $side eq "rhs")){
    print "WARNING: $start,$strand,$cigar. No valid side argument provided\n"; 
    return(0);
  }

  if($cigar eq "*" or $cigar eq ""){
    return($start);
  }

  my @len1 = split (/\D+/,$cigar); # storing the length per operation
  my @ops1 = split (/\d+/,$cigar); # storing the operation
  shift @ops1; # remove the empty first element

  my $x=0;
  my $len1 =0;
  my $i=0;
  for($i=0; $i<scalar@ops1; $i++){
      if($ops1[$i] eq "M" or $ops1[$i] eq "D" ){
        $len1 += $len1[$i]; 
      }
  }

  if($ops1[0] eq "S" && $ops1[-1] ne "S"){
    $x=$start;
  }
  elsif($ops1[-1] eq "S" && $ops1[0] ne "S"){
    $x=$start+$len1;
  }
  elsif($ops1[0] eq "S" &&  $ops1[-1] eq "S" && $side eq "lhs"){  #$len1[0] > $len1[-1]
    $x=$start; 
  }
  elsif($ops1[0] eq "S" &&  $ops1[-1] eq "S" && $side eq "rhs"){ #$len1[-1] > $len1[0]
    $x=$start+$len1;
  }
  else{
    $x=$start;
  }
  return($x);
}

sub get_read_mapping_span {
  my ($start,$cigar) =@_;
  if($cigar eq "*"){
    return($start);
  }

  my @len1 = $cigar =~ m/(\d+)I/g;
  my @len2 = $cigar =~ m/(\d+)M/g;
  
  my $length=0;
  foreach (@len1) { $length += $_; }
  foreach (@len2) { $length += $_; }
  return(($length+$start));
}

sub unique{
  my ($string,$delim,$number) = @_;
  if(!defined $delim){
    $delim = " ";
  }
  if(!defined $number){
    $number = "false";
  }

  $string =~ s/$delim$//;
  my @temp = split($delim,$string);
  my %hash;
  foreach my $i(@temp){
    next if($i eq "" or !defined $i);
    $hash{$i}++;
  }
  if($number eq "true"){
    return(scalar(keys%hash));
  }elsif($number eq "false"){
     return(join($delim,keys%hash));
  }else{
    return("ERROR");
  }  
}

sub round {
  my ($n, $places) = @_;
  my $abs = abs $n;
  my $val = substr($abs + ('0.' . '0' x $places . '5'),
                   0,
                   length(int($abs)) +
                     (($places > 0) ? $places + 1 : 0)
                  );
  ($n < 0) ? "-" . $val : $val;
}

sub get_softclipseq {
  my ($seq,$cigar) = @_;
    my ($a) = $cigar=~ /^(\d+)S\S+/;
    my ($b) = $cigar=~ /\D(\d+)S$/;
    $a=0 if(!defined($a) || $a eq"");  ## left hand side
    $b=0 if(!defined($b) || $b eq"");  ## right hand side 
    if($a < 5 and $b < 5){
      return("");
    }

    if($a > $b){
      return(substr($seq,0,$a));
    }elsif($b > $a){
      return(substr($seq,-$b));
    }else{
      return("");
    }
}

sub max{
  my @i = @_;
  @i = grep { $_ ne'NA' } @i;
  @i = grep { $_ ne'-' } @i;
  @i = grep { $_ ne'*' } @i;
  #my @i = split(",",$i);
  @i = sort{$b<=>$a} @i if(scalar(@i)>1);
  return($i[0]);
}

sub min{
  my @i = @_;
  @i = grep { $_ ne 'NA' } @i;
  @i = grep { $_ ne'-' } @i;
  @i = grep { $_ ne'*' } @i;
  #my @i = split(",",$i);
  @i = sort{$a<=>$b} @i if(scalar(@i)>1);
  return($i[0]);
}

sub get_fasta{
  my %seqs;
  my $header; 
  my $first = 0;
  
  open(FILE, "<@_") or die("Cannot open FASTA file\n");
  my @lines = <FILE>;
  foreach my $line(@lines){
    chomp($line);
    next if($line eq "");
    if ($line =~ /^>/){
      $header = $line;
      $header =~ s/^>//;
      $header =~ s/s.*//;
      $header =~ s/\n//g;
      $header =~ s/\r//;
      $header =~ s/\t.*$//;
      if ($first == 0){
        $first = 1;
      }
      next;
    }
    if ($first == 0){ die("Not a standard FASTA file\n"); }
    $seqs{$header} .= $line;
  }
  close(FILE);
  while(my($i,$j)=each%seqs){
    next if($i eq "");
    $seqs{$i} = length($j);
  }
  return %seqs;
}

# Input: #(0)seq, (1)cigar, (2)strand, (3)readstart
sub getstrand{
  my($r1,$r2) = @_;
  $r2=",,,0" if(!($r2));
  my @r1 = split(",",$r1);
  my @r2 = split(",",$r2);
    
    #(0)seq, (1)cigar, (2)strand, (3)readstart
    my ($a1) = $r1[1]=~ /^(\d+)S\S+/;
    my ($b1) = $r1[1]=~ /\D(\d+)S$/;
    my ($a2) = $r2[1]=~ /^(\d+)S\S+/;
    my ($b2) = $r2[1]=~ /\D(\d+)S$/;
    $a1=0 if(!defined($a1) || $a1 eq"");  ## left hand side
    $b1=0 if(!defined($b1) || $b1 eq"");  ## right hand side 
    $a2=0 if(!defined($a2) || $a2 eq"");  ## left hand side
    $b2=0 if(!defined($b2) || $b2 eq"");  ## right hand side 
    
    my $s1="";
    my $s2="";
    #my $s1c="";
    #my $s2c="";
    my $x=10; ## require that 7/10 should be either A or T

    if($a1>$b1){
      $x=$a1 if($a1<10);
      $s1 = substr($r1[0],($a1-$x),$x);
      #$s1c = substr($r1[0],0,$a1) if($a1>10);
    }else{
      $x=$b1 if($b1<10);     
      $s1 = substr($r1[0],(length($r1[0])-$b1),$x);
      #$s1c = substr($r1[0],-$b1) if($b1>10);
    }
    if($a2>$b2){
      $x=$a2 if($a2<10);
      $s2 = substr($r2[0],($a2-$x),$x);
      #$s2c = substr($r2[0],0,$a2) if($a2>10);
    }else{
      $x=$b2 if($b2<10);     
      $s2 = substr($r2[0],(length($r2[0])-$b2),$x);
      #$s2c = substr($r2[0],-$b2) if($b2>10);
    }

    if($s1 eq "" and $s2 eq ""){
      return("ERROR-Strand");
    }
    if(length($s1)==0 and length($s2)==0){
      print join("\t",@r1),"\n",join("\t",@r2),"\n";
      exit 1;
    }

    my %freq;
    $freq{"PolyA"} = Utilitiesround((($s1=~ tr/A//)/length($s1)),2);
    $freq{"PolyT"} = Utilitiesround((($s1=~ tr/T//)/length($s1)),2);
    
    if(length($s2)>0){
     if($freq{"PolyA"} < Utilitiesround((($s2=~ tr/A//)/length($s2)),2)){
         $freq{"PolyA"} = Utilitiesround((($s2=~ tr/A//)/length($s2)),2);
     }
     if($freq{"PolyT"} < Utilitiesround((($s2=~ tr/T//)/length($s2)),2)){
         $freq{"PolyT"} = Utilitiesround((($s2=~ tr/T//)/length($s2)),2);
     }
    }    
    #print "$r1[4]\t$r1[5]\t$r2[4]\t$r2[5]\n"; 
    #print "$s1\n$s2\n";
    #print $freq{"PolyA"},"\t",$freq{"PolyT"},"\n";  

    if($freq{"PolyA"} > $freq{"PolyT"} and $freq{"PolyA"} >= 0.7){
      return("PolyA");
    }elsif($freq{"PolyT"} >= 0.7){
      return("PolyT");
    }else{
        
      return("NA");
    }
}

sub get_teMap_clusterfreq {
  my ($j, $tr) = @_;
  $j =~ s/\n//g;
  $j =~ s/;$//; 
  my %transp = %{$tr};
  my $gdist=10;  # some threshold to merge TE mapping coordinate
  my %freq;
  my %pos;
  
  foreach my $te (sort keys %transp){
    $te = $te.":";
    my @y = grep(/$te/,split(";",$j));
    next if(scalar @y ==0);
    for (@y) {s/$te//g;}
    ## for PolyA tails, all mapping can be attributed to start.
    if($te eq "PolyA:"){
      @y = map { ($_ * 0) +1}  @y;
    }
    ## for SVAs, since there is a 30bp repeat expansion in the consensus, any mapping within this region should be awarded a start of 30
    ## Highly customized for SVAs based on the provided consensus
    if($te eq "SVA:"){
      @y = map { if($_ <30){ ($_ * 0) +30} else{ $_}}  @y;
    }

    my $num =1;
    if(scalar(@y)==1){
      $pos{$te.$num} = $y[0];
      $freq{$te.$num}=1;
    }else{
      @y = sort{$a <=> $b} @y;
      my $offset=0;
      for(my $i=0;$i<=(scalar(@y)-2);$i++){
        if(($y[($i+1)]-$y[$i])<$gdist){
          $offset++;
          $freq{$te.$num}++;
          if($offset==1){
            $freq{$te.$num}++;          
            $pos{$te.$num} .= $y[$i].",";
          }
          $pos{$te.$num} .= $y[$i+1].",";
        }else{
          $num++;
        }
      }
    }    
  } 

  ## Following part accomodates PolyA as separate TE class (needed when TEA based transposon assemlbies are used)
  my $out ="";
  my $i=0;
  my @sortedkeys;
  if(scalar (keys %freq)>1){
    @sortedkeys = sort{$freq{$b} <=> $freq{$a}} keys %freq;
  }else{
    @sortedkeys = keys %freq;
  }
  foreach my $o (@sortedkeys){
    $i++;
    if($i==1){
      $pos{$o} =~ s/,$//;
      my @x = split(",",$pos{$o});
      @x = sort {$a <=> $b} @x if(scalar@x >1);
      my $coord = $x[0];      
      my %x;
      for (@x) {$x{$_}++;}     
      if(scalar keys %x >1){
         @x = sort{ $x{$b} <=> $x{$a}} keys %x;
         $coord= $x[0];
      } 
      my $t = $freq{$o};
      $o =~ s/:\d*$//;
      if($out eq ""){
        $out .= $o.":".$t.":".$coord.";";
      }elsif($out ne "" and $o ne "PolyA"){
        $out .= $o.":".$t.":".$coord.";";
      }
      if($o eq "PolyA"){
        $i=0;
      }else{
        last;  
      }    
    }
  }
  $out =~ s/\n//g;
  $out =~ s/;$//;
     
  ## output   
  my @res = ("-","-","-","FALSE"); ## (0)RecrpClustLoc, (1)RecrpClustReads, (2)RecprTE, (3)RecprPolyA
  if($out eq ""){
    return(\@res); 
  }

  my @t = split(";",$out);
  if(scalar @t == 1){
    my ($ele,$numread, $clustloc) = $t[0] =~ m/(\S*):(\S*):(\S*)/;
    $res[0] = $clustloc;
    $res[1] = $numread;
    $res[2] = $ele;
    if($ele eq "PolyA"){
        $res[3] = "TRUE"; 
    }
  }elsif( scalar @t ==2){
    my ($ele,$numread, $clustloc) = $t[1] =~ m/(\S*):(\S*):(\S*)/;
    $res[0]  = $clustloc;
    $res[1] = $numread;
    $res[2] = $ele;
    $res[3] = "TRUE"; 
  }else{
    print "  error in finding reciprocal clusters while completing the annotation object\n";
  }
  return(\@res);
}

sub get_longestseq{
  my $temp = shift;
  if(!defined($temp) or $temp eq ""){
    return("");
  }
  $temp =~ s/,$//;
  my @seq = split(",",$temp);
  my %s;
  foreach my $sq(@seq){
    $s{$sq} =length($sq);
  }
  @seq = sort {$s{$b}<=>$s{$a}} keys%s;
  return($seq[0]);
}

sub getOriRatio{
  my ($i) =shift;
    $i =~ s/,$//;
    my @i = Utilities::unique($i,",","false");
    my %res;
    foreach my $l (@i){
      #print $l,"\n";
       my @m = $l =~ /\d+([+-])\d+/g;  # chr2156259243+101Mchr5126259981+11S90M_2
       my $Ori="";
       if(scalar@m==2){
         if($m[0] ne $m[1]){
            $Ori="FR";
          }else{
            $Ori="FF";
          }
        }
      $res{$Ori}++;
    }
    if(exists($res{"FF"}) and exists($res{"FR"}) ){
      return(Utilitiesround( ($res{"FR"})/($res{"FF"}+$res{"FR"})  ,2));
    }elsif(exists($res{"FF"}) and !exists($res{"FR"})){
      return("0.00");
    }elsif(!exists($res{"FF"}) and exists($res{"FR"})){
      return("1.00");
    }else{
      return("NA");
    }
}

## Input: # (0)seq, (1)cigar, (2)strand, (3)$read.start,(4)clust.start,(5)clust.end
sub tsd {
  my ($a) =@_;
  my @read = @{$a};
  #print join("\t",@read),"\t";
  if($read[4] == $read[5]){
    return("-");
  }

  if($read[4]<$read[3] or $read[5]<$read[3]){
    return("ERROR-TSD");
  }

  # (0)seq, (1)cigar, (2)strand, (3)$read.start,(4)clust.start,(5)clust.end
  my$s =$read[4]-$read[3];
  my$l = $read[5]-$read[4]+$s;

  if($read[1]=~ m/^(\d+)S\S+/){ ## correct the coords when softclipped bases present at hte lefthand mapping
    $s+=$1;
    $l+=$1;
  }

  my @len = split (/\D+/,$read[1]); # storing the length per operation
  my @ops = split (/\d+/,$read[1]); # storing the operation
  shift @ops; # remove the empty first element
  unless (scalar @len == scalar @ops){
    return("");
  }
  
  my $ops_string="";
  foreach my $i (1 .. scalar@ops){
    $ops_string .=$ops[$i-1] x $len[$i-1];
  }

  my @ops_string = split(//, $ops_string);
  my @seq = split(//, $read[0]);
  my $outseq="";
  for(my $i=$s; $i<=$l;$i++){
    next if(!$seq[$i]);
     $outseq .=$seq[$i];
     if($ops_string[$i] eq "I"){
        $l++;
      }
  }
  #print $outseq,"\n";
  return($outseq);
}


1;