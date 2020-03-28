#!/usr/bin/perl
# HiTEA

package src::Utilities;
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
  
  $cigar=~ s/H/S/g;  ## replace hardclip by softclip in CIGAR

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
  @temp = grep { $_ ne $delim } @temp;
  my %hash;
  foreach my $i(@temp){
    next if($i eq "" or !defined $i or $i eq $delim);
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
  
  open(FILE, "<@_") or die("Cannot open FASTA file @_\n");
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

sub get_fasta_seqs{
  my %seqs;
  my $header; 
  my $first = 0;
  
  open(FILE, "<@_") or die("Cannot open FASTA file @_\n");
  my @lines = <FILE>;
  foreach my $line(@lines){
    chomp($line);
    $line =~ s/\n//g;
    $line =~ s/\r//;
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
    $freq{"PolyA"} = round((($s1=~ tr/A//)/length($s1)),2);
    $freq{"PolyT"} = round((($s1=~ tr/T//)/length($s1)),2);
    
    if(length($s2)>0){
     if($freq{"PolyA"} < round((($s2=~ tr/A//)/length($s2)),2)){
         $freq{"PolyA"} = round((($s2=~ tr/A//)/length($s2)),2);
     }
     if($freq{"PolyT"} < round((($s2=~ tr/T//)/length($s2)),2)){
         $freq{"PolyT"} = round((($s2=~ tr/T//)/length($s2)),2);
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

sub common{
  my ($i,$j,$delim,$number) = @_;
  if($i eq "" or $i eq "-" or $i eq "*" or $j eq "" or $j eq "-" or $j eq "*"){
    return(0);
  }
  if(!defined $i or!defined $j){
    return(0);
  }

  $i =~ s/$delim$//;
  $j =~ s/$delim$//;
  my @temp = split($delim,$i);
  my %hash = map { $_ => 1 } @temp;
  my @temp2 = split($delim,$j);
  
  my $cnt=0;
  my @out;
  foreach my $l (@temp2){
    $cnt++ if($hash{$l});
    push(@out,$l) if($hash{$l});
  }

  if($number eq 'true'){
    return($cnt);
  }elsif($number eq "false"){
    return(join($delim,@out));
  }else{
    return(0);
  }  
}

sub get_reciprocalClust {
  my ($jj, $te,$ids) = @_;
  $jj =~ s/\n//g;
  $jj =~ s/,$//; 
  $ids =~ s/\n//g;
  $ids =~ s/,$//; 
  
  my $gdist=10;  # some threshold to merge TE mapping coordinate
  my @out = ("-",0,0, $ids); ## (0)pos, (1)freq, (2)fraction of total reads (3) clustered read-ids

  my %res;
  my @ids=split(",",$ids);
  my @y = split(",",$jj);
  my @idx = sort { $y[$a] <=> $y[$b] } 0 .. $#y;
  @y = @y[@idx];
  @ids = @ids[@idx];

  ## for PolyA tails, all mapping can be attributed to start.
  if($te eq "PolyA"){
    @y = map { ($_ * 0) +1}  @y;
  }
  ## for SVAs, since there is a 30bp repeat expansion in the consensus, any mapping within this region should be awarded a start of 30
    ## Highly customized for SVAs based on the provided consensus
  if($te eq "SVA"){
    @y = map { if($_ <30){ ($_ * 0) +30} else{ $_}}  @y;
  }

  if(scalar @y==1){
    @out =  ($y[0],1,1,$ids);
  }else{
    #159 160 174 182 182 183 187 188 191 192 194 194 197
    for(my $i=0; $i<=(scalar(@y)-1);$i++){ 
      $res{$i}{pos} .= $y[$i]." ";
      $res{$i}{ids} .= $ids[$i]." ";
      $res{$i}{freq}++;
      for(my $j=($i+1); $j<=(scalar(@y)-1);$j++){
        if( ($y[$j] - $y[$i]) <= $gdist){
          $res{$i}{pos} .= $y[$j]." ";
          $res{$i}{ids} .= $ids[$j]." ";
          $res{$i}{freq}++;
        }
      }
      #print $res{$i}{pos},"\n";
    }

    my @sortedkeys = sort{ no warnings;
                       my ($side1) = $res{$a}{freq};
                       my ($side2) = $res{$b}{freq};
                       $side2 <=> $side1; } keys %res;
    my @keys = grep {$res{$_}{freq} == $res{$sortedkeys[0]}{freq}} keys %res;
    @keys = sort{$a <=> $b } @keys;
    $out[1] = src::Utilities::unique($res{$keys[0]}{ids}," ","true");
    $out[3] = src::Utilities::unique($res{$keys[0]}{ids}," ","false");
    $out[3] =~ s/ /,/g;

    my $xx = src::Utilities::unique($ids,",","false");
    $xx = scalar(split(",",$xx));
    $out[2] = src::Utilities::round($out[1]/$xx,2);
    my %p;
    my @p = split(" ",$res{$keys[0]}{pos});

    if($out[1]>1){
      foreach(@p){ $p{$_}++;}
      if(scalar (keys %p)>1){
        my @sortedp = sort{$p{$b} <=> $p{$a} } keys %p;
        my @keysp = grep {$p{$_} == $p{$sortedp[0]}} keys %p;
        @keysp = sort{$a <=> $b} @keysp if(scalar @keysp>1);
        $out[0] = $keysp[0]; 
      }else{
        $out[0] = $p[0];
      }  
    }else{
      @out =("-",0,0,"");
    }
  }
  return(\@out);
}

sub array_sort{
  my $j = shift;
  $j=~ s/,$//;
  my @ttt = split(",",$j);
  @ttt = sort {$b <=> $a} @ttt if(scalar @ttt>1);
  return($ttt[0]);           
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
    my @i = unique($i,",","false");
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
      return(round( ($res{"FR"})/($res{"FF"}+$res{"FR"})  ,2));
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

sub checkREatClip{
  my ($seq,$motif,$cl) = @_;
  my $dist = 500;
  ## finds the nearest distance between the clip coordinate and ligation motif
  ## 0 if clip-coordinate overlaps ligation motif
  ## $cl ca be set to -1to check presence of motif in the sequence
  
  if($motif eq ""){
    print " Motif not defined!! exiting! \n";
    exit 1;
  }
  
  if($cl == -1){
    if($seq =~ m/$motif/g){
      $dist = 0;
    }
  }else{
     my @dist;
     while ($seq =~ /$motif/g){
       push @dist, abs($cl - 1 - $-[0]);
       push @dist, abs($cl - $+[0]);
       push (@dist, 0) if(($cl-1)>= $-[0] and $cl<= $+[0]);
       #print $cl,"\t",$-[0],"\t",$+[0],"\n";
    }
    @dist = sort {$a <=> $b} @dist;
    #print join("\t",@dist),"\n";
    $dist = $dist[0] if(scalar @dist>0); 
  }
  return($dist);
}

sub getATFreq{
  my ($rr1, $rr2) =@_;
  $rr2 = "." if(!defined $rr2 or length($rr2)==0);
  $rr1 = "." if(!defined $rr1 or length($rr1)==0);
  my @t;
  $t[0] = () = $rr1 =~ /A|a/g;
  $t[1] = () = $rr1 =~ /T|t/g;
  $t[2] = () = $rr2 =~ /A|a/g;
  $t[3] = () = $rr2 =~ /T|t/g;
  
  @t = ($t[0]/length($rr1), $t[1]/length($rr1), $t[2]/length($rr2), $t[3]/length($rr2));
  @t = sort{$b <=> $a}@t;
  my $x = src::Utilities::round($t[0],2);
  return($x);
}

sub binarysearch { #binarysearch($pos, \@locations)
  my ($x,$zz, $a) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);          # lower, upper end of search interval
  my $i;                               # index of probe
  my $q1=0;
  my $q2=0;

  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
    } 
    else {
      # return $i+1; #original 
      $q1 = abs($a->[$i] - $x);
      $q2 = $q1;  
      if($a->[$i+1]){
        $q2 = abs($a->[$i+1] - $x);
      }
      if($q1 > $q2){
        return($a->[$i+1]); # found
      }else{
        return($a->[$i]); # found
      }
    }
  }
  
  $q2 = abs($a->[$l-1] - $x);
  $q1 = $q2;
  if($a->[$l]){
    $q1 = abs($a->[$l] - $x);  
  }
  if($q1 < $q2){
    return($a->[$l]); # found
  }else{
    return($a->[$l-1]); # found
  }
}

sub hompolymer_check{
  my ($seq,$aa,$cut) =@_;
  my $start = 0;
  my $end = length($seq);
  if($aa>0){
    $start= ($aa-10);
    $end = $aa+10;
  }
  $start = 0 if($start <0);
  $end = length($seq) if($end>length$seq);
  my $sseq = substr($seq,$start, ($end-$start));
  my %fq;
  $fq{A} = $sseq =~ tr/A|a//;
  $fq{T} = $sseq =~ tr/T|t//;
  $fq{G} = $sseq =~ tr/G|g//;
  $fq{C} = $sseq =~ tr/C|c//;
  my @keys = keys %fq;
  @keys = sort{ $fq{$b} <=> $fq{$a} or $b cmp $a} @keys if(scalar @keys >1);
  
  my$rat = src::Utilities::round( ($fq{$keys[0]}*100)/length($sseq),2);
  #print $sseq,"\t",length $sseq,,"\t",$rat,"\n";

  my $out ="*";
  if($rat>=$cut){
     $out = "poly".$keys[0]."(".$rat."%)";
  }
  return($out);
}

sub getClonalClipPercent {
  my ($ty,$delim) = @_;
  $ty =~ s/$delim$//;
  $ty =~ s/$delim{2,}/$delim/g;  
  $ty = $ty.$delim;
  $ty =~ s/=/chr/g; $ty =~ s/\*/chr/g; $ty =~ s/chr/ /g;
  my @ty = split(" ",$ty);
  @ty = grep { $_ !~ $delim } @ty;
  @ty = grep { $_ ne ' ' } @ty;
  #print join("\n",@ty);
  my %ty; $ty{$_}++ for(@ty);
  my @keys = sort { $ty{$b} <=> $ty{$a} } keys %ty;
  my $RT=0;
  $RT = src::Utilities::round($ty{$keys[0]}/scalar(@ty),2) if(scalar @keys>0);          
  return($RT);
}

sub getCorrectClipSeq{
  my ($rr,$side)=@_;
  $rr =~ s/,$//;
  $rr =~ s/,{2,}/,/g;
  my @rr = split(",",$rr);
  if($side eq "a" or $side eq "b"){
    my %res;
    foreach my $jj(@rr){
      my $sq;
      $sq = substr($jj,length($jj)-10,10) if($side eq "a");
      $sq = substr($jj,0,10) if($side eq "b");
      $res{$jj}{A} = () = $sq=~ /A|a/g;
      $res{$jj}{T} = () = $sq=~ /T|t/g;
      $res{$jj}{F} = 0;
      if($res{$jj}{A}>$res{$jj}{T}){
        $res{$jj}{F}=$res{$jj}{A};
      }else{
        $res{$jj}{F}=$res{$jj}{T};
      }
    }
    my @keys = keys %res;
    @keys = sort{$res{$b}{F} <=> $res{$a}{F} or length($b) <=> length($a) or $b cmp $a} keys %res if(scalar @keys >1);
    return($keys[0]);
  }else{
    @rr = sort{length($b) <=> length($a) or $b cmp $a} @rr if(scalar @rr > 1);
    return($rr[0]);
  }
}

sub gridCheck_PolyA{
  my ($CLP,$side,$locus,$file,$x,$y) =@_;
  
  open FH,"samtools view $file '$locus'|" or next "Can't open bamfile $file";
  #open FH, $file or die $!;
  my @tmp = <FH>;
  close FH;

  my %out;
  foreach my $line (@tmp){
   my @sam = split("\t",$line);
   my ($sidedness) = $line =~ m/.*;side=(\S*);clip.*/;  
   my ($aa) = $sam[5]=~ /^(\d+)S\S+/;
   my ($ab) = $sam[5]=~ /\D(\d+)S$/;
   $aa=0 if(!defined($aa) || $aa eq"");  ## left hand side
   $ab=0 if(!defined($ab) || $ab eq"");  ## right hand side   
   my $id = $sam[2].$sam[3].$sam[5];
      
   if($aa>=10 and $side eq "a"){
     my @res;
     my $clp = src::Utilities::get_clip_coordV1($sam[3],"*",$sam[5],"lhs");
     if( $clp<= ($CLP+$y) and $clp>= ($CLP+$x) ){
       my $sq = substr($sam[9],0,$aa);
       $res[0]= () = $sq =~ /A|a/g;
       $res[0] = $res[0]/length($sq);
       $res[1]= () = $sq =~ /T|t/g;
       $res[1] = $res[1]/length($sq);
       $res[2]= () = $sq =~ m/A{7,}|a{7,}|T{7,}|t{7,}/g;
       $sq = substr($sq,length($sq)-10,10);
       $res[3]= () = $sq =~ /A|a/g;
       $res[4]= () = $sq =~ /T|t/g;
       $out{$clp}{ids} .=$id."," if($res[3]>=7 or $res[4]>=7 or $res[2]>0 or $res[0]>0.67 or $res[1]>0.67);
     }
   }elsif($ab>=10 and $side eq "b"){
     my @res;
     my $clp = src::Utilities::get_clip_coordV1($sam[3],"*",$sam[5],"rhs");
     if(  $clp<= ($CLP-$x) and $clp>= ($CLP-$y) ){
       my $sq = substr($sam[9],$ab*(-1));
       $res[0]= () = $sq =~ /A|a/g;
       $res[0] = $res[0]/length($sq);
       $res[1]= () = $sq =~ /T|t/g;
       $res[1] = $res[1]/length($sq);
       $res[2]= () = $sq =~ m/A{7,}|a{7,}|T{7,}|t{7,}/g;
       $sq = substr($sq,0,10);
       $res[3]= () = $sq =~ /A|a/g;
       $res[4]= () = $sq =~ /T|t/g;
       $out{$clp}{ids} .=$id."," if($res[3]>=7 or $res[4]>=7 or $res[2]>0 or $res[0]>0.67 or $res[1]>0.67);  
     }
   }
  }

  my $out = "n/a";
  if(scalar keys %out >0){
    my @keys = keys %out;
    @keys = sort{ abs($a - $CLP) <=> abs($b - $CLP) or src::Utilities::unique($out{$b}{ids},",","true") <=> src::Utilities::unique($out{$b}{ids},",","true"); 
              } @keys if(scalar @keys >1);
    $out = $keys[0]."|".src::Utilities::unique($out{$keys[0]}{ids},",","true");              
  }

  #print $out,"\n\n";
  #print Dumper %out;
  return($out);
}

1;


