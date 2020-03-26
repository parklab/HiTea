#!/usr/bin/perl
# HiTEA
use strict;
use warnings;

our %redb;
## MboI
$redb{"MboI"}{motif} = "GATC|gatc";
$redb{"MboI"}{motifa} = "GATC|gatc";
$redb{"MboI"}{ligmotif} = "GATCGATC|CTAGCTAG|gatcgatc|ctagctag";
$redb{"MboI"}{mlen} = "4";
$redb{"MboI"}{ligmotifln} = "8";
## DpnII
$redb{"DpnII"}{motif} = "GATC|gatc";
$redb{"DpnII"}{motifa} = "GATC|gatc";
$redb{"DpnII"}{ligmotif} = "GATCGATC|CTAGCTAG|gatcgatc|ctagctag";
$redb{"DpnII"}{mlen} = "4";
$redb{"DpnII"}{ligmotifln} = "8";
## HindIII
$redb{"HindIII"}{motif} = "[A]{0,1}AGCTT|AAGCT[T]{0,1}|[a]{0,1}agctt|aagct[t]{0,1}";
$redb{"HindIII"}{ligmotif} = "AAGCTAGCTT|TTCGATCGAA|aagctagctt|ttcgatcgaa";
$redb{"HindIII"}{mlen} = "6";
$redb{"HindIII"}{ligmotifln} = "10";
## Arima
$redb{"Arima"}{motif} = "GATC|GA[ATGC]{1}T[C]{0,1}|[G]{0,1}A[ATGC]{1}TC|gatc|ga[atgc]{1}t[c]{0,1}|[g]{0,1}a[atgc]{1}tc";
$redb{"Arima"}{ligmotif} = "GATCGATC|CTAGCTAG|GA[ATGC]{1}TA[ATGC]{1}TC|GATCA[ATGC]{1}TC|GA[ATGC]{1}TGATC|CT[GATC]{1}AT[AGTC]{1}AG|CT[ATGC]{1}ACTAG|CTAGT[ATGC]{1}AG|gatcgatc|ctagctag|ga[atgc]{1}ta[atgc]{1}tc|gatca[atgc]{1}tc|ga[atgc]{1}tgatc|ct[gatc]{1}at[agtc]{1}ag|ct[atgc]{1}actag|ctagt[atgc]{1}ag";
$redb{"Arima"}{mlen} = "4";
$redb{"Arima"}{ligmotifln} = "8";
## NcoI
$redb{"NcoI"}{motif} = "[C]{0,1}CATGG|CCATG[G]{0,1}|[c]{0,1}catgg|ccatg[g]{0,1}";
$redb{"NcoI"}{ligmotif} = "CCATGCATGG|GGTACGTACC|ccatgcatgg|ggtacgtacc";
$redb{"NcoI"}{mlen} = "6";
$redb{"NcoI"}{ligmotifln} = "10";

our (%chrs);
%chrs =(
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


our $PVAL_CUTOFF = 0.05;
our $RAMCountCUTOFF = 3;
our $PURITY_CUTOFF = 0.25;

