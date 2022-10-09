#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 2){
  die "Usage: perl $0 <expr_matrix.txt> <outPrefix>\n";
}

my $celltype = "/path_to/GSE135133_clusterAssignments.txt"; #the meta data file
my $exprMa = $ARGV[0];
my $outprefix = $ARGV[1];
my @output = ("Amacrine","Astrocyte","Bipolar","Cones","Horizontal",
              "Muller","Myeloid","RGC","Rods","RPE","Vascular"); #these were got by analysis in R
my %match;

open(IN,$celltype) or die "OpenError: $celltype, $!\n";
my $n=0;
while(<IN>){
  $_=~s/[\n\r]//g;
  if($n ==0){
    $n++;
    next;
  }
  my @data=split(/\t/);
  $match{$data[0]}=$data[4]; #the column "Cell_type_full"
}
close IN;

open(IN,$exprMa) or die "OpenError: $exprMa, $!\n";
open(OUT0,">expr_$outprefix\_$output[0]\.csv")  or die "OpenError: expr_$outprefix\_$output[0]\.csv, $!\n";
open(OUT1,">expr_$outprefix\_$output[1]\.csv")  or die "OpenError: expr_$outprefix\_$output[1]\.csv, $!\n";
open(OUT2,">expr_$outprefix\_$output[2]\.csv")  or die "OpenError: expr_$outprefix\_$output[2]\.csv, $!\n";
open(OUT3,">expr_$outprefix\_$output[3]\.csv")  or die "OpenError: expr_$outprefix\_$output[3]\.csv, $!\n";
open(OUT4,">expr_$outprefix\_$output[4]\.csv")  or die "OpenError: expr_$outprefix\_$output[4]\.csv, $!\n";
open(OUT5,">expr_$outprefix\_$output[5]\.csv")  or die "OpenError: expr_$outprefix\_$output[5]\.csv, $!\n";
open(OUT6,">expr_$outprefix\_$output[6]\.csv")  or die "OpenError: expr_$outprefix\_$output[6]\.csv, $!\n";
open(OUT7,">expr_$outprefix\_$output[7]\.csv")  or die "OpenError: expr_$outprefix\_$output[7]\.csv, $!\n";
open(OUT8,">expr_$outprefix\_$output[8]\.csv")  or die "OpenError: expr_$outprefix\_$output[8]\.csv, $!\n";
open(OUT9,">expr_$outprefix\_$output[9]\.csv")  or die "OpenError: expr_$outprefix\_$output[9]\.csv, $!\n";
open(OUT10,">expr_$outprefix\_$output[10]\.csv")  or die "OpenError: expr_$outprefix\_$output[10]\.csv, $!\n";

$n=1;
while(<IN>){
  $_=~s/[\n\r]//g;
  if($n == 1){ #print the header
    my @data=split(/\t/);
    my $out = join(",",@data);
    print OUT0 "$out\n";
    print OUT1 "$out\n";
    print OUT2 "$out\n";
    print OUT3 "$out\n";
    print OUT4 "$out\n";
    print OUT5 "$out\n";
    print OUT6 "$out\n";
    print OUT7 "$out\n";
    print OUT8 "$out\n";
    print OUT9 "$out\n";
    print OUT10 "$out\n";
    $n++;
    next;
  }
  my @data=split(/\t/);
  my $out=join(",",@data);
  if(exists $match{$data[0]}){
    if($match{$data[0]} eq $output[0]){
      print OUT0 "$out\n";
    }elsif($match{$data[0]} eq $output[1]){
      print OUT1 "$out\n";
    }elsif($match{$data[0]} eq $output[2]){
      print OUT2 "$out\n";
    }elsif($match{$data[0]} eq $output[3]){
      print OUT3 "$out\n";
    }elsif($match{$data[0]} eq $output[4]){
      print OUT4 "$out\n";
    }elsif($match{$data[0]} eq $output[5]){
      print OUT5 "$out\n";
    }elsif($match{$data[0]} eq $output[6]){
      print OUT6 "$out\n";
    }elsif($match{$data[0]} eq $output[7]){
      print OUT7 "$out\n";
    }elsif($match{$data[0]} eq $output[8]){
      print OUT8 "$out\n";
    }elsif($match{$data[0]} eq $output[9]){
      print OUT9 "$out\n";
    }elsif($match{$data[0]} eq $output[10]){
      print OUT10 "$out\n";
    }elsif($match{$data[0]} eq "NA"){ #skip the cells annotated as "NA"
      next;
    }else{
      print "Warning: Unexpected cell type: $match{$data[0]}\n";
    }
  }else{
    print "Warning: unexpected cell id: $data[0]\n";
  }
}
close IN;
close OUT0;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;
