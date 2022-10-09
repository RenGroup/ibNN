#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 4){
  die "Usage: perl $0 <input_type> <id_type> <row_or_column> <input>\nSee comments in the script for id_type";
}
#example: perl /Users/chix/neuralNetwork/BIB/reviewerComments/newData/0.convert_to_geneID.pl csv ENSG row expr_GSM3988006_Astrocyte.csv
#note: all of the quotes "\"" will be removed, since the output csv files from R may contain quotes

my ($input_type,$id_type,$direct,$input) = @ARGV;
#input_type is either "csv" or "tsv", meaning expression values are separated by "," or by tab ("\t")
#id_type is the type of gene's id in the input file. Options are:
#"ENSG", "NM","ENST","NP","ENSP","symbol","synonyms","sysy" ##"sysy" means using both symbols and synonyms
#$direct specifies if the first row or the first column contains the gene's identifier for convertion. Options: "row","column"
#the output is in csv format

#load the matching table according to $id_type
my $matching_table = "/Users/chix/neuralNetwork/BIB/reviewerComments/matching_table/geneID_match_to_otherIDs.tsv";
open(IN,$matching_table) or die "OpenError: $matching_table, $!\n";
my %match;
while(<IN>){
  $_=~s/[\n\r]//g;
  my @data = split(/\t/);
  if($id_type eq "ENSG"){
    if($data[2] ne "-"){
      $match{$data[2]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
  }elsif($id_type eq "NM"){
    if($data[3] ne "-"){
      $match{$data[3]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
  }elsif($id_type eq "ENST"){
    if($data[4] ne "-"){
      $match{$data[4]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
  }elsif($id_type eq "NP"){
    if($data[5] ne "-"){
      $match{$data[5]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
  }elsif($id_type eq "ENSP"){
    if($data[6] ne "-"){
      $match{$data[6]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
  }elsif($id_type eq "symbol"){
    if($data[7] ne "-"){
      $match{$data[7]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
  }elsif($id_type eq "synonyms"){
    if($data[8] ne "-"){
      my @tmp = split(/\|/,$data[8]); #multiple synonyms are separated by "|"
      foreach my $tmpGene (@tmp){
        $match{$tmpGene}=$data[1];
      }
    }
  }elsif($id_type eq "sysy"){
    if($data[7] ne "-"){
      $match{$data[7]}=$data[1]; #$data[1] is the ncbi's gene ID
    }
    if($data[8] ne "-"){
      my @tmp = split(/\|/,$data[8]); #multiple synonyms are separated by "|"
      foreach my $tmpGene (@tmp){
        $match{$tmpGene}=$data[1];
      }
    }
  }
}
close IN;
my $num = keys %match;
print "Loading match table finished. Selected input id type: $id_type.";
print "Loaded $num matching relations. Start converting.\n";

if($direct eq "row"){
  print "The first row will be converted\n";
  open(IN, $input) or die "OpenError: $input, $!\n";
  $num=0; #reset $num. Here it calculates the number of matched ids.
  my $line = <IN>; # the first line is stored in $line now
  $line=~s/[\n\r]//g;
  $line=~s/\"//g; #remove all the quotes "\""
  my @data;
  if($input_type eq "csv"){
    @data=split(/\,/,$line);
  }elsif($input_type eq "tsv"){
    @data=split(/\t/,$line);
  }
  for(my $i=1; $i<=$#data;$i++){ #loop for all the identifiers in @data; skip $i=0. the first element of @data should be the header of the cell id's column
    if(exists $match{$data[$i]}){
      $data[$i]=$match{$data[$i]};
      $num++;
    }else{
      $data[$i]="-";
    }
  }
  my $header = join(",",@data);
  print "Converted $num/$#data identifiers. Now output the results.\n";
  $input=~s/(\.csv$|\.tsv|.txt)//; #remove the suffix
  my $output="$input\_geneID.csv";
  open(OUT,">$output") or die "OpenError: $output, $!\n"; #didn't use "sed" since the syntaxes of sed in GNU and non-GNU (like macOS) are different
  print OUT "$header\n";
  while(<IN>){ #this start from the second line now
    $_=~s/\"//g; #remove all the quotes "\""
    print OUT "$_";
  }
  close IN;
  close OUT;
  print "Convertion finished\n";
}else{
  die "The \"column\" option is still in development\n";
}
