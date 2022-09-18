#!/use/bin/perl
use strict;
use warnings;

my $sym = "Homo_sapiens.gene_info";
my $gene2ensembl = "hs_gene2ensembl.tsv"; #only do this for human
my $output = "geneID_match_to_otherIDs.tsv";

open(IN,$sym) or die "OpenError: $sym, $!\n";
my %hash;
while(<IN>){
  next if (/^\#/);
  $_=~s/[\n\r]//g;
  my @data =split(/\t/);
  $hash{$data[1]}{"symbol"}=$data[2];
  $hash{$data[1]}{"synonyms"}=$data[4]; # guess "-" for no synonyms
}
close IN;

open(IN,$gene2ensembl) or die "OpenError: $gene2ensembl, $!\n";
open(OUT,">$output") or die "OpenError: $output, $!\n";
while(<IN>){
  $_=~s/[\n\r]//g;
  my @data=split(/\t/);
  if(exists $hash{$data[1]}){
    print OUT "$_\t";
    if(exists $hash{$data[1]}{"symbol"}){
      print OUT "$hash{$data[1]}{'symbol'}\t";
    }else{
      print "Warning: synonyms exist but symbol does not\n";
      print OUT "-\t";
    }
    if(exists $hash{$data[1]}{"synonyms"}){
      print OUT "$hash{$data[1]}{'synonyms'}\n";
    }else{
      print OUT "\-\n";
    }
  }else{
    print OUT "$_\t\-\t\-\n"; #"-" for no data
  }
}
close IN;
close OUT;
