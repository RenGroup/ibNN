#!/usr/bin/perl
use strict;
use warnings;

my @output = ("Amacrine","Astrocyte","Bipolar","Cones","Horizontal",
              "Muller","Myeloid","RGC","Rods","RPE","Vascular");
for(my $j=0;$j<=$#output;$j++){
  print "Start processing $output[$j]\n";
  my $filelist = `ls expr_*$output[$j]\_geneID.csv`;
  #print "$filelist";
  my @filearray=split(/\n/,$filelist);
  my $n=@filearray;
  for(my $i=1;$i<=$#filearray;$i++){ #start from the second file; the first file keeps the header
    system("awk '(NR>1)' $filearray[$i] >tmp_$filearray[$i]");
  }
  my $firstFile = $filearray[0];
  $filelist = `ls tmp_*`;
  @filearray=split(/\n/,$filelist);
  my $files=join(" ",@filearray); #all the file names, including the first file
  my $cmd = "cat $firstFile $files > merged_expr_$output[$j]\_geneID.csv";
  print "The cmd will be: $cmd\n";
  print "Merging started for $output[$j]\n";
  system("$cmd");
  print "Removing temporary files\n";
  system("rm $files");
  print "Done\n";

}
