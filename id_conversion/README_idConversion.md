# Tips for ID conversion
## Step 1. Generate matching table
ID conversion is always a problematic process with issues of multi-mapping or no-matched terms. There are already number of online tools and packages to deal with this.
Here we present our own solution for the conversion between gene symbols, NCBI gene IDs, and ensembl gene ids.</br>
First, we download the most updated id match table from official resources: NCBI's [gene info](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)
and [gene2ensembl](https://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz). Unzip the two files in the same folder, then run the following in the same folder:
```
grep "^9606" gene2ensembl > hs_gene2ensembl.tsv
```
This will generate a matching file with ensembl ID. Then copy the perlscript "1.append_symbol_synonyms.pl" in the same folder and run:
```
perl 1.append_symbol_synonyms.pl
```
The script will automatically look for the file "hs_gene2ensembl.tsv" and "Homo_sapiens.gene_info" and generate a matching table "geneID_match_to_otherIDs.tsv". 
This file can be downloaded from this folder, however, as NCBI update these files very frequently, users can get the most updated matching table following the above tips

## Step 2. Convert IDs to NCBI's gene IDs
Edit line 19 of the perl script 2.convert_to_geneID.pl where users should change $matching_table's value to their own path to the file "geneID_match_to_otherIDs.tsv" or the
downloaded file "geneID_match_to_otherIDs_220904.tsv". Then run
```
perl 2.convert_to_geneID.pl
```
to check the help info. Example code:
```
perl /path_to/2.convert_to_geneID.pl csv ENSG row expr_GSM3988006_Astrocyte.csv
```
Note that current script can only convert the first row's ID since our data are all arranged in the way that each row is a cell, and each column is a gene
