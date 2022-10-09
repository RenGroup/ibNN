# The records of processing data from GSE135133
- Processed data (raw counts) "GSE135133_RAW.tar" can be downloaded [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135133) together with a file "GSE135133_clusterAssignments.txt.gz" containing the meta-data of each cell. Unzip the two files. 
- Use vi or vim to add a "\t" at the start of the first line for all of the 11 files</br>
```
vi GSM3988006_SAM24362284.txt
```
- Transpose the matrices:</br>
```
for i in `ls GSM39880*`; do echo $i;cat $i|datamash transpose > tr_$i; done
```
- Now, each row is a cell. Split the matrices by cell type: (do this for all GSM* files)
```
perl 1.split_cellType.pl tr_GSM3988016_SAM24362294.txt GSM3988016
```
- 
