# The records of data processing for GSE135133
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
- Convert the gene identifiers into gene id:
```
for i in `ls expr_*`; do echo $i; perl 0.convert_to_geneID.pl csv ENSG row $i; done
```
- Now merge all the files of the same cell type:
```
perl 2.merge_samples.pl
```
This script will automatically search for all the files of the same cell type in the folder, keep one header and remove all the other headers, then merge them into one matrix starting with "merged_".
- Generate masked and testing data:
```
for i in `ls merged_expr_*`; do echo $i; python 3.generate_testingData.py /path_to_merged_files/ $i; done
```
This python script will do the following things:</br>
(1) It will remove the genes which failed in the conversion step</br>
(2) It will filter the cells by the number of genes per UMI, suggested by the QC [here](https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/seurat-qc-cell-level-filtering.html)</br>
(3) It will randomly mask one tenth of the non-zero values in one cell into zeros, producing the masked files. The original values of the masked genes will be stored in the testing file</br>
(4) It will create two new folders to hold the masked and testing files</br>
- Now we can do the training and imputation:
```
for i in ls `masked_oneTenth*`; do echo $i; python 2.5.2.train_impute_ibNN.py -d /path/to/ -i $i; done
```
This script will create three folders: ../log/log_...txt which records the information during training; ../imputed/imputed_...csv which is the imputed matrix; ../wma/who_...txt and ../wma/wih...txt, which are the weight matrices of input to hidden (wih) and hidden to output (who).

