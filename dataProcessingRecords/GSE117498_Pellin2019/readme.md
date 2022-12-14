# Processing of the scRNA-seq data from Pellin et al. 2019
* Unzip the raw count files of HSC and MPP, then use the codes in preprocess_PellinData.R to convert the gene symbols to gene ids, and remove the redundant gene ids. Remember to modify the file path accordingly.
* The output files should then be: exprMa_hsc_rmDup_pellin_github.txt and exprMa_mpp_rmDup_pellin_github.txt
* Follow the codes in the Jupyter notebook files to train ibNN and perform imputations.
### Notes
The codes in the Jupyter notebook files are previous versions of ibNN. The core algorithms are the same with current version. The imputation result might be slightly different from the paper, since the id conversion was done by R package biomaRt, which uses online resources that may have been updated.
