# ibNN
Interpretable bionic neural network

## General introduction</br>
ibNN is a simple fully connected neural network with only three-layers (including the input layer). However, its structure simulates human signaling and gene regulatory network. By incorporating existing biological knowledge, it can learn the molecular relations from single cell/nuclei RNA-seq data. The core of the network is built upon the conversion of the adjacency matrix of the directed graph of signaling network and TF-target relations to the initial weight matrices of ibNN. Therefore, each node in the network has been assigned an explicit name of the genes, and the trained weight matrices can be converted back to the relations between molecules, which provides clear intepretation of the meaning of the trained network.

## Before you start
ibNN is specially designed with several assumptions, so please check whether ibNN suits your data before going on:</br>
- The input of ibNN should be the raw counts or other values without log-transformation (raw counts are prefered)
- The input of ibNN should be a group of cells with clear cell types/subtypes. This is because ibNN learns the inner mechanisms of signaling and transcription factor(TF)-target regulations. Input of multiple cell types may confuses ibNN since its weight matrices can represent only one set of mechanisms.
- The current version of ibNN does not read the whole input file into memory, therefore is not memory intensive. It can process large datasets (n_cells > 5000), even with 400k cells, on a laptop (it took about 27 hours, tested on M2 macBook Air).
- The output of ibNN is not limited to the imputation results. The weight matrices are also informative since they represent the molecular interactions in signaling and gene-regulatory networks. Protocols are still under development, but examples will be available in our paper.

## Contents
- [Installation](https://github.com/RenGroup/ibNN/blob/main/README.md#installation)
- [Dependencies](https://github.com/RenGroup/ibNN/blob/main/README.md#dependencies)
- [The format of input file](https://github.com/RenGroup/ibNN/blob/main/README.md#the-format-of-input-file)
- [Run and optimize parameters](https://github.com/RenGroup/ibNN/blob/main/README.md#run-and-optimize-parameters)
- [The outputs](https://github.com/RenGroup/ibNN/blob/main/README.md#the-outputs)
- [Limitations](https://github.com/RenGroup/ibNN/blob/main/README.md#limitations)
- [Design of ibNN](https://github.com/RenGroup/ibNN/blob/main/README.md#design-of-ibnn)

## Installation</br>
Users can simply download the python script https://raw.githubusercontent.com/RenGroup/ibNN/main/ibNN_main/2.5.train_impute_ibNN.py to the local computer. 
```
wget https://raw.githubusercontent.com/RenGroup/ibNN/main/ibNN_main/2.5.train_impute_ibNN.py
```
Then download the PPrel and TF-target matrix file: [wMa.pprel.txt](https://github.com/RenGroup/ibNN/blob/main/initial_weight_matrices/wMa.pprel.txt.zip) and the four split files of wMa.tf.txt: [wMa.tf.txt.01.zip](https://github.com/RenGroup/ibNN/blob/main/initial_weight_matrices/wMa.tf.txt.01.zip),[wMa.tf.txt.02.zip](https://github.com/RenGroup/ibNN/blob/main/initial_weight_matrices/wMa.tf.txt.02.zip), [wMa.tf.txt.03.zip](https://github.com/RenGroup/ibNN/blob/main/initial_weight_matrices/wMa.tf.txt.03.zip), and [wMa.tf.txt.04.zip](https://github.com/RenGroup/ibNN/blob/main/initial_weight_matrices/wMa.tf.txt.04.zip) </br>
After unzip the files, users should concatenate the files in the correct order to restore the file of wMa.tf.txt:
```
cat wMa.tf.txt.01 wMa.tf.txt.02 wMa.tf.txt.03 wMa.tf.txt.04 > wMa.tf.txt
```
Before running ibNN, there are still two things to do:</br>
Modify *line 13* and *line 14* in 2.3.train_impute_ibNN.py, change *"/path_to/"* to the dir of where you put the matrix files. Then check the dependencies:
## Dependencies</br>
ibNN was built upon commonly-used packages, and tested on both intel- and M1/M2-based macOS. If observed error messages or any other unexpected behavior, please first check the versions of the packages.
```
Python 3.8.5
numpy 1.19.5
scipy 1.5.0
re 2.2.1
argparse 1.1
pandas 1.0.5
```
If the errors persist, please report the bugs to chix@big.ac.cn.</br>
If users want to check the version of these packages (assume they are already installed), the following codes may be helpful:
```
import numpy
import scipy.special
import re
import argparse
import pandas as pd #for batch output
import subprocess

cmd = "python3 --version"
cmdOut = subprocess.check_output(cmd,shell=True).decode("utf-8")
cmdOut = re.sub("[\n\r]","",cmdOut)
print(cmdOut)
print("numpy "+numpy.__version__)
print("scipy "+scipy.__version__)
print("re "+re.__version__)
print("argparse "+argparse.__version__)
print("pandas "+pd.__version__)
```
Copy and paste the above codes into jupyter notebook and directly run, or into an empty python script named for example, "check_version.py", and run in cmd:
```
python3 check_version.py
```
## The format of input file</br>
### General format requirements
ibNN expects single cell/nuclei RNA-seq data file in csv format(data separated by ","), ending in ".csv" or ".txt" suffix. The data should be the raw counts of UMIs, or other values which have not been log-transformed. Each row should be a cell, and each column should be a gene. The gene identifier should be converted to NCBI's gene ID. We understand that ID conversion is always problematic with multi-mapping issues, so we made extra [tips](https://github.com/RenGroup/ibNN/blob/main/id_conversion/README_idConversion.md) for converting the ID using the most updated official ID mapping files. Users can build the ID mapping files by their own following the instructions, and then write their own script to format the input data file, or use our scripts if the data format matches our examples.</br>
Run the following codes in cmd to check the 10 lines, 5 columns of the example input data:
```
cut -f1-5 -d ',' masked_oneTenth_merged_expr_Vascular_geneID.csv |head
```
We'll get the following output:
```
,1,10,100,1000
CELL3898054,0,0,0,1
CELL3898448,0,0,1,4
CELL3898560,0,0,0,1
CELL3899207,0,0,0,0
CELL3899320,0,0,0,0
CELL3899527,0,0,0,0
CELL3899866,0,0,0,0
CELL3900024,0,0,0,1
CELL3900050,0,0,0,1
```
Note that the first line is started with a ",". 
### Transpose the expression matrix
If your data is in the transposed format, you may consider "datamash" to transpose the large dataset with low memory and time consumption. "datamash" can be installed like this:
```
#download source code
curl -0 http://ftp.gnu.org/gnu/datamash/datamash-1.8.tar.gz --output datamash-1.8.tar.gz
#unzip and make
tar -xzf datamash-1.8.tar.gz
cd datamash-1.8
./configure
make
make check
sudo make install
datamash --help
```
To transpose a csv file, this example code may be helpful:
```
cat Cells.UMI.collapsed.csv|datamash transpose -t ',' > tr_Cells.UMI.collapsed.csv
```
### Handling large dataset in .loom format
Loom format is an efficient way of accessing large datasets without loading the entire matrix into the memory. For large datasets like GSE156793 which contains more than four million cells, only a .loom file was provided to store the count matrix. Currently ibNN is not capable of directly access the loom files, therefore cannot benefit from the conveniences of loom format, but we provide alternative solutions to handle the data, at the expenses of time consumption. We will implement the ability of accessing loom files in ibNN in future.
</br>
To access loom files in R, one can use SeuratDisk:
```
library(SeuratDisk)
sc.793.loom <- Connect(filename = "/path_to/GSE156793_S3_gene_count.loom", mode = "r")
sc.793.loom[["matrix"]][1:5,1:5]
```
The output is like this:
</br>
![image](https://user-images.githubusercontent.com/109563761/190891770-3b4b8278-5a94-4400-a523-ac7dbba2c94c.png)
</br>
Unlike python package "Dask", SeuratDisk allows positional indexing.
</br>
To output the expression data of desired cells with large numbers, we recommend the R package "vroom". For example, if you want to export the entire matrix of "GSE156793_S3_gene_count.loom", you may try the following R code:
```
install.packages("vroom")
library(vroom)
n=0
for(i in seq(1,4065000,15000)){
  n=n+1
  end.this<-i+14999
  tmp.1<-sc.793.loom[["matrix"]][i:end.this,]
  tmp.1<-as.data.frame(tmp.1)
  vroom_write(tmp.1, paste0("/path_to_output/file",n,".csv.gz"), 
              delim = ",", quote = "none", col_names = F, num_threads = 6)
  rm(list="tmp.1")
}
```

## Run and optimize parameters</br>
Two parameters are required by ibNN to specify the path to the input .csv file (-d), and the file name of the .csv file (-i). Additional parameters can be found using:
```
python3 2.3.train_impute_ibNN.py -h
```
The result of this cmd:
```
usage: 2.3.train_impute_ibNN.py [-h] [-d path_to_csv] [-i input_file_name] [-e Max_num_training] [-l Learning_rate] [-c Num_cell_for_training] [-b Num_cell_in_one_batch]

Train ibNN and impute gene expressions

optional arguments:
  -h, --help            show this help message and exit
  -d path_to_csv, --dir path_to_csv
                        The path to the input csv.
  -i input_file_name, --filename input_file_name
                        The name of the input csv file.
  -e Max_num_training, --nRounds Max_num_training
                        Training will end if flag turns to 3 or rounds of training reach this number
  -l Learning_rate, --learningRate Learning_rate
                        set by experiences.
  -c Num_cell_for_training, --nTrain Num_cell_for_training
                        Num of cells for training, default '-1' for auto detect (2/3 of total cell number, may take time to count)
  -b Num_cell_in_one_batch, --batchSize Num_cell_in_one_batch
                        The mean of the errors of the cells in one batch will be used for back propagation
```
A quick start of the script is like this:
```
python3 2.3.train_impute_ibNN.py -d /path_to_file/ -i masked_oneTenth_merged_expr_Vascular_geneID.csv
```
ibNN has been optimized for small (<100 cells) or large (>5000 cells) datasets. It will check the number of cells if "-c" option (number of cells for training) is not specified. Automatic adjustment of the parameters may be triggered when:</br>
- **For small dataset**  
If (n_train (number of cells for training) < 100) & (n_rounds (max number of epochs) < 40), then n_rounds will be adjusted to 40. The reason for this adjustment is that when the number of cell used for training is low, there's higher possibility that the internal-calculated MSE (i.e. MSE calculated using the 1/3 cells that were not used in the training step) will be higher than normal (> 0.4), which indicate the failure of learning the general features of the data but instead, ibNN were over-fitted to the trained cells.</br>
- **For large dataset**  
If (n_train > 5000) & (n_batch == 1) (i.e. the number of cells for training is more than 5000, and the number of cells used in each time of training is set to 1), then n_batch will be reset to int(n_train/1000). In this way, when n_train is > 5000, ibNN will only be trained for 1000 times in each round (i.e. epoch). For each time ibNN is trained, n_batch cells will be used to query the network and calculate the loss (i.e. target - output), then the mean of the loss of each gene was calculated, then propagated back to the network to update the weight matrices. We tested mean, median and max, and found that "mean" gave the best MSE (the lowest) for all the datasets we tested. Users are welcome to test the methods. Simply change "numpy.mean" at line 204 to other functions.
### Notes about the parameters</br>
#### The -c and -b option for large dataset
The default -c option takes 2/3 of the total number of cells for training, the rest 1/3 for testing (i.e. calculating the MSE). The script will print the median MSE for each round (i.e. epoch). From our experiences, when n_training is larger than 500, the benefit of large number of cells for training drops quickly. This is why we control the batch size to be n_train/1000; each batch only updates the weight matrices once, therefore 1000 batches will update the weight matrices for 1000 times, just like when n_train is 1000 and n_batch is 1. 
#### Further optimization for small dataset
For small groups of cells, the MSE is usually higher. To lower down the MSE, we can extend the rounds (i.e. epochs) of training, just like what the auto-adjustment do. However, if the MSE is still high, users may try modifying an inner parameter at line 340 named "epochs_inner" to further increase the number of updates. This parameter is originally 1, meaning each cell is used once before training with another cell's data. Increasing to 2 or 3 may increase the number of updates by two or three times, but the MSE may only drop marginally. This is usually not recommended since it may result in the over-fitting of the model.
#### The learning rate
The default is 0.02, which has been tested and performed best for most of the datasets so far. Other tested learning rates are: 0.1, 0.01, 0.05, and decay. We did't provide an option of sophisticated decaying method of learning rate, which will be implemented in future versions.
#### When to end the training
The training will be stopped if the max rounds of training (the -e option) is reached or ibNN detected the following situation:

> abs(median(MSE<sub>thisRound</sub>) - median(MSE<sub>previousRound</sub>))/median(MSE<sub>previousRound</sub>) < 0.001

and 
> median(MSE<sub>thisRound</sub>) < 0.4

for two sequential rounds. When running ibNN, the value of a "flag" will be printed, which is usually 1. If met the above situation, the value of flag will be added by 1; if not, flag will be set to 1 again. The training will be stopped when flag == 3.
#### When the training failed
When the training is done, the median MSE will be checked. If the median MSE is > 0.4, then the imputation is risky since the model might not learned the general regulation rules, and ibNN will stop. Usually this happens when the n_train is small. If n_train > 100 and MSE is > 0.4, then the input cells might be too heterogenous (consists of two or more cell types) or the structure of ibNN does not suit the cell type. For the former one, users may try to use sub-groups of the cells to reduce the heterogeneity. For the later one, we are developing more comprehensive and reasonable neural network structures to incorporate gene expression regulatory mechanisms in addition to current signaling and TF-target networks.

## The outputs</br>
The outputs of ibNN are consisted of: an imputation result file started with "imputed_", two weight matrix files start with "wih_" (weight matrix of input to hidden) and "who_" (weight matrix of hidden to output), and one log file start with "log_". The messages printed to screen contains where to find these files.
## Limitations</br>
There are several limitations of ibNN:</br>
- ibNN needs input cells to be of the same cell type, or with similar gene regulatory mechanisms. This means unless the scRNA-seq library is prepared with clear cell-surface markers, ibNN has to wait for the cell identity assignments to be done.
- ibNN only incorporates two of the mechanisms that control gene expressions, however, many other well-known mechanisms have been demonstrated to greatly affect gene expressions, i.e. epigenetics. Although we have successfully trained ibNN for many cell types from several tissues, we cannot rule out the possibility that under certain conditions, there are cell types whose gene expressions are not reflected by cell signaling networks.
- ibNN currently is only able to process human data. Although KEGG also collected PPrel data for other animals such as mouse and monkey, the TF-target relations seems to be also important, since using data other than GeneHancer's estimations failed in mouse and monkey (GeneHancer provides cis-element estimations only for human). The training of mouse and monkey data showed signs of learning (the MSE drops), but failed in the accuracy control (MSE > 2)
- ibNN's code have not taken the advantage of GPU processing. It can be further optimized to increase the speed of training.

## Design of ibNN</br>
### The core of ibNN</br>
The core of ibNN is the conversion of the directed graph of biological networks to the initial weight matrices of neural network. The original idea of the conversion is as the following figure:
<img width="1273" alt="adjacencyMatrixConvertion" src="https://user-images.githubusercontent.com/109563761/189895896-2aee0246-b5b4-49f5-99da-e72b0e2a000a.png">

**Figure 1.** The conversion of graphical knowledges of regulation to weight matrices of neural network. (A) A simple directed graph of three proteins. Protein A activates protein B, and protein B activates protein C. This type of graph is often seen in the signaling cascade, e.g. the PPrel of KEGG database. (B) The adjacency matrix of A. (C) A two-layer directed graph rearranged from A. (D) The adjacency matrix of the directed graph C with three vertices in a general format. (E) The layout of the input and middle layers of ibNN. For demonstration purpose, only three genes: A, B and C were shown. (F) The formula to calculate the input signal to the middle layer. The W<sub>input_middle</sub> was derived from the transpose of the adjacency matrix in (D).

### Building the global networks of human intra-cellular signaling and GRN</br>
One of the important gene expression regulations is the transcription factor(TF)-target regulations. The TFs are in turn, regulated through various mechanisms, one of which is the signaling networks. The signaling cascades converge on TFs, which then regulate gene expressions[ref]. To construct the basic structure of ibNN, we first compiled the signaling networks from KEGG's protein-protein relations (PPrel) data, then the GRNs by merging the TF-binding sites from GTRD and the estimations of cis-element to neighbor genes from GeneHancer database. The adjacency matrices were properly adjusted (assigned normally-distributed small numbers to the zeros and scaled according to the conventions of neural network construction). The final matrices which could be loaded into the neural network can be found in the folder "initial_weight_matrices".

to be continued...
