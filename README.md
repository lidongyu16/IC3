# IC3
Infer Cell-Cell Communication with single cell resolution

## Installation

IC3 is a R package that inferring cell-cell communication based on the spatial transcriptomic data. It can be installed by following:

```R
install.packages("devtools")
library(devtools)
install_github("lidongyu16/IC3/IC3")
library(IC3)
```

## Tutorial

### Import packages

```R
install.packages("progress")
library(progress)
library(stringr)
```
### Download the data

Download a spatial transcriptome dataset (e.g. seqFISH+ olfactory bulb dataset) and ligand-receptor database(e.g. CellTalkDB). The original spatial transcritome data is from https://github.com/CaiGroup/seqFISH-PLUS. The ligand-receptor database is from https://github.com/ZJUFanLab/CellTalkDB

We extracted example data from seqFISH+ and CellTalkDB and placed it in the following folder https://github.com/lidongyu16/IC3/tree/master/IC3/data/obdataforIC3.

Load the data after downloading the data to the local path:

```R
A=read.table("~/obdataforIC3/A.txt",header=TRUE);
cellinfo=read.table("~/obdataforIC3/cellinfo.txt",header=TRUE);
lrinfo=read.table(~/obdataforIC3/lrinfo.txt",header=TRUE);
```

### Input of the main function IC3

A is the gene expression count matrix, each row represents a cell and each colume represents a gene.

cellinfo is the cell information matrix, each row represents a cell.The first column represents the x-axis coordinates of the cell, and the second column represents the y-axis coordinates of the cell. The third column represents the cell type of the cell. The row name is the name of the cell, and its order should be consistent with the order of gene expression matrix A;

lrinfo represents the ligand-receptor database, each row represents a ligand-receptor pair, ligand is in the first column.

Then we use IC3 function to construct the cell-cell communication with single-cell resolution.

```R
result=IC3(A,cellinfo,lrinfo)
```
This line of code will take approximately 40 minutes.

### Output of the main function IC3

result contains 5 parts. The first part is the communication probability matrix at cell type level. It is a T by T upper triangular matrix where T is the number of cell type. The second part is the communication probability matrix at single cell level. The last 3 parts are the parameter estimation of lambda, beta and r.

## Verification of results

### Benchmark 

For the communication between different cell types in the Olfactory Bulb and SVZ (Subventurarian Zone) of mouse brain, we conducted an extensive review of the literature, with data sourced from the citeDB database (link: https://github.com/shanny01/benchmark/) and PubMed. We have put the summarized results into:


We substituted the same data into the other five methods and obtained the interaction matrix at the cell type level. We compared the results of different methods with the results of existing papers to obtain the ROC curve. The specific process is as follows 


```R
typeresult=result[[1]];
cellresult=result[[2]];
lambda=result[[3]]
beta=result[[4]];
r=result[[5]]

alldata=read.table("/home/lidy/obOMresult.txt")
typepairnum=dim(alldata)[1]
IC3result=rep(0,typepairnum);
typename=rownames(typeresult)
for (i in 1:typepairnum)
{
    typeindexone=which(typename==alldata[i,1]);
    typeindextwo=which(typename==alldata[i,2]);
    indexone=min(typeindexone,typeindextwo);
    indextwo=max(typeindexone,typeindextwo);
    IC3result[i]=typeresult[indexone,indextwo]
}
```
