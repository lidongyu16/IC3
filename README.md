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
lrinfo=read.table("~/obdataforIC3/lrinfo.txt",header=TRUE);
```

### Input of IC3
Input contains 3 parts: A, cellinfo and lrinfo.

A is the gene expression count matrix, each row represents a cell and each colume represents a gene.

cellinfo is the cell information matrix, each row represents a cell.The first column represents the x-axis coordinates of the cell, and the second column represents the y-axis coordinates of the cell. The third column represents the cell type of the cell. The row name is the name of the cell, and its order should be consistent with the order of gene expression matrix A;

lrinfo represents the ligand-receptor database, each row represents a ligand-receptor pair, ligand is in the first column.

Then we use IC3 function to construct the cell-cell communication with single-cell resolution.

```R
result=IC3(A, cellinfo, lrinfo, alpha = 0.01,maxitr = 5)
```
This line of code will take approximately 75 minutes.

### Output of IC3

Output result contains 5 parts. The first part is the communication probability matrix at cell type level. It is a T by T upper triangular matrix where T is the number of cell type. The second part is the communication probability matrix at single cell level. The last 3 parts are the parameter estimation of lambda, beta and r.

## Verification of communication results

### Benchmark 

For the communication between different cell types in the Olfactory Bulb and SVZ (Subventurarian Zone) of mouse brain, we conducted an extensive review of the literature, with data sourced from the citeDB database (link: https://github.com/shanny01/benchmark/) and PubMed. We have put the summarized results into benchmark file : https://github.com/lidongyu16/IC3/tree/master/IC3/data/obbenchmark.txt and https://github.com/lidongyu16/IC3/tree/master/IC3/data/svzbenchmark.txt 

The benchmark file contains 3 columes, "Pubmed ID" means the Pubmed ID of reference. "Type1" "Type2" means the two interacting cell types' name. 

Load the benchmark file after downloading the benchmark data to the local path:

```R
realob=read.table("~/obbenchmark.txt",header=TRUE) 
```

### Other Method Result
We substituted the same data into the other five methods (STANN, CellphoneDB, Giotto, SpaOTsc, COMMOT). We collect the results of these methods and organize them into  https://github.com/lidongyu16/IC3/tree/master/IC3/data/obOMresult.txt
and  https://github.com/lidongyu16/IC3/tree/master/IC3/data/obsvzresult.txt. 

The Other Method result contains 7 columes. First 2 colume represent 2 cell types, last 5 columes represent the communication strength estimated by 5 method. 

Load the Other Method result file after downloading data to the local path:

```R
OMdata=read.table("~/obOMresult.txt")
```

### Plot ROC curve


```R
typeresult=result[[1]];  ## get the IC3 result
typepairnum=dim(OMdata)[1];   
internum=dim(realob)[1];

IC3result=rep(0,typepairnum);   
typename=rownames(typeresult)
for (i in 1:typepairnum)
{
    typeindexone=which(typename==OMdata[i,1]);
    typeindextwo=which(typename==OMdata[i,2]);
    indexone=min(typeindexone,typeindextwo);
    indextwo=max(typeindexone,typeindextwo);
    IC3result[i]=typeresult[indexone,indextwo];     ## change IC3 result matrix to a line that consistent with OMdata
}
wzreal=rep(0,typepairnum);
for (i in 1:typepairnum)
{
   for (j in 1:internum)
   {
     if(realob[j,1]==OMdata[i,1] &&  realob[j,2]==OMdata[i,2])
      {
              wzreal[i]=1                         ## change benchmark matrix to a line that consistent with OMdata
      }
     if(realob[j,1]==OMdata[i,2] &&  realob[j,2]==OMdata[i,1])
      {
              wzreal[i]=1
      }
   }
}
library(pROC)
library(ggplot2)
stanndata=roc(wzreal,OMdata[,3],direction="<")
cellphonedata=roc(wzreal,OMdata[,4],direction="<")
giottodata=roc(wzreal,OMdata[,5],direction="<")
spaotscdata=roc(wzreal,OMdata[,6],direction="<")
commotdata=roc(wzreal,OMdata[,7],direction="<")
cellchatdata=roc(wzreal,OMdata[,8],direction="<")

IC3data=roc(wzreal,IC3result,direction="<")                     ## calculate the ROC between benchmark and all the methods.
library(stringr)
IC3str=paste("IC3 AUC =",round(IC3data[["auc"]],3),sep="")
STANNstr=paste("STANN AUC =",round(stanndata[["auc"]],3),sep="")
CellphoneDBstr=paste("CellphoneDB AUC =",round(cellphonedata[["auc"]],3),sep="")
Giottostr=paste("Giotto AUC =",round(giottodata[["auc"]],3),sep="")
SpaOTscstr=paste("SpaOTsc AUC =",round(spaotscdata[["auc"]],3),sep="")
COMMOTstr=paste("COMMOT AUC =",round(commotdata[["auc"]],3),sep="")
CellChatstr=paste("CellChat AUC =",round(cellchatdata[["auc"]],3),sep="")

p=ggroc(list(STANN=stanndata,CellphoneDB=cellphonedata,Giotto=giottodata,SpaOTsc=spaotscdata,COMMOT=commotdata,CellChat=cellchatdata,IC3=IC3data), legacy.axes = TRUE)

p=p+  annotate("text", x = 0.75, y = 0.5,label =IC3str,color="red")
p=p+  annotate("text", x = 0.75, y = 0.45,label =STANNstr)
p=p+  annotate("text", x = 0.75, y = 0.4,label =CellphoneDBstr)
p=p+  annotate("text", x = 0.75, y = 0.35,label =Giottostr)
p=p+  annotate("text", x = 0.75, y = 0.3,label =SpaOTscstr)
p=p+  annotate("text", x = 0.75, y = 0.25,label =COMMOTstr)
p=p+  annotate("text", x = 0.75, y = 0.2,label =CellChatstr)
p <- p + scale_color_manual(values = c("IC3" = "red", "STANN" = "blue", "CellphoneDB" = "green", "Giotto" = "purple", "SpaOTsc" = "orange", "COMMOT" = "yellow","CellChat"="black"))
p    
```

So we can get the following plot:


![image](https://github.com/lidongyu16/IC3/blob/master/IC3/data/obPRcurve.png)

## Find communication hotspot with scan statistic using IC3 cell-level result

Our results include not only cell type level communication, but also cell level. For all cells of two cell types, we can analyze whether the communications between the two cell types are clustered through the cell-level communication network. For this purpose we consider all possible cell pairs consisting of two cell types, and select the midpoint as a representative for each cell pair. Input three columns of information matrix: the midpoint horizontal axis coordinate, the vmidpoint ertical axis coordinate and the communication status of the cell pair. We can use the IC3 function scan_p_value to get the hotspot and significance level of the communication hotspot between the two cell types.

We use the interaction between pericytes and endothelial cells in the MERFISH data mouse 1 slice 5 as an example. First, we download the example data from https://github.com/lidongyu16/IC3/tree/master/IC3/data/scanexampledata.

Load the data after downloading the data to the local path:

```R
set.seed(1)
A=read.table("~/obdataforIC3/scanexampledata.txt",header=TRUE);
result <- scan_p_value(A, num_permutations = 1000)
```
The running process takes three minutes. The result is divided into five parts. The first two parts are the scan statistics and the corresponding p value of this data. The third and fourth parts are the center and side length of the corresponding hot spot area. The last part is the statistic result obtained by shuffling 1000 times.

```R
print(result$observed_statistic)
print(result$p_value)
print(result$best_center)
print(result$best_size)
```

We can get the result:

7.448827;0.002;-2793.993 -4050.203;125.3125

To visualize the results, we can first plot the null distribution of the scan statistic and the location of the sample data using density estimation.


```R
library(ggplot2)
data <- data.frame(null_statistics = result$null_statistics)
ggplot(data, aes(x = null_statistics)) +
  geom_density(fill = "lightblue", alpha = 0.5) +  
  geom_vline(aes(xintercept = result$observed_statistic), color = "red", linetype = "dashed", size = 1) +  
  annotate("text", x = result$observed_statistic, y = 0.02, label = paste0("Observed\nScan\nStatistic: ", round(result$observed_statistic, 2)),
           color = "red", angle = 90, vjust = -0.5, hjust = -0.2) +  
  ggtitle("Null Distribution of Scan Statistic") +
  theme_minimal() +
  labs(x = "Statistic", y = "Density") +
  annotate("text", x = result$observed_statistic -2 ,y = 0.2, label = paste0("p-value: ", round(result$p_value, 4)),
           color = "black", hjust = -0.2, vjust = -0.5)  
```

![image](https://github.com/lidongyu16/IC3/blob/master/IC3/data/scanstatnull.png)
We can also use ggplot to draw a dot plot of the data and the location of the hotspots:

```R
A_df <- data.frame(x = A[,1], y = A[,2], communication_status = as.factor(A[,3]))
center_x <- result$best_center[1];center_y <- result$best_center[2]
size <- result$best_size / 2  
bottom_right_x <- center_x + size;bottom_right_y <- center_y - size
top_right_x <- center_x + size;top_right_y <- center_y + size;bottom_left_x <- center_x - size;bottom_left_y <- center_y - size
top_left_x <- center_x - size;top_left_y <- center_y + size;offset <- size
label_x <- top_right_x + offset;label_y <- top_right_y + offset
ggplot(A_df, aes(x = x, y = y)) +
  geom_point(aes(color = communication_status), size = 0.2)+
  scale_color_manual(name = "Communication Status", 
                     values = c("1" = "purple", "0" = "#FFD700"), 
                     labels = c("1" = "Yes", "0" = "No")) +   theme_minimal() +
  ggtitle("Scatter Plot of Points with Highlighted Hotspot Region") +
  annotate("segment", x = center_x - size, xend = center_x + size, y = center_y - size, yend = center_y - size, color = "red") + 
  annotate("segment", x = center_x - size, xend = center_x + size, y = center_y + size, yend = center_y + size, color = "red") +  
  annotate("segment", x = center_x - size, xend = center_x - size, y = center_y - size, yend = center_y + size, color = "red") + 
  annotate("segment", x = center_x + size, xend = center_x + size, y = center_y - size, yend = center_y + size, color = "red") + 
  geom_text(aes(x = label_x, y = label_y), 
            label = paste0("p-value: ", round(result$p_value, 4)), color = "blue", size = 5) +   coord_fixed()  
```
Then we get:

![image](https://github.com/lidongyu16/IC3/blob/master/IC3/data/exampleregion.png)

