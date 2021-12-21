# ILnc
ILnc is an R-package based on the individual identification of lncRNA signature sets of 14 immune cell types to estimate its infiltration.
# Install
ILnc is a part of R/Bioconductor and is availble on Linux, macOS and Windows platforms.<br>
The latest version of ILnc can be installed from GitHub using devtools package, which can take up to a few minutes to install all the dependencies:<br>
```r
library(devtools)
install_github("hamijia/ILnc/ILnc")
```
# Usage
```r
library(ILnc)
data(td)
lnc_sig
sig_label
exprMatrix = read.table(mixture_file,header=TRUE,row.names=1, as.is=TRUE)
ILnc(exprMatrix,c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg"))
```
#### Load data
lnc_sig:<br>
sig_label:


# Data input
| **Arguments** | **Detail** |
| --- | --- |
| **x** | The expression matrix should be a matrix with genes in rows and samples in columns. The rownames should be Ensembl ID.|
| **y** | 1 or more names of the 14 cell type:("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg").|
