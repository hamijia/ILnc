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
#Load the td file under the data folder in the ILnc package.
data(td)
#lnc_sig：lnc_sig is the matrix of lncRNA signatures of 14 immune cells.
lnc_sig
#sig_label：sig_label is a set of 14 characteristic lncRNAs of immune cells.
sig_label
#Load the example_data file under the data folder in the ILnc package.
data(example_Data)
#ILnc_example：An example of an expression matrix for ILnc analysis.
ILnc_example
#Download the data and save it to your work path.
setwd("C:/Users/xumengjia/Desktop/aaa")
write.table(ILnc_example,"ILnc_example.txt",col.names = T,row.names = T,quote=F,sep = "\t")
ILnc("ILnc_example.txt",c("NK","Treg"))
```
# Data input
### ILnc(x,y)
| **Arguments** | **Detail** |
| --- | --- |
| **x** | The expression matrix should be a matrix with genes in rows and samples in columns. The rownames should be Ensembl ID.|
| **y** | 1 or more names of the 14 cell type:("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg").|
