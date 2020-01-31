# MiMIC
Manifold based algorithm for integrative clustering of multi-omics data.


Inorder to execute the R code for the **Kidney Carcinoma (KIDNEY) data set **,  within the R environment execute:
>source("KIDNEYdemo.R")



Low-rank subspace corresponding to the best modality Xbest is written to file : *UbestStar.txt*
UbestStar.txt contains a *(n x r)* matrix.
Here *n* is the number of samples in the data set and *r* is the optimal/required rank of the joint subspace.


*k*-means clustering can be performed on the rows of UbestStar matrix to get the clusters in the data set. 
The cluster assignments are written to the file * KIDNEY-ClusterAssignment.txt* for the KIDNEY data set.

The file *MiMIC.R* contains the R implementation of the MiMIC algorithm as a function `ManifoldBestMinimize`. 
Details of the fuctions is as follows:

Function Name: `ManifoldBestMinimize`

###### #Usage 
`ManifoldBestMinimize(Data,K,rank=NULL,simFromFile=FALSE,mod=NULL,lambda=NULL,Beta=NULL,dfEps=NULL,etaval=NULL,delta_decrease=NULL)
`


Arguments
*Data*:  A list object containing *M* data matrices representing *M* different omic data types measured in a set of *n* samples. 
For each matrix, the rows represent samples, and the columns represent genomic features.
The matrices in the list can have variable numbe of columns(features), but all must have the same number of *n* rows(samples).

*K*: The number of clusters in the data set.

*rank*: The rank of the individual and joint Laplacian. 
Default value: NULL.
if *rank=NULL*, the algorithm varies the rank between *K* to 20 and selects the optimal rank of the subspace.
M is the number of modalities in the data set.

*etaval*: The step length for optimizaing over *k*-Means and Stiefel manifolds.
Default etaval: NULL.
if *etaval=NULL*, the algorithm sets etaval to *0.5*, as stated in paper.

*Beta*: Value of Damping factor for the step length.
Default value: NULL.
if *Beta=NULL*, value of beta is set to *0.5* as in paper.

*lambda*: Value of lambda parameter for reducing the negative entries in Ubest while optimizaing over the *k*-Means Manifold.
Default value: NULL.
if *lambda=NULL*, value of beta is set to *0.01* as in paper.

*dfEps*: Value of convergence parameter eplison for convergence criterion of MiMIC algorithm.
Default value: NULL.
if *dfEps=NULL*, value of beta is set to *0.06* as in paper.

*delta_decrease*: Value of Armijo parameter for choice of step length of MiMIC algorithm.
Default value: NULL.
if *delta_decrease=NULL*, value of beta is set to *1e-04* as in paper.

*simFromFile*: Boolean value (TRUE/FALSE)
if FALSE, algorithm considers the matrices in the 'Data' list as feature-based representation, i.e., as *(n x d_m)* data matrices,
and computes the graph Laplacian from data matrices.
if TRUE, algorithm considers the matrices in the 'Data' list as graph-based representation, i.e., as *(n x n)* similarity matrices,
and computes the graph Laplacian from similarity matrices.
Default value: FALSE.

*mod*: Array containing names of modalities
Default value: NULL.
if mod=NULL, the algorithm names the modalities as* 1,2,...M*.




# Example call:

```r
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/KIDNEY/mDNA",n),sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/KIDNEY/RNA",n),sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/KIDNEY/miRNA",n),sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/KIDNEY/RPPA",n),sep=" ",header=TRUE,row.names=1))
K=3

#Log Transform of Sequence based Gene and miRNA modality
LogData=Data
LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)


source("MiMIC.R")
ManifoldBestMinimize(Data=LogData,K=K,rank=6,mod=modalities)
```
