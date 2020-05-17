# MiMIC
Multi-Manifold Integrative Clustering

The MiMIC algorithm permforms integrative clustering on high-dimensional multimodal data sets. A multimodal data set consists of ``M`` modalities X<sub>1</sub>, ..., X<sub>m</sub>, ..., X<sub>M</sub>. Each modality X<sub>m</sub> represents the observations for same set of ``n`` samples from the ``m``-th data source.

Inorder to execute the R code for the lower grade glioma (LGG) data set,  within the ``R`` environment execute:
>source("LGGdemo.R")



Optimal low-rank joint subspace corresponding is written to file : ``UjointStar.txt``   
UjointStar.txt contains a ``(n x r)`` matrix.   
Here ``n`` is the number of samples in the data set and ``r`` is the optimal/required rank of the joint subspace.   

``k``-means clustering can be performed on the rows of UbestStar matrix to get the clusters in the data set.   
The cluster assignments are written to the file ``LGG-ClusterAssignment.txt`` for the LGG data set.  

The file ``MiMICjoint.R`` contains the ``R`` implementation of the MiMIC algorithm as a function `ManifoldBestMinimize`. 
Details of the fuctions is as follows:

Function Name: `ManifoldBestMinimize`

###### #Usage 
`ManifoldJointMinimize<-function(Data,K,rank=NULL,modname="RNA",simFromFile=0)
`


Arguments
``Data``:  A list object containing ``M`` data matrices representing ``M`` different omic data types measured in a set of ``n`` samples.    
For each matrix, the rows represent samples, and the columns represent genomic features.
The matrices in the list can have variable numbe of columns(features), but all must have the same number of *n* rows(samples).

``K``: The number of clusters in the data set.

``rank``: The rank of the individual and joint Laplacian. 
Default value: ``NULL``.
if ``rank=NULL``, the algorithm varies the rank between ``K`` to 50 and selects the optimal rank of the subspace.

``simFromFile``: Boolean value (``TRUE``/``FALSE``)   
if `FALSE`, algorithm considers the matrices in the `Data` list as feature-based representation, i.e., as `(n x d_m)` data matrices, and computes the graph Laplacian from data matrices.   
if ``TRUE``, algorithm considers the matrices in the `Data` list as graph-based representation, i.e., as `(n x n)` similarity matrices, and computes the graph Laplacian from similarity matrices.
Default value: ``FALSE``.

`mod`: Array containing names of modalities
Default value: `NULL`.
if `mod=NULL`, the algorithm names the modalities as  `1,2,...M`.




# Example call:

```r
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/LGG/mDNA",n),sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/LGG/RNA",n),sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/LGG/miRNA",n),sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/LGG/RPPA",n),sep=" ",header=TRUE,row.names=1))
K=3

M=length(Data)  #The number of modalities
modalities=c("mDNA","RNA","miRNA","RPPA")

#Log Transformation of sequence based RNA and miRNA modality
LogData=Data
#Log Transform Sequence based Gene(RNA) expression and miRNA expression modalities
#Replace the 0 conuts by 1 before log transformation
#Skip this step if not required for your data set

LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)

#Pass data set to the data integration algorithm CoALa
source("MiMICjoint.R")
Algo="MiMIC"
mimic=ManifoldJointMinimize(Data=LogData,K=K,rank=6,mod=modalities)


#Perform K-means clustering on joint subspace
UbestSub=as.matrix(read.table("UjointStar.txt",sep=" ",header=FALSE))
cat("\n First few rows of Ujoint* subspace:\n")
print(UbestSub[1:5,1:K])
cat("\n Subspace Dimension: ",dim(UbestSub)[1]," rows",dim(UbestSub)[2]," columns")
cat("\nClustering on First k columns")
UbestSubK=UbestSub[,1:K]
km=kmeans(UbestSubK,K)$cluster
df=data.frame(cbind(samples,km))
write.table(df,quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste0(DataSet,"-ClusterAssignment.txt"))
cat("\n\nFinal cluster assignments written to file:",paste0(DataSet,"-ClusterAssignment.txt\n\n"))
```
