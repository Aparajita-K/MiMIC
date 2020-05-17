rm(list=ls())
DataSet="LGG"  	#Data set name (if any)
n=267			#Number of samples in the data set
nfeat=2000		#Maximum number of features to consider from a modality
K=3				#Number of clusters in the data set


#Extract sample names
truedata=read.table(paste0("DataSets/",DataSet,"/TCGA Subtype Labels",n),stringsAsFactors=FALSE)
true=as.vector(truedata[,2],mode="numeric")
samples=as.vector(truedata[,1])
#Load the data into a list
Data<-list()

# DNA Methylation and Gene Expression have very large number of features, so we take the first 2000 features
#The features in all modalities are sorted in decreasing order of variance
Data[[1]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/mDNA",n), sep=" ",header=TRUE,row.names=1))[,1:nfeat]
Data[[2]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/RNA",n), sep=" ",header=TRUE,row.names=1))[,1:nfeat]

#miRNA and Protein modalities have relatively fewer features (<1000), so we take all features
Data[[3]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/miRNA",n), sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/RPPA",n), sep=" ",header=TRUE,row.names=1))

M=length(Data)  #The number of modalities
modalities=c("mDNA","RNA","miRNA","RPPA")

#Log Transformation of sequence based RNA and miRNA modality
LogData=Data
#Log Transform Sequence based Gene(RNA) expression and miRNA expression modalities
#Replace the 0 counts by 1 before log transformation
#Skip this step if not required for your data set

LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)

#Pass data set to the data integration algorithm CoALa
source("MiMICjoint.R")
Algo="MiMIC"
mimic=ManifoldJointMinimize(Data=LogData,K=K,rank=NULL,mod=modalities)


#Perform K-means clustering on joint subspace
UjointSub=as.matrix(read.table("UjointStar.txt",sep=" ",header=FALSE))
cat("\n First few rows of Ujoint* subspace:\n")
print(UjointSub[1:5,1:K])
cat("\n Subspace Dimension: ",dim(UjointSub)[1]," rows",dim(UjointSub)[2]," columns")
cat("\nClustering on First k columns")
UjointSubK=UjointSub[,1:K]
km=kmeans(UjointSubK,K)$cluster
df=data.frame(cbind(samples,km))
write.table(df,quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste0(DataSet,"-ClusterAssignment.txt"))
cat("\n\nFinal cluster assignments written to file:",paste0(DataSet,"-ClusterAssignment.txt\n\n"))
