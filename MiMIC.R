LaplacianFromFile<-function(W)
{
     n=dim(W)[1]
     Dg=rowSums(W)
     DH=diag(Dg^(-0.5))
     AN=DH%*%W%*%DH
     I=diag(n)
     L=I+AN
     return(list(L=L,W=W,Dg=Dg))
}
GraphLaplacian<-function(Dmat,rbfsg=1000,mod=1)
{
     show=5
     n=dim(Dmat)[1]
     W=matrix(0,n,n)
     for(i in 1:n)
     {
        for(j in 1:i)
        {
             W[i,j]=W[j,i]=exp(-sum((Dmat[i,]-Dmat[j,])^2)/(2*rbfsg*rbfsg))
        }
     }
     Wun=W
     Dg=rowSums(W)
     DH=diag(Dg^(-0.5))
     W=DH%*%W%*%DH
     I=diag(n)
     L=I+W
     return(list(L=L,W=Wun,Dg=Dg))
#	  L=WNormalize(W)
#	  return(list(L=W,W=W))
}
WNormalize<-function(W)
{
	n=nrow(W)
     for(i in 1:n)
     {
		sm=sum(W[i,])-W[i,i]
     	for(j in 1:n)
     		{
     			if(i!=j)
     				W[i,j]=W[i,j]/(2*sm)
     			else
     				W[i,j]=0.5
     		}
     }
     return(W)
}
rowNormalize<-function(Data)
{
	n=dim(Data)[1]
	for(i in 1:n)
	{	if(sum(Data[i,]^2)>0)
			Data[i,]=Data[i,]/sqrt(sum(Data[i,]^2))
	}
	return(Data)
}
negativeEntries<-function(Y)
{
	r=nrow(Y)
	c=ncol(Y)
	Yneg=matrix(0,r,c)
	for(i in 1:r)
	{
		for(j in 1:c)
		{
			if(Y[i,j]<0)
				Yneg[i,j]=(Y[i,j])
		}
	}
	return(Yneg)
}
check_Majority<-function(Y)
{
	r=nrow(Y)
	c=ncol(Y)
	for(j in 1:c)
	{
		countneg=length(which(Y[,j]<0))
		if(countneg>(r-countneg))
			Y[,j]=-Y[,j]

	}
	return(Y)
}
minimize_KmeansManifold<-function(L,Ut,lambda=0.01, eta=0.01)
{
	r=nrow(Ut)
	c=ncol(Ut)
#	cat("\nUt:\n")
#	print(Ut[1:5,])
	Utneg=negativeEntries(Ut)
#	cat("\nUtneg:\n")
#	print(Utneg[1:10,])
	Qt=L%*%Ut-lambda*Utneg
	one=matrix(rep(1,r),r,1)
	x=(1/r)*(Qt%*%t(Ut)%*%one)
#	cat("\n x:\n")
#	print(x[1:10,])
#	cat("\n dim x=",dim(x))
	x1=(x%*%t(one)+one%*%t(x))
	WY=t(Qt)%*%Ut+t(Ut)%*%Qt
	Omg=(1/4)*(WY-2*(t(Ut)%*%x1%*%Ut))
#	cat("\nOmg:\n")
#	print(Omg)
	Zt=Qt-2*Ut%*%Omg-x1%*%Ut
#	cat("\nUt:\n")
#	print(Ut[1:5,])
#	cat("\nZt:\n")
#	print(Zt[1:5,])
	Zt1=Ut+(eta*Zt)
#	cat("\nZt1:\n")
#	print(Zt1[1:5,])
	library("Matrix")
	A=t(Ut)%*%Zt1
	A1=Ut%*%A%*%t(Ut)
	B=Zt1%*%t(Ut)-Ut%*%t(Zt1)-2*A1
#	cat("\nDim B=",dim(B),"Dim A1=",dim(A1))
	Ut1=as.matrix(expm(B)%*%expm(A1)%*%Ut)
#	cat("\nUt1\n")
#	print(Ut1[1:5,])
	return(list(Ut1=Zt1,Qt=Qt))
}
maximize_SteifelManifold<-function(L,Ut,eta=0.01)
{
	n=nrow(Ut)
	c=ncol(Ut)
	Qt=L%*%Ut
	Zt=Qt-0.5*Ut%*%(t(Ut)%*%Qt+t(Qt)%*%Ut)
#	Zt=(diag(n)-Ut%*%t(Ut))%*%Qt
	Zt1=Ut+(eta*Zt)
	sv=svd(Zt1)
	Ut1R=sv$u[,1:c]%*%t(sv$v[,1:c])
	Ut1=Ut1R
	return(Ut1)
}
checkEigen<-function(prevEig,nextEig)
{
	nc=ncol(prevEig)
	for(i in 1:nc)
	{
		dorg=sum((prevEig[,i]-nextEig[,i])^2)
		dminus=sum((prevEig[,i]+nextEig[,i])^2)
		if(dminus<dorg)
			nextEig[,i]=-nextEig[,i]
	}
	return(nextEig)
}
ManifoldBestMinimize<-function(Data,K,rank=NULL,simFromFile=FALSE,mod=NULL,lambda=NULL,Beta=NULL,dfEps=NULL,etaval=NULL,delta_decrease=NULL)
{

	#Setting parameters
	if(is.null(rank))
    {
		startRk=K
		maxRk=20 	#Set the maximum rank r to vary upto
    }
	else
		startRk=maxRk=rank
	if(is.null(lambda))
		lambda=0.01
	if(is.null(Beta))	
	  	Beta=0.5
	if(is.null(dfEps))
		dfEps=0.04
	if(is.null(etaval))
  		etaval=0.05			
	if(is.null(delta_decrease))
		delta_decrease=1e-04
    if(is.null(mod))
    	mod=seq(1,M)
	maxiter=100
		
		cat("\nStarting Run")
		fOpt=Inf
		UbestOpt=NULL
		rkSeq=seq(startRk,maxRk)
		show=5
		nstart=100
		ShiLm=list()
		Datas=list()
		Rel=vector(mode="numeric",length=M)
		order=vector(mode="numeric",length=M)
		Dist=vector(mode="numeric",length=M)
		sigma=vector(mode="numeric",length=M)
		n=dim(Data[[1]])[1]
		for(rkFrac in rkSeq)
		{
			cat("\nOptimizing At Rank ",rkFrac)
			nsc=rkFrac
		    evind=seq(1,nsc)
		    for(m in 1:M)
		    {
			    if(simFromFile==TRUE)
			       Lp=LaplacianFromFile(Data[[m]])
			    else
			    {
			        Datas[[m]]=scale(Data[[m]],center=T,scale=F)
			        sgFrac=0.5
			        Dist[m]=max(as.numeric(dist(Datas[[m]])))
			        sigma[m]=sgFrac*Dist[m]
		            Lp=GraphLaplacian(Datas[[m]],rbfsg=sigma[m],mod=m)
		        }
		        Lm=Lp$L
		        evi=eigen(Lm)
		        evi$values[ abs(evi$values)<1e-10 ] <- 0
		        evi$values<-round(evi$values,digits=10)
		        ln=length(evi$values)
		        ind=2
		        evind=which(evi$values!=2)[1:nsc]
		        Um=evi$vectors[,evind,drop=FALSE]
		        Dm=evi$values[evind]
		        Dmat=evi$vectors[,ind,drop=FALSE]
		        km=kmeans(Dmat,2,iter.max=100,nstart=nstart)
		        Rel[m]=(evi$values[ind]*(1+silhouette(Dmat,km$cluster)))*0.25
		        cat("\n Modality: ",mod[m]," Relevance=",Rel[m])
			    ShiLm[[m]]=list(L=Lm,D=Dm,U=Um)
		    }			
		    order=sort(Rel,index.return=TRUE,decreasing=TRUE)$ix
			Best=order[1]
			Rest=setdiff(order,Best)
			if(nsc==K)
			{
				cat("\n Order=",modname[order],file=fileLbest,append=TRUE)
				cat("\n Best=",modname[Best],"\n Rest=",modname[Rest],file=fileLbest,append=TRUE)
			}
######################### Gradient Descent Manifold Optimization #################################################
			eta=etaval
			Ubest=ShiLm[[Best]]$U
			Lbest=ShiLm[[Best]]$L
			URestList=list()
			LRestList=list()
			URestListTemp=list()
			fprev=0
			for(i in 1:(M-1))
			{
				URestList[[i]]=ShiLm[[Rest[i]]]$U
				LRestList[[i]]=ShiLm[[Rest[i]]]$L
			}
			BestKernel=Ubest%*%t(Ubest)
		    SumKernel=matrix(0,n,n)
		    for(i in 1:(M-1))
		        SumKernel=SumKernel+URestList[[i]]%*%t(URestList[[i]])
		    Ubest=check_Majority(Ubest)
		    cat("\nStep_Length=",etaval,"\nKmeans-Manifold-Lambda=",lambda,"\nStepReduce_Beta=",Beta,"\nConvergence epsilon=",dfEps)
	
			cat("\nMiMIC-iteration-",0)
		    fclust=-sum(diag(t(Ubest)%*%(Lbest)%*%Ubest))
			Ubestneg=negativeEntries(Ubest)
			fclust=fclust+sum(Ubestneg^2)
			fclust=fclust/nsc
			fdag=0
			for(i in 1:(M-1))
				fdag=fdag-sum(diag(t(URestList[[i]])%*%(BestKernel)%*%URestList[[i]]))
			fdag=fdag/nsc
			f=fclust+fdag
			fprev=f
	    	cat("\tf=",f,"fclust=",fclust,"fdag=",fdag)
#### Manifold Optimization Initialization Ends Here
		    t=0
		    while(eta>1e-06)
		    {
		      	#K-means Manifold Optimization
		      	KmOpt=minimize_KmeansManifold(SumKernel+Lbest,Ubest,eta=eta,lambda=lambda)
				UbestTemp=KmOpt$Ut1
			    fclust=-sum(diag(t(UbestTemp)%*%(Lbest)%*%UbestTemp))
			    Ubestneg=negativeEntries(UbestTemp)
			    fclust=fclust+sum(Ubestneg^2)
			    fclust=fclust/nsc
				BestKernel=UbestTemp%*%t(UbestTemp)
				fdag=0
				for(i in 1:(M-1))
				{
					RestiPrev=URestList[[i]]
					URestListTemp[[i]]=maximize_SteifelManifold((BestKernel),URestList[[i]],eta=eta)
					fdag=fdag-sum(diag(t(URestListTemp[[i]])%*%((BestKernel))%*%URestListTemp[[i]]))
				}
				fdag=fdag/nsc
				f=fclust+fdag
				df=fprev-f
				Ca=fprev-f-delta_decrease*eta*sum(diag(t(KmOpt$Qt)%*%KmOpt$Qt))
				params=paste("\tf=",f,"fclust=",fclust,"fdag=",fdag,"Diff-f=",df,"Armijo criterion Ca=",Ca)
				if(df>=dfEps && Ca>=0)
				{
					t=t+1
					Ubest=UbestTemp
					for(i in 1:(M-1))
						URestList[[i]]=URestListTemp[[i]]
				  	cat("\nMiMIC-iteration-",t,params)
					SumKernel=matrix(0,n,n)
					for(i in 1:(M-1))
						SumKernel=SumKernel+URestList[[i]]%*%t(URestList[[i]])
					fprev=f
				}
				else
				{
					eta=Beta*eta
					cat("\nf=",f,"Difference in f=",df,"Objective not improving, reducing eta=",eta)
				}
    		}  #End While
            cat("\nAt rank r=",nsc,"  Objective f=",f)
            if(f<fOpt)
            {
            	UbestOpt=Ubest
            	fOpt=f
            	rStar=nsc
            }		
    	}
		cat("\n\nOptimal rank or Given Rank r*=",rStar)
		cat("\nSubspace Ubest* corresponding to Best modality written to file: UbestStar.txt")
		write.table(UbestOpt,row.names=FALSE,col.names=FALSE,quote=FALSE,file="UbestStar.txt")
		return(list(Ubest=UbestOpt,f=fOpt)) 
}




#Internal Indices
#Silhouette
sumsq<-function(v1,v2)
{
	eqd=sum((v1-v2)^2)
	return(eqd)
}
proximity<-function(Dmat)
{
	n=dim(Dmat)[1]
	proxim=matrix(0,n,n)
	for(i in 1:n)
	{
		for(j in 1:(i-1))
		{
			proxim[i,j]=sqrt(sumsq(Dmat[i,],Dmat[j,]))
			proxim[j,i]=proxim[i,j]
		}
	}
	return(proxim)
}
silhouette<-function(Dmat,kmclust)
{
    proxim=proximity(Dmat)
    K=max(kmclust)
	cl=seq(1:K)
	n=dim(Dmat)[1]
	s<-vector(mode="numeric",length=n)
	bvec<-vector(mode="numeric",length=K-1)
	for(i in 1:n)
	{
		own=kmclust[i]
		temp=proxim[i,]
		restown=setdiff(which(kmclust==own),i)
		if(length(restown)>0)
		{
			owndist=temp[restown]
			ai=mean(owndist)
		}
		else
			ai=0
		oc=setdiff(cl,own)
		for(j in 1:length(oc))
		{
			oth=oc[j]
			otherclust=which(kmclust==oth)
			othdist=temp[otherclust]
			bvec[j]=mean(othdist)
		}
		bi=min(bvec)
		s[i]=(bi-ai)/max(ai,bi)
	}
	clustsil=rep(0,K)
	for(i in 1:K)
	{
		temp=which(kmclust==i)
		clsil=s[temp]
		clustsil[i]=mean(clsil)
	}
    sil=mean(s)
	return(sil)
}
