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
minimize_KmeansManifold<-function(L,Ut,lambda=0.01, alpha=0.01)
{
	r=nrow(Ut)
	c=ncol(Ut)
	Utneg=negativeEntries(Ut)
	Qt=L%*%Ut-lambda*Utneg
	one=matrix(rep(1,r),r,1)
	x=(1/r)*(Qt%*%t(Ut)%*%one)
	x1=(x%*%t(one)+one%*%t(x))
	WY=t(Qt)%*%Ut+t(Ut)%*%Qt
	Omg=(1/4)*(WY-2*(t(Ut)%*%x1%*%Ut))
	Zt=Qt-2*Ut%*%Omg-x1%*%Ut
	Zt1=Ut+(alpha*Zt)
	library("Matrix")
	A=t(Ut)%*%Zt1
	A1=Ut%*%A%*%t(Ut)
	B=Zt1%*%t(Ut)-Ut%*%t(Zt1)-2*A1
	Ut1R=as.matrix(expm(B)%*%expm(A1)%*%Ut)
	return(list(Ut1=Zt1,Qt=Qt))
}
maximize_SteifelManifold<-function(L,Ut,alpha=0.01)
{
	n=nrow(Ut)
	c=ncol(Ut)
	Qt=L%*%Ut
	# Zt=Qt-0.5*Ut%*%(t(Ut)%*%Qt+t(Qt)%*%Ut)
	Zt=(diag(n)-Ut%*%t(Ut))%*%Qt
	Zt1=Ut+(alpha*Zt)
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
ManifoldJointMinimize<-function(Data,K,rank=NULL,modname="GENE",simFromFile=0)
{
    set.seed(2)
	err=1e-03
	maxiter=100
	lambda=0.01
  	Beta=0.5
  	dfEps=0.005          #Set dfEps=0.001 for benchmark data sets, dfEps=0.005 for Omics data sets
  	alphaval=0.05
  	sigma_decrease=1e-05
    	dampF=2              #Set dampF=1 for benchmark data sets, dampF=2 for Omics data sets
    	startRk=K
    	maxRk=50
	if(is.null(mod))
    		mod=seq(1,M)
    	cat("\n Step length=",alphaval)
    	M=length(Data)
    	nstart=100
        silOpt=-Inf
        UjointOpt=NULL
		rkSeq=seq(startRk,maxRk,1)
		show=5
		ShiLm=list()
		Datas=list()
		Rel=vector(mode="numeric",length=M)
		Alpha=vector(mode="numeric",length=M)
		order=vector(mode="numeric",length=M)
		Dist=vector(mode="numeric",length=M)
		sigma=vector(mode="numeric",length=M)
		n=dim(Data[[1]])[1]
		for(rkFrac in rkSeq)
		{
			cat("\n\n",DataSet,"At Rank ",rkFrac)
			nsc=rkFrac
            dK=min(K,nsc)
		    for(m in 1:M)
		    {
			    if(simFromFile==1)
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
#		        cat("\n Modality ",m," lambda2=",evi$values[2],"sil=",silhouette(Dmat,km$cluster))
		        Rel[m]=(evi$values[ind]*(1+silhouette(Dmat,km$cluster)))*0.25
			    Lmr=Um%*%diag(Dm)%*%t(Um)
			    ShiLm[[m]]=list(L=Lm,U=Um,Lr=Lmr)
		    }
		    tstart=proc.time()[3]
		    order=sort(Rel,index.return=TRUE,decreasing=TRUE)$ix
            alphaW=rep(0,M)
	    	for(m in 1:M)
	        	alphaW[order[m]]=Rel[order[m]]*(1/(dampF^m))
	    	Alpha=alphaW/sum(alphaW)
            if(rkFrac == rkSeq[1])
            {
                cat("\nRelevance=",Rel)
                cat("\norder=",order)
                cat("\nAlpha=",Alpha)
            }
	    	LJoint=matrix(0,n,n)
            for(m in 1:M)
                LJoint=LJoint+Alpha[m]*ShiLm[[m]]$Lr
			Rest=seq(1,M)
			Best=order[1]
######################### Gradient Descent Manifold Optimization #################################################
			alpha=alphaval
            eigLJoint=eigen(LJoint)
            # cat("\n eigen vals of LJoint= ",eigLJoint$values[1:nsc])
            UJoint=eigLJoint$vectors[,1:nsc,drop=FALSE]
			URestList=list()
			LRestList=list()
			URestListTemp=list()
			fprev=0
			for(i in 1:M)
			{
				URestList[[i]]=ShiLm[[Rest[i]]]$U
				LRestList[[i]]=ShiLm[[Rest[i]]]$Lr
			}
			JointKernel=UJoint%*%t(UJoint)
		    SumKernel=matrix(0,n,n)
		    for(i in 1:M)
		        SumKernel=SumKernel+URestList[[i]]%*%t(URestList[[i]])
		    UJoint=check_Majority(UJoint)
		    fclust=-sum(diag(t(UJoint)%*%(LJoint)%*%UJoint))
			UJointneg=negativeEntries(UJoint)/2
			fclust=fclust+lambda*sum(UJointneg^2)
			fclust=fclust/(2*nsc)
			fdag=0
			for(i in 1:M)
				fdag=fdag-sum(diag(t(URestList[[i]])%*%(LRestList[[i]]+JointKernel)%*%URestList[[i]]))
			fdag=fdag/(2*nsc*M)
			f=fclust+fdag
			fprev=f
#### Manifold Optimization Initialization Ends Here
		    t=0
		    while(alpha>1e-03)
		    {
                cat("\niteration ",t, "f=",f)
		      	#K-means Manifold Optimization
		      	KmOpt=minimize_KmeansManifold(SumKernel+LJoint,UJoint,alpha=alpha,lambda=lambda)
				UJointTemp=KmOpt$Ut1
			    fclust=-sum(diag(t(UJointTemp)%*%(LJoint)%*%UJointTemp))
			    UJointneg=negativeEntries(UJointTemp)/2
			    fclust=fclust+lambda*sum(UJointneg^2)
			    fclust=fclust/(2*nsc)
				JointKernel=UJointTemp%*%t(UJointTemp)
				fdag=0
				for(i in 1:M)
				{
                    URestListTemp[[i]]=maximize_SteifelManifold((LRestList[[i]]+JointKernel),URestList[[i]],alpha=alpha)
					fdag=fdag-sum(diag(t(URestListTemp[[i]])%*%((LRestList[[i]]+JointKernel))%*%URestListTemp[[i]]))

				}
				fdag=fdag/(2*nsc*M)
				f=fclust+fdag
				df=fprev-f
				Ca=fprev-f-sigma_decrease*alpha*sum(diag(t(KmOpt$Qt)%*%KmOpt$Qt))

				if(df>=dfEps && Ca>=0)
				{
					t=t+1
					UJoint=UJointTemp
					for(i in 1:M)
						URestList[[i]]=URestListTemp[[i]]
					SumKernel=matrix(0,n,n)
					for(i in 1:M)
						SumKernel=SumKernel+URestList[[i]]%*%t(URestList[[i]])
					fprev=f
				}
				else
				{
					alpha=Beta*alpha
					cat("\nf=",f,"Objective not improving significantly, reducing step length to",alpha)
				}
    		}  #End While
            km=kmeans(UJoint[,1:dK],K,iter.max=100,nstart=nstart)$cluster
            sil=silhouette(UJoint[,1:dK],km)
            cat("\nAt rank r=",nsc,"  Objective f=",f)
            if(sil>silOpt)
            {
                UjointOpt=UJoint
                silOpt=sil
                rStar=nsc
            }
		} #Rank Loop End
        cat("\n\nOptimal rank or Given Rank r*=",rStar)
        cat("\nOptimal Subspace Ujoint* corresponding to Joint view written to file: UjointStar.txt")
        write.table(UjointOpt,row.names=FALSE,col.names=FALSE,quote=FALSE,file="UjointStar.txt")
        return(list(Ujoint=UjointOpt))
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
