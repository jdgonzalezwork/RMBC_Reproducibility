compute_estimator_and_performance_measures=function(seed,method,case,alphaTrim=0.05){
  set.seed(seed)
  returnDataGeratingProcces=dataGeneratingProcess(sem=seed,case=case)
  Y=returnDataGeratingProcces$Y
  indicesDGP=returnDataGeratingProcces$indicesDGP
  actualMu=returnDataGeratingProcces$actualMu
  actualSigma=returnDataGeratingProcces$actualSigma
  actualAlpha=returnDataGeratingProcces$actualAlpha
  alphaCont=returnDataGeratingProcces$alphaCont
  sem=returnDataGeratingProcces$sem
  case=returnDataGeratingProcces$case
  actualk=returnDataGeratingProcces$actualk
  n=nrow(Y)
  p=ncol(Y)
  if(method=="RMBC"){
    salRMBC= RMBC(Y,actualk)
    thetaNew.mu=salRMBC$mu;
    thetaNew.sigma=salRMBC$Sigma;  
    thetaNew.alpha=salRMBC$alpha
    estimated_clusters=salRMBC$cluster
    estimated_clusters[salRMBC$outliers]=0
  }
  
  if(method=="otrimle"){
    salot=otrimle(Y,G=actualk,ncores=1)
    #salot=otrimle(Y,G=actualk)
    tamOtrimle=salot$pi[2:(actualk+1)]
    alphaOtrimle=tamOtrimle/sum(tamOtrimle)
    thetaNew.alpha=rep(0,actualk); 
    thetaNew.mu=vector(mode="list", length=actualk)    
    thetaNew.sigma=vector(mode="list", length=actualk)
    for (j in 1:actualk){
      thetaNew.alpha[[j]]=alphaOtrimle[j]
      thetaNew.mu[[j]]=salot$mean[,j]
      thetaNew.sigma[[j]]=salot$cov[,,j]
    }
    estimated_clusters=salot$cluster
  }
  
  
  ###################################################
  # Comparaciones con otros:  Caso Clasico. mclust
  #Chris Fraley , # Adrian E. Raftery , Scruca Et al. 
  ######################################################
  if(method=="mclust"){
    mod1= densityMclust(Y,G=actualk,plot = F)
    thetaNew.alpha=rep(0,actualk); 
    thetaNew.mu=vector(mode="list", length=actualk)    
    thetaNew.sigma=vector(mode="list", length=actualk)
    for (j in 1:actualk){
      thetaNew.alpha[[j]]=mod1$parameters$pro[j]
      thetaNew.mu[[j]]=mod1$parameters$mean[,j]
      thetaNew.sigma[[j]]=mod1$parameters$variance$sigma[,,j]
    }
    estimated_clusters=mod1$classification
    
  }
  
  
  
  ###################################################
  if(method=="tclust"){
    saltclust=tclust(Y,k=actualk,alpha = alphaTrim)
    thetaNew.alpha=rep(0,actualk); 
    thetaNew.mu=vector(mode="list", length=actualk)    
    thetaNew.sigma=vector(mode="list", length=actualk)
    casoEspecial=!(saltclust$size==actualk); 
    for (j in 1:actualk){
      # some times tclust gives empty groups,
      if (j<=length(saltclust$size)){
        #thetaNew.alpha[[j]]=saltclust$size[j]/sum(saltclust$size)
        thetaNew.alpha[[j]]=saltclust$weights[j]
        thetaNew.mu[[j]]=saltclust$centers[,j]
        thetaNew.sigma[[j]]=saltclust$cov[,,j]
      }
      # some times tclust gives empty groups,
      if (j>length(saltclust$size)){
        thetaNew.alpha[[j]]=0
        thetaNew.mu[[j]]=saltclust$centers[,1]*0 
        thetaNew.sigma[[j]]=saltclust$cov[,,1]
      }
    }
    estimated_clusters=saltclust$cluster
    
    # the trimmed set is assigned to the cluster that maximizes
    # the probability according the computed mixture parameters.
    trimmedset= estimated_clusters==0
    YTrimmed=Y[trimmedset,]
    # There is a conflict if Y[trimmedset,] is not a matrix (when Size of trimmedset == 1.) 
    if( sum(trimmedset)==1){
      YTrimmed=matrix(YTrimmed,nrow=sum(trimmedset))
    }
    #test
    #plot(Y);points(Y[trimmedset, ],pch=19,col="red")
    #points(Y[estimated_clusters==3,], col=2,pch=19)
    #points(Y[estimated_clusters==1,], col=4,pch=19)
    #points(Y[estimated_clusters==2,], col=5,pch=19)
    
    valoresDisc2=RMBC::quad_disc(Y=YTrimmed, theta.alpha=thetaNew.alpha, 
                                 theta.mu=thetaNew.mu,theta.sigma=thetaNew.sigma)
    
    indices2=apply(valoresDisc2,1,function(x) which(x==max(x))[1])
    estimated_clusters[trimmedset]=indices2
    estan2=is_in_gr(Y=Y, cutoff=1-1e-3,theta.mu=thetaNew.mu,
                    theta.sigma=thetaNew.sigma)
    estimated_clusters[!estan2]=0
    
    
  }
  
  ## ADDING ORACLE FOR tclust
  
  if(method=="tclustoracle"){
    alphaTrim=alphaCont
    saltclust=tclust(Y,k=actualk,alpha = alphaTrim)
    thetaNew.alpha=rep(0,actualk); 
    thetaNew.mu=vector(mode="list", length=actualk)    
    thetaNew.sigma=vector(mode="list", length=actualk)
    casoEspecial=!(saltclust$size==actualk); 
    for (j in 1:actualk){
      # some times tclust gives empty groups,
      if (j<=length(saltclust$size)){
        #thetaNew.alpha[[j]]=saltclust$size[j]/sum(saltclust$size)
        thetaNew.alpha[[j]]=saltclust$weights[j]
        thetaNew.mu[[j]]=saltclust$centers[,j]
        thetaNew.sigma[[j]]=saltclust$cov[,,j]
      }
      # some times tclust gives empty groups,
      if (j>length(saltclust$size)){
        thetaNew.alpha[[j]]=0
        thetaNew.mu[[j]]=saltclust$centers[,1]*0 
        thetaNew.sigma[[j]]=saltclust$cov[,,1]
      }
    }
    
    estimated_clusters=saltclust$cluster
    
    # the trimmed set is assigned to the cluster that maximizes
    # the probability according the computed mixture parameters.
    trimmedset= estimated_clusters==0
    YTrimmed=Y[trimmedset,]
    # There is a conflict if Y[trimmedset,] is not a matrix (when Size of trimmedset == 1.)
    if( sum(trimmedset)==1){
      YTrimmed=matrix(YTrimmed,nrow=sum(trimmedset))
    }
    #test
    #plot(Y);points(Y[trimmedset,1 ],Y[trimmedset,2 ],pch=19,col="red")
    #points(Y[estimated_clusters==3,], col=2,pch=19)
    #points(Y[estimated_clusters==1,], col=4,pch=19)
    #points(Y[estimated_clusters==2,], col=5,pch=19)
    
    valoresDisc2=RMBC::quad_disc(Y=YTrimmed, theta.alpha=thetaNew.alpha, 
                                 theta.mu=thetaNew.mu,theta.sigma=thetaNew.sigma)

    indices2=apply(valoresDisc2,1,function(x) which(x==max(x))[1])
    estimated_clusters[trimmedset]=indices2
    estan2=is_in_gr(Y=Y, cutoff=1-1e-3,theta.mu=thetaNew.mu,
                    theta.sigma=thetaNew.sigma)
    estimated_clusters[!estan2]=0
    }
  
  ## Performance measures computation

  ret1=performance_measures(Y,indicesDGP,thetaNew.alpha,
                            thetaNew.mu,thetaNew.sigma,actualAlpha,
                            actualSigma,actualMu,estimated_clusters,
                            nkl=1000,cutoff=1-1e-3)
  

}
