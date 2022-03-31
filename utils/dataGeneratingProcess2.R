# ATENCION CAMBIOS EN EL ARXIV
# sunspot 5 alpha cont 0.005 y alphas redondeados. c(0.15,0.30,0.10,0.15,0.30)
# Matriz Sigma de Side Noise 2. 
# Dimension 20 en Side Noise 2 H 
# df igual a 3 en Sidenoise3T
dataGeneratingProcess=function(sem,case){
  # generate synthetic data, 
  #' sem indicates the seed value 
  #' case must be one of the following 
  #' c("SideNoise2M","SideNoise2HM",  "SideNoise3TM", "SunSpot5M","RandomScatter",  "RandomScatterH")
  #' Return a list 
  #' 
  #'Y: the generated data indicesDGP=indicesDGP,actualMu=actualMu,actualSigma=actualSigma,
  #'  actualAlpha, the probability of each mixture component
  #' sem: the random seed with the data are generated 
  #' case: the case's name 
  #' actualk: number of clusters
  #' alphaCont: contamination level. 
   
  n=1000;
  set.seed(sem)

  if (case=="SunSpot5"){
    actualk=5;
    actualMu=vector(mode="list", length=actualk); 
    actualSigma=vector(mode="list", length=actualk);
    thetaOld.mu=vector(mode="list", length=actualk);
    thetaOld.sigma=vector(mode="list", length=actualk);
    thetaNew.mu=vector(mode="list", length=actualk);
    thetaNew.sigma=vector(mode="list", length=actualk);
    alphaCont=0.005
    ReportedActualAlpha=c(0.15,0.30,0.10,0.15,0.30) # checked
    actualAlpha=(1-alphaCont)*ReportedActualAlpha
    a=c(30,40)# checked
    b=c(30,40)# checked
    p=2; # checked
    n=1000;# checked
    
    actualMu[[1]]=c(0,3) # checked
    actualMu[[2]]=c(7,1) # checked
    actualMu[[3]]=c(5,9) # checked
    actualMu[[4]]=c(-13,5)# checked
    actualMu[[5]]=c(-9,5)# checked
    
    actualSigma[[1]]=diag(p);actualSigma[[1]][1,2]=0.5;actualSigma[[1]][2,1]=0.5# checked
    actualSigma[[2]]=2*diag(p);actualSigma[[2]][1,2]=-1.5;actualSigma[[2]][2,1]=-1.5# checked
    actualSigma[[3]]=2*diag(p);actualSigma[[3]][1,2]=1.3;actualSigma[[3]][2,1]=1.3# checked
    actualSigma[[4]]=diag(p)*.5# checked
    actualSigma[[5]]=diag(p)*2.5# checked
    
    datos=generate_mixture_cont(n=n,theta.alpha=actualAlpha,
                                theta.mu=actualMu,theta.sigma=actualSigma,
                                alphaCont =alphaCont,aunif = a,bunif=b)
    
    Y=datos$Y
    indicesDGP=datos$grupos
  }
  
  
  
  
  if (case=="SideNoise3"){# CHECKED 
    n=1000;
    actualk=3;
    alphaCont=0.1#checked 
    ReportedActualAlpha=c(0.28, 0.33, 0.39) #checked 
    # actual Alpha that considers alphaCont is:
    actualAlpha=(1-alphaCont)*ReportedActualAlpha 
    a = c(-20,15)#checked
    b = c(-50,5)#checked
    p = 2; #checked
    actualMu=vector(mode="list", length=actualk); 
    actualSigma=vector(mode="list", length=actualk);
    thetaOld.mu=vector(mode="list", length=actualk);
    thetaOld.sigma=vector(mode="list", length=actualk);
    thetaNew.mu=vector(mode="list", length=actualk);
    thetaNew.sigma=vector(mode="list", length=actualk);
    
    
    actualMu[[1]]=c(-2,-2)  # CHECKED 
    actualMu[[2]]=c(7,1)    # CHECKED 
    actualMu[[3]]=c(15,19)  # CHECKED 
    
    actualSigma[[1]]=diag(p);actualSigma[[1]][1,2]=0.5;actualSigma[[1]][2,1]=0.5 # CHECKED
    actualSigma[[2]]=2*diag(p);actualSigma[[2]][1,2]=-1.5;actualSigma[[2]][2,1]=-1.5 # CHECKED
    actualSigma[[3]]=2*diag(p);actualSigma[[3]][1,2]=1.3;actualSigma[[3]][2,1]=1.3 # CHECKED
    datos=generate_mixture_t_cont(n=n,theta.alpha=actualAlpha,
                                  theta.mu=actualMu,theta.sigma=actualSigma,
                                  alphaCont =alphaCont,aunif = a,bunif=b,df1=3) # CHECKED
    datos=generate_mixture_cont(n=n,theta.alpha=actualAlpha,
                                theta.mu=actualMu,theta.sigma=actualSigma,
                                alphaCont =alphaCont,aunif = a,bunif=b) # CHECKED  
    
    ## Removing oultiers in the 99% confidence region
    Yout=datos$Y[datos$grupos==0,]
    nCont=nrow(Yout)
    p=ncol(Yout)
    alphaGR1=1e-2;
    for (i in 1:nCont){
      # outlier that are in the confidence region are replaced by a new outlier
      while (
        RMBC::is_in_gr(t(as.matrix(Yout[i,])), 
                       cutoff=1-alphaGR1,
                       theta.mu=actualMu,
                       theta.sigma=actualSigma)
      ){
          Yout[i,1]=runif(n = 1,min=a[1],max=a[2])
          Yout[i,2]=runif(n = 1,min=b[1],max=b[2])
        }  
      }
  
    ##########################################################
    #########################################################
  datos$Y[datos$grupos==0]=Yout
  Y = datos$Y
  indicesDGP = datos$grupos
  }
  
  if (case=="SideNoise2"){ 
    n=1000;     
    actualk=2;
    alphaCont=0.1
    ReportedActualAlpha=c(0.75,0.25)
    actualAlpha=(1-alphaCont)*ReportedActualAlpha #CHECKED 
    a=c(-50,5)#CHECKED
    b=c(-50,5)#CHECKED
    p=2; #CHECKED
    
    actualMu=vector(mode="list", length=actualk); 
    actualSigma=vector(mode="list", length=actualk);
    thetaOld.mu=vector(mode="list", length=actualk);
    thetaOld.sigma=vector(mode="list", length=actualk);
    thetaNew.mu=vector(mode="list", length=actualk);
    thetaNew.sigma=vector(mode="list", length=actualk);
    
    
    actualMu[[1]]=c(-10,5) #CHECKED
    actualMu[[2]]=c(3,13) #CHECKED
    actualSigma[[1]]=.4*diag(p); #CHECKED
    actualSigma[[2]]=1.5*diag(p);actualSigma[[2]][1,2]=-1.1;actualSigma[[2]][2,1]=-1.1; #CHECKED
    
    
    
    datos=generate_mixture_cont(n=n,theta.alpha=actualAlpha,
                                theta.mu=actualMu,theta.sigma=actualSigma,
                                alphaCont =alphaCont,aunif = a,bunif=b)
    Y=datos$Y
    indicesDGP=datos$grupos
  }
  
  
  
  if (case=="SideNoise2H"){ 
    n=2000;     
    actualk=2;
    alphaCont=0.1
    # ReportedActualAlpha=c(0.85,0.15) # this was a bug. 
    ReportedActualAlpha=c(0.75,0.25)
    # CHECKED
    actualAlpha=(1-alphaCont)*ReportedActualAlpha #CHECKED 
    # CHECKED
    
    a=c(-50,5)#CHECKED
    b=c(-50,5)#CHECKED
    p=2; #CHECKED
    
    actualMu=vector(mode="list", length=actualk); 
    actualSigma=vector(mode="list", length=actualk);
    thetaOld.mu=vector(mode="list", length=actualk);
    thetaOld.sigma=vector(mode="list", length=actualk);
    thetaNew.mu=vector(mode="list", length=actualk);
    thetaNew.sigma=vector(mode="list", length=actualk);
    
    
    actualMu[[1]]=c(-10,5) #CHECKED
    actualMu[[2]]=c(3,13) #CHECKED
    actualSigma[[1]]=.4*diag(p); #CHECKED
    actualSigma[[2]]=1.5*diag(p);actualSigma[[2]][1,2]=-1.1;actualSigma[[2]][2,1]=-1.1; #CHECKED
    
    
    datos=generate_mixture_cont(n=n,theta.alpha=actualAlpha,
                                theta.mu=actualMu,theta.sigma=actualSigma,
                                alphaCont =alphaCont,aunif = a,bunif=b)
    
    YhighDimension=mvtnorm::rmvnorm(n=n,mean = rep(0,18),diag(18))
    Y1=datos$Y
    indicesDGP=datos$grupos
    Y=cbind(Y1,YhighDimension)#Checked 
    
    actualMu[[1]]=c(actualMu[[1]],rep(0,18))#Checked 
    actualMu[[2]]=c(actualMu[[2]],rep(0,18))#Checked 
    actualSigmaAux1=actualSigma[[1]]#Checked 
    actualSigmaAux2=actualSigma[[2]]
    actualSigma[[1]] =diag(20)
    actualSigma[[1]][1:2,1:2] =actualSigmaAux1#Checked 
    actualSigma[[2]] =diag(20)
    actualSigma[[2]][1:2,1:2] =actualSigmaAux2#Checked 
  }
  
  
  if (case=="RandomScatter"){ #Checked 
    p=2;#Checked 
    a=3;#Checked 
    actualk=6#Checked 
    actualMu=vector(mode="list", length=actualk); 
    actualSigma=vector(mode="list", length=actualk);
    thetaOld.mu=vector(mode="list", length=actualk);
    thetaOld.sigma=vector(mode="list", length=actualk);
    thetaNew.mu=vector(mode="list", length=actualk);
    thetaNew.sigma=vector(mode="list", length=actualk);
    
    n=1200
    ReportedActualAlpha=c(1/11, 2/11, 2/11, 2/11, 2/11, 2/11)
    alphaCont=0.05
    actualAlpha=(1-alphaCont)*ReportedActualAlpha  
    ncont=floor(0.1*n)
    as=3*(1:actualk - 3)
    for (j in 1:actualk){
        actualMu[[j]]=rep(1,p)*as[j]
    }
    
    ################################################
    ################### BUILDING SIGMA_J  ##########
    ################################################
    # Sigma Random Scatter setting  
    for (j in 1:actualk){
      UUU=2*runif(p^2)-1
      SigmaAux=matrix(UUU,nrow=p,ncol=p)
      SigmaAux=t(SigmaAux)%*%(SigmaAux)
      actualSigma[[j]]=SigmaAux;
    }
    
    ################################################################
    ########### CLEAN DATA GENERATION  ##################################
    ################################################################
    Ysal=generate_mixture_clean(n=n,theta.alpha=actualAlpha,
                                theta.mu=actualMu,
                                theta.sigma=actualSigma,
                                devolverIndices = TRUE)
    
    ###############################################################
    ########### OUTLIERS GENERATION  ##############################
    ###############################################################
    
    maxBound=apply(Ysal$Y,2,max) 
    minBound=apply(Ysal$Y,2,min)
    nCont=floor(alphaCont*n)+1
    Yout=matrix(0,nrow=nCont,ncol=p)
    factorAmplitud=2
    for (j in 1:p){ 
      pmedio=.5*(minBound[j]+maxBound[j]); 
      delta=.5*(maxBound[j]-minBound[j])
      a=pmedio - factorAmplitud*delta 
      b=pmedio + factorAmplitud*delta
      Yout[,j]= runif(nCont,a,b)
    }
    ################################################################
    ################################################################
    ###############################################################
    ################################################################
      
    alphaGR1=1e-2;
    for (i in 1:nCont){
      # outlier that are in the confidence ellipse are replaced by a new outlier
      while (RMBC::is_in_gr(t(as.matrix(Yout[i,])), cutoff=1-alphaGR1,theta.mu=actualMu,theta.sigma=actualSigma)){
        for (j in 1:p){ 
          pmedio=.5*(minBound[j]+maxBound[j]); 
          delta=.5*(maxBound[j]-minBound[j])
          a=pmedio - factorAmplitud*delta 
          b=pmedio + factorAmplitud*delta
          Yout[i,j]= runif(1,a,b)
        }  
      }
    }
    Zcont=sample(1:n,nCont)
    Y=Ysal$Y
    indicesDGP=Ysal$indicesDGP
    Y[Zcont, ]=Yout
    indicesDGP[Zcont]=0
    Y=as.matrix(Y)
    ##########################################################
    #########################################################
  }
  

  if (case=="RandomScatterH"){ #CHECKED
    p=10;#Checked 
    a=3;#Checked 
    actualk=6#Checked 
    actualMu=vector(mode="list", length=actualk); 
    actualSigma=vector(mode="list", length=actualk);
    thetaOld.mu=vector(mode="list", length=actualk);
    thetaOld.sigma=vector(mode="list", length=actualk);
    thetaNew.mu=vector(mode="list", length=actualk);
    thetaNew.sigma=vector(mode="list", length=actualk);
    
    n=1200
    ReportedActualAlpha=c(1/11, 2/11, 2/11, 2/11, 2/11, 2/11)
    alphaCont=0.05
    actualAlpha=(1-alphaCont)*ReportedActualAlpha  
    ncont=floor(0.1*n)
    as=3*(1:actualk - 3)
    for (j in 1:actualk){
      actualMu[[j]]=rep(1,p)*as[j]
    }
    
    ################################################
    ################### BUILDING SIGMA_J  ##########
    ################################################
    # Sigma Random Scatter setting  
    for (j in 1:actualk){
      UUU=2*runif(p^2)-1
      SigmaAux=matrix(UUU,nrow=p,ncol=p)
      SigmaAux=t(SigmaAux)%*%(SigmaAux)
      actualSigma[[j]]=SigmaAux;
    }
    
    ################################################################
    ########### CLEAN DATA GENERATION  ##################################
    ################################################################
    Ysal=generate_mixture_clean(n=n,theta.alpha=actualAlpha,
                                theta.mu=actualMu,
                                theta.sigma=actualSigma,
                                devolverIndices = TRUE)
    
    ###############################################################
    ########### OUTLIERS GENERATION  ##############################
    ###############################################################
    
    maxBound=apply(Ysal$Y,2,max) 
    minBound=apply(Ysal$Y,2,min)
    nCont=floor(alphaCont*n)+1
    Yout=matrix(0,nrow=nCont,ncol=p)
    factorAmplitud=2
    for (j in 1:p){ 
      pmedio=.5*(minBound[j]+maxBound[j]); 
      delta=.5*(maxBound[j]-minBound[j])
      a=pmedio - factorAmplitud*delta 
      b=pmedio + factorAmplitud*delta
      Yout[,j]= runif(nCont,a,b)
    }
    ################################################################
    ################################################################
    ###############################################################
    ################################################################
    
    alphaGR1=1e-2;
    for (i in 1:nCont){
      # outlier that are in the confidence are replaced by a new outlier
      while (
        RMBC::is_in_gr(t(as.matrix(Yout[i,])), 
                       cutoff=1-alphaGR1,
                       theta.mu=actualMu,
                       theta.sigma=actualSigma)
        ){
        for (j in 1:p){ 
          pmedio=.5*(minBound[j]+maxBound[j]); 
          delta=.5*(maxBound[j]-minBound[j])
          a=pmedio - factorAmplitud*delta 
          b=pmedio + factorAmplitud*delta
          Yout[i,j]= runif(1,a,b)
        }  
      }
    }
    Zcont=sample(1:n,nCont)
    Y=Ysal$Y
    indicesDGP=Ysal$indicesDGP
    Y[Zcont, ]=Yout
    indicesDGP[Zcont]=0
    Y=as.matrix(Y)
    ##########################################################
    #########################################################
  }

  ret=list(Y=Y,indicesDGP=indicesDGP,actualMu=actualMu,actualSigma=actualSigma,
           actualAlpha=actualAlpha,sem=sem,case=case,actualk=actualk,alphaCont=alphaCont)
  ret
}











############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
