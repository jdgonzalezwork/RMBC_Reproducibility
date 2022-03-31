library(combinat)

mcrComun=function(Y,theta.alpha,theta.mu,theta.sigma,actualMu,indicesTrue){
  # mcrComun: computes the misclassification Rate (in percentage %) by estimating the labels according 
  # the closest centers of each group. The procedure is to find the labeling that achieve the Hausdorff  
  # distance between the sets theta.mu and actualMu. 
  
  # Pros: is fast. 
  # Cons: When the estimated centers are poorly determined  this function 
  # can overestimate the "true" misclassification rate
  actualk=length(theta.alpha); 
  p = length(actualMu[[1]]); 
  # computes misclassificationRate 
  estan2=is_in_gr(Y, cutoff=1-1e-3,theta.mu=theta.mu,theta.sigma=theta.sigma)
  valoresDisc2=quad_disc(Y,theta.alpha=theta.alpha,theta.mu=theta.mu,theta.sigma=theta.sigma)
  indices2=apply(valoresDisc2,1,function(x) which(x==max(x))[1])
  indices2[!estan2]=0
  
  # computes the "right" permutation
  theta1= matrix(unlist(actualMu),nrow=actualk,ncol=p,byrow =T)
  theta2= matrix(unlist(theta.mu),nrow=actualk,ncol=p,byrow =T)
  permut=HausdorffPermutation(theta1matrix =theta1, theta2matrix = theta2)$permutation
  indices2Acomodado=indices2
  for (j in 1:actualk){
    indices2Acomodado[indices2==permut[j,2]]=permut[j,1]
  }
  mcrNew=(1-mean((indicesTrue[indicesTrue!=0]-indices2Acomodado[indicesTrue!=0])==0))*100
  mcrNew
}

## mcr, including outlier  as cluster. 
mcr_agregando0=function(Y,theta.alpha,theta.mu,theta.sigma,actualMu,indicesTrue){
  
  # mcragregando00: computes the misclassification Rate (in percentage %) by permuting 
  # all the labels, and then takes 
  # the minumum among all case. The sufix agregando00  means that outliers are 
  # added as a new cluster with label equal to zero 
  # Pros:  Is exact. 
  # Cons: could be very slow  or impracticable 
  # when the number of clusters K is intermediate or high (say, K=10)
  
  
  actualk=length(theta.alpha); p = length(actualMu[[1]]); 
  #  Computes estimated index
  estan2=is_in_gr(Y, cutoff=1-1e-3,theta.mu=theta.mu,theta.sigma=theta.sigma)
  valoresDisc2=quad_disc(Y,theta.alpha=theta.alpha,theta.mu=theta.mu,theta.sigma=theta.sigma)
  indices2=apply(valoresDisc2,1,function(x) which(x==max(x))[1])
  indices2[!estan2]=0
  mcr_agregando00=MCRpermn(cluster1=indices2,clusterTrue=indicesTrue,actualk=actualk+1)$MCR*100
  mcr_agregando00
}

mcr_sinOutliers=function(Y,theta.alpha,theta.mu,theta.sigma,actualMu,indicesTrue){
  #  computes the misclassification Rate (in percentage %) by permuting 
  # all the labels, and then takes the minumum among all case. 
  # The procedure is computed for the true observations, that is the algorithm
  # only is computed for observations that have indicesTrue != 0
  # Pros:  Is exact. 
  # Cons: could be very slow  or impracticable 
  # when the number of clusters K is intermediate or high (say, K=10)
   actualk=length(theta.alpha); 
   p = length(actualMu[[1]]); 
   
  #  Computes estimated outliers
  estan2=is_in_gr(Y, cutoff=1-1e-3,theta.mu=theta.mu,theta.sigma=theta.sigma)
  
  #  Computes the quadratic discriminant  
  valoresDisc2=quad_disc(Y,theta.alpha=theta.alpha,theta.mu=theta.mu,theta.sigma=theta.sigma)
  
  #  Computes the indices by taking the maximum by files 
  indices2=apply(valoresDisc2,1,function(x) which(x==max(x))[1])
  # sets the outliers to zero 
  indices2[!estan2]=0
  
  mcr_sinOutliers=MCRpermn(cluster1=indices2[indicesTrue!=0],clusterTrue=indicesTrue[indicesTrue!=0],actualk=actualk)$MCR*100
  mcr_sinOutliers
}

############################
## Performance measures##### 
############################
performance_measures_OLD=function(Y,indicesDGP,thetaNew.alpha,
                              thetaNew.mu,thetaNew.sigma,
                              actualAlpha,actualSigma,
                              actualMu,nkl=1000,cutoff=1-1e-3,
                              estimated_clusters=estimated_clusters){
  nonOutliers=indicesDGP>0;
  TrueOutliers=indicesDGP==0;
  actualk=length(thetaNew.alpha); 

  # Kullback- Leibler distance
  kl1=kl(thetaNew.alpha,thetaNew.mu,thetaNew.sigma,
         theta0.alpha=actualAlpha/sum(actualAlpha),theta0.mu=actualMu,theta0.sigma=actualSigma,n=nkl)
  
  ####################################################### #####
  ## MCR case I: based on parameters estimation \mu,\sigma,\alpha #####
  #############################################################
  
  # mcrComun: computes the misclassification Rate by estimating the labels according 
  # the closest centers of each group. 
  
  mcr1c=mcrComun(Y,thetaNew.alpha,thetaNew.mu,thetaNew.sigma,actualMu=actualMu,indicesTrue =indicesDGP)
  
  # mcragregando00: computes the misclassification Rate by permuting 
  # all the labels, and then takes 
  # the minumun among each case. The sufix agregando00  means that outliers are 
  # added as a new cluster with label equal to zero 
  
  mcr1_agregando00=mcr_agregando0(Y,thetaNew.alpha,thetaNew.mu,thetaNew.sigma,
                                  actualMu=actualMu,indicesTrue =indicesDGP)

  # mcragregando00: computes the misclassification Rate by estimating permuting 
  # all the labels, and then takes 
  # the minumun among each case. The sufix sinOutliers2021 means that outliers are not 
  # considered as a group, so that the comparison is taken only in the regular observations.
  
  mcr_sinOutliers2021=mcr_sinOutliers(Y,thetaNew.alpha,thetaNew.mu,thetaNew.sigma,
                                      actualMu=actualMu,indicesTrue =indicesDGP)
  
  
  
  ####################################################### #####
  ## MCR case II: based on estimation of clusters label as ### 
  ####it is given by the packages          ####################
  #############################################################
  
  
  mcrOutputsClusters=MCRpermn(cluster1=estimated_clusters[indicesDGP!=0],
                              clusterTrue=indicesDGP[indicesDGP!=0],
                              actualk=actualk)$MCR*100
  ########################
  ########################
  RandIndex=mclust::adjustedRandIndex(estimated_clusters[indicesDGP!=0],indicesDGP[indicesDGP!=0])
  
  ######################
  ######################
  
  
  ##############################################
  ####### Specificity and Sensitivity ##########
  #############################################
  
  outliers1=!is_in_gr(Y, cutoff=cutoff,theta.mu=thetaNew.mu,
                      theta.sigma=thetaNew.sigma);
  
  specificity1=sum(which((!outliers1)) %in% which(nonOutliers))/sum(nonOutliers)
  
  sensitivity1=sum(which((outliers1))%in% which(TrueOutliers))/sum(TrueOutliers)
  ### Add some performance measures for the referee. 
  wrongly_detected_outliers=sum(which((outliers1))%in% which(nonOutliers))
  
  wrongly_detected_outliers_over_detected_outliers=wrongly_detected_outliers/sum(outliers1)
  
  wrongly_detected_outliers_over_true_outliers=wrongly_detected_outliers/sum(TrueOutliers)
  
  ret1=list(kl=kl1, mcrSinCero2021=mcr_sinOutliers2021,
            mcrSinCero=mcr1c, # this is computed from parameters Sigma, mu alpha, only over regular true observation 
            mcr=mcr1_agregando00,# this is computed from parameters Sigma, mu alpha, including true outliers as clusters
            mcrOutputsClusters=mcrOutputsClusters, #this is computed from outputs programs only over regular true observation 
            specificity=specificity1,
            sensitivity=sensitivity1,
            thetaNew.alpha=thetaNew.alpha, 
            RandIndex=RandIndex,
            thetaNew.mu=thetaNew.mu,thetaNew.sigma=thetaNew.sigma,
            wrongly_detected_outliers_over_detected_outliers=wrongly_detected_outliers_over_detected_outliers,
            wrongly_detected_outliers_over_true_outliers=wrongly_detected_outliers_over_true_outliers 
            )
  ret1
}



performance_measures=function(Y,indicesDGP,thetaNew.alpha,
                                  thetaNew.mu,thetaNew.sigma,
                                  actualAlpha,actualSigma,
                                  actualMu,nkl=1000,cutoff=1-1e-3,
                                  estimated_clusters=estimated_clusters){
  nonOutliers=indicesDGP>0;
  TrueOutliers=indicesDGP==0;
  actualk=length(thetaNew.alpha); 
  
  # Kullback- Leibler distance
  kl1=kl(thetaNew.alpha,thetaNew.mu,thetaNew.sigma,
         theta0.alpha=actualAlpha/sum(actualAlpha),theta0.mu=actualMu,theta0.sigma=actualSigma,n=nkl)
  
  ####################################################### #####
  ## MCR case I: based on parameters estimation \mu,\sigma,\alpha #####
  #############################################################
  
  # mcrComun: computes the misclassification Rate by estimating the labels according 
  # the closest centers of each group. 
  
  mcr1c=mcrComun(Y,thetaNew.alpha,thetaNew.mu,thetaNew.sigma,actualMu=actualMu,indicesTrue =indicesDGP)
  
  # mcragregando00: computes the misclassification Rate by permuting 
  # all the labels, and then takes 
  # the minumun among each case. The sufix agregando00  means that outliers are 
  # added as a new cluster with label equal to zero 
  
  mcr1_agregando00=mcr_agregando0(Y,thetaNew.alpha,thetaNew.mu,thetaNew.sigma,
                                  actualMu=actualMu,indicesTrue =indicesDGP)
  
  # mcragregando00: computes the misclassification Rate by estimating permuting 
  # all the labels, and then takes 
  # the minumun among each case. The sufix sinOutliers2021 means that outliers are not 
  # considered as a group, so that the comparison is taken only in the regular observations.
  
  mcr_sinOutliers2021=mcr_sinOutliers(Y,thetaNew.alpha,thetaNew.mu,thetaNew.sigma,
                                      actualMu=actualMu,indicesTrue =indicesDGP)
  
  
  
  ####################################################### #####
  ## MCR case II: based on estimation of clusters label as ### 
  ####it is given by the packages          ####################
  #############################################################
  
  
  mcrOutputsClusters=MCRpermn(cluster1=estimated_clusters[indicesDGP!=0],
                              clusterTrue=indicesDGP[indicesDGP!=0],
                              actualk=actualk)$MCR*100
  ########################
  ########################
  RandIndex=mclust::adjustedRandIndex(estimated_clusters[indicesDGP!=0],indicesDGP[indicesDGP!=0])
  
  ######################
  ######################
  
  
  ##############################################
  ####### Specificity and Sensitivity ##########
  #############################################
  
  outliers1=!is_in_gr(Y, cutoff=cutoff,theta.mu=thetaNew.mu,theta.sigma=thetaNew.sigma);
  specificity1=sum(which((!outliers1)) %in% which(nonOutliers))/sum(nonOutliers)
  
  sensitivity1=sum(which((outliers1))%in% which(TrueOutliers))/sum(TrueOutliers)
  ### Add some performance measures for the referee. 
  wrongly_detected_outliers=sum(which((outliers1))%in% which(nonOutliers))
  
  wrongly_detected_outliers_over_detected_outliers=wrongly_detected_outliers/sum(outliers1)
  
  wrongly_detected_outliers_over_true_outliers=wrongly_detected_outliers/sum(TrueOutliers)
  
  ret1=list(kl=kl1, mcrSinCero2021=mcr_sinOutliers2021,
            mcrSinCero=mcr1c, # this is computed from parameters Sigma, mu alpha, only over regular true observation 
            mcr=mcr1_agregando00,# this is computed from parameters Sigma, mu alpha, including true outliers as clusters
            mcrOutputsClusters=mcrOutputsClusters, #this is computed from outputs programs only over regular true observation 
            specificity=specificity1,
            sensitivity=sensitivity1,
            thetaNew.alpha=thetaNew.alpha, 
            RandIndex=RandIndex,
            thetaNew.mu=thetaNew.mu,thetaNew.sigma=thetaNew.sigma,
            wrongly_detected_outliers_over_detected_outliers=wrongly_detected_outliers_over_detected_outliers,
            wrongly_detected_outliers_over_true_outliers=wrongly_detected_outliers_over_true_outliers 
  )
  ret1
}




######################################
######################################
### KULL BACK LEIBLER DIVERGENCE######
######################################
######################################

gaux_kl= function(X, theta.alpha,theta.mu,theta.sigma,theta0.alpha,theta0.mu,theta0.sigma){
  FX=0; GX=0;  
  # Consider 2 cases, just to avoid computational errors when the number of mixture 
  # components are different. 
  for (j in 1:length(theta.alpha)){
    GX= GX + theta.alpha[j]*mvtnorm::dmvnorm(as.matrix(X),mean=theta.mu[[j]],sigma=as.matrix(theta.sigma[[j]]))
  }
  for (j in 1:length(theta0.alpha)){
    FX= FX + theta0.alpha[j]*mvtnorm::dmvnorm(as.matrix(X),mean=theta0.mu[[j]], sigma=as.matrix(theta0.sigma[[j]]))
  }
  ret = log(FX/GX)
}

kl= function(theta.alpha,theta.mu,theta.sigma,
             theta0.alpha,theta0.mu,theta0.sigma,n=1000){
  # As there are no exact formula for computing K-L when K>2, this function 
  # computes the K-L divergence by a Monte Carlo simulation, where n is the size of 
  # the Monte Carlo sample 
  #   Kl= \int f(x) * ln(f(x)/g(x)) dx = E_{f(x)}{ ln(g(x)/f(x))}
  Xr=generate_mixture_clean(n,theta0.alpha,theta0.mu,theta0.sigma); 
  logg=gaux_kl(Xr, theta.alpha,theta.mu,theta.sigma,theta0.alpha,theta0.mu,theta0.sigma)
  kl=mean(logg); vv=var(logg)
  c(kl,kl-2*sqrt(vv/n),kl+2*sqrt(vv/n))
}


HausdorffDistance=function(theta1matrix,theta2matrix){
  distij=dist2(theta1matrix,theta2matrix)
  D1= max(apply(distij, 1,  min ))
  D2= max(apply(distij, 2, min ))
  DH=max(D1,D2)
  return(DH)
}

dist2 <- function(x, y){ 
  # source from flexClust. See that package to get different distnance 
  # options.
  #  Copyright (C) 2005 Friedrich Leisch
  #  $Id: distances.R 3 2013-06-12 10:06:43Z leisch $
  if(ncol(x)!=ncol(y))
    stop(sQuote("x")," and ",sQuote("y"),
         " must have the same number of columns")
  z <- matrix(0, nrow=nrow(x), ncol=nrow(y))
  for(k in 1:nrow(y)){
    z[,k] <- sqrt( colSums((t(x) - y[k,])^2) )
  }
  z
}


MCRpermn=function(cluster1,clusterTrue,actualk){    
  # computes the misclassification Rate (a number between 0 and 1, Not percentage) by considering all 
  #permutations for the cluster-labels, and then takes 
  # the minumun among all cases. 
  permutaciones= permn(1:actualk); 
  h=1;
  valorMCR=rep(0,length(permutaciones))
  permutAux=matrix(0,nrow=actualk,ncol=2)
  permutAux[,1]=1:actualk; 
  cluster1Aux=0*cluster1; 
  for (h in 1:length(permutaciones)){ 
    permutAux[,2]=permutaciones[[h]]
    for(j in 1:actualk){
      cluster1Aux[cluster1==permutAux[j,2]]=permutAux[j,1]
    } 
    valorMCR[h]=mean(abs(cluster1Aux-clusterTrue)>0)
  }
  
  MCRmin=min(valorMCR)
  indargmin= which(valorMCR==MCRmin)[1]
  permut=permutAux; 
  permut[,2]=permutaciones[[indargmin]]
  return(list(MCR=MCRmin,permut=permut,valorMCR=valorMCR))
}




HausdorffPermutation=function(theta1matrix,theta2matrix){
  # para conjuntos con mismo numero de elementos esta fucion,
  # devuelve la permutacion optima donde se "alcanza" la 
  # distancia de Hausdorff
  drawing=FALSE;
  distij=dist2(theta1matrix,theta2matrix)
  # la idea es ir tomando maxmaxmin por filas/ columnas y suprimiendo 
  # filas y colunas de la matriz distij. 
  # la distancia de Hausdorff  debiera conseguirse ejecutando la siguente linea. 
  ## DH = dist2(theta1matrix,theta2matrix)[permutation[1,1],permutation[1,2]]
  distij0=distij
  actualk=dim(theta1matrix)[1]; 
  permutation=matrix(0,nrow=actualk, ncol=2)
  minimoLocal=rep(0, actualk)
  for (l in 1:(actualk-1)){
    D1= max(apply(distij, 1,  min )) # 1 toma el maximo por fila.
    D2= max(apply(distij, 2, min ))  # 2 toma el maximo por columna.
    indmax=which.max(c(D1,D2)) 
    # indmax coincide con el criterio de busqueda fila columna donde se alcanza el minimo.   
    val1=which.max((apply(distij, indmax,  min )))
    if (indmax == 1) {
      minLocal=min(distij[val1,])
      val2=which.min(distij[val1,]);
      # Tachando filas y columnas: 
      ## esto luce complicado porque r convierte matrices fila en 
      # columna. Y lo arreglo con lo siguiente 
      distij=matrix(distij[-val1,],nrow=dim(distij)[1]-1)[,-val2]
    }  
    if (indmax == 2) {
      minLocal=min(distij[,val1]);
      val2=which.min(distij[,val1]);
      # tachando filas y columnas:  
      distij= as.matrix(distij[,-val1])[-val2,]
    }
    minimoLocal[l]=minLocal
  }
  minimoLocal[actualk]=distij;
  # Defino la permutacion. 
  for (ll in 1:actualk){
    permutation[ll,]=which(distij0==minimoLocal[ll], arr.ind = T)[1,]
  }
  
  if (drawing==TRUE){
    plot(c(theta1matrix[,1],theta2matrix[,1]), c(theta1matrix[,2],theta2matrix[,2]))
    points(theta1matrix[,1],theta1matrix[,2],pch=19,col="red")
    for (l in 1:actualk){
      points(theta1matrix[permutation[l,1],1],theta1matrix[permutation[l,1],2],pch=1,col=palette()[l],cex=3)
      points(theta2matrix[permutation[l,2],1],theta2matrix[permutation[l,2],2],pch=1,col=palette()[l],cex=3)
    }
  }
  return(list(permutation=permutation,distanciasMinimas=minimoLocal))
}

# # testHausdorffPermutaion:
# thet1=rbind(1:8,1:8)
# ruido=0.1*runif(16)
# thet2=thet1 + ruido;
# thet2=thet2[,sample(8)]
# uu=HausdorffPermutation(theta1matrix = t(thet1),theta2matrix = t(thet2))
# plot(thet1[1,],thet1[2,],col=1:8,pch=19)
# points(thet2[1,],thet2[2,],col=1:8,cex=2)
# thet11=thet1[,uu$permutation[,1]]
# thet22=thet2[,uu$permutation[,2]]
# plot(thet11[1,],thet11[2,],col=1:8,pch=19)
# points(thet22[1,],thet22[2,],col=1:8,cex=2)

klfor2normals<-function(theta1.mu,theta1.sigma,theta2.mu,theta2.sigma){
  sigma2inv=ginv(theta2.sigma)
  sigma2inv_mult_sigma1=sigma2inv%*%theta1.sigma;
  
  ter1=sum(diag(sigma2inv_mult_sigma1)) - 
       log(abs(det(sigma2inv_mult_sigma1)))-
        dim(sigma2inv_mult_sigma1)[1];
  
  ter2=t((theta1.mu -theta2.mu ))%*% sigma2inv%*%(theta1.mu -theta2.mu );
  ret=0.5*(ter1+ter2) 
}


sumkl=function(thetaNew.mu,thetaNew.sigma,thetaOld.mu,thetaOld.sigma){
  KK=length(thetaNew.mu)
  aux=rep(0,KK)
  for (j1 in 1:KK){
    aux[j1]=klfor2normals(theta1.mu=thetaNew.mu[[j1]],
                         theta1.sigma=thetaNew.sigma[[j1]],
                         theta2.mu=thetaOld.mu[[j1]],
                         theta2.sigma=thetaOld.sigma[[j1]])
  }
  ret=sum(aux)
  ret 
}


##########################################################################
##########################################################################
#################### test kl vs klfor2normals ############################ 
##########################################################################


# value1=klfor2normals(theta1.mu=salGenerarDatos$actualMu[[j1]],
#                      theta1.sigma=salGenerarDatos$actualSigma[[j1]],
#                      theta2.mu=salGenerarDatos$actualMu[[j1]]+4,
#                      theta2.sigma=salGenerarDatos$actualSigma[[j1]]*0.15)
# 
# 
# value2=kl(theta0.mu=list(salGenerarDatos$actualMu[[j1]]),
#           theta0.sigma=list(salGenerarDatos$actualSigma[[j1]]),
#           theta.mu=list(salGenerarDatos$actualMu[[j1]]+4),
#           theta.sigma=list(salGenerarDatos$actualSigma[[j1]]*0.15),
#           theta.alpha = 1,theta0.alpha = 1,n = 1000)
# 
# #test1: 
# abs(round(value1,2)-round(value2[1],2))==0
# 
# abs(value2[1]-value1)< (value2[3]-value2[1])
# 
##########################################################################
##########################################################################
##########################################################################

generarMixture=function(n,theta.alpha,theta.mu,theta.sigma,devolverIndices=FALSE){
  theta.alpha=theta.alpha/sum(theta.alpha); 
  Y=c();ns=n*theta.alpha; ns[ns!=floor(ns)]=ns[ns!=floor(ns)] +1;indicesDGP=c()
  for (j in 1:length(theta.alpha)){
    #rmvnorm(x, mean = rep(0, p), sigma = diag(p), log = FALSE)
    YN=mvtnorm::rmvnorm(floor(ns[j]),mean=theta.mu[[j]],sigma=theta.sigma[[j]])
    Y= rbind(Y,YN)
    indicesDGP=c(indicesDGP,rep(j,ns[j]))
  }
  Ysal=Y[1:n,]
  indicesDGP=indicesDGP[1:n]
  if (devolverIndices){
    # this is to solve an technical error 
    Ysal=list(Y=Y,indicesDGP=indicesDGP,grupos=indicesDGP)  
  }
  Ysal
}
