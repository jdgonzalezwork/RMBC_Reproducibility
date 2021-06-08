generate_mixture_cont = function(n,theta.alpha,theta.mu,theta.sigma,alphaCont,aunif,bunif){
  Y=matrix(0,nrow=n,ncol=length(theta.mu[[1]]));
  theta.alpha=  c(alphaCont,theta.alpha)
  theta.alpha=  theta.alpha/sum(theta.alpha)
  grupos = generate_contaminated_labels(n,alpha=theta.alpha)
  for (j in 1:(length(theta.alpha)-1)){
    naux=sum(grupos==j)
    Y[grupos==j,]=mvtnorm::rmvnorm(naux,mean=theta.mu[[j]],sigma=theta.sigma[[j]])
  }
  # contaminacion
  naux=sum(grupos==0)
  if(naux>0){
    Y[grupos==0,1]=runif(n = naux,min=aunif[1],max=aunif[2])
    Y[grupos==0,2]=runif(n = naux,min=bunif[1],max=bunif[2])
  }
  
  return(list(Y=Y,grupos=grupos))
}





generate_mixture_t_cont=function(n,theta.alpha,theta.mu,theta.sigma,alphaCont,aunif,bunif,df1=df1){
  Y=matrix(0,nrow=n,ncol=length(theta.mu[[1]]));
  theta.alpha=  c(alphaCont,theta.alpha)
  grupos = generate_contaminated_labels(n,alpha=theta.alpha)
  for (j in 1:(length(theta.alpha)-1)){
    naux=sum(grupos==j)
    Y[grupos==j,]=mvtnorm::rmvt(n=naux,df=df1, delta=theta.mu[[j]],sigma=theta.sigma[[j]])
  }
  # contaminacion
  naux=sum(grupos==0)
  if(naux>0){
    Y[grupos==0,1]=runif(n = naux,min=aunif[1],max=aunif[2])
    Y[grupos==0,2]=runif(n = naux,min=bunif[1],max=bunif[2])
  }
  
  return(list(Y=Y,indicesDGP=grupos,grupos=grupos))
}




generate_mixture_clean=function(n,theta.alpha,theta.mu,theta.sigma,devolverIndices=FALSE){
  theta.alpha=theta.alpha/sum(theta.alpha); 
  p=length(theta.mu[[1]])
  Y=matrix(0,nrow=n,ncol=p);
  
  indicesDGP=generate_clean_labels(n,alpha=theta.alpha)
  
  for (j in 1:length(theta.alpha)){
    Z = indicesDGP==j
    ns = sum(Z)
    Y[Z,] = mvtnorm::rmvnorm(ns,mean=theta.mu[[j]],sigma=theta.sigma[[j]])
  }
  
  Ysal=Y
  if (devolverIndices){
    # a pesar de que sea ineficiente, 
    #devuelvo los indices con 2 variables... 
    Ysal=list(Y=Y,indicesDGP=indicesDGP,grupos=indicesDGP)  
  }
  Ysal
}



generate_clean_labels=function(n,alpha){
 #' generates labels labels between 1:K, where K is the length
 #' of alpha 
 #' INPUTS: 
 #' n: numbers of observations.
 #' alpha: vector of probabilities (sum(alpha=1)) 
 #' VALUE
 #' An array of length(n) with labels according to alpha Prob.
  A=rmultinom(1,size=n,prob = alpha)
  B=rbind(t(A),1:length(A))
  ret=sample(unlist(apply( B,MARGIN = 2, FUN = function(x){rep(x[2],x[1])})))
  ret
}

generate_contaminated_labels=function(n,alpha){
  #' generates labels labels between 0:(K-1), where K is the length
  #' of alpha, labels equal to zero are considered outliers labels 
  #' (zero-labels will be replaced by outliers further)
  #' INPUTS: 
  #' n: numbers of observations.
  #' alpha: vector of probabilities (sum(alpha=1)) 
  #' VALUE
  #' An array of length n with labels according to alpha Prob.
  #' The first position in alpha[1] indicates the contamination 
  #' level
  #' Example 
  #' alpha=c(0.1,0.3,0.4,0.2)
  #' labels=  generate_contaminated_labels(100)
  A=rmultinom(1,size=n,prob = alpha)
  B=rbind(t(A),0:(length(A)-1))
  ret=sample(unlist(apply( B,MARGIN = 2, FUN = function(x){rep(x[2],x[1])})))
  ret
}
