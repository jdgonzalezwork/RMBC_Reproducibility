auxComputeClusters<-function(Y, theta.mu, theta.sigma,theta.alpha){
  
  salqs=quad_disc(Y = Y,theta.alpha = theta.alpha,
                  theta.mu = theta.mu,theta.sigma = theta.sigma)
  
  nonOutliersEstimadosRMBC=is_in_gr( Y = Y,theta.mu =theta.mu,theta.sigma = theta.sigma,
                                     cutoff=0.999)
  indices=apply(salqs,1,function(x) which(x==max(x))[1])
  indices 
}
  
