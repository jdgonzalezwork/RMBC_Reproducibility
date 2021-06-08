auxPlotGroupsModelBasedCluster=function(){
  coordinates=c(1,2)
  plot(YGraf[,coordinates],
       main=title1,cex.main=3,xlab="",ylab="",lwd=1,pch=19, cex=.7, type="n",
       xlim=c(0,1.1),ylim=c(0,43))
  actualk=K
  uu=c(2,3)
  mu1Y= theta.mu[[1]][2]
  mu2Y= theta.mu[[2]][2]
  light = which.max(c(mu1Y,mu2Y))
  dark = which.min(c(mu1Y,mu2Y))
  colors=c("","")
  symbols=c(0,0)
  colors[light]=uu[1]
  colors[dark]=uu[2]
  symbols[light]=1
  symbols[dark]=2
  
  for (j in 1:length(colors)){
    Ygrupo= YGraf[(indicesGraf==j) & nonOutliersEstimadosGraf,]; 
    ngrupo=dim(Ygrupo)[1]; 
    points(Ygrupo,col=colors[j],lwd=1.5,pch=symbols[j],cex=1)
    points(theta.mu[[j]][1],theta.mu[[j]][2],
           pch=22, bg=colors[j],cex=4,lwd=3)
  }
  
  points(YGraf[!nonOutliersEstimadosGraf,1],
         YGraf[!nonOutliersEstimadosGraf,2],
         lwd=4,cex=2,pch=4)

}

