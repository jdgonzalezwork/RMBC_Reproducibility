rm(list=ls())
library(RMBC)
source("utils/dataGeneratingProcess2_small_n_balanced.R")
source("utils/aux_generate_mixture.R")
source("utils/performance_measures.R")
source("utils/compute_estimator_and_performance_measures.R")
packages=c("tclust","RSKC","GSE","otrimle","mclust","mvtnorm",
           "ktaucenters","combinat","testthat"); 

for (i in 1:length(packages)){
  if (packages[i] %in% rownames(installed.packages())==FALSE){ 
    install.packages(packages[i])}
  require(packages[i], character.only = TRUE)
}
options(warn=-1)
nrep=37; 
#cases2D=c( "SideNoise3", "SunSpot5","SideNoise2",
#         "RandomScatter", "SideNoise2H","RandomScatterH")
cases2D=c( "SideNoise3", "SunSpot5","SideNoise2",
         "RandomScatter")
dat=c()
  seed=1; 
  j=1; 
  for (j in 1:length(cases2D)){
    data=dataGeneratingProcess(sem=seed,case=cases2D[j])
    Aux=tibble(x=data$Y[,1],y=data$Y[,2], cluster=as.factor(data$indicesDGP),case=cases2D[j])
    dat=rbind(dat,Aux)
  }
  
  library(ggplot2)
  p=ggplot(data=dat,aes(x=x,y=y,
                      col=(cluster),
                      shape=(cluster)))+
                      geom_point()+
                      facet_wrap(~case,scales="free")
  print(p)
  
