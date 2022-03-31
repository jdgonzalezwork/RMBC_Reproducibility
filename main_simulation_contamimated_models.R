rm(list=ls())
library(RMBC)
library(parallel)
source("utils/dataGeneratingProcess2.R")
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

nrep=5; 
cases=c( "SideNoise3", "SunSpot5","SideNoise2",
         "RandomScatter", "SideNoise2H","RandomScatterH")
metodos=c("RMBC","otrimle","mclust","tclust","tclustoracle")
INIC=array(NaN,dim=c(nrep,length(cases),length(metodos)))
sensitivityGrilla=INIC
specificityGrilla=INIC
MCRGrilla=INIC
klGrilla=INIC
MCRGrilla2019=INIC
MCRGrilla2021=INIC
MCRGrillaOutputsClusters=INIC
RandIndexGrilla=INIC
seedGrid=1:nrep;
t1=date()
options(warn=-1)
library(parallel)
PARALELO=T
meth=1;
l=1;
jj=498
#A=compute_estimator_and_performance_measures(seed = 1,method=metodos[2],case=cases[5])
#A$wrongly_detected_outliers_over_true_outliers
nom="simulationResults2021n_rep_500_partial.RDATA"
for (l in 1:length(cases)){
  for (meth in 1:length(metodos)){
    method=metodos[meth]
    case=cases[l]
    status1 = paste( "doing", "scenario",l, "from a total of", length(cases),"with the method number",
                     meth, "from a total of", length(metodos) ,"methods")
    print(status1)
    status2 = paste( "running scenario", case, ", with method ",
                     method)
    print(status2)
    
    
    if(!PARALELO){
      salidaPar=lapply(seedGrid,compute_estimator_and_performance_measures,method=method,case=case)
    }

    if(PARALELO){
      
     no_cores=detectCores()-1
     cl<-makeCluster(no_cores,type="FORK")
     salidaPar=parLapply(cl,seedGrid,compute_estimator_and_performance_measures,method=method,case=case)
     stopCluster(cl)
    }
    for (jj in 1:length(seedGrid)){
      MCRGrilla2019[jj,l,meth]=salidaPar[[jj]]$mcrSinCero
      MCRGrilla2021[jj,l,meth]=salidaPar[[jj]]$mcrSinCero2021
      MCRGrilla[jj,l,meth]=salidaPar[[jj]]$mcr # adding zero as a new cluster
      MCRGrillaOutputsClusters[jj,l,meth]=salidaPar[[jj]]$mcrOutputsClusters
      klGrilla[jj,l,meth]=salidaPar[[jj]]$kl[1]
      sensitivityGrilla[jj,l,meth]=salidaPar[[jj]]$sensitivity
      specificityGrilla[jj,l,meth]=salidaPar[[jj]]$specificity
      RandIndexGrilla[jj,l,meth]=salidaPar[[jj]]$RandIndex
    } 
  }
  nom=paste("simulationResults",format(Sys.time(), "%B_%d_%Y"),"_rep_", nrep,"_partial.RDATA",sep="")
  save.image(nom)
}
t2=date()

nom=paste("simulationResults",format(Sys.time(), "%B_%d_%Y"), nrep,"_final.RDATA",sep="")
save.image(nom)




######## SCRIPT TO BUILD THE TABLE #####
# test.file
# load("RData/simulationResults2021n_rep_500_final.RDATA")

big=T
contaminated=T
source("table_building_aux.R")
gtres2


