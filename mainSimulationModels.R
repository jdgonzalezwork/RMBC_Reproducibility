rm(list=ls())
library(RMBC)
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

nrep=500; 
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

#library(parallel)
PARALELO=F
meth=1;
l=1;
jj=498

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
    
 #   if(PARALELO){
  #    no_cores=detectCores()-1
   #   cl<-makeCluster(no_cores,type="FORK")
    #  salidaPar=parLapply(cl,seedGrid,compute_estimator_and_performance_measures,method=method,caso=caso,ORACLE=estadoOracle)
     # stopCluster(cl)
  #  }
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
  nom=paste("simulationResults2021n_rep_", nrep,"_partial.RDATA",sep="")
  save.image(nom)
}
t2=date()

nom=paste("simulationResults2021n_rep_", nrep,"_final.RDATA",sep="")
save.image(nom)



##########################################################
################## TABLE BUILDING  ######################
###########################################################

packages=c("dplyr","gt","tidyr");
for (i in 1:length(packages)){
  if (packages[i] %in% rownames(installed.packages())==FALSE){
    install.packages(packages[i])}
  require(packages[i], character.only = TRUE)
}

library(dplyr)

simulation=c()
for (l in 1: length(cases)){
  for (meth in 1:length(metodos)){
    case=cases[l]
    method=metodos[meth]
    aux=bind_rows(
      tibble(sem=1:nrep, perf="mcr", value=MCRGrilla2021[,l,meth],method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="mcr0", value=MCRGrilla[,l,meth],method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="kl", value=klGrilla[,l,meth],method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="sensitivty", value=100*sensitivityGrilla[,l,meth],method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="randindex", value=100*RandIndexGrilla[,l,meth],method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="mcr*", value=MCRGrillaOutputsClusters[,l,meth],method=metodos[meth],scenario=case)
    )
    simulation=bind_rows(simulation,aux)
  }
}
table(MCRGrilla2021[,l,meth])


library(dplyr)
meaninf=function(x){mean(x[!is.infinite(x)],na.rm=T)}
res1=simulation%>%group_by(perf,method,scenario)%>%summarize(valorProm=meaninf(value))
unique(res1$perf)
res2=tidyr::pivot_wider(data = res1,names_from=scenario,values_from=valorProm)%>%
  mutate( method= factor(method,levels=c("RMBC","otrimle","tclust","tclustoracle","mclust")))%>%arrange(perf,method)

res2=res2%>%filter(perf %in% c("kl","mcr*","sensitivty"))%>%
  mutate(perf=factor(perf,levels=c("mcr*","kl","sensitivty")))

res2=res2[,c(1,2,8,5,6,7,3,4)]
gtres2=res2%>%group_by(perf)%>%arrange(perf)%>%gt::gt()%>%gt::fmt_number(
  columns=c( "SideNoise3", "SunSpot5","SideNoise2",
             "RandomScatter", "SideNoise2H","RandomScatterH"),
  decimals = 2)

tab_latex <- gt::as_latex(gtres2) 
tab_latex %>%  as.character() %>%  cat()



