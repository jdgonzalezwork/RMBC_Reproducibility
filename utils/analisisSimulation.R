rm(list=ls())
library(RMBC)
library(dplyr)
source("dataGeneratingProcess.R")
source("aux_generate_mixture.R")
source("../performance_measures.R")
source("compute_estimator_and_performance_measures.R")
packages=c("tclust","RSKC","GSE","otrimle","mclust","mvtnorm",
           "ktaucenters","combinat","testthat"); 
for (i in 1:length(packages)){
  if (packages[i] %in% rownames(installed.packages())==FALSE){ 
    install.packages(packages[i])}
  require(packages[i], character.only = TRUE)
}
nom="resultadoDeSimulacion2021n_rep_508_finalSoloTCLUSToracle.RDATA"
load(nom)
cases
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

unique(simulation$method)
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