rm(list=ls())
###################### 
# Set the simulation acoording the sample size and 
# in the presence of outliers 
#(variable big can be T or False, contaminated True or False)
big=T
contaminated=F

# The raw simulations results of the paper were saved in the directory 
# "RData"
# If you desired run your own simulation, just run 

# For big samples: 
# main_contaminated_models.R
# main_clean_models.R

# For small samples: 
# main_simulation_clean_and_small_n_models_balanced.R
# main_simulation_small_n_contaminated_models_balanced.R

# Take into account that this can be time consuming. 


# This code is for load data of simulation results.
# You should change this part, if you want to 
# run your own simulation over raw data. 
rds=list.files("RData")
clean_small_N=paste0("RData/",rds[1])
contaminated_small_N=paste0("RData/",rds[2])
contaminated_big_N=paste0("RData/",rds[4])
clean_big_N=paste0("RData/",rds[3])

if(big&contaminated){ load(contaminated_big_N);
  cption= "Contaminated Big N"
}
if(!big&!contaminated){load(clean_small_N)
  cption= "Clean small N"
}
if(big&!contaminated){load(clean_big_N)
  cption= "Clean Big N"}
if(!big&contaminated){load(contaminated_small_N);
  cption= "contaminated Small N";
}
nomarchivo=paste0(gsub(" ","",cption),".tex");


###########################################
#### Setting the enviroment ###############
###########################################
library(RMBC)
library(parallel)
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

cases=c( "SideNoise3", "SunSpot5","SideNoise2",
         "RandomScatter", "SideNoise2H","RandomScatterH")
metodos=c("RMBC","otrimle","mclust","tclust","tclustoracle")



#########################################################
################## TABLE BUILDING  ######################
#########################################################

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
      tibble(sem=1:nrep, perf="mcr*", value=MCRGrillaOutputsClusters[,l,meth],method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="WDoverdetected",
             value=100*wrongly_detected_over_detectedGrilla[,l,meth],
             method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="WDoveroutliers",
             value=100*wrongly_detected_over_actualOutliersGrilla[,l,meth],
             method=metodos[meth],scenario=case),
      tibble(sem=1:nrep, perf="1-specificity",
             value=100*(1-specificityGrilla[,l,meth]),
             method=metodos[meth],scenario=case)
      )
    simulation=bind_rows(simulation,aux)
  }
}

table(MCRGrilla2021[,l,meth])
sensitivityGrilla[,3,4]

if(F){
datgg=simulation%>%filter(perf %in% c("kl","mcr*","mcr0","mcr","sensitivty","1-specificity"))%>%
  mutate(perf=factor(perf,levels=c("mcr*","mcr0","mcr","kl","sensitivty","1-specificity")))
library(ggplot2)
  ggplot(datgg%>%filter(scenario==cases[1]),aes(x=value))+
    geom_histogram()+
    facet_wrap(perf~method,scale = "free_y")+
    xlim(0,5)
#@?facet_wrap
}

library(dplyr)
meaninf=function(x){mean(x[!is.infinite(x)],na.rm=T)}
sdinf=function(x){sd(x[!is.infinite(x)],na.rm=T)/sqrt(length(x[!is.infinite(x)]))}
mean_sd=function(x){
  ret=paste(NA," (",NA,")",sep="")
  prom=round(meaninf(x),2)
  sd1=format(round(sdinf(x),3), nsmall = 3)
  xvalor=(!is.na(prom)) & (!is.na(sd1))
  if (xvalor){
     #  prom=round(mean(x),2)
     #  sd1=format(round(sd(x),3), nsmall = 3)
      if(prom==0){
        prom=format(meaninf(x), nsmall = 2,scientific=T)
      }
      if (sd1=="0.000"){
        sd1=format(sdinf(x), nsmall = 2,scientific=T)
        }
      ret=paste(prom," (",sd1,")",sep="")
  }
  ret
}

unique(simulation$perf)
simulation%>%filter(perf=="sensitivty" ,method=="mclust",scenario=="RandomScatter")
res1=simulation%>%group_by(perf,method,scenario)%>%summarize(valorProm=mean_sd(value))
res2=tidyr::pivot_wider(data = res1,names_from=scenario,values_from=valorProm)%>%
  mutate( method= factor(method,levels=c("RMBC","otrimle","tclust","tclustoracle","mclust")))%>%arrange(perf,method)

unique(res2$perf)

res2$perf[res2$perf=="mcr0"]="emcr"

res2=res2%>%filter(perf %in% c("kl","emcr","mcr","sensitivty","1-specificity"))%>%
mutate(perf=factor(perf,levels=c("emcr","mcr","kl","sensitivty","1-specificity")))
#res2=res2%>%filter(method=="mclust")
res2=res2[,c(1,2,8,5,6,7,3,4)]
if(!big){res2=res2[,-c(5,8)]}
#gtres2=res2%>%group_by(perf)%>%arrange(perf)%>%gt::gt()%>%gt::fmt_number(
#  columns=c( "SideNoise3", "SunSpot5","SideNoise2",
#            "RandomScatter", "SideNoise2H","RandomScatterH"),
#  decimals = 2)
gtres2=res2%>%group_by(perf)%>%arrange(perf)%>%gt::gt(caption=cption)
tab_latex <- gt::as_latex(gtres2)
tab_latex %>%  as.character() %>%  cat()
gtsave(gtres2, filename =nomarchivo)
gt::as_latex(gtres2)
gtres2
