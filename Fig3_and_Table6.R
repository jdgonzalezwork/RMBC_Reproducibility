#########################################################
## This code generates a pdf file with the figure 3 ##### 
## and the table 6 of  the manuscript ###### 
#########################################################

### Firstly, you should install the RMBC package 
## following the instructions in 
## https://github.com/jdgonzalezwork/RMBC
### Install and load other Packages
packages=c("tclust","RSKC","GSE","otrimle","mclust","mvtnorm",
           "ktaucenters","combinat","testthat","xtable"); 
for (i in 1:length(packages)){
  if (packages[i] %in% rownames(installed.packages())==FALSE){ 
    install.packages(packages[i])}
  require(packages[i], character.only = TRUE)
}
dirutils= "utils"
source(paste(dirutils,"/auxPlotGroupsModelBasedCluster.R",sep=""))
source(paste(dirutils,"/auxComputeClusters.R",sep=""))
source(paste(dirutils,"/performance_measures.R",sep=""))
source(paste(dirutils,"/aux_generate_mixture.R",sep=""))
options(warn=-1)
# Load the Workspace. 
library(RMBC)
set.seed(1)
Y=phytoplankton_acoustic_data$Y
outliers_index=phytoplankton_acoustic_data$outliers_index
n=dim(Y)[1]
Yclean=Y[-outliers_index,]
trueOutliers=Y[outliers_index,]
K=2


############################################
##  1 ROBUST ESTIMACION: RMBC ##########
############################################
resultRMBC=RMBC(Y = Y,K = K)
estimated_clustersRMBC=resultRMBC$cluster
estimated_clustersRMBC[resultRMBC$outliers]=0
##  resultRMBC keep the parameters estimation.

##############################################################
# 2 Classic case # ###########################################
#Chris Fraley , # Adrian E. Raftery , Scruca Et al.        ###
##############################################################
mod1= densityMclust(Y,G=K)
estimated_clustersMclust=mod1$classification
############################
######## 3.  OTRIMLE ############
###########################
salot=otrimle(Y,G=K)
estimated_clustersotrimle=salot$cluster

######################################################
## 4. TCLUST_20 (ORACLE)
######################################################
saltclustOracle=tclust(Y,k=K,alpha=0.2)
estimated_clusterstoracle=saltclustOracle$cluster



#############################################
## All the estimators have been computed ####
#############################################



##############################################################
# Compute clusters for Classic Case
##############################################################
thetaNewClas.alpha=rep(0,K);
thetaNewClas.mu=vector(mode="list", length=K)
thetaNewClas.sigma=vector(mode="list", length=K)
for (j in 1:K){
  thetaNewClas.alpha[[j]]=mod1$parameters$pro[j]
  thetaNewClas.mu[[j]]=mod1$parameters$mean[,j]
  thetaNewClas.sigma[[j]]=mod1$parameters$variance$sigma[,,j]
}

# We estimate Allocation of Points to Clusters. As is stated in Sec. 5.1
# of the manuscript

salqs=quad_disc(Y = Y,theta.alpha =   thetaNewClas.alpha,
                theta.mu =thetaNewClas.mu,theta.sigma = thetaNewClas.sigma)

indices2MCLUST=apply(salqs,1,function(x) which(x==max(x))[1])

##############################################################
# Compute clusters for RMBC
##############################################################
indices2RMBC=auxComputeClusters(Y=Y,
                                theta.mu=resultRMBC$mu,
                                theta.sigma=resultRMBC$Sigma,
                                theta.alpha=resultRMBC$alpha
)

nonOutliersEstimadosRMBC=is_in_gr(Y = Y,theta.mu =  resultRMBC$mu,
                                  theta.sigma =  resultRMBC$Sigma,
                                  cutoff=0.999)

##########################################
# Compute clusters for Otrimle  ##########
##########################################

tamOtrimle=salot$size[2:(K+1)]
alphaOtrimle=tamOtrimle/sum(tamOtrimle)
#alphaOtrimleOLD=salot$pi[2:(K+1)]
thetaNewOtrimle.alpha=rep(0,K);
thetaNewOtrimle.mu=vector(mode="list", length=K)
thetaNewOtrimle.sigma=vector(mode="list", length=K)
for (j in 1:K){
  thetaNewOtrimle.alpha[[j]]=alphaOtrimle[j]
  thetaNewOtrimle.mu[[j]]=salot$mean[,j]
  thetaNewOtrimle.sigma[[j]]=salot$cov[,,j]
}
indices2Otrimle=auxComputeClusters(Y, theta.alpha = thetaNewOtrimle.alpha,
                                   theta.mu = thetaNewOtrimle.mu,
                                   theta.sigma = thetaNewOtrimle.sigma)

nonOutliersEstimadosOtrimle=is_in_gr(Y = Y,theta.mu =  thetaNewOtrimle.mu,
                                     theta.sigma =  thetaNewOtrimle.sigma,
                                     cutoff=0.999)

##############################################################
# Compute clusters for TCLUST_20  (ORACLE)
##############################################################
thetaNewTclustOracle.alpha=rep(0,K);
thetaNewTclustOracle.mu=vector(mode="list", length=K)
thetaNewTclustOracle.sigma=vector(mode="list", length=K)
for (j in 1:K){
  thetaNewTclustOracle.alpha[[j]]= saltclustOracle$weights[j]
  thetaNewTclustOracle.mu[[j]]= saltclustOracle$centers[,j]
  thetaNewTclustOracle.sigma[[j]]= saltclustOracle$cov[,,j]
}

## Observation with zero labels are corrected 
ZZZ=estimated_clusterstoracle==0
indices2TclustOracleAUX=auxComputeClusters(Y = Y[ZZZ,],
                                        theta.alpha =  thetaNewTclustOracle.alpha,
                                        theta.mu =  thetaNewTclustOracle.mu,
                                        theta.sigma =  thetaNewTclustOracle.sigma)
estimated_clusterstoracle[ZZZ]=indices2TclustOracleAUX;  
nonOutliersEstimadosTclustOracle=is_in_gr(Y = Y,theta.mu =  thetaNewTclustOracle.mu,
                                          theta.sigma =  thetaNewTclustOracle.sigma,
                                          cutoff=0.999)
estimated_clusterstoracle[!nonOutliersEstimadosTclustOracle] =0
indices2TclustOracle=estimated_clusterstoracle
###############################################################
## Reference Estimator:   Classic Estimator in Clean Data   ###
## Chris Fraley , # Adrian E. Raftery , Scruca Et al.        ###
###############################################################
mod1= densityMclust(Yclean,G=K)
thetaNewClasRef.alpha=rep(0,K);
thetaNewClasRef.mu=vector(mode="list", length=K)
thetaNewClasRef.sigma=vector(mode="list", length=K)
for (j in 1:K){
  thetaNewClasRef.alpha[[j]]=mod1$parameters$pro[j]
  thetaNewClasRef.mu[[j]]=mod1$parameters$mean[,j]
  thetaNewClasRef.sigma[[j]]=mod1$parameters$variance$sigma[,,j]
}

salqs=quad_disc(Y = Yclean,theta.alpha =   thetaNewClasRef.alpha,theta.mu =   thetaNewClasRef.mu,theta.sigma =   thetaNewClasRef.sigma)
indicesReference=apply(salqs,1,function(x) which(x==max(x))[1])
indicesDGP=0*as.vector(Y[,1])
indicesDGP[-outliers_index]=indicesReference


######################################
## The results are plotted ############
######################################
######################################
par(mfrow=c(1,1)); plot(0,0)
par(mfrow=c(2,3))

# Reference solution
title1=paste("Reference  (clean Data)")
YGraf=Yclean; indicesGraf=indicesReference; theta.mu=thetaNewClasRef.mu
nonOutliersEstimadosGraf=rep(T,length(indicesGraf))
theta.sigma=thetaNewClasRef.sigma
auxPlotGroupsModelBasedCluster()

title1=paste("MCLUST")
YGraf=Y; 
#indicesGraf=indices2MCLUST;
indicesGraf=estimated_clustersMclust
theta.mu=thetaNewClas.mu
nonOutliersEstimadosGraf=rep(T,length(indicesGraf))
theta.sigma=thetaNewClas.sigma
auxPlotGroupsModelBasedCluster()



# This work (RMBC)
title1=paste("RMBC (This work)")
nonOutliersEstimadosGraf=nonOutliersEstimadosRMBC
YGraf=Y; indicesGraf=indices2RMBC;
theta.mu=resultRMBC$mu
theta.sigma=resultRMBC$Sigma
auxPlotGroupsModelBasedCluster()


# TCLUST ORACLE
title1=paste("TCLUST (Oracle)")
nonOutliersEstimadosGraf=nonOutliersEstimadosTclustOracle
indicesGraf=indices2TclustOracle;
#nonOutliersEstimadosGraf=estimated_clusterstoracle>0
#indicesGraf=estimated_clusterstoracle;
YGraf=Y; 
theta.mu= thetaNewTclustOracle.mu
theta.sigma= thetaNewTclustOracle.sigma
auxPlotGroupsModelBasedCluster()




# OTRIMLE
title1=paste("OTRIMLE")
YGraf=Y; 
nonOutliersEstimadosGraf=nonOutliersEstimadosOtrimle
indicesGraf=indices2Otrimle; 
#nonOutliersEstimadosGraf= estimated_clustersotrimle>0
#indicesGraf=estimated_clustersotrimle; 
theta.mu=thetaNewOtrimle.mu
theta.sigma=thetaNewOtrimle.sigma
auxPlotGroupsModelBasedCluster()

###############################################################
#####PERFORMANCE MEASURES    ########################
###############################################################
perf=data.frame(RMBC=rep(NaN,6),TCLUSTOracle=rep(NaN,6),MCLUST=rep(NaN,6),Otrimle=rep(NaN,6))
thetaNew.mu=resultRMBC$mu;
thetaNew.sigma=resultRMBC$Sigma
thetaNew.alpha=resultRMBC$alpha
medRMBC=performance_measures(Y,indicesDGP,thetaNew.alpha,
                             thetaNew.mu,thetaNew.sigma,actualAlpha=thetaNewClasRef.alpha,
                             actualSigma=thetaNewClasRef.sigma,actualMu=thetaNewClasRef.mu,
                             nkl=100000,estimated_clusters = estimated_clustersRMBC)
MCRs=medRMBC$mcrSinCero
MCR0s=medRMBC$mcr
Specificitys=medRMBC$specificity
Sensitivitys=medRMBC$sensitivity
KLs=medRMBC$kl[1]
MCRstar=medRMBC$mcrOutputsClusters
todo=c(MCRs, MCR0s,MCRstar,Specificitys,Sensitivitys,KLs);
perf$RMBC<-todo




###############################################################
###############################################################



medTCLUSTOracle=performance_measures(Y,indicesDGP,thetaNewTclustOracle.alpha,
                                     thetaNewTclustOracle.mu,thetaNewTclustOracle.sigma,
                                     actualAlpha=thetaNewClasRef.alpha,
                                     actualSigma=thetaNewClasRef.sigma,
                                     actualMu=thetaNewClasRef.mu,
                                     nkl=10000,
                                     estimated_clusters = estimated_clusterstoracle)

Specificitys=  medTCLUSTOracle$specificity
Sensitivitys=  medTCLUSTOracle$sensitivity
MCRs=medTCLUSTOracle$mcrSinCero
MCR0s=medTCLUSTOracle$mcr
KLs=medTCLUSTOracle$kl[1]
MCRstar=medTCLUSTOracle$mcrOutputsClusters
todo=c(MCRs, MCR0s,MCRstar,Specificitys,Sensitivitys,KLs);
perf$TCLUSTOracle<-todo


medMCLUST=performance_measures(Y,indicesDGP,thetaNewClas.alpha,
                               thetaNewClas.mu,thetaNewClas.sigma,
                               actualAlpha=thetaNewClasRef.alpha,
                               actualSigma=thetaNewClasRef.sigma,
                               actualMu=thetaNewClasRef.mu,
                               nkl=10000,
                               estimated_clusters = estimated_clustersMclust)

Specificitys=medMCLUST$specificity
Sensitivitys=medMCLUST$sensitivity
MCRs=medMCLUST$mcrSinCero
MCR0s=medMCLUST$mcr
KLs=medMCLUST$kl[1]
# Fix an error  (2020)
MCRstar=medMCLUST$mcrSinCero
todo=c(MCRs, MCR0s,MCRstar,Specificitys,Sensitivitys,KLs);
perf$MCLUST<-todo

medOtrimle=performance_measures(Y,indicesDGP,thetaNewOtrimle.alpha,
                                thetaNewOtrimle.mu,thetaNewOtrimle.sigma,
                                actualAlpha=thetaNewClasRef.alpha,
                                actualSigma=thetaNewClasRef.sigma,
                                actualMu=thetaNewClasRef.mu,
                                nkl=10000,
                                estimated_clusters = estimated_clustersotrimle)
Specificitys=medOtrimle$specificity
Sensitivitys=medOtrimle$sensitivity
MCRs=medOtrimle$mcrSinCero
MCR0s=medOtrimle$mcr
KLs=medOtrimle$kl[1]
# Fix an error
#MCRstar=medOtrimle$mcrOutputsClusters
MCRstar=medOtrimle$mcrSinCero
todo=c(MCRs, MCR0s,MCRstar,Specificitys,Sensitivitys,KLs);
perf$Otrimle<-todo



######################################
## The results are ploted #############
#######################################
#######################################
FinalPlot=TRUE
if(FinalPlot){
  
  
  
  #  par(mfrow=c(1,1)); plot(0,0)
  ancho=4
  pdf(file = "Figure2.pdf",width = 3*ancho,height =2*ancho)
  #  par(mfrow=c(2,3))
  nf <- layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow=TRUE),
               respect=FALSE) 
  # Reference solution
  title1=paste("Reference  (clean Data)")
  YGraf=Yclean; indicesGraf=indicesReference; theta.mu=thetaNewClasRef.mu
  nonOutliersEstimadosGraf=rep(T,length(indicesGraf))
  theta.sigma=thetaNewClasRef.sigma
  #source("auxPlotGroupsModelBasedCluster.R")
  auxPlotGroupsModelBasedCluster()
  
  title1=paste("MCLUST")
  YGraf=Y; indicesGraf=indices2MCLUST; theta.mu=thetaNewClas.mu
  nonOutliersEstimadosGraf=rep(T,length(indicesGraf))
  theta.sigma=thetaNewClas.sigma
  #source("auxPlotGroupsModelBasedCluster.R")
  auxPlotGroupsModelBasedCluster()
  
  
  
  # este trabajo (RMBC)
  title1=paste("RMBC (This work)")
  nonOutliersEstimadosGraf=nonOutliersEstimadosRMBC
  YGraf=Y; indicesGraf=indices2RMBC;
  theta.mu=thetaNew.mu
  theta.sigma=thetaNew.sigma
  #source("auxPlotGroupsModelBasedCluster.R")
  auxPlotGroupsModelBasedCluster()
  
  
  # TCLUST ORACLE
  title1=paste("TCLUST (Oracle)")
  nonOutliersEstimadosGraf=nonOutliersEstimadosTclustOracle
  indicesGraf=indices2TclustOracle;
  YGraf=Y; indicesGraf=indices2TclustOracle;
  theta.mu= thetaNewTclustOracle.mu
  theta.sigma= thetaNewTclustOracle.sigma
  #source("auxPlotGroupsModelBasedCluster.R")
  auxPlotGroupsModelBasedCluster()
  
  
  
  # OTRIMLE
  title1=paste("OTRIMLE")
  nonOutliersEstimadosGraf=nonOutliersEstimadosOtrimle
  YGraf=Y; indicesGraf=indices2Otrimle; theta.mu=thetaNewOtrimle.mu
  theta.sigma=thetaNewOtrimle.sigma
  auxPlotGroupsModelBasedCluster()
  
  dev.off()
}
#################################################
# the table is generated by library xtable ######
#################################################

library(xtable)
u=colnames(perf)
u[2]="TCLUST"
colnames(perf)<-u;
rownames(perf)<-c("MCR", "MCR_0","MCR*", "Specificity","Sensitivity","KLD")
#only MCR*, Sensitivity and KL are reported in the paper 
perf <- perf[c(3,5,6),]
perf <- perf[c(1,3,2),]
perf[3,]=100*perf[3,]
uu = rownames(perf)
uu[1]="MCR"
rownames(perf)<-uu
captionTable="Performance of the compared clustering procedures when applied to the phy-
toplankton data. The subscripts on TCLUST are the values of the trimming parameter alpha. The choice 0.20 corresponds to the actual known fraction of outliers in the given data
(ORACLE). textbf{MCR}:Misclassification Rate. textbf{KL} Kullback-Leibler divergence."
write.csv(perf,"Table6.csv",row.names= T)
xtable(perf)


