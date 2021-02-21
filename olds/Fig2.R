#########################################################
## This code generates a pdf file with the figure 2 ##### 
## and a .tex file whith table 3 of  the manuscript###### 
#########################################################

### Install and load Packages
packages=c("tclust","RSKC","GSE","otrimle","mclust","mvtnorm",
           "ktaucenters","xtable","ellipse","dplyr","combinat"); 
for (i in 1:length(packages)){
  if (packages[i] %in% rownames(installed.packages())==FALSE){ 
         install.packages(packages[i])}
  require(packages[i], character.only = TRUE)
}
source("auxPlotGroupsModelBasedCluster.R")
source("auxComputeClusters.R")

# Load the Workspace. 
library(RMBC)
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
##  resultRMBC keep the parameters estimation.

##############################################################
# 2 Classic case # ###########################################
#Chris Fraley , # Adrian E. Raftery , Scruca Et al.        ###
##############################################################
mod1= densityMclust(Y,G=K)

############################
######## 3.  OTRIMLE ############
###########################
salot=otrimle(Y,G=K)

######################################################
## 4. TCLUST_20 (ORACLE)
######################################################
saltclustOracle=tclust(Y,k=K,alpha=0.2)

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
    thetaNewTclustOracle.alpha[[j]]= saltclustOracle$size[j]/sum(saltclustOracle$size)
  thetaNewTclustOracle.mu[[j]]= saltclustOracle$centers[,j]
  thetaNewTclustOracle.sigma[[j]]= saltclustOracle$cov[,,j]
}

indices2TclustOracle=auxComputeClusters(Y = Y,
                                        theta.alpha =  thetaNewTclustOracle.alpha,
                                        theta.mu =  thetaNewTclustOracle.mu,
                                        theta.sigma =  thetaNewTclustOracle.sigma)

nonOutliersEstimadosTclustOracle=is_in_gr(Y = Y,theta.mu =  thetaNewTclustOracle.mu,
                                    theta.sigma =  thetaNewTclustOracle.sigma,
                                    cutoff=0.999)


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
  YGraf=Y; indicesGraf=indices2MCLUST; theta.mu=thetaNewClas.mu
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
  YGraf=Y; indicesGraf=indices2TclustOracle;
  theta.mu= thetaNewTclustOracle.mu
  theta.sigma= thetaNewTclustOracle.sigma
  auxPlotGroupsModelBasedCluster()
  

  # OTRIMLE
  title1=paste("OTRIMLE")
  nonOutliersEstimadosGraf=nonOutliersEstimadosOtrimle
  YGraf=Y; 
  indicesGraf=indices2Otrimle; 
  theta.mu=thetaNewOtrimle.mu
  theta.sigma=thetaNewOtrimle.sigma
  auxPlotGroupsModelBasedCluster()

  
