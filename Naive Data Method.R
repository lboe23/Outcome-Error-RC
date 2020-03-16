library(truncnorm)
library(pracma)
library(tidyr)
library(MASS)
library(gdata)

#Functions for this program:
######################################################################################################################
# 1
#Get rid of any 1's after first true positive
after_first_pos <- function(x){
  npos<-cumsum(x==1)
  (npos==0) | (npos==1 & x==1)
}

# 2
# Output for model for the truth
Data_output_true_naive <- function(dataset,beta_1,beta_2,beta_3){

  ### avg across simulations
  mean_delta<-mean(dataset$delta)

  #Mean of regression coefficients
  mean_beta1<-mean(dataset$beta1)
  mean_beta2<-mean(dataset$beta2)
  mean_beta3<-mean(dataset$beta3)

  #ASE and ESE
  ASE_beta1<-mean(dataset$se_beta1)
  ASE_beta2<-mean(dataset$se_beta2)
  ASE_beta3<-mean(dataset$se_beta3)
  ESE_beta1<-sd(dataset$beta1)
  ESE_beta2<-sd(dataset$beta2)
  ESE_beta3<-sd(dataset$beta3)

  #Mean percent bias corrected = (estimated-target)/target
  bias_beta1<-(dataset$beta1-beta_1)/beta_1
  bias_beta2<-(dataset$beta2-beta_2)/beta_2
  bias_beta3<-(dataset$beta3-beta_3)/beta_3

  mean_bias_beta1<-mean(bias_beta1)*100
  mean_bias_beta2<-mean(bias_beta2)*100
  mean_bias_beta3<-mean(bias_beta3)*100

  #Coverage Probability Calculation
  beta_1_l95<-dataset$beta1-qnorm(0.975)*(dataset$se_beta1)
  beta_1_r95<-dataset$beta1+qnorm(0.975)*(dataset$se_beta1)

  beta_2_l95<-dataset$beta2-qnorm(0.975)*(dataset$se_beta2)
  beta_2_r95<-dataset$beta2+qnorm(0.975)*(dataset$se_beta2)

  beta_3_l95<-dataset$beta3-qnorm(0.975)*(dataset$se_beta3)
  beta_3_r95<-dataset$beta3+qnorm(0.975)*(dataset$se_beta3)

  beta_1_ci<-as.data.frame(cbind(beta_1_l95,beta_1_r95))
  beta_2_ci<-as.data.frame(cbind(beta_2_l95,beta_2_r95))
  beta_3_ci<-as.data.frame(cbind(beta_3_l95,beta_3_r95))

  # check the proportion of intervals containing the parameter
  CP1<-mean(apply(beta_1_ci, 1, findInterval, x = beta_1) == 1)
  CP2<-mean(apply(beta_2_ci, 1, findInterval, x = beta_2) == 1)
  CP3<-mean(apply(beta_3_ci, 1, findInterval, x = beta_3) == 1)
  Coverage<-cbind(CP1,CP2,CP3)
  #row.names(Coverage) <- (c("Coverage"))

  beta1_results<-cbind(mean_bias_beta1,ASE_beta1,ESE_beta1,CP1)
  beta2_results<-cbind(mean_bias_beta2,ASE_beta2,ESE_beta2,CP2)
  beta3_results<-cbind(mean_bias_beta3,ASE_beta3,ESE_beta3,CP3)
  results<-rbind(beta1_results,beta2_results,beta3_results)
  results<-round(results,4)
  colnames(results) <- (c("Bias","ASE","ESE","Coverage"))

  return(list(results,mean_delta))
}


#Simulation for the naive discrete PH model: use binomial glm with cloglog

######################################################################
################Begin data generating mechanism#######################
######################################################################
#Since Gu assumes time is continuous but with distinct points of follow-up
#we will generate data from an exponential distribution.

####################################################
###### Begin Data Generation and Simulation  #######
####################################################

#test times is vector of pre-scheduled test times
set.seed(1512)

N<-1000
nvalid<-500

#Varying Sens/Spec:
sensitivity<-0.80; specificity<-0.90
#sensitivity<-0.90; specificity<-0.80

#varying Negative Predictive Value (Web Appendix Table A1):
negpred<-1
#negpred<-0.90
#negpred<-0.98

########Event Rate

#want testtimes to be integers for visit times to be clinically meaningful
#Assumes exponential distribution for time to event of interest
#Give a baseline hazard rate

#TABLE 1 - 2 settings
#testtimes<-c(1,3,4,6) #censoring rate 0.55
#blambda<-0.094 #censoring rate 0.55

testtimes<-c(2,5,7,8) #censoring rate 0.90
blambda<-0.012 #censoring rate 0.90

#Table 2 - 2 settings
#testtimes<-c(1,3,4,6) #censoring rate 0.55
#blambda<-0.076 #censoring rate 0.55

#testtimes<-c(2,5,7,8) #censoring rate 0.90
#blambda<-0.008 #censoring rate 0.90

#Data is going to be in long form: creates 4 IDs per subject 1-1000
ntest <- length(testtimes)
ID <- rep(1:N, each = ntest)
time <- rep(testtimes, times = N)

#Set a reasonable value for our true beta_1 and true beta_2

#Table 1 beta_1
beta_1<-log(1.5)

#Table 2 beta_1
#beta_1<-log(3)

beta_2<-log(.7)
beta_3<-log(1.3)
betas<-c(beta_1,beta_2,beta_3)
nbeta=length(betas)

#Initializing values for parameters gamma_1,...,gamma_4 and for the beta that will be used in optim
beta1_0<-beta2_0<-beta3_0<-0.5
gamma1_0<-gamma2_0<-gamma3_0<-gamma4_0<-.1
param0b<-c(beta1_0,beta2_0,beta3_0)

NSIM<-1000

betamat_naive<-NULL
se_mat_naive<-NULL
answermat_naive<-NULL
deltamat<-NULL

for(iter in 1:NSIM){
  time <- rep(testtimes, times = N)
  mu <- rep(0,3)
  Sigma1 <- matrix(.3, nrow=3, ncol=3) + diag(3)*.7
  x_z_data <- (mvrnorm(n=N, mu=mu, Sigma=Sigma1))

  #If using stratified approach, create strata:
  #strat <- sample(c('cat_1', 'cat_2','cat_3','cat_4'), N, replace=TRUE)

  colnames(x_z_data) <- c("X_1","Z_1","Z_2")

  x_1<-(x_z_data[,1])
  z_1<-(x_z_data[,2])
  z_2<-(x_z_data[,3])

  #Observed Measurement X* includes random error (mean 0 independent of x) and systematic error (alllowed to depend on true X)
  #X*=alpha0+alphax*X+U
  alpha0<-1
  alphaX<-0.8
  alphaZ1<-0.3
  alphaZ2<-0.5

  #Create error-prone x_star: 2 settings
  x_star<-alpha0+alphaX*x_z_data[,1]+alphaZ1*x_z_data[,2]+alphaZ1*x_z_data[,3]+ rnorm(N,0,1.31) #delta = 0.3
  #x_star<-alpha0+alphaX*x_z_data[,1]+alphaZ1*x_z_data[,2]+alphaZ1*x_z_data[,3]+ rnorm(N,0,0.77) #delta = 0.6

  #Create error-prone x_star: Non-Normal Error (T dist w/ 4 DF or Mixture of Normals, Table 3)
  #x_star<-alpha0+alphaX*x_z_data[,1]+alphaZ1*x_z_data[,2]+alphaZ1*x_z_data[,3]+ rt(N,4) #T dist
  #x_star<-alpha0+alphaX*x_z_data[,1]+alphaZ1*x_z_data[,2]+alphaZ1*x_z_data[,3]+ rnormmix(N, lambda=c(0.4,0.6), mu=c(0,2), sigma=c(1,1.5)) #Mixture of Normals

  x_star_star<-x_z_data[,1]+ rnorm(N,0,.25)

  #################################################################
  #Subset the data - first nvalid observations
  valid_subset_data<-as.data.frame(cbind(x_star,x_star_star,z_1,z_2))
  names(valid_subset_data)<-c("x_star_v","x_starstar_v","z_1_v","z_2_v")
  index<-sample(1:nrow(valid_subset_data),nvalid,replace=FALSE)
  valid_subset_data<-valid_subset_data[index,]

  correction<-with(valid_subset_data,lm(x_starstar_v~z_1_v+z_2_v+x_star_v))
  lm<-summary(correction)
  alp0<-lm$coefficients[1,1]
  la1<-lm$coefficients[4,1]
  la21<-lm$coefficients[2,1]
  la22<-lm$coefficients[3,1]

  covla<-lm$cov.unscaled
  covla_fin<-covla[-c(1),-c(1)]

  zero_k1_k2<-matrix(0, nrow=3, ncol=6)
  zero_k2_k1<-matrix(0, nrow=6, ncol=3)
  zero_k2_k2<-matrix(0, nrow=6, ncol=6)

  rbind1 <- rbind(covla_fin, zero_k2_k1)
  rbind2 <- rbind(zero_k1_k2, zero_k2_k2)
  covmatlam<-as.matrix(cbind(rbind1,rbind2))

  eta1<-matrix(c(la1,0,0,la21,1,0,la22,0,1), nrow=3, ncol=3)

  corx1z1<-cor(x_1,z_1)
  corx1z2<-cor(x_1,z_2)
  corz1z2<-cor(z_1,z_2)
  corrmat<-rbind(corrmat,cbind(corx1z1,corx1z2,corz1z2))

  delta<-eta1[1,1]
  deltamat<-rbind(deltamat,delta)

  lambda<-blambda*exp(c(x_1*beta_1+z_1*beta_2+z_2*beta_3))
  ET <- rexp(N, lambda)
  ET <- ET[ID]
  x_star <- x_star[ID]

  #put data in long form
  x_1<- x_1[ID]
  z_1<- z_1[ID]
  z_2<- z_2[ID]

  #Indicator for whether visit time >= actual event time
  occur <- time > ET
  true_result<-as.numeric(occur)
  probs <- ifelse(occur, sensitivity, 1 - specificity)
  result <- rbinom(length(occur), 1, probs)

  mydata_naive<-data.frame(ID,time,ET,result,true_result,x_1,z_1,z_2,x_star,strat)

  #If using stratified approach:
  datas1_s <- mydata_naive[which(mydata_naive$strat=="cat_1"),]
  datas2_s <- mydata_naive[which(mydata_naive$strat=="cat_2"),]
  datas3_s <- mydata_naive[which(mydata_naive$strat=="cat_3"),]
  datas4_s <- mydata_naive[which(mydata_naive$strat=="cat_4"),]


  keep<-unlist(tapply(mydata_naive$result,mydata_naive$ID,after_first_pos))
  #unlist - unlist a list of vectors into a single vector
  #tapply applies my "after_first_pos" function to the "result" of data vector
  datafinal_naive<-mydata_naive[keep,]

  #Non Error-Prone x
  datafinal_naive$time<-as.factor(datafinal_naive$time)

  fit_naive<-glm(result~time+x_star+z_1+z_2,family=binomial(link="cloglog"),data=datafinal_naive)

  #Note: for stratified approach, add interaction between strata and time
  #fit_naive<-glm(result~time+x_star+z_1+z_2+time*strat,family=binomial(link="cloglog"),data=datafinal_naive)

  fitsum_naive<-summary(fit_naive)
  beta1est_naive<-fitsum_naive$coefficients[5,1]
  beta2est_naive<-fitsum_naive$coefficients[6,1]
  beta3est_naive<-fitsum_naive$coefficients[7,1]

  beta1se_naive<-fitsum_naive$coefficients[5,2]
  beta2se_naive<-fitsum_naive$coefficients[6,2]
  beta3se_naive<-fitsum_naive$coefficients[7,2]

  betamat_naive<-rbind(betamat_naive,cbind(beta1est_naive,beta2est_naive,beta3est_naive))
  se_mat_naive<-rbind(se_mat_naive,cbind(beta1se_naive,beta2se_naive,beta3se_naive))
  answermat_naive<-cbind(betamat_naive,se_mat_naive,deltamat)

}
answermat_naive<-data.frame(answermat_naive)
names(answermat_naive)<-c("beta1","beta2","beta3","se_beta1","se_beta2","se_beta3","delta")
outnaive<-Data_output_true_naive(answermat_naive,beta_1,beta_2,beta_3)

outnaive
