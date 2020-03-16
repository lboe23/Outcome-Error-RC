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
Data_output_truth <- function(dataset,beta_1,beta_2,beta_3){

  ### avg across simulations
  mean_beta1<-mean(dataset$beta1)
  mean_beta2<-mean(dataset$beta2)
  mean_beta3<-mean(dataset$beta3)

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

  return((results))
}


#Simulation for the truth: use binomial glm with cloglog

####################################################
###### Begin Data Generation and Simulation  #######
####################################################

#test times is vector of pre-scheduled test times
set.seed(1512)

N<-1000
nvalid<-500


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

NSIM<-1000

betamat_truth<-se_mat_truth<-answermat_truth<-deltamat<-NULL

for(iter in 1:NSIM){
  time <- rep(testtimes, times = N)
  mu <- rep(0,3)
  Sigma1 <- matrix(.3, nrow=3, ncol=3) + diag(3)*.7
  x_z_data <- (mvrnorm(n=N, mu=mu, Sigma=Sigma1))
  #covm <- matrix(x_z_data)
  #x_z_data <- matrix(rnorm(nbeta * N), ncol = nbeta)

  colnames(x_z_data) <- c("X_1","Z_1","Z_2")

  x_1<-(x_z_data[,1])
  z_1<-(x_z_data[,2])
  z_2<-(x_z_data[,3])

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

  mydata_truth<-data.frame(ID,time,ET,true_result,x_1,z_1,z_2)
  head(mydata_truth)

  keep<-unlist(tapply(mydata_truth$true_result,mydata_truth$ID,after_first_pos))
  #unlist - unlist a list of vectors into a single vector
  #tapply applies my "after_first_pos" function to the "result" of data vector
  datafinal_truth<-mydata_truth[keep,]

  #Non Error-Prone x
  datafinal_truth$time<-as.factor(datafinal_truth$time)
  fit_truth=glm(true_result~time+x_1+z_1+z_2,family=binomial(link="cloglog"),data=datafinal_truth)

  fitsum_truth<-summary(fit_truth)
  beta1est_truth<-fitsum_truth$coefficients[5,1]
  beta2est_truth<-fitsum_truth$coefficients[6,1]
  beta3est_truth<-fitsum_truth$coefficients[7,1]

  beta1se_truth<-fitsum_truth$coefficients[5,2]
  beta2se_truth<-fitsum_truth$coefficients[6,2]
  beta3se_truth<-fitsum_truth$coefficients[7,2]

  betamat_truth<-rbind(betamat_truth,cbind(beta1est_truth,beta2est_truth,beta3est_truth))
  se_mat_truth<-rbind(se_mat_truth,cbind(beta1se_truth,beta2se_truth,beta3se_truth))
  answermat_truth<-cbind(betamat_truth,se_mat_truth)

}
answermat_truth<-data.frame(answermat_truth)
names(answermat_truth)<-c("beta1","beta2","beta3","se_beta1","se_beta2","se_beta3")
outtruth<-Data_output_truth(answermat_truth,beta_1,beta_2,beta_3)

outtruth

