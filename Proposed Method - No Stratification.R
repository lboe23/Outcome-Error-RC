library(truncnorm)
library(pracma)
library(tidyr)
library(MASS)
library(mixtools)
library(reshape2)

#Note that the file and directory below will need to be changed in order for another user to get the code to run.
Rcpp::sourceCpp('U:/Paper 1/Final Code/RcppFunc.cpp')

#Get rid of any 1's after first true positive
after_first_pos <- function(x){
  npos<-cumsum(x==1)
  (npos==0) | (npos==1 & x==1)
}

#Function to Assess Method Performance
coverCI<-function(betahat,SE,betatrue){
  UL<-betahat+qnorm(.975)*SE
  LL<-betahat+qnorm(.025)*SE
  cover<-ifelse( (UL > betatrue) & (LL < betatrue),1,0)
  return(cover)
}

Data_output_all <- function(dataset,beta_1,beta_2,beta_3){
  dataset<-na.omit(dataset)
  ### avg across simulations
  mean_delta1<-mean(dataset$delta1)
  mean_delta2<-mean(dataset$delta2)
  mean_delta3<-mean(dataset$delta3)
  delta<-rbind(mean_delta1,mean_delta2,mean_delta3)

  mean_truecensrate<-mean(dataset$truecensrate)


  #Average regression coefficients
  mean_beta1<-mean(dataset$beta1)
  mean_beta2<-mean(dataset$beta2)
  mean_beta3<-mean(dataset$beta3)

  #Apply VarBeta fix and get ASE and ESE
  ASE_beta1<-mean(dataset$se_beta1)
  ASE_beta2<-mean(dataset$se_beta2)
  ASE_beta3<-mean(dataset$se_beta3)
  ESE_beta1<-sd(dataset$beta1)
  ESE_beta2<-sd(dataset$beta2)
  ESE_beta3<-sd(dataset$beta3)

  ASE_beta1X<-mean(dataset$se_betaX1)
  ASE_beta2X<-mean(dataset$se_betaX2)
  ASE_beta3X<-mean(dataset$se_betaX3)
  ESE_beta1X<-sd(dataset$betaX1)
  ESE_beta2X<-sd(dataset$betaX2)
  ESE_beta3X<-sd(dataset$betaX3)

  ASE_beta1Y<-mean(dataset$se_betaY1)
  ASE_beta2Y<-mean(dataset$se_betaY2)
  ASE_beta3Y<-mean(dataset$se_betaY3)
  ESE_beta1Y<-sd(dataset$betaY1)
  ESE_beta2Y<-sd(dataset$betaY2)
  ESE_beta3Y<-sd(dataset$betaY3)

  ASE_beta1T<-mean(sqrt(dataset$se_betaT1))
  ASE_beta2T<-mean(sqrt(dataset$se_betaT2))
  ASE_beta3T<-mean(sqrt(dataset$se_betaT3))
  ESE_beta1T<-sd(dataset$betaT1)
  ESE_beta2T<-sd(dataset$betaT2)
  ESE_beta3T<-sd(dataset$betaT3)

  ASE_beta1N<-mean(sqrt(dataset$se_betaN1))
  ASE_beta2N<-mean(sqrt(dataset$se_betaN2))
  ASE_beta3N<-mean(sqrt(dataset$se_betaN3))
  ESE_beta1N<-sd(dataset$betaN1)
  ESE_beta2N<-sd(dataset$betaN2)
  ESE_beta3N<-sd(dataset$betaN3)

  #Mean percent bias corrected = (estimated-target)/target
  bias_beta1<-(dataset$beta1-beta_1)/beta_1
  bias_beta2<-(dataset$beta2-beta_2)/beta_2
  bias_beta3<-(dataset$beta3-beta_3)/beta_3

  bias_beta1X<-(dataset$betaX1-beta_1)/beta_1
  bias_beta2X<-(dataset$betaX2-beta_2)/beta_2
  bias_beta3X<-(dataset$betaX3-beta_3)/beta_3

  bias_beta1Y<-(dataset$betaY1-beta_1)/beta_1
  bias_beta2Y<-(dataset$betaY2-beta_2)/beta_2
  bias_beta3Y<-(dataset$betaY3-beta_3)/beta_3

  bias_beta1N<-(dataset$betaN1-beta_1)/beta_1
  bias_beta2N<-(dataset$betaN2-beta_2)/beta_2
  bias_beta3N<-(dataset$betaN3-beta_3)/beta_3

  bias_beta1T<-(dataset$betaT1-beta_1)/beta_1
  bias_beta2T<-(dataset$betaT2-beta_2)/beta_2
  bias_beta3T<-(dataset$betaT3-beta_3)/beta_3

  mean_bias_beta1<-mean(bias_beta1)*100
  mean_bias_beta2<-mean(bias_beta2)*100
  mean_bias_beta3<-mean(bias_beta3)*100

  mean_bias_beta1X<-mean(bias_beta1X)*100
  mean_bias_beta2X<-mean(bias_beta2X)*100
  mean_bias_beta3X<-mean(bias_beta3X)*100

  mean_bias_beta1Y<-mean(bias_beta1Y)*100
  mean_bias_beta2Y<-mean(bias_beta2Y)*100
  mean_bias_beta3Y<-mean(bias_beta3Y)*100

  mean_bias_beta1N<-mean(bias_beta1N)*100
  mean_bias_beta2N<-mean(bias_beta2N)*100
  mean_bias_beta3N<-mean(bias_beta3N)*100

  mean_bias_beta1T<-mean(bias_beta1T)*100
  mean_bias_beta2T<-mean(bias_beta2T)*100
  mean_bias_beta3T<-mean(bias_beta3T)*100

  #Coverage Probability Calculation
  CP1<-mean(coverCI(dataset$beta1,dataset$se_beta1,beta_1),na.rm=TRUE)
  CP2<-mean(coverCI(dataset$beta2,dataset$se_beta2,beta_2),na.rm=TRUE)
  CP3<-mean(coverCI(dataset$beta3,dataset$se_beta3,beta_3),na.rm=TRUE)

  CP1Y<-mean(coverCI(dataset$betaY1,dataset$se_betaY1,beta_1),na.rm=TRUE)
  CP2Y<-mean(coverCI(dataset$betaY2,dataset$se_betaY2,beta_2),na.rm=TRUE)
  CP3Y<-mean(coverCI(dataset$betaY3,dataset$se_betaY3,beta_3),na.rm=TRUE)

  CP1X<-mean(coverCI(dataset$betaX1,dataset$se_betaX1,beta_1),na.rm=TRUE)
  CP2X<-mean(coverCI(dataset$betaX2,dataset$se_betaX2,beta_2),na.rm=TRUE)
  CP3X<-mean(coverCI(dataset$betaX3,dataset$se_betaX3,beta_3),na.rm=TRUE)

  CP1T<-mean(coverCI(dataset$betaT1,sqrt(dataset$se_betaT1),beta_1),na.rm=TRUE)
  CP2T<-mean(coverCI(dataset$betaT2,sqrt(dataset$se_betaT2),beta_2),na.rm=TRUE)
  CP3T<-mean(coverCI(dataset$betaT3,sqrt(dataset$se_betaT3),beta_3),na.rm=TRUE)

  CP1N<-mean(coverCI(dataset$betaN1,sqrt(dataset$se_betaN1),beta_1),na.rm=TRUE)
  CP2N<-mean(coverCI(dataset$betaN2,sqrt(dataset$se_betaN2),beta_2),na.rm=TRUE)
  CP3N<-mean(coverCI(dataset$betaN3,sqrt(dataset$se_betaN3),beta_3),na.rm=TRUE)

  beta1_results<-cbind(mean_bias_beta1,ASE_beta1,ESE_beta1,CP1)
  beta2_results<-cbind(mean_bias_beta2,ASE_beta2,ESE_beta2,CP2)
  beta3_results<-cbind(mean_bias_beta3,ASE_beta3,ESE_beta3,CP3)
  results<-rbind(beta1_results,beta2_results,beta3_results)
  results<-round(results,4)
  colnames(results) <- (c("Bias","ASE","ESE","Coverage"))

  beta1_resultsX<-cbind(mean_bias_beta1X,ASE_beta1X,ESE_beta1X,CP1X)
  beta2_resultsX<-cbind(mean_bias_beta2X,ASE_beta2X,ESE_beta2X,CP2X)
  beta3_resultsX<-cbind(mean_bias_beta3X,ASE_beta3X,ESE_beta3X,CP3X)
  resultsX<-rbind(beta1_resultsX,beta2_resultsX,beta3_resultsX)
  resultsX<-round(resultsX,4)
  colnames(resultsX) <- (c("Bias","ASE","ESE","Coverage"))

  beta1_resultsY<-cbind(mean_bias_beta1Y,ASE_beta1Y,ESE_beta1Y,CP1Y)
  beta2_resultsY<-cbind(mean_bias_beta2Y,ASE_beta2Y,ESE_beta2Y,CP2Y)
  beta3_resultsY<-cbind(mean_bias_beta3Y,ASE_beta3Y,ESE_beta3Y,CP3Y)
  resultsY<-rbind(beta1_resultsY,beta2_resultsY,beta3_resultsY)
  resultsY<-round(resultsY,4)
  colnames(resultsY) <- (c("Bias","ASE","ESE","Coverage"))

  beta1_resultsN<-cbind(mean_bias_beta1N,ASE_beta1N,ESE_beta1N,CP1N)
  beta2_resultsN<-cbind(mean_bias_beta2N,ASE_beta2N,ESE_beta2N,CP2N)
  beta3_resultsN<-cbind(mean_bias_beta3N,ASE_beta3N,ESE_beta3N,CP3N)
  resultsN<-rbind(beta1_resultsN,beta2_resultsN,beta3_resultsN)
  resultsN<-round(resultsN,4)
  colnames(resultsN) <- (c("Bias","ASE","ESE","Coverage"))

  beta1_resultsT<-cbind(mean_bias_beta1T,ASE_beta1T,ESE_beta1T,CP1T)
  beta2_resultsT<-cbind(mean_bias_beta2T,ASE_beta2T,ESE_beta2T,CP2T)
  beta3_resultsT<-cbind(mean_bias_beta3T,ASE_beta3T,ESE_beta3T,CP3T)
  resultsT<-rbind(beta1_resultsT,beta2_resultsT,beta3_resultsT)
  resultsT<-round(resultsT,4)
  colnames(resultsT) <- (c("Bias","ASE","ESE","Coverage"))

  return(list(Proposed=results,CorrX=resultsX,CorrY=resultsY,True=resultsT,Naive=resultsN,delta=delta,CR=mean_truecensrate))
}

#Functions needed for proposed method variance calculation
CovAfunc_helper<-function(A,k,CovarLam,j_1,j_2,n_r,i_1,i_2){
      covA_val<-A[i_1,1]*A[1,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[1,1])+
        A[i_1,1]*A[2,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[2,1])+
        A[i_1,1]*A[3,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[3,1])+

        A[i_1,1]*A[1,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[1,2])+
        A[i_1,1]*A[2,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[2,2])+
        A[i_1,1]*A[3,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[3,2])+

        A[i_1,1]*A[1,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[1,3])+
        A[i_1,1]*A[2,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[2,3])+
        A[i_1,1]*A[3,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[3,3])


  return(covA_val)
}

CovAfunction<-function(A,k,CovarLam,j_1,j_2,n_r){
  covA_val<-matrix(c(CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=1),CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=1),
            CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=1),CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=2),
            CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=2),CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=2),
            CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=3),CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=3),
            CovAfunc_helper(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=3)),nrow=k,ncol=k)
  return(covA_val)
}

VarB_valid<-function(SigmaBeta,CovarLam,CorrectA,Beta_Star,n_r){
  k<-nrow(SigmaBeta)
  Ainv<-solve(CorrectA)
  ASigBR<-t(Ainv)%*%SigmaBeta%*%Ainv
  varbeta<-matrix(NA,nrow=k,ncol=k)
  for (j_1 in 1:k){
    for (j_2 in 1:k){
      varbeta[j_1,j_2]<-ASigBR[j_1,j_2]+t(Beta_Star)%*%CovAfunction(Ainv,k,CovarLam,j_1,j_2,n_r)%*%Beta_Star
    }
  }
  return(varbeta)
}

####################################################
###### Begin Data Generation and Simulation  #######
####################################################

#test times is vector of pre-scheduled test times
set.seed(1512)

N<-1000
#N<-65000
nvalid<-500

#Varying Sens/Spec:
#sensitivity<-0.61; specificity<-0.995
#sensitivity<-0.80; specificity<-0.90
sensitivity<-0.90; specificity<-0.80

#varying Negative Predictive Value (Web Appendix Table A1):
negpred<-1
#negpred<-0.90
#negpred<-0.98
#negpred<-0.96

########Event Rate

#want testtimes to be integers for visit times to be clinically meaningful
#Assumes exponential distribution for time to event of interest
#Give a baseline hazard rate

#TABLES 1,3,4,S1,S2 - 2 settings
  #testtimes<-c(1,3,4,6) #censoring rate 0.55
  #blambda<-0.094 #censoring rate 0.55

  testtimes<-c(2,5,7,8) #censoring rate 0.90
  blambda<-0.012 #censoring rate 0.90
#Table 2 - 2 settings
  #testtimes<-c(1,3,4,6) #censoring rate 0.55
  #blambda<-0.076 #censoring rate 0.55

  #testtimes<-c(2,5,7,8) #censoring rate 0.90
  #blambda<-0.008 #censoring rate 0.90

#TABLE S2 - WHI Censoring Rate

  #testtimes<-c(2,5,7,8) #Error-prone censoring rate 0.95
  #blambda<-0.004 #Error-prone censoring rate 0.95

#Data is going to be in long form: creates 4 IDs per subject 1-1000
ntest <- length(testtimes)
ID <- rep(1:N, each = ntest)
time <- rep(testtimes, times = N)

#Set a reasonable value for our true beta_1 and true beta_2

#Table 1 beta_1
  beta_1<-log(1.5)

#Table 2 beta_1
  #beta_1<-log(3)

#Table 5 beta_1 (Type I Error)
  #beta_1<-0

  beta_2<-log(.7)
  beta_3<-log(1.3)

betas<-c(beta_1,beta_2,beta_3)
nbeta<-length(betas)

NSIM<-1000
betamat<-betasdmat<-totalmat<-Etamat<-deltamat<-SDBeta<-VarBeta<-betaHatStar<-TrueCensRate<-Pmat<-Zmat<-ErrXMat<-SDBetaX<-VarBetaX<-SDBetaY<-ErrCensRate<-NaiveBetaMat<-NaiveSDMat<-TrueBetaMat<-TrueSDMat<-NULL

for(iter in 1:NSIM){


mu <- rep(0,3)
Sigma1 <- matrix(.3, nrow=3, ncol=3) + diag(3)*.7
x_z_data <- (mvrnorm(n=N, mu=mu, Sigma=Sigma1))
  #covm <- matrix(x_z_data)
#x_z_data <- matrix(rnorm(nbeta * N), ncol = nbeta)
colnames(x_z_data) <- paste0("cov", 1:nbeta)

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
#Subset the data
valid_subset_data<-as.data.frame(cbind(x_star,x_star_star,z_1,z_2))
names(valid_subset_data)<-c("x_star_v","x_starstar_v","z_1_v","z_2_v")
index<-sample(1:nrow(valid_subset_data),nvalid,replace=FALSE)
valid_subset_data<-valid_subset_data[index,]

correction<-with(valid_subset_data,lm(x_starstar_v~x_star_v+z_1_v+z_2_v))
lm<-summary(correction)
alp0<-lm$coefficients[1,1]
la1<-lm$coefficients[2,1]
la21<-lm$coefficients[3,1]
la22<-lm$coefficients[4,1]
covla<-vcov(lm)
covla_fin<-covla[-c(1),-c(1)]


eta1<-matrix(c(la1,0,0,la21,1,0,la22,0,1), nrow=nbeta, ncol=nbeta)

delta1<-eta1[1,1]
delta2<-eta1[1,2]
delta3<-eta1[1,3]

deltamat<-rbind(deltamat,cbind(delta1,delta2,delta3))

#################################################################

lambda1 <- blambda * exp(c(x_z_data %*% betas))

####
ET <- rexp(N, lambda1)
x_z_data <- x_z_data[ID, , drop = F]

ET_true<-ET[ID]
ET[rbinom(N, 1, 1 - negpred) == 1] <- 0

ET <- ET[ID]
x_star <- x_star[ID]

occur <- time > ET
occur_TRUE<- time > ET_true
true_result<-as.numeric(occur_TRUE)
####

probs <- ifelse(occur, sensitivity, 1 - specificity)
result <- rbinom(length(occur), 1, probs)

data <- data.frame(ID, x_z_data,x_star, testtime = time, result = result,true_result=true_result)

data_wide_trueresult <- dcast(data, ID ~ testtime, value.var="true_result")
colnames(data_wide_trueresult)<-c("ID","Vis1T","Vis2T","Vis3T","Vis4T")

data_wide_errresult <- dcast(data, ID ~ testtime, value.var="result")
colnames(data_wide_errresult)<-c("ID","Vis1","Vis2","Vis3","Vis4")

TrueCensRate<-rbind(TrueCensRate,mean(data_wide_trueresult$Vis4==0))
ErrCensRate<-rbind(ErrCensRate,mean(data_wide_errresult$Vis4==0))

#Fit Model for Error in X Only
keep_groupedsurv<-unlist(tapply(data$result,data$ID,after_first_pos))
#unlist - unlist a list of vectors into a single vector
#tapply applies my "after_first_pos" function to the "result" of data vector
datafinal_errX<-data[keep_groupedsurv,]

#Non Error-Prone x - Correct error in X but not Y
datafinal_errX$testtime<-as.factor(datafinal_errX$testtime)
fit_errX<-glm(result~testtime+x_star+cov2+cov3,family=binomial(link="cloglog"),data=datafinal_errX)

fitsum_errX<-summary(fit_errX)
beta_correrrX<-fitsum_errX$coefficients[c("x_star","cov2","cov3"),1]
betaSE_correrrX<-vcov(fitsum_errX)[c("x_star","cov2","cov3"),c("x_star","cov2","cov3")]

corrected_errX_beta<-t(as.matrix(beta_correrrX))%*%solve(eta1)
ErrXMat<-rbind(ErrXMat,corrected_errX_beta)

VarBetaCorrectedX<-VarB_valid(betaSE_correrrX,covla_fin,eta1,beta_correrrX,nvalid)
SDBetaX<-rbind(SDBetaX,c(sqrt(VarBetaCorrectedX[1,1]),sqrt(VarBetaCorrectedX[2,2]),sqrt(VarBetaCorrectedX[3,3])))
VarBetaX<-rbind(VarBetaX,c(VarBetaCorrectedX[1,1],VarBetaCorrectedX[2,2],VarBetaCorrectedX[3,3]))

#True Model
#Fit Model for Error in X Only
keep_groupedsurv_truth<-unlist(tapply(data$true_result,data$ID,after_first_pos))
#unlist - unlist a list of vectors into a single vector
#tapply applies my "after_first_pos" function to the "result" of data vector
datafinal_true<-data[keep_groupedsurv_truth,]
datafinal_true$testtime<-as.factor(datafinal_true$testtime)

fit_truth<-glm(true_result~testtime+cov1+cov2+cov3,family=binomial(link="cloglog"),data=datafinal_true)

fitsum__true<-summary(fit_truth)
beta_truth<-fitsum__true$coefficients[c("cov1","cov2","cov3"),1]
betaSE_truth<-vcov(fitsum__true)[c("cov1","cov2","cov3"),c("cov1","cov2","cov3")]

NaiveBetaMat<-rbind(NaiveBetaMat,beta_correrrX)
NaiveSDMat<-rbind(NaiveSDMat,diag(betaSE_correrrX))

TrueBetaMat<-rbind(TrueBetaMat,beta_truth)
TrueSDMat<-rbind(TrueSDMat,diag(betaSE_truth))


subject<-ID
testtime<-testtimes
result<-result

formula=result~x_star+cov2+cov3

initsurv = 0.1


id <- eval(substitute(subject), data, parent.frame())
time <- eval(substitute(testtime), data, parent.frame())
result <- eval(substitute(result), data, parent.frame())
ord <- order(id, time)

if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    result <- result[ord]
    data <- data[ord, ]}
utime <- sort(unique(time))
timen0 <- (time != 0)
Dm <- dmat(id[timen0], time[timen0], result[timen0], sensitivity,
           specificity, negpred)
J <- ncol(Dm) - 1
nsub <- nrow(Dm)

lami <- log(-log(seq(1, initsurv, length.out = J + +1)[-1]))
lami <- c(lami[1], diff(lami))
tosurv <- function(x) exp(-exp(cumsum(x)))
lowlam <- c(-Inf, rep(0, J - 1))

Xmat <- model.matrix(formula, data = data)[, -1, drop = F]
beta.nm <- colnames(Xmat)
nbeta <- ncol(Xmat)
parmi <- c(lami, beta_correrrX)
uid <- getrids(id, nsub)

Xmat <- Xmat[uid, , drop = F]

q <- optim(parmi, loglikC, gradlikC, lower = c(lowlam, rep(-Inf, nbeta)), Dm = Dm, Xmat = Xmat, method = "L-BFGS-B", hessian = T)

loglik <- -q$value
lam <- q$par[1:J]
surv <- tosurv(lam)
survival <- data.frame(time = utime[utime != 0], surv = surv)
cov <- as.matrix(solve(q$hessian)[-(1:J), -(1:J)])
rownames(cov) <- colnames(cov) <- beta.nm
beta.fit <- q$par[-(1:J)]
beta.sd <- sqrt(diag(cov))
SDBetaY<-rbind(SDBetaY,beta.sd)

betaHatStar<-rbind(betaHatStar,beta.fit)
corrected_eta_beta<-t(as.matrix(beta.fit))%*%solve(eta1)
Etamat<-rbind(Etamat,corrected_eta_beta)

EtaInverse<-solve(eta1)

VarBetaCorrected<-VarB_valid(cov,covla_fin,eta1,beta.fit,nvalid)
SDBeta<-rbind(SDBeta,c(sqrt(VarBetaCorrected[1,1]),sqrt(VarBetaCorrected[2,2]),sqrt(VarBetaCorrected[3,3])))
VarBeta<-rbind(VarBeta,c(VarBetaCorrected[1,1],VarBetaCorrected[2,2],VarBetaCorrected[3,3]))

betamat<-rbind(betamat,beta.fit)
betasdmat<-rbind(betasdmat,beta.sd)

beta_z <- c(corrected_eta_beta[1]/sqrt(VarBetaCorrected[1,1]),corrected_eta_beta[2]/sqrt(VarBetaCorrected[2,2]),corrected_eta_beta[3]/sqrt(VarBetaCorrected[3,3]))
pvalues <- 2 * (1 - pnorm(abs(beta_z)))
Zmat<-rbind(Zmat,beta_z)
Pmat<-rbind(Pmat,pvalues)

totalmat<-cbind(betaHatStar,Etamat,ErrXMat,SDBetaY,SDBeta,SDBetaX,deltamat,TrueCensRate,Zmat,Pmat,ErrCensRate,NaiveBetaMat,NaiveSDMat,TrueBetaMat,TrueSDMat)

print(iter)
}

totalmat<-data.frame(totalmat)

names(totalmat)<-c("betaY1","betaY2","betaY3","beta1","beta2","beta3","betaX1","betaX2","betaX3","se_betaY1","se_betaY2","se_betaY3","se_beta1","se_beta2","se_beta3","se_betaX1","se_betaX2","se_betaX3",
                   "delta1","delta2","delta3","truecensrate","z_beta1","z_beta2","z_beta3","p_beta1","p_beta2","p_beta3","errcensrate",
                   "betaN1","betaN2","betaN3","se_betaN1","se_betaN2","se_betaN3","betaT1","betaT2","betaT3","se_betaT1","se_betaT2","se_betaT3")
Data_output_all(totalmat,beta_1,beta_2,beta_3)
