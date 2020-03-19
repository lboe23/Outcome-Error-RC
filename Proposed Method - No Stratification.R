library(tables)
library(truncnorm)
library(pracma)
library(tidyr)
library(MASS)
library(mixtools)
library(reshape2)

#Function to Assess Method Performance
Data_output <- function(dataset,beta_1,beta_2,beta_3){

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

  return(list(results,delta,CR=mean_truecensrate))
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
  testtimes<-c(1,3,4,6) #censoring rate 0.55
  blambda<-0.094 #censoring rate 0.55

  #testtimes<-c(2,5,7,8) #censoring rate 0.90
  #blambda<-0.012 #censoring rate 0.90

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
betamat<-betasdmat<-totalmat<-Etamat<-deltamat<-SDBeta<-VarBeta<-betaHatStar<-TrueCensRate<-Pmat<-Zmat<-NULL

for(iter in 1:NSIM){


mu <- rep(0,3)
Sigma1 <- matrix(.3, nrow=3, ncol=3) + diag(3)*.7
x_z_data <- (mvrnorm(n=N, mu=mu, Sigma=Sigma1))
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
#Subset the data - first 200 observations
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
ET <- rexp(N, lambda1)
x_z_data <- x_z_data[ID, , drop = F]

ET <- ET[ID]
x_star <- x_star[ID]

occur <- time > ET
true_result<-as.numeric(occur)


probs <- ifelse(occur, sensitivity, 1 - specificity)
result <- rbinom(length(occur), 1, probs)

data <- data.frame(ID, x_z_data,x_star, testtime = time, result = result,true_result=true_result)

data_wide_trueresult <- dcast(data, ID ~ testtime, value.var="true_result")
colnames(data_wide_trueresult)<-c("ID","Vis1","Vis2","Vis3","Vis4")

TrueCensRate<-rbind(TrueCensRate,mean(data_wide_trueresult$Vis4==0))


subject<-ID
testtime<-testtimes
result<-result

formula=result~x_star+cov2+cov3

betai = param0b
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
parmi <- c(lami, betai)
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

totalmat<-cbind(betaHatStar,Etamat,SDBeta,deltamat,TrueCensRate,Zmat,Pmat)


}

totalmat<-data.frame(totalmat)

names(totalmat)<-c("betastar1","betastar2","betastar3","beta1","beta2","beta3","se_beta1","se_beta2","se_beta3","delta1","delta2","delta3","truecensrate","z_beta1","z_beta2","z_beta3","p_beta1","p_beta2","p_beta3")
output<-Data_output(totalmat,beta_1,beta_2,beta_3)
output


