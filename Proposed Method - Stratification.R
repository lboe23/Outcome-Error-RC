library(tables)
library(truncnorm)
library(pracma)
library(tidyr)
library(MASS)

#Note that the file and directory below will need to be changed in order for another user to get the code to run.
Rcpp::sourceCpp('U:/Paper 1/Final Code/RcppFunc.cpp')

Data_output_strat <- function(dataset,beta_1,beta_2,beta_3){

  ### avg across simulations
  mean_delta1<-mean(dataset$delta1)
  mean_delta2<-mean(dataset$delta2)
  mean_delta3<-mean(dataset$delta3)
  delta<-rbind(mean_delta1,mean_delta2,mean_delta3)

  mean_strat1<-mean(dataset$strat_1)
  mean_strat2<-mean(dataset$strat_2)
  mean_strat3<-mean(dataset$strat_3)
  mean_strat4<-mean(dataset$strat_4)
  strat_all<-rbind(mean_strat1,mean_strat2,mean_strat3,mean_strat4)

  mean_truecensrate<-mean(dataset$TrueCensRate)

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

  beta1_results<-cbind(mean_bias_beta1,ASE_beta1,ESE_beta1,CP1)
  beta2_results<-cbind(mean_bias_beta2,ASE_beta2,ESE_beta2,CP2)
  beta3_results<-cbind(mean_bias_beta3,ASE_beta3,ESE_beta3,CP3)
  results<-rbind(beta1_results,beta2_results,beta3_results)
  results<-round(results,4)
  colnames(results) <- (c("Bias","ASE","ESE","Coverage"))

  return(list(results,delta,mean_truecensrate,strat_all))
}


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

betamat<-betasdmat<-totalmat<-TrueCensRate<-stat_table<-Etamat<-deltamat<-SDBeta<-VarBeta<-betaHatStar<-Zmat<-Pmat<-NULL

mu <- rep(0,3)
Sigma1 <- matrix(.3, nrow=3, ncol=3) + diag(3)*.7

for(iter in 1:NSIM){

  x_z_data <- (mvrnorm(n=N, mu=mu, Sigma=Sigma1))
  strat <- sample(c('cat_1', 'cat_2','cat_3','cat_4'), N, replace=TRUE)
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
  #x_star<-alpha0+alphaX*x_z_data[,1]+ rnorm(N,0,1.31) #delta = 0.3
  x_star<-alpha0+alphaX*x_z_data[,1]+alphaZ1*x_z_data[,2]+alphaZ1*x_z_data[,3]+ rnorm(N,0,1.31) #delta = 0.3
  #x_star<-alpha0+alphaX*x_z_data[,1]+alphaZ1*x_z_data[,2]+alphaZ1*x_z_data[,3]+ rnorm(N,0,0.77) #delta = 0.6

  x_star_star<-x_z_data[,1]+ rnorm(N,0,.25)

  #################################################################
  #Subset the data - first 200 observations
  valid_subset_data<-as.data.frame(cbind(x_star,x_star_star,z_1,z_2))
  names(valid_subset_data)<-c("x_star_v","x_starstar_v","z_1_v","z_2_v")
  index<-sample(1:nrow(valid_subset_data),nvalid,replace=FALSE)
  valid_subset_data<-valid_subset_data[index,]

  correction<-with(valid_subset_data,lm(x_starstar_v~x_star_v+z_1_v+z_2_v))
  summary_correction<-summary(correction)
  alp0<-summary_correction$coefficients[1,1]
  la1<-summary_correction$coefficients[2,1]
  la21<-summary_correction$coefficients[3,1]
  la22<-summary_correction$coefficients[4,1]
  covla<-vcov(summary_correction)
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
  strat <- strat[ID]
  x_star <- x_star[ID]

  occur <- time > ET
  true_result<-as.numeric(occur)

  EventRate<-rbind(EventRate,cbind(mean(true_result)))

  probs <- ifelse(occur, sensitivity, 1 - specificity)
  result <- rbinom(length(occur), 1, probs)

  thisdata <- data.frame(ID, x_z_data,x_star,strat, testtime = time, result = result,true_result=true_result)

  data_wide_trueresult <- dcast(thisdata, ID ~ testtime, value.var="true_result")
  colnames(data_wide_trueresult)<-c("ID","Vis1","Vis2","Vis3","Vis4")

  TrueCensRate<-rbind(TrueCensRate,mean(data_wide_trueresult$Vis4==0))


  datas1 <- thisdata[which(thisdata$strat=="cat_1"),]
  datas2 <- thisdata[which(thisdata$strat=="cat_2"),]
  datas3 <- thisdata[which(thisdata$strat=="cat_3"),]
  datas4 <- thisdata[which(thisdata$strat=="cat_4"),]

  prop_strat<-prop.table(table(thisdata$strat))
  stat_table<-rbind(stat_table,prop_strat)


  subject1<-datas1$ID
  testtime1<-datas1$testtime
  result1<-datas1$result

  subject2<-datas2$ID
  testtime2<-datas2$testtime
  result2<-datas2$result

  subject3<-datas3$ID
  testtime3<-datas3$testtime
  result3<-datas3$result

  subject4<-datas4$ID
  testtime4<-datas4$testtime
  result4<-datas4$result

  formula=result~x_star+cov2+cov3

  negpred = 1
  betai = param0b
  initsurv = 0.1

  id1 <- eval(substitute(subject1), datas1, parent.frame())
  time1 <- eval(substitute(testtime1), datas1, parent.frame())
  result1 <- eval(substitute(result1), datas1, parent.frame())
  ord1 <- order(id1, time1)

  id2 <- eval(substitute(subject2), datas2, parent.frame())
  time2 <- eval(substitute(testtime2), datas2, parent.frame())
  result2 <- eval(substitute(result2), datas2, parent.frame())
  ord2 <- order(id2, time2)

  id3 <- eval(substitute(subject3), datas3, parent.frame())
  time3 <- eval(substitute(testtime3), datas3, parent.frame())
  result3 <- eval(substitute(result3), datas3, parent.frame())
  ord3 <- order(id3, time3)

  id4 <- eval(substitute(subject4), datas4, parent.frame())
  time4 <- eval(substitute(testtime4), datas4, parent.frame())
  result4 <- eval(substitute(result4), datas4, parent.frame())
  ord4 <- order(id4, time4)

  if (is.unsorted(ord1)) {
    id1 <- id1[ord1]
    time1 <- time1[ord1]
    result1 <- result1[ord1]
    data1 <- data1[ord1, ]}

  if (is.unsorted(ord2)) {
    id2 <- id2[ord2]
    time2 <- time2[ord2]
    result2 <- result2[ord2]
    data2 <- data2[ord2, ]}

  if (is.unsorted(ord3)) {
    id3 <- id3[ord3]
    time3 <- time3[ord3]
    result3 <- result3[ord3]
    data3 <- data3[ord3, ]}

  if (is.unsorted(ord4)) {
    id4 <- id4[ord4]
    time4 <- time4[ord4]
    result4 <- result4[ord4]
    data4 <- data4[ord4, ]}

  utime1 <- sort(unique(time1))
  utime2 <- sort(unique(time2))
  utime3 <- sort(unique(time3))
  utime4 <- sort(unique(time4))

  timen01 <- (time1 != 0)
  timen02 <- (time2 != 0)
  timen03 <- (time3 != 0)
  timen04 <- (time4 != 0)

  Dm1 <- dmat(id1[timen01], time1[timen01], result1[timen01], sensitivity,
              specificity, negpred)
  Dm2 <- dmat(id2[timen02], time2[timen02], result2[timen02], sensitivity,
              specificity, negpred)
  Dm3 <- dmat(id3[timen03], time3[timen03], result3[timen03], sensitivity,
              specificity, negpred)
  Dm4 <- dmat(id4[timen04], time4[timen04], result4[timen04], sensitivity,
              specificity, negpred)
  J1 <- ncol(Dm1) - 1
  J2 <- ncol(Dm2) - 1
  J3 <- ncol(Dm3) - 1
  J4 <- ncol(Dm4) - 1

  nsub1 <- nrow(Dm1)
  nsub2 <- nrow(Dm2)
  nsub3 <- nrow(Dm3)
  nsub4 <- nrow(Dm4)

  #param 3
  lami1 <- log(-log(seq(1, initsurv, length.out = J1 + +1)[-1]))
  lami1 <- c(lami1[1], diff(lami1))
  lami2 <- log(-log(seq(1, initsurv, length.out = J2 + +1)[-1]))
  lami2 <- c(lami2[1], diff(lami2))
  lami3 <- log(-log(seq(1, initsurv, length.out = J3 + +1)[-1]))
  lami3 <- c(lami3[1], diff(lami3))
  lami4 <- log(-log(seq(1, initsurv, length.out = J4 + +1)[-1]))
  lami4 <- c(lami4[1], diff(lami4))

  tosurv <- function(x) exp(-exp(cumsum(x)))

  lowlam1 <- c(-Inf, rep(0, J1 - 1))
  lowlam2 <- c(-Inf, rep(0, J2 - 1))
  lowlam3 <- c(-Inf, rep(0, J3 - 1))
  lowlam4 <- c(-Inf, rep(0, J4 - 1))



  Xmat1 <- model.matrix(formula, data = datas1)[, -1, drop = F]
  Xmat2 <- model.matrix(formula, data = datas2)[, -1, drop = F]
  Xmat3 <- model.matrix(formula, data = datas3)[, -1, drop = F]
  Xmat4 <- model.matrix(formula, data = datas4)[, -1, drop = F]

  beta.nm1 <- colnames(Xmat1)
  nbeta1 <- ncol(Xmat1)
  uid1<- getrids(id1, nsub1)
  Xmat1 <- Xmat1[uid1, , drop = F]

  beta.nm2 <- colnames(Xmat2)
  nbeta2 <- ncol(Xmat2)
  uid2<- getrids(id2, nsub2)
  Xmat2 <- Xmat2[uid2, , drop = F]

  beta.nm3 <- colnames(Xmat3)
  nbeta3 <- ncol(Xmat3)
  uid3<- getrids(id3, nsub3)
  Xmat3 <- Xmat3[uid3, , drop = F]

  beta.nm4 <- colnames(Xmat4)
  nbeta_4 <- ncol(Xmat4)
  uid4<- getrids(id4, nsub4)
  Xmat4 <- Xmat4[uid4, , drop = F]

  parmall <- c(lami1,lami2,lami3,lami4,betai)

  loglikStrat <- function(parmsS,Dm1,Dm2,Dm3,Dm4,Xmat1,Xmat2,Xmat3,Xmat4){
    parmi1<-parmsS[c(1:4,17:19)]
    parmi2<-parmsS[c(5:8,17:19)]
    parmi3<-parmsS[c(9:12,17:19)]
    parmi4<-parmsS[c(13:16,17:19)]
    strataLL<-loglikC(parm=parmi1,Dm=Dm1,Xmat=Xmat1)+loglikC(parm=parmi2,Dm=Dm2,Xmat=Xmat2)+loglikC(parm=parmi3,Dm=Dm3,Xmat=Xmat3)+loglikC(parm=parmi4,Dm=Dm4,Xmat=Xmat4)
    return(strataLL)
  }

  q <- optim(parmall, loglikStrat, lower = c(rep(lowlam1,4), rep(-Inf, nbeta)), Dm1 = Dm1, Dm2 = Dm2, Dm3 = Dm3, Dm4 = Dm4, Xmat1 = Xmat1, Xmat2 = Xmat2, Xmat3 = Xmat3, Xmat4 = Xmat4, method = "L-BFGS-B", hessian = T)


  loglik <- -q$value
  totJ<-J1+J2+J3+J4
  lam <- q$par[1:totJ]
  surv <- tosurv(lam)
  survival <- data.frame(time = utime1[utime1 != 0], surv = surv)
  cov <- as.matrix(solve(q$hessian)[-(1:totJ), -(1:totJ)])
  rownames(cov) <- colnames(cov) <- beta.nm1
  beta.fit <- q$par[-(1:totJ)]
  beta.sd <- sqrt(diag(cov))


  betaHatStar<-rbind(betaHatStar,beta.fit)
  corrected_eta_beta<-t(as.matrix(beta.fit))%*%solve(eta1)
  Etamat<-rbind(Etamat,corrected_eta_beta)

  EtaInverse<-solve(eta1)

  VarBetaCorrected<-VarB_valid(cov,covla_fin,eta1,beta.fit,nvalid)
  beta_z <- c(corrected_eta_beta[1]/sqrt(VarBetaCorrected[1,1]),corrected_eta_beta[2]/sqrt(VarBetaCorrected[2,2]),corrected_eta_beta[3]/sqrt(VarBetaCorrected[3,3]))
  pvalues <- 2 * (1 - pnorm(abs(beta_z)))
  Zmat<-rbind(Zmat,beta_z)
  Pmat<-rbind(Pmat,pvalues)

  SDBeta<-rbind(SDBeta,c(sqrt(VarBetaCorrected[1,1]),sqrt(VarBetaCorrected[2,2]),sqrt(VarBetaCorrected[3,3])))
  VarBeta<-rbind(VarBeta,c(VarBetaCorrected[1,1],VarBetaCorrected[2,2],VarBetaCorrected[3,3]))

  betamat<-rbind(betamat,beta.fit)
  betasdmat<-rbind(betasdmat,beta.sd)

  totalmat<-cbind(betaHatStar,Etamat,SDBeta,deltamat,TrueCensRate,Zmat,Pmat,stat_table)


}

totalmat<-data.frame(totalmat)

names(totalmat)<-c("betastar1","betastar2","betastar3","beta1","beta2","beta3","se_beta1","se_beta2","se_beta3","delta1","delta2","delta3","TrueCensRate","z_beta1","z_beta2","z_beta3","p_beta1","p_beta2","p_beta3","strat_1","strat_2","strat_3","strat_4")

beta_1<-log(1.5)
beta_2<-log(.7)
beta_3<-log(1.3)
output_strat<-Data_output_strat(totalmat,beta_1,beta_2,beta_3)
output_strat
