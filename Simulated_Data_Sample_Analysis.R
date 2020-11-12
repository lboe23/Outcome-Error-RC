#Rcpp::sourceCpp('RcppFunc.cpp')
#source('Variance_functions.R')

Rcpp::sourceCpp('U:/Paper 1/Final Code/GitHub October 2020/RcppFunc.cpp')
source('U:/Paper 1/Final Code/GitHub October 2020/Variance_functions.R')


#As suggested by the name of our data set, the input data is in long form, where each row represents one time point and each subject has multiple rows. Below we read our data into R and then present the first 6 rows of the data.

#Read in the simulated data
  data_long<-read.csv("U:/Paper 1/Final Code/GitHub October 2020/Simulated_Data_Example_5Cov_Long.csv",header=TRUE,
                      sep=",",check.names = FALSE)
  head(data_long)




#To apply the proposed method, we will begin by fitting the calibration equations. First, we create a new dataset that only has one row per subject and only includes the members of the calibration subset.

    unique_data<-data_long[!duplicated(data_long$ID),]
  calibration_data<-unique_data[which(unique_data$subset_ind==1),]


#Next, we will fit the calibration model by regressing the covariate measure with classical measurement error, $X_1^{**}$ on the covariate with prone to more extreme error, $X^*$, and other covariates, $Z_1$, $Z_2$, $Z_3$, and $Z_4$. We note that the model below corresponds to equation 6 of in main manuscript: $X^{**}_i=\delta_{(0)}+\delta_{(1)}X_i^{*}+\delta_{(2)}Z_i+V_i$.

    #Perform calibration - get elements needed to correct for exposure error
    calibration_eq<-lm(x_1_starstar~x_1_star+z_1+z_2+z_3+z_4,data=calibration_data)


#We will now save the summary data from the calibration equation and use this to create our multivariate correction factor from equation 8 of the main manuscript, which recall has the following form:

  calibration_summary<-summary(calibration_eq)
  nbeta<-length(calibration_summary$coefficients[-1,1])
  Delta_Mat<-rbind(as.matrix(t(calibration_summary$coefficients[-1,1])),
                   cbind(matrix(0,nrow=nbeta-1,ncol=1),diag(nbeta-1)))


#Next, we want to save the elements of the variance-covariance matrix from the calibration equation, as this will be used later in the computation of the variance-covariance matrix $\Sigma$ for $\hat{\beta}$. Note that we do not need the elements of the variance-covariance matrix that correspond to the intercept term for this approach.

  covla<-vcov(calibration_summary)
  covla_fin<-covla[-c(1),-c(1)]


  #We will now begin the process of fitting our outcome models. As we did for the WHI data example in the text, we will consider 3 approaches: (1) the naive method ignoring error in the outcome and covariate, (2) the regression calibration method that corrects for error in the covariate only, and (3) the proposed method. First, let's assign sensitivity, specificty, and negative predictive value.

  sensitivity<-0.60
  specificity<-0.98
  negpred<-0.95


  #We will now fit our first outcome model, which corresponds to the naive approach. To estimate regression coefficients for the naive grouped continuous time Cox proportional hazards model, we will fit a generalized linear model with a binomial outcome and assume a complementary log-log link.

  data_long$t<-as.factor(data_long$t)
  fit_naive<-glm(y~t+x_1_star+z_1+z_2+z_3+z_4,family=binomial(link="cloglog"),
  data=data_long)
  fitsum_naive<-summary(fit_naive)
  ntest<-length(unique(data_long$t))
  param0b<-fitsum_naive$coefficients[-c(1:ntest),1]


  #Now we will get our data in the format required to use the proposed method. First, we will define the formula that we want to use for our outcome model.

  formula=y~x_1_star+z_1+z_2+z_3+z_4


  #Now, let's make sure our data is ordered properly before we begin calculating sum of the likelihood components.

  initsurv = 0.1
  id <- eval(substitute(data_long$ID), data_long, parent.frame())
  time <- eval(substitute(t), data_long, parent.frame())
  y <- eval(substitute(y), data_long, parent.frame())
  ord <- order(id, time)

  if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    y <- y[ord]
    data <- data[ord, ]}
  utime <- sort(unique(time))
  timen0 <- (time != 0)


  #Now that our data is in an appropriate form, we can calculate the $D$ matrix, defined in section 2.1 of the main manuscript. Additionally, we calculate $J$ and the number of rows in $D$.

  Dm <- dmat(id[timen0], time[timen0], y[timen0], sensitivity,
               specificity, negpred)
  J <- ncol(Dm) - 1
  nsub <- nrow(Dm)


  # As we get ready to maximize our log-likelihood, we want to think of starting values for our survival parameters. To avoid maximization problems due to the ordered constraint of the survival parameters $1=S_1>S_2>...>S_{J+1}>0$, we re-parameterize these terms for optimization. The re-parameterization that we use is a log-log transformation of
  #survival function for $S_2$, and a change in log-log of the survival function for all other parameters. We consider initial values of 0.1 for our survival parameters, then transform these based on this re-parameterization. Additionally, we define a lower bound of $-\infty$ for the first re-parameterized survival function and 0 for the subsequent $J-1$ terms.

  initsurv <- 0.1
  lami <- log(-log(seq(1, initsurv, length.out = J +1)[-1]))
  lami <- c(lami[1], diff(lami))
  tosurv <- function(x) exp(-exp(cumsum(x)))
  lowlam <- c(-Inf, rep(0, J - 1))


  #  Next, we want to create a matrix version of our covariate data which will be used in the maximization of the log-likelihood.

  Xmat <- model.matrix(formula, data = data_long)[, -1, drop = F]
  beta.nm <- colnames(Xmat)
  uid <- getrids(id, nsub)
  Xmat <- Xmat[uid, , drop = F]


  #  We will now maximize our log-likelihood function that corrects for outcome error only using the ``L-BFGS-B" method in the optim function. We will give the lower bound $lowlam$ defined above for our survival function parameter estimates and a lower bound of $\infty$ for our regression coefficient estimates. We will use the $lami$ values defined above as our initial values for our baseline survival functions. We will use the estimated regression parameters from the naive method as our starting values for $\beta_{X1}$, $\beta_{Z1}$, $\beta_{Z2}$, $\beta_{Z3}$, and $\beta_{Z4}$ in the proposed method.

  parmi <- c(lami, param0b)
  loglikeoptimize <- optim(parmi, loglikC, gradlikC,
  lower = c(lowlam, rep(-Inf, nbeta)), Dm = Dm, Xmat = Xmat,
  method = "L-BFGS-B", hessian = T)

  #We can now  invert the Hessian matrix to calculate $\hat{\Sigma}_{\beta^*}$.

  cov <- as.matrix(solve(loglikeoptimize$hessian)[-(1:J), -(1:J)])
  rownames(cov) <- colnames(cov) <- beta.nm
  beta_fit <- loglikeoptimize$par[-(1:J)]


  #It is finally time to apply the proposed method. Below, we calculate our corrected vector of estimated regression coefficients of interest, using equation (7) from the main manuscript: $\hat{\beta}=\hat{\beta}^*\hat{\Delta}^{-1}$. Recall that we computed $\hat{\Delta}$ above using regression calibration.

  corrected_beta<-t(as.matrix(beta_fit))%*%solve(Delta_Mat)
  corrected_beta<-as.data.frame(t(corrected_beta))
  rownames(corrected_beta) <-beta.nm


  #Lastly, we compute the variance for the proposed approach. To do this, we use the function ``Proposed\_Var" from the Variance\_functions.R file that we imported above. This code for the variance calculation can accommodate 1 error-prone covariate and up to 19 precisely-measured covariates, for a total of 20 covariates in the calibration and outcome models. The input values for this function, in order, are the following: (1) $\hat{\Sigma}_{\beta^*}$, the variance-covariance matrix from the method that corrects for outcome error only; (2) the variance-covariance matrix from the calibration model; (3) the estimated multivariate correction factor from regression calibration, $\hat{\Delta}$; and (4) the estimated regression parameters obtained by fitting the model that corrects for outcome error only.

  corrected_beta_var<-Proposed_Var(cov,covla_fin,Delta_Mat,beta_fit)
  SDBeta<-sqrt(corrected_beta_var[1,1])


  #Now, to complete our results table, we will use regression calibration to find the results for the method that corrects for covariate error only:

  corrected_beta_x<-t(as.matrix(param0b))%*%solve(Delta_Mat)
  corrected_beta_x<-as.data.frame(t(corrected_beta_x))
  cov_cloglog<-(vcov(fitsum_naive)[-c(1:ntest),-c(1:ntest)])
  corrected_beta_var_x<-Proposed_Var(cov_cloglog,covla_fin,Delta_Mat,param0b)
  SDBeta_x<-sqrt(corrected_beta_var_x[1,1])


  #The last step is to exponentiate our regression parameters and corresponding confidence interval bounds and put them into a table so that we can present the results for all three methods simultaneously. The final results are presented below:

    #Put results in table
    corrected_HR_myMethod<-cbind(HR=exp(corrected_beta[1,]),
                                 Lower=exp(corrected_beta[1,]-qnorm(0.975)*SDBeta),
                                 Upper=exp(corrected_beta[1,]+qnorm(0.975)*SDBeta))
  naive_HR<-cbind(HR=exp(fitsum_naive$coefficients[c("x_1_star"),1]),Lower=exp(fitsum_naive$coefficients[c("x_1_star"),1]-qnorm(0.975)*fitsum_naive$coefficients[c("x_1_star"),2]),Upper=exp(fitsum_naive$coefficients[c("x_1_star"),1]+qnorm(0.975)*fitsum_naive$coefficients[c("x_1_star"),2]))
  corrected_HR_x<-cbind(HR=exp(corrected_beta_x[1,]),Lower=exp(corrected_beta_x[1,]-qnorm(0.975)*SDBeta_x),Upper=exp(corrected_beta_x[1,]+qnorm(0.975)*SDBeta_x))
  all_results<-rbind(naive_HR,corrected_HR_x,corrected_HR_myMethod)
  rownames(all_results)<-c("Naive","Regression Calibration","Proposed")
  round(all_results,3)
