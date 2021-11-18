# Outcome-Error-RC
Code for simulation studies in paper, "An approximate quasi-likelihood approach for error-prone failure time outcomes and exposures" (Boe et al. 2021, https://doi.org/10.1002/sim.9108) 

Code files include:

(1) R code that runs simulations to assess the performance of the proposed method, the method that corrects for outcome error only, the method that corrects for covariate error only, the naive method that ignores error in the failure time outcome and the exposure, and the method based on the true data in which we have no measurement error in the failure time outcome and the exposure (no stratification): Proposed Method - No Stratification.R

(2) R code that runs simulations to assess the performance of the proposed method, the method that corrects for outcome error only, the method that corrects for covariate error only, the naive method that ignores error in the failure time outcome and the exposure, and the method based on the true data in which we have no measurement error in the failure time outcome and the exposure, all accomodating a baseline hazard that varies across strata: Proposed Method - Stratification.R

(3) A simulated dataset representing hypothetical data with one error prone covariate and four precisely measured covariates that will be used to illustrate our method: Simulated_Data_Example_5Cov_Long.csv

(4) R code that performs a sample data analysis using the simulated data above: Simulated_Data_Sample_Analysis.R. This code illustrates the analysis methods we used in our manuscript to compute hazard ratio (HR) and 95% confidence interval (CI) estimates of a health outcome based on (a) the naive method ignoring error in the outcome and covariate, (b) the method corrected for error in the covariate only, and (c) the proposed method. This is identical to the code that appears in the first section of the Supplementary Materials of Boe et al. 

(5) An R script, Variance_functions.R, containing a single function called "Proposed_Var" which calculates the variance-covariance matrix for the proposed method. This function for the variance calculation can accommodate 1 error-prone covariate and up to 19 precisely-measured covariates, for a total of 20 covariates in the calibration and outcome models. The function from this script will need to be sourced in order to apply the proposed method in files (1), (2), and (4).

Files (1), (2), and (4) rely on functions "loglikC" and "gradlikC" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/loglikC.cpp and "dmat" and "getrids" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/dataproc.cpp in order to run. These Rcpp functions were developed by Gu, Ma, and Balasubramanian (2015) for the following paper: "Semiparametric time to event models in the presence of error-prone, self-reported outcomes—With application to the women’s health initiative."
