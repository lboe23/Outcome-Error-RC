# Outcome-Error-RC
Code for simulation studies in paper, "An Approximate Quasi-Likelihood Approach to Analyzing Error-Prone Failure Time Outcomes and Exposures" by Boe et al.

Code files include:

(1) R code that runs simulations to assess the performance of the proposed method (no stratification): Proposed Method - No Stratification.R

(2) R code that runs simulations to assess the performance of the proposed method that accommodates a baseline hazard that varies across strata: Proposed Method - Stratification.R

(3) R code that runs simulations to assess the performance of the naive method, which ignores error in the failure time outcome and the exposure: Naive Data Method.R

(4) R code that runs simulations to assess the performance of the "truth," or the scenario in which we have no measurement error in the failure time outcome and the exposure: True Data Method.R

(5) A simulated dataset representing hypothetical data with similar features to the Women's Health Intiative (WHI) that will be used to illustrate our method: Simulated_Data_Example_Wide_Form.csv

(6) R code that performs a sample data analysis using the simulated data above: Simulated_Data_Sample_Analysis.R. This code illustrates the analysis methods we used in our manuscript to compute hazard ratio (HR) and 95% confidence interval (CI) estimates of a health outcome for a 20% increase in consumption of a nutrient based on (a) the naive method ignoring error in the outcome and covariate, (b) the method corrected for error in the covariate only, and (c) the proposed method. 

Files (1) and (2) rely on functions "loglikC" and "gradlikC" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/loglikC.cpp and "dmat" and "getrids" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/dataproc.cpp in order to run. These Rcpp functions were developed by Gu, Ma, and Balasubramanian (2015) for the following paper: "Semiparametric time to event models in the presence of error-prone, self-reported outcomes—With application to the women’s health initiative."
