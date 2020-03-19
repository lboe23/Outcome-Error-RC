# Outcome-Error-RC
Code for simulation studies in paper, "An Approximate Quasi-Likelihood Approach to Analyzing Error-Prone Failure Time Outcomes and Exposures" by Boe et al.

Code files include:

(1) A .R file that runs simulations to assess the performance of the proposed method (no stratification).

(2) A .R file that runs simulations to assess the performance of the proposed method that accommodates a baseline hazard that varies across strata.

(3) A .R file that runs simulations to assess the performance of the naive method, which ignores error in the failure time outcome and the exposure.

(4) A .R file that runs simulation to assess the performance of the "truth," or the scenario in which we have no measurement error in the failure time outcome and the exposure.

Files (1) and (2) use functions "loglikC" and "gradlikC" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/loglikC.cpp and "dmat" and "getrids" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/dataproc.cpp. These Rcpp functions were developed by Gu, Ma, and Balasubramanian (2015) for the following paper: "Semiparametric time to event models in the presence of error-prone, self-reported outcomes—With application to the women’s health initiative."
