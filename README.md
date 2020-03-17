# Outcome-Error-RC
Code for simulation studies in paper, "An Approximate Quasi-Likelihood Approach to Analyzing Error-Prone Failure Time Outcomes and Exposures" by Boe et al.

Code files include:

(1) A .R file that runs simulations to assess the performance of the proposed method (no stratification).

(2) A .R file that runs simulations to assess the performance of the proposed method that accommodates a baseline hazard that varies across strata.

(3) A .R file that runs simulations to assess the performance of the naive method, which ignores error in the failure time outcome and the exposure.

(4) A .R file that runs simulation to assess the performance of the "truth," or the scenario in which we have no measurement error in the failure time outcome and the exposure.

(5) A .Rcpp file that contains Rcpp functions that are used to run the code in files (1) and (2) above.
