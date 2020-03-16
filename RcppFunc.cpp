#include <Rcpp.h>
using namespace Rcpp;

// parameterization type C
// parm[0] is log(-log(S[0])) corresponding to first visit time
// parm[j] is log(-log(S[j])) - log(-log(S[j-1])) corresponding to change
//   in log(-log(S)) in that time period
// This is inspired by representing survival function as exp(-exp(lambda + beta*Z))

// [[Rcpp::export]]
double loglikfunction(NumericVector parm, NumericMatrix Dm, NumericMatrix Xmat) {
  int nsub = Dm.nrow(), J = Dm.ncol() - 1, nbeta = Xmat.ncol(), i, j, k;
  double result = 0, temp, templik, b;
  NumericVector lamb(J), beta(nbeta);
  for (i = 0; i < J; i++) lamb[i] = parm[i];
  for (i = 0; i < nbeta; i++) beta[i] = parm[J + i];
  for (i = 0; i < nsub; i++) {
    templik = Dm(i, 0);
    b = 0;
    for (k = 0; k < nbeta; k++) b += Xmat(i, k)*beta[k];
    temp = b;
    for (j = 0; j < J; j++) {
      temp += parm[j];
      templik += Dm(i, j+1)*exp(-exp(temp));
    }
    result += log(templik);
  }
  return -result;
}

// [[Rcpp::export]]
NumericVector gradientlikfunction(NumericVector parm, NumericMatrix Dm, NumericMatrix Xmat) {
  int nsub = Dm.nrow(), J = Dm.ncol() - 1, nbeta = Xmat.ncol(), i, j, k;
  double temp, templik, b, likj;
  NumericVector lamb(J), beta(nbeta), Dlamb(J), Dbeta(nbeta), result(J + nbeta);
  for (i = 0; i < J; i++) lamb[i] = parm[i];
  for (i = 0; i < nbeta; i++) beta[i] = parm[J + i];
  for (i = 0; i < nsub; i++) {
    templik = Dm(i, 0);
    b = 0;
    for (k = 0; k < nbeta; k++) b += Xmat(i, k)*beta[k];
    temp = b;
    Dlamb.fill(0);
    Dbeta.fill(0);
    for (j = 0; j < J; j++) {
      temp += parm[j];
      likj = Dm(i, j+1)*exp(-exp(temp));
      templik += likj;
      for (k = 0; k <= j; k++) Dlamb[k] += likj*exp(temp);
      for (k = 0; k < nbeta; k++) Dbeta[k] += likj*exp(temp)*Xmat(i, k);
    }
    for (j = 0; j < J; j++) result[j] += Dlamb[j]/templik;
    for (j = 0; j < nbeta; j++) result[J + j] += Dbeta[j]/templik;
  }
  return result;
}

// Functions for processing data for likelihood functions

// [[Rcpp::export]]
NumericMatrix dmat(NumericVector id, NumericVector time, IntegerVector result, double phi1, double phi0, double negpred) {
  NumericVector utime = unique(time);
  utime.sort();
  int J = utime.size(), nsub = unique(id).size(), nobs = id.size(), i, j;
  NumericVector rid(nobs);
  NumericMatrix Cm(nsub, J + 1), Dm(nsub, J + 1);
  utime.push_back(utime[J-1] + 1);
  //initiate with 1 for Cm
  std::fill(Cm.begin(), Cm.end(), 1);
  //Calculate Cm
  j = 0;
  for (i = 1; i < nobs; i++) {
    if (id[i] != id[i-1]) j++;
    rid[i] = j;
  }
  for (j = 0; j <= J; j++) {
    for (i = 0; i < nobs; i++) {
      if (result[i]==0) {
        if (time[i] >= utime[j]) {
          Cm(rid[i], j) *= 1 - phi1;
        } else {
          Cm(rid[i], j) *= phi0;
        }
      } else {
        if (time[i] >= utime[j]) {
          Cm(rid[i], j) *= phi1;
        } else {
          Cm(rid[i], j) *= 1 - phi0;
        }
      }
    }
  }
  //Calculate Dm
  for (i = 0; i < nsub; i++) {
    Dm(i, 0) = Cm(i, 0);
    for (j = 1; j < J+1; j++) {
      Dm(i, j) = negpred*(Cm(i, j) - Cm(i, j-1));
    }
  }
  return Dm;
}

// [[Rcpp::export]]
IntegerVector getrids(NumericVector id, int nsub){
  IntegerVector rid(nsub);
  int j = 0;
  rid[0] = 1;
  for (int i = 1; i < id.size(); i++) {
    if (id[i] != id[i-1]) {
      j++;
      rid[j] = i + 1;
    }
  }
  return rid;
}

