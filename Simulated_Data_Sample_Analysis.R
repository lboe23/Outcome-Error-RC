library(pastecs)
library(survival)
library(survminer)
library(dplyr)
library(splitstackshape)
library(rms)
library(psych)

#16 covariates
CovAfunc_helper_16<-function(A,k,CovarLam,j_1,j_2,n_r,i_1,i_2){
  covA_val<-A[i_1,1]*A[1,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[1,1])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[2,1])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[3,1])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[4,1])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[5,1])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[6,1])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[7,1])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[8,1])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[9,1])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[10,1])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[11,1])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[12,1])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[13,1])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[14,1])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[15,1])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[1,j_2]*(CovarLam[16,1])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[1,2])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[2,2])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[3,2])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[4,2])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[5,2])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[6,2])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[7,2])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[8,2])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[9,2])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[10,2])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[11,2])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[12,2])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[13,2])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[14,2])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[15,2])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[2,j_2]*(CovarLam[16,2])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[1,3])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[2,3])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[3,3])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[4,3])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[5,3])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[6,3])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[7,3])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[8,3])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[9,3])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[10,3])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[11,3])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[12,3])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[13,3])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[14,3])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[15,3])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[3,j_2]*(CovarLam[16,3])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[1,4])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[2,4])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[3,4])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[4,4])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[5,4])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[6,4])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[7,4])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[8,4])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[9,4])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[10,4])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[11,4])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[12,4])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[13,4])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[14,4])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[15,4])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[4,j_2]*(CovarLam[16,4])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[1,5])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[2,5])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[3,5])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[4,5])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[5,5])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[6,5])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[7,5])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[8,5])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[9,5])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[10,5])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[11,5])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[12,5])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[13,5])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[14,5])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[15,5])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[5,j_2]*(CovarLam[16,5])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[1,6])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[2,6])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[3,6])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[4,6])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[5,6])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[6,6])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[7,6])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[8,6])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[9,6])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[10,6])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[11,6])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[12,6])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[13,6])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[14,6])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[15,6])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[6,j_2]*(CovarLam[16,6])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[1,7])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[2,7])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[3,7])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[4,7])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[5,7])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[6,7])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[7,7])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[8,7])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[9,7])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[10,7])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[11,7])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[12,7])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[13,7])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[14,7])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[15,7])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[7,j_2]*(CovarLam[16,7])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[1,8])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[2,8])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[3,8])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[4,8])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[5,8])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[6,8])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[7,8])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[8,8])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[9,8])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[10,8])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[11,8])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[12,8])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[13,8])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[14,8])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[15,8])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[8,j_2]*(CovarLam[16,8])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[1,9])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[2,9])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[3,9])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[4,9])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[5,9])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[6,9])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[7,9])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[8,9])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[9,9])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[10,9])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[11,9])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[12,9])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[13,9])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[14,9])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[15,9])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[9,j_2]*(CovarLam[16,9])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[1,10])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[2,10])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[3,10])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[4,10])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[5,10])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[6,10])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[7,10])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[8,10])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[9,10])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[10,10])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[11,10])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[12,10])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[13,10])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[14,10])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[15,10])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[10,j_2]*(CovarLam[16,10])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[1,11])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[2,11])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[3,11])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[4,11])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[5,11])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[6,11])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[7,11])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[8,11])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[9,11])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[10,11])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[11,11])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[12,11])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[13,11])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[14,11])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[15,11])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[11,j_2]*(CovarLam[16,11])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[1,12])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[2,12])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[3,12])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[4,12])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[5,12])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[6,12])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[7,12])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[8,12])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[9,12])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[10,12])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[11,12])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[12,12])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[13,12])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[14,12])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[15,12])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[12,j_2]*(CovarLam[16,12])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[1,13])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[2,13])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[3,13])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[4,13])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[5,13])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[6,13])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[7,13])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[8,13])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[9,13])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[10,13])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[11,13])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[12,13])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[13,13])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[14,13])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[15,13])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[13,j_2]*(CovarLam[16,13])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[1,14])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[2,14])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[3,14])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[4,14])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[5,14])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[6,14])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[7,14])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[8,14])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[9,14])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[10,14])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[11,14])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[12,14])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[13,14])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[14,14])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[15,14])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[14,j_2]*(CovarLam[16,14])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[1,15])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[2,15])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[3,15])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[4,15])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[5,15])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[6,15])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[7,15])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[8,15])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[9,15])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[10,15])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[11,15])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[12,15])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[13,15])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[14,15])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[15,15])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[15,j_2]*(CovarLam[16,15])+

    A[i_1,1]*A[1,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[1,16])+
    A[i_1,1]*A[2,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[2,16])+
    A[i_1,1]*A[3,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[3,16])+
    A[i_1,1]*A[4,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[4,16])+
    A[i_1,1]*A[5,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[5,16])+
    A[i_1,1]*A[6,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[6,16])+
    A[i_1,1]*A[7,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[7,16])+
    A[i_1,1]*A[8,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[8,16])+
    A[i_1,1]*A[9,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[9,16])+
    A[i_1,1]*A[10,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[10,16])+
    A[i_1,1]*A[11,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[11,16])+
    A[i_1,1]*A[12,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[12,16])+
    A[i_1,1]*A[13,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[13,16])+
    A[i_1,1]*A[14,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[14,16])+
    A[i_1,1]*A[15,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[15,16])+
    A[i_1,1]*A[16,j_1]*A[i_2,1]*A[16,j_2]*(CovarLam[16,16])


  return(covA_val)
}



CovAfunction_16<-function(A,k,CovarLam,j_1,j_2,n_r){
  covA_val<-matrix(c(CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=1),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=1),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=2),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=2),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=3),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=3),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=4),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=4),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=5),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=5),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=6),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=6),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=7),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=7),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=8),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=8),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=9),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=9),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=10),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=10),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=11),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=11),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=12),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=12),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=13),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=13),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=14),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=14),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=15),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=15),

                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=1,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=2,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=3,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=4,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=5,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=6,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=7,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=8,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=9,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=10,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=11,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=12,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=13,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=14,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=15,i_2=16),
                     CovAfunc_helper_16(A,k,CovarLam,j_1,j_2,n_r,i_1=16,i_2=16)),
                   nrow=k,ncol=k)
  return(covA_val)
}



VarB_valid_16<-function(SigmaBeta,CovarLam,CorrectA,Beta_Star,n_r){
  k<-nrow(SigmaBeta)
  Ainv<-solve(CorrectA)
  ASigBR<-t(Ainv)%*%SigmaBeta%*%Ainv
  varbeta<-matrix(NA,nrow=k,ncol=k)
  for (j_1 in 1:k){
    for (j_2 in 1:k){
      varbeta[j_1,j_2]<-ASigBR[j_1,j_2]+t(Beta_Star)%*%CovAfunction_16(Ainv,k,CovarLam,j_1,j_2,n_r)%*%Beta_Star
    }
  }
  return(varbeta)
}

#Read in the simulated data (users will have to change this path to get code to run)
sim_data_wide<-read.csv("U://Paper 1//Final Code//Simulated_Data_Example_Wide_Form.csv",header=TRUE,sep=",",check.names = FALSE)

#Perform calibration - get elements needed to correct for exposure error
correction_energy<-lm(Biom_Energy_m~log_FFQ_energy_c+bmi_c+age_c+black+hispanic+other+income2+income3+income5+income6+texpwk+educ3+educ5+htn+alcohol1_nondrinker+alcohol3_moderateheavy,data=sim_data_wide)

lm_energy<-summary(correction_energy)
nbeta_e<-length(lm_energy$coefficients[-1,1])
eta1_energy<-rbind(as.matrix(t(lm_energy$coefficients[-1,1])),cbind(matrix(0,nrow=nbeta_e-1,ncol=1),diag(nbeta_e-1)))

covla_e<-vcov(lm_energy)
covla_fin_e<-covla_e[-c(1),-c(1)]

#Assign sensitivity and specificity
sensitivity<-0.61
specificity<-0.995
negpred<-0.96

#Convert data from wide to long form
Sim_data_long<-melt(sim_data_wide, value.name = "result",id.vars=c("log_FFQ_energy_c","bmi_c","age_c", "black", "hispanic", "other","income2", "income3","income5","income6","texpwk",
                                             "educ3", "educ5", "Biom_Energy_m","htn", "alcohol1_nondrinker","alcohol3_moderateheavy","ID","os_dm_age"))
Sim_data_long$testtime<-as.numeric(Sim_data_long$variable)

#Sort by ID so we can see all subject's data organized together for each visit time
Sim_data_long<-Sim_data_long[order(Sim_data_long$ID),]

#Discrete - naive - ignore error in Y and X
#Function that deletes all positives after first self-report
after_first_pos <- function(x){
  npos<-cumsum(x==1)
  (npos==0) | (npos==1 & x==1)
}

keep_naive<-unlist(tapply(Sim_data_long$result,Sim_data_long$ID,after_first_pos))
datafinal_naive<-Sim_data_long[keep_naive,]

datafinal_naive$testtime<-as.factor(datafinal_naive$testtime)
fit_naive_energy_strat=glm(result~testtime+testtime*os_dm_age+log_FFQ_energy_c+bmi_c+age_c+black+hispanic+other+income2+income3+income5+income6+
                             texpwk+educ3+educ5+htn+alcohol1_nondrinker+alcohol3_moderateheavy,family=binomial(link="cloglog"),data=datafinal_naive)

cbind(HR=exp(fitsum_naive_strat_e$coefficients[14,1])^log(1.2),Lower=exp(fitsum_naive_strat_e$coefficients[14,1]-1.96*fitsum_naive_strat_e$coefficients[14,2])^log(1.2),Upper=exp(fitsum_naive_strat_e$coefficients[14,1]+1.96*fitsum_naive_strat_e$coefficients[14,2])^log(1.2))

fitsum_naive_strat_e<-summary(fit_naive_energy_strat)

nend<-ntest+(nstrat-1)
coeffnum<-nend+1
nstart_e<-(length(fitsum_naive_strat_e$coefficients[,1])-((ntest-1)*(nstrat-1)))+1
param0b_e<-fitsum_naive_strat_e$coefficients[-c(1:nend,nstart_e:length(fitsum_naive_strat_e$coefficients[,1])),1]

datas1 <- Sim_data_long[which(Sim_data_long$os_dm_age=="OS_age50s"),]
datas2 <- Sim_data_long[which(Sim_data_long$os_dm_age=="OS_age60s"),]
datas3 <- Sim_data_long[which(Sim_data_long$os_dm_age=="OS_age70s"),]
datas4 <- Sim_data_long[which(Sim_data_long$os_dm_age=="DM_age50s"),]
datas5 <- Sim_data_long[which(Sim_data_long$os_dm_age=="DM_age60s"),]
datas6 <- Sim_data_long[which(Sim_data_long$os_dm_age=="DM_age70s"),]

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

subject5<-datas5$ID
testtime5<-datas5$testtime
result5<-datas5$result

subject6<-datas6$ID
testtime6<-datas6$testtime
result6<-datas6$result

formula_energy=result~log_FFQ_energy_c+bmi_c+age_c+black+hispanic+other+income2+income3+income5+income6+texpwk+educ3+educ5+htn+alcohol1_nondrinker+alcohol3_moderateheavy

negpred = 0.96
betai_e = param0b_e

initsurv = 0.1
id1 <- eval(substitute(subject1), datas1, parent.frame())
time1 <- eval(substitute(testtime), datas1, parent.frame())
result1 <- eval(substitute(result1), datas1, parent.frame())
ord1 <- order(id1, time1)

id2 <- eval(substitute(subject2), datas2, parent.frame())
time2 <- eval(substitute(testtime), datas2, parent.frame())
result2 <- eval(substitute(result2), datas2, parent.frame())
ord2 <- order(id2, time2)

id3 <- eval(substitute(subject3), datas3, parent.frame())
time3 <- eval(substitute(testtime), datas3, parent.frame())
result3 <- eval(substitute(result3), datas3, parent.frame())
ord3 <- order(id3, time3)

id4 <- eval(substitute(subject4), datas4, parent.frame())
time4 <- eval(substitute(testtime), datas4, parent.frame())
result4 <- eval(substitute(result4), datas4, parent.frame())
ord4 <- order(id4, time4)

id5 <- eval(substitute(subject5), datas5, parent.frame())
time5 <- eval(substitute(testtime), datas5, parent.frame())
result5 <- eval(substitute(result5), datas5, parent.frame())
ord5 <- order(id5, time5)

id6 <- eval(substitute(subject6), datas6, parent.frame())
time6 <- eval(substitute(testtime), datas6, parent.frame())
result6 <- eval(substitute(result6), datas6, parent.frame())
ord6 <- order(id6, time6)

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

if (is.unsorted(ord5)) {
  id5 <- id5[ord5]
  time5 <- time5[ord5]
  result5 <- result5[ord5]
  data5 <- data5[ord5, ]}

if (is.unsorted(ord6)) {
  id6 <- id6[ord6]
  time6 <- time6[ord6]
  result6 <- result6[ord6]
  data6 <- data6[ord6, ]}

utime1 <- sort(unique(time1))
utime2 <- sort(unique(time2))
utime3 <- sort(unique(time3))
utime4 <- sort(unique(time4))
utime5 <- sort(unique(time5))
utime6 <- sort(unique(time6))

timen01 <- (time1 != 0)
timen02 <- (time2 != 0)
timen03 <- (time3 != 0)
timen04 <- (time4 != 0)
timen05 <- (time5 != 0)
timen06 <- (time6 != 0)


Dm1 <- dmat(id1[timen01], time1[timen01], result1[timen01], sensitivity,
            specificity, negpred)
Dm2 <- dmat(id2[timen02], time2[timen02], result2[timen02], sensitivity,
            specificity, negpred)
Dm3 <- dmat(id3[timen03], time3[timen03], result3[timen03], sensitivity,
            specificity, negpred)
Dm4 <- dmat(id4[timen04], time4[timen04], result4[timen04], sensitivity,
            specificity, negpred)
Dm5 <- dmat(id5[timen05], time5[timen05], result5[timen05], sensitivity,
            specificity, negpred)
Dm6 <- dmat(id6[timen06], time6[timen06], result6[timen06], sensitivity,
            specificity, negpred)

J1 <- ncol(Dm1) - 1
J2 <- ncol(Dm2) - 1
J3 <- ncol(Dm3) - 1
J4 <- ncol(Dm4) - 1
J5 <- ncol(Dm5) - 1
J6 <- ncol(Dm6) - 1


nsub1 <- nrow(Dm1)
nsub2 <- nrow(Dm2)
nsub3 <- nrow(Dm3)
nsub4 <- nrow(Dm4)
nsub5 <- nrow(Dm5)
nsub6 <- nrow(Dm6)

lami1 <- log(-log(seq(1, initsurv, length.out = J1 + +1)[-1]))
lami1 <- c(lami1[1], diff(lami1))
lami2 <- log(-log(seq(1, initsurv, length.out = J2 + +1)[-1]))
lami2 <- c(lami2[1], diff(lami2))
lami3 <- log(-log(seq(1, initsurv, length.out = J3 + +1)[-1]))
lami3 <- c(lami3[1], diff(lami3))
lami4 <- log(-log(seq(1, initsurv, length.out = J4 + +1)[-1]))
lami4 <- c(lami4[1], diff(lami4))
lami5 <- log(-log(seq(1, initsurv, length.out = J5 + +1)[-1]))
lami5 <- c(lami5[1], diff(lami5))
lami6 <- log(-log(seq(1, initsurv, length.out = J6 + +1)[-1]))
lami6 <- c(lami6[1], diff(lami6))

tosurv <- function(x) exp(-exp(cumsum(x)))

lowlam1 <- c(-Inf, rep(0, J1 - 1))
lowlam2 <- c(-Inf, rep(0, J2 - 1))
lowlam3 <- c(-Inf, rep(0, J3 - 1))
lowlam4 <- c(-Inf, rep(0, J4 - 1))
lowlam5 <- c(-Inf, rep(0, J5 - 1))
lowlam6 <- c(-Inf, rep(0, J6 - 1))

Xmat1_e <- model.matrix(formula_energy, data = datas1)[, -1, drop = F]
Xmat2_e <- model.matrix(formula_energy, data = datas2)[, -1, drop = F]
Xmat3_e <- model.matrix(formula_energy, data = datas3)[, -1, drop = F]
Xmat4_e <- model.matrix(formula_energy, data = datas4)[, -1, drop = F]
Xmat5_e <- model.matrix(formula_energy, data = datas5)[, -1, drop = F]
Xmat6_e <- model.matrix(formula_energy, data = datas6)[, -1, drop = F]

beta.nm1_e <- colnames(Xmat1_e)
nbeta1_e <- ncol(Xmat1_e)
uid1_e<- getrids(id1, nsub1)
Xmat1_e <- Xmat1_e[uid1_e, , drop = F]

beta.nm2_e <- colnames(Xmat2_e)
nbeta2_e <- ncol(Xmat2_e)
uid2_e<- getrids(id2, nsub2)
Xmat2_e <- Xmat2_e[uid2_e, , drop = F]

beta.nm3_e <- colnames(Xmat3_e)
nbeta3_e <- ncol(Xmat3_e)
uid3_e<- getrids(id3, nsub3)
Xmat3_e <- Xmat3_e[uid3_e, , drop = F]

beta.nm4_e <- colnames(Xmat4_e)
nbeta4_e <- ncol(Xmat4_e)
uid4_e<- getrids(id4, nsub4)
Xmat4_e <- Xmat4_e[uid4_e, , drop = F]

beta.nm5_e <- colnames(Xmat5_e)
nbeta5_e <- ncol(Xmat5_e)
uid5_e<- getrids(id5, nsub5)
Xmat5_e <- Xmat5_e[uid5_e, , drop = F]

beta.nm6_e <- colnames(Xmat6_e)
nbeta6_e <- ncol(Xmat6_e)
uid6_e<- getrids(id6, nsub6)
Xmat6_e <- Xmat6_e[uid6_e, , drop = F]

parmall_e <- c(lami1,lami2,lami3,lami4,lami5,lami6,param0b_e)
nbeta_e<-length(param0b_e)

loglikStrat <- function(parmsS,Dm1,Dm2,Dm3,Dm4,Dm5,Dm6,Xmat1,Xmat2,Xmat3,Xmat4,Xmat5,Xmat6){
  parmi1<-parmsS[c(1:length(testtimes),(length(testtimes)*6+1):length(parmsS))]
  parmi2<-parmsS[c((length(testtimes)+1):(length(testtimes)*2),(length(testtimes)*6+1):length(parmsS))]
  parmi3<-parmsS[c((length(testtimes)*2+1):(length(testtimes)*3),(length(testtimes)*6+1):length(parmsS))]
  parmi4<-parmsS[c((length(testtimes)*3+1):(length(testtimes)*4),(length(testtimes)*6+1):length(parmsS))]
  parmi5<-parmsS[c((length(testtimes)*4+1):(length(testtimes)*5),(length(testtimes)*6+1):length(parmsS))]
  parmi6<-parmsS[c((length(testtimes)*5+1):(length(testtimes)*6),(length(testtimes)*6+1):length(parmsS))]

  strataLL<-loglikC(parm=parmi1,Dm=Dm1,Xmat=Xmat1)+loglikC(parm=parmi2,Dm=Dm2,Xmat=Xmat2)+loglikC(parm=parmi3,Dm=Dm3,Xmat=Xmat3)+
    loglikC(parm=parmi4,Dm=Dm4,Xmat=Xmat4)+loglikC(parm=parmi5,Dm=Dm5,Xmat=Xmat5)+loglikC(parm=parmi6,Dm=Dm6,Xmat=Xmat6)
  return(strataLL)
}

q_energy <- optim(parmall_e, loglikStrat, lower = c(rep(lowlam1,6), rep(-Inf, nbeta_e)), Dm1 = Dm1, Dm2 = Dm2, Dm3 = Dm3, Dm4 = Dm4, Dm5 = Dm5,Dm6 = Dm6,
                  Xmat1 = Xmat1_e, Xmat2 = Xmat2_e, Xmat3 = Xmat3_e, Xmat4 = Xmat4_e,Xmat5 = Xmat5_e,Xmat6 = Xmat6_e, method = "L-BFGS-B", hessian = T)

totJ<-J1+J2+J3+J4+J5+J6

cov_e <- as.matrix(solve(q_energy$hessian)[-(1:totJ), -(1:totJ)])
rownames(cov_e) <- colnames(cov_e) <- beta.nm1_e
beta_fit_e <- q_energy$par[-(1:totJ)]

corrected_eta_beta_e<-t(as.matrix(beta_fit_e))%*%solve(eta1_energy)
corrected_eta_beta_e<-as.data.frame(t(corrected_eta_beta_e))

rownames(corrected_eta_beta_e) <-beta.nm1_e

#SE and CI for Energy
VarBetaCorrected_e<-VarB_valid_16(cov_e,covla_fin_e,eta1_energy,beta_fit_e,nvalid)
SDBeta_e<-sqrt(VarBetaCorrected_e[1,1])

#Naive in  Y Only - Correct for error in x
corrected_eta_beta_e_x<-t(as.matrix(param0b_e))%*%solve(eta1_energy)
corrected_eta_beta_e_x<-as.data.frame(t(corrected_eta_beta_e_x))

cov_cloglog_e<-(vcov(fitsum_naive_strat_e)[-c(1:nend,nstart_e:length(fitsum_naive_strat_e$coefficients[,1])),-c(1:nend,nstart_e:length(fitsum_naive_strat_e$coefficients[,1]))])

VarBetaCorrected_e_x<-VarB_valid_16(cov_cloglog_e,covla_fin_e,eta1_energy,param0b_e,nvalid)
SDBeta_e_x<-sqrt(VarBetaCorrected_e_x[1,1])

#Put results in table
corrected_HR_myMethod<-cbind(HR=exp(corrected_eta_beta_e[1,])^log(1.2),Lower=exp(corrected_eta_beta_e[1,]-1.96*SDBeta_e)^log(1.2),Upper=exp(corrected_eta_beta_e[1,]+1.96*SDBeta_e)^log(1.2))
naive_HR<-cbind(HR=exp(fitsum_naive_strat_e$coefficients[14,1])^log(1.2),Lower=exp(fitsum_naive_strat_e$coefficients[14,1]-1.96*fitsum_naive_strat_e$coefficients[14,2])^log(1.2),Upper=exp(fitsum_naive_strat_e$coefficients[14,1]+1.96*fitsum_naive_strat_e$coefficients[14,2])^log(1.2))
corrected_HR_x<-cbind(HR=exp(corrected_eta_beta_e_x[1,])^log(1.2),Lower=exp(corrected_eta_beta_e_x[1,]-1.96*SDBeta_e_x)^log(1.2),Upper=exp(corrected_eta_beta_e_x[1,]+1.96*SDBeta_e_x)^log(1.2))

energy_results<-rbind(naive_HR,corrected_HR_x,corrected_HR_myMethod)
rownames(energy_results)<-c("Naive HR (95% CI)","Corrected HR for x only (95% CI)","Proposed Method HR (95% CI)")

round(energy_results,3)


