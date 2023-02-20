//
//  main.c
//  test
//
//  Created by YinYiyi on 4/25/18.
//  Copyright © 2018 Yiyi Yin. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

// sign function
int sign (double v){
    if (v > 0){
        return 1;
    }
    else if (v == 0){
        return 0;
    }
    else{
        return -1;
    }
}

// soft thresholding function
// lambda >= 0
double softT (double v, double lambda){
    double posi_part = fabs(v) - lambda;
    if (posi_part > 0){
        return sign(v)*posi_part;
    }
    else{
        return 0;
    }
}

// crossprod of two vectors
double crossprod (double * vec1, double * vec2, int p){
    double ans = 0;
    for (int i = 0; i < p; i++){
        ans += vec1[i] * vec2[i];
    }
    return ans;
}

// maximum of three variables
double max (double a, double b, double c){
    double ans = (a > b ? a : b);
    ans = (ans > c ? ans : c);
    return ans;
}

//// create sequence with equal difference
//// min_v <= max_v
//// if increase == true, start with min_v and end with max_v
//// if increase == false, start with max_v and end with min_v
//void seq (double min_v, double max_v, int length, bool increase,
//          double * sequence){
//    double diff = (max_v - min_v)/(length - 1), curr = (increase? min_v: max_v);
//    for (int i = 0; i < length; i++){
//        sequence[i] = curr;
//        curr = (increase? curr + diff: curr - diff);
//    }
//}

// matrix multiplication
// mat1 is a symmetric matrix
void multiply (double * mat1, double * mat2, int p,
               double * mat){
    int i, j;
    for (i = 0; i < p; i++){
        for (j = 0; j < p; j++){
            mat[j*p+i] = crossprod(& mat1[i*p], & mat2[j*p], p);
        }
    }
}

//// D-trace loss function
//// 1/2*tr(Theta^T*Sigma*Theta)-tr(Theta)
//double d_loss (double * Sigma, int p, double * Theta){
//    double ans = 0, * tmp = malloc(p*p*sizeof(double));
//    multiply(Sigma, Theta, p, tmp);
//    for (int i = 0; i < p; i++){
//        ans += crossprod(& Theta[i*p], & tmp[i*p], p)/2 - Theta[i*(p+1)];
//    }
//    free(tmp);
//    return ans;
//}

// Frobenius norm of the difference of two vectors
double f_norm2 (double * v1, double * v2, int p){
    double ans = 0;
    for (int i = 0; i < p; i++){
        ans += (v1[i]-v2[i])*(v1[i]-v2[i]);
    }
    return sqrt(ans);
}

// Frobenius norm of the vector
double f_norm1 (double * v1, int p){
    double ans = 0;
    for (int i = 0; i < p; i++){
        ans += v1[i]*v1[i];
    }
    return sqrt(ans);
}

// get matrix diagonal elements
void getdiag (double * mat, int p,
              double * sig){
    for (int i = 0; i < p; i++){
        sig[i] = mat[i*(p+1)];
    }
}

//// find the lambda max
//// for step 2 estimator
//// weight is symmetric
//double lamMax (double * Sigma, int p, double alpha, double * weight){
//    double * sig = malloc(p*sizeof(double));
//    getdiag (Sigma, p, sig);
//    double tmp, ans = 0;
//    int i, j;
//    for (j = 0; j < p-1; j++){
//        for (i = j+1; i < p; i++){
//            tmp = fabs(Sigma[j*p+i]*(1/(sig[i]+alpha)+1/(sig[j]+alpha))/2/weight[j*p+i]);
//            if (tmp > ans){
//                ans = tmp;
//            }
//        }
//    }
//    free(sig);
//    return ans;
//}

// check gradient of inactive elements
// for ith column of step 1 estimator
// inactive elements that violate gradient will be added into active set
// num of active set will also be changed if violations exist
void check_inactive (double * Sigma, int p, int i, double * theta_vec,
                     double lambda, double error,
                     int * act_set, int * num_act){
    int * tmp_set = malloc(p*sizeof(int)), j;
    double grad, tmp;
    memset(tmp_set, 0, p*sizeof(int));
    for (j = 0; j < * num_act; j++){
        tmp_set[act_set[j]] = 1;
    }
    for (j = 0; j < p; j++){
        if (tmp_set[j] == 1) continue;
        tmp = crossprod (& Sigma[j*p], theta_vec, p);
        if (j == i){
            tmp -= 1;
        }
        grad = fabs(tmp) - lambda;
        if (grad > error){
            act_set[* num_act] = j;
            * num_act += 1;
        }
    }
    free(tmp_set);
}

// check gradient of active elements
// for ith column of step 1 estimator
// return whether or not all active elements satisfy KKT condition
bool check_active (double * Sigma, int p, int i, double * theta_vec,
                  double lambda, double alpha, double error,
                  int * act_set, int num_act){
    double grad, tmp;
    int ind;
    for (int j = 0; j < num_act; j++){
        ind = act_set[j];
        tmp = crossprod (& Sigma[ind*p], theta_vec, p) + alpha * theta_vec[ind];
        if (ind == i){
            tmp -= 1;
        }
        if (theta_vec[ind] > 0){
            grad = fabs(tmp + lambda);
        }
        else if (theta_vec[ind] < 0){
            grad = fabs(tmp - lambda);
        }
        else{
            grad = fabs(tmp) - lambda;
        }
        if (grad > error){
            return false;
        }
    }
    return true;
}

// solution path
// for ith column of step 1 estimator
// lamList in decreasing order
void step1column (double * Sigma, int p, int i,
                  double * lamList, int lamLength, double alpha,
                  int maxit, double error,
                  double * theta_vec_list, int * niterList){
    double * theta_vec = malloc(p*sizeof(double)), * old_theta_vec = malloc(p*sizeof(double));
    double lam, tmp_st;
    int t, j, ind, time;
    int num_act = 0, old_num_act, * act_set = malloc(p*sizeof(int));
    memset(theta_vec, 0, p*sizeof(double));
    for (t = 0; t < lamLength; t++){
        lam = lamList[t];
        time = 0;
        while (true){
            old_num_act = num_act;
            check_inactive(Sigma, p, i, theta_vec, lam, error, act_set, &num_act);
            if (num_act == old_num_act){
                if (check_active(Sigma, p, i, theta_vec, lam, alpha, error, act_set, num_act)){
                    memcpy(& theta_vec_list[t*p], theta_vec, p*sizeof(double));
                    niterList[t] = time;
                    break;
                }
            }
            while (time <= maxit){
                memcpy(old_theta_vec, theta_vec, p*sizeof(double));
                time++;
                for (j = 0; j < num_act; j++){
                    ind = act_set[j];
                    tmp_st = Sigma[ind*(p+1)]*theta_vec[ind] - crossprod(& Sigma[ind*p], theta_vec, p);
                    if (ind == i){
                        tmp_st += 1;
                    }
                    theta_vec[ind] = softT(tmp_st, lam)/(Sigma[ind*(p+1)] + alpha);
                }
                if (f_norm2(theta_vec, old_theta_vec, p) < max(1, f_norm1(theta_vec, p), f_norm1(old_theta_vec, p)) * error){
                    break;
                }
            }
            if (time > maxit){
                memcpy(& theta_vec_list[t*p], theta_vec, p*sizeof(double));
                niterList[t] = time;
                break;
            }
        }
    }
    free(theta_vec);
    free(old_theta_vec);
    free(act_set);
}

// solution path
// of step 1 estimator
// lamList in decreasing order
void step1matrix (double * Sigma, int * p, double * lamList, int * lamLength, double * alpha,
                  int * maxit, double * error,
                  double * precisionList, int * niterList){
    int p0 = p[0], lamLength0 = lamLength[0], maxit0 = maxit[0];
    double alpha0 = alpha[0], error0 = error[0];
//    memset (niterList, 0, p0 * sizeof(int));
    double * column_list = malloc(lamLength0 * p0 * sizeof(double));
    int * column_niter_list = malloc(lamLength0 * sizeof(int));
    int i, j;
    for (i = 0; i < p0; i++){
        step1column (Sigma, p0, i, lamList, lamLength0, alpha0, maxit0, error0, column_list, column_niter_list);
        for (j = 0; j < lamLength0; j++){
            memcpy (& precisionList[j*p0*p0 + i*p0], & column_list[j*p0], p0 * sizeof(double));
            if (column_niter_list[j] > niterList[j]) {
                niterList[j] = column_niter_list[j];
            }
        }
    }
    free(column_list);
    free(column_niter_list);
}



// check gradient of inactive elements
// for step 1 estimator matrix
// indexs of inactive elements that violate gradient will be added into rowind, colind
// num of active set will also be changed if violations exist
void check_inactive2 (double * Sigma, double * theta, int p,
                      double lambda, double alpha, double * weight, double error,
                      int *rowind, int *colind, int *num_act){
    int * ma = malloc(p*p*sizeof(int)), i, j;
    double grad = 0, tmp;
    memset (ma, 0, p*p*sizeof(int));
    for (i = 0; i < *num_act; i++){
        ma[colind[i]*p+rowind[i]] = 1;
    }
    for (j = 0; j < p-1; j++){
        for (i = j+1; i < p; i++){
            if (ma[j*p+i] == 1) continue;
            tmp = (crossprod (&Sigma[i*p], &theta[j*p], p) + crossprod (&Sigma[j*p], &theta[i*p], p))/2;
            grad = fabs(tmp+alpha*theta[j*p+i])-lambda*weight[j*p+i];
            if (grad > error){
                rowind[*num_act] = i;
                colind[*num_act] = j;
                *num_act += 1;
            }
        }
    }
    free(ma);
}

// check gradient of active elements
// for step 1 estimator matrix
// return whether or not all active elements satisfy KKT condition
bool check_active2 (double * Sigma, double * theta, int p,
                    double lambda, double alpha, double * weight, double error,
                    int *rowind, int *colind, int num_act){
    double grad = 0, tmp;
    int i, j;
    for (int ind = 0; ind < num_act; ind++){
        i = rowind[ind];
        j = colind[ind];
        tmp = (crossprod(&Sigma[i*p], &theta[j*p], p) + crossprod(&Sigma[j*p], &theta[i*p], p))/2;
        tmp += alpha*theta[j*p+i];
        if (i == j) {
            tmp -= 1;
        }
        if (theta[j*p+i] > 0){
            grad = fabs(tmp+lambda*weight[j*p+i]);
        }
        else if (theta[j*p+i] < 0){
            grad = fabs(tmp-lambda*weight[j*p+i]);
        }
        else {
            grad = fabs(tmp)-lambda*weight[j*p+i];
        }
        if (grad > error){
            return false;
        }
    }
    return true;
}

//// single solution for step 2 estimator
//// input: Sigma, lambda, alpha, weight, maxit, error
//// output(change): precision, niter
//void precision_single (double * Sigma, int * p, double * lambda, double * alpha, double * weight,
//                       int * maxit, double * error,
//                       double * precision, int * niter){
//    int p0 = p[0], maxit0 = maxit[0];
//    double error0 = error[0], alpha0 = alpha[0];
//    double *sig = malloc(p0*sizeof(double));
//    double lam = lambda[0];
//    double * theta = malloc(p0*p0*sizeof(double)), * old_theta = malloc(p0*p0*sizeof(double));
//    int *actrowset = malloc(p0*p0*sizeof(int)), *actcolset = malloc(p0*p0*sizeof(int));
//    int i, j, time = 0, idx, idy, num_act = p0, old_num_act = p0;
//    getdiag(Sigma, p0, sig);
//    memset(theta, 0, p0*p0*sizeof(double));
//    for (i = 0; i < p0; i++){
//        theta[i*(p0+1)] = 1/(sig[i]+alpha0);
//    }
//    for (i = 0; i < p0 ; i++) {
//        actrowset[i] = i;
//        actcolset[i] = i;
//    }
//    while (true){
//        old_num_act = num_act;
//        check_inactive2 (Sigma, theta, p0, lam, alpha0, weight, error0, actrowset, actcolset, &num_act);
//        if (num_act == old_num_act){
//            if (check_active2 (Sigma, theta, p0, lam, alpha0, weight, error0, actrowset, actcolset, num_act)){
//                memcpy(& precision[0], theta, p0*p0*sizeof(double));
//                niter[0] = time;
//                break;
//            }
//        }
//        while (time <= maxit0){
//            time++;
//            memcpy(old_theta, theta, p0*p0*sizeof(double));
//            for (j = 0; j < num_act; j++){
//                idx = actrowset[j];
//                idy = actcolset[j];
//                if (idx == idy){
//                    theta[idx*(p0+1)] = (1 - crossprod(&theta[idx*p0], &Sigma[idx*p0], p0) + sig[idx]*theta[idx*(p0+1)]) / (sig[idx]+alpha0);
//                }
//                else{
//                    theta[idy*p0+idx] = softT (-crossprod(&Sigma[idx*p0], &theta[idy*p0], p0) - crossprod(&Sigma[idy*p0], &theta[idx*p0], p0) + theta[idy*p0+idx]*(sig[idx]+sig[idy]), 2*lam*weight[idy*p0+idx]) / (sig[idx]+sig[idy]+2*alpha0);
//                    theta[idx*p0+idy] = theta[idy*p0+idx];
//                }
//            }
//            if (f_norm2(theta, old_theta, p0*p0)/max(1, f_norm1(theta, p0*p0), f_norm1(old_theta, p0*p0)) < error0){
//                break;
//            }
//        }
//        if (time > maxit0){
//            memcpy(& precision[0], theta, p0*p0*sizeof(double));
//            niter[0] = time;
//            break;
//        }
//    }
//    free(sig);
//    free(theta);
//    free(old_theta);
//    free(actrowset);
//    free(actcolset);
//}

// solution path for step 2 estimator
// input: Sigma, lamList, N, alpha, weight, maxit, error
// output(change): precisionList, niterList
void precision (double * Sigma, int * p, double * lamList, int * lamLength, double * alpha, double * weight,
                int * maxit, double * error,
                double * precisionList, int * niterList){
    int p0 = p[0], lamLength0 = lamLength[0], maxit0 = maxit[0];
    double error0 = error[0], alpha0 = alpha[0];
    double * sig = malloc(p0*sizeof(double));
    double * theta = malloc(p0*p0*sizeof(double)), * old_theta = malloc(p0*p0*sizeof(double));
    int *actrowset = malloc(p0*p0*sizeof(int)), *actcolset = malloc(p0*p0*sizeof(int));
    int i, j, time, idx, idy, num_act = p0, old_num_act;
    double lam;
//    if (lamList[0] == 0){
//        double laMax = lamMax(Sigma, p0, alpha0, weight);
//        seq(laMax/N0, laMax, N0, false, lamList);
//    }
    getdiag (Sigma, p0, sig);
    memset(theta, 0, p0*p0*sizeof(double));
    for (i = 0; i < p0; i++){
        theta[i*(p0+1)] = 1/(sig[i]+alpha0);
    }
    for (i = 0; i < p0 ; i++) {
        actrowset[i] = i;
        actcolset[i] = i;
    }
    for (i = 0; i < lamLength0; i++){
        lam = lamList[i];
        time = 0;
        while (true){
            old_num_act = num_act;
            check_inactive2(Sigma, theta, p0, lam, alpha0, weight, error0, actrowset, actcolset, &num_act);
            if (num_act == old_num_act){
                if (check_active2(Sigma, theta, p0, lam, alpha0, weight, error0, actrowset, actcolset, num_act)){
                    memcpy(& precisionList[i*p0*p0], theta, p0*p0*sizeof(double));
                    niterList[i] = time;
                    break;
                }
            }
            while (time <= maxit0){
                time++;
                memcpy(old_theta, theta, p0*p0*sizeof(double));
                for (j = 0; j < num_act; j++){
                    idx = actrowset[j];
                    idy = actcolset[j];
                    if (idx == idy){
                        theta[idx*(p0+1)] = (1 - crossprod (&theta[idx*p0], &Sigma[idx*p0], p0) + sig[idx]*theta[idx*(p0+1)]) / (sig[idx]+alpha0);
                    }
                    else{
                        theta[idy*p0+idx] = softT (-crossprod (&Sigma[idx*p0], &theta[idy*p0], p0) - crossprod(&Sigma[idy*p0], &theta[idx*p0], p0) + theta[idy*p0+idx]*(sig[idx]+sig[idy]), 2*lam*weight[idy*p0+idx]) / (sig[idx]+sig[idy]+2*alpha0);
                        theta[idx*p0+idy] = theta[idy*p0+idx];
                    }
                }
                if (f_norm2(theta, old_theta, p0*p0)/max(1, f_norm1(theta, p0*p0), f_norm1(old_theta, p0*p0)) < error0){
                    break;
                }
            }
            if (time > maxit0){
                memcpy(& precisionList[i*p0*p0], theta, p0*p0*sizeof(double));
                niterList[i] = time;
                break;
            }
        }
    }
    free(sig);
    free(theta);
    free(old_theta);
    free(actrowset);
    free(actcolset);
}



int main(){
    
    return 0;
}
