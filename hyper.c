/*
    hyper.c
*/
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_sf_psi.h>

static double update_hyper(double param, int *vec, int **mat, int nrow, int ncol){
    double numer = 0.0;
    double denom = 0.0;
    int i, j;
    
    for(i = 0;i < nrow;i++)
        for(j = 0;j < ncol;j++)
            numer += gsl_sf_psi((double)mat[i][j] + param);
    numer -= (double)nrow * (double)ncol * gsl_sf_psi(param);
    
    for(i = 0;i < nrow;i++)
        denom += (double)ncol * gsl_sf_psi((double)vec[i] + param * (double)ncol);
    denom -= (double)nrow * (double)ncol * gsl_sf_psi(param * (double)ncol);
    
    return param * numer / denom;
}

void update_alpha(double *alpha, int *n_m, int **n_mz, int ndocs, int nclass){
    double sum_alpha = 0.0;
    double numer = 0.0;
    double denom = 0.0;
    int d, k;
    for(k = 0;k < nclass;k++)
        sum_alpha += alpha[k];
    for(d = 0;d < ndocs;d++)
        denom += gsl_sf_psi((double)n_m[d] + sum_alpha);
    denom -= (double)ndocs * gsl_sf_psi(sum_alpha);
    
    for(k = 0;k < nclass;k++){
        numer = 0.0;
        for(d = 0;d < ndocs;d++)
            numer += gsl_sf_psi((double)n_mz[d][k] + alpha[k]);
        numer -= (double)ndocs * gsl_sf_psi(alpha[k]);
        alpha[k] = alpha[k] * numer / denom;
    }
}

double update_beta(double beta, int *n_z, int **n_zw, int nclass, int nlex){
    double new_beta = update_hyper(beta, n_z, n_zw, nclass, nlex);
    return new_beta;
}

double update_gamma(double gamma, int *m_z, int **m_zt, int nclass, int ntlex){
    double new_gamma = update_hyper(gamma, m_z, m_zt, nclass, ntlex);
    return new_gamma;
}
