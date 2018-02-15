/*
    likelihood.h
*/
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

extern double loglikelihood(int **n_mz, int **n_zw, int *n_m,
        int **m_mz, int **m_zt, int *m_z,
        int nclass, int nlex, int ntlex, int ndocs, int Mall,
        double *alpha, double beta, double gamma);

#endif
