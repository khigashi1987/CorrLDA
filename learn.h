/*
    learn.h
*/
#ifndef LEARN_H
#define LEARN_H
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "feature.h"

//#define RANDOM ((double)rand()/(double)RAND_MAX)

int sampling_multinomial(gsl_rng *r, double *p, double *cum_sum_p, int len_p);

extern void corrlda_learn(document *data,
        document *tag,
        double *alpha,
        double beta, double gamma,
        int nclass, int nlex, int ntlex, int dlenmax, int tlenmax,
        int maxiter,
        double **phi, double **theta, double **psi,
        int **n_mz, int **n_zw, int **m_mz, int **m_zt,
        FILE *likp, FILE *hyperp,
        unsigned long int random_seed);

#endif
