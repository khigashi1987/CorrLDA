/*
    learn.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "learn.h"
#include "feature.h"
#include "imatrix.h"
#include "dmatrix.h"
#include "util.h"
#include "likelihood.h"
#include "hyper.h"

void corrlda_learn(document *data, document *tag, double *alpha, double beta, double gamma, int nclass, int nlex, int ntlex, int dlenmax, int tlenmax, int maxiter, double **phi, double **theta, double **psi, int **n_mz, int **n_zw, int **m_mz, int **m_zt, FILE *likp, FILE *hyperp, unsigned long int random_seed){
    document *dp;
    document *tp;
    int ndocs;
    int ntdocs;
    int *n_m;
    int *n_z;
    int *m_m;
    int *m_z;
    int *m_0t;
    int ***topics;
    int ***tagtopics;
    int word_index;
    int word_num;
    int tag_index;
    int tag_num;
    double sum_alpha;
    double *left;
    double *center;
    double *right;
    double *p_z;
    double *cum_sum_p_z;
    double sum_p_z, sum_p_r;
    double p_r0, p_r1;
    double lik;
    double **temp_phi;
    double **temp_theta;
    double **temp_psi;
    bool infinite_prob;
    int Mall;
    int z;
    int it;
    int m, w, t, i, j, k;
    const gsl_rng_type *T;
    gsl_rng *r;
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    // count data length
    for(dp = data, ndocs = 0;(dp->len) != -1;dp++, ndocs++)
        ;
    for(dp = tag, ntdocs = 0;(dp->len) != -1;dp++, ntdocs++)
        ;
    if(ndocs != ntdocs){
        fprintf(stderr,"corrlda_learn:: size of documents and tagsets are different.\n");
        return;
    }
    
    // initialize buffers
    if((n_m = calloc(ndocs,sizeof(int))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate n_m.\n");
        return;
    }
    if((n_z = calloc(nclass,sizeof(int))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate n_z.\n");
        return;
    }
    if((m_m = calloc(ndocs,sizeof(int))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate m_m.\n");
        return;
    }
    if((m_z = calloc(nclass,sizeof(int))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate m_z.\n");
        return;
    }
    if((left = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate left.\n");
        return;
    }
    if((center = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate center.\n");
        return;
    }
    if((right = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate right.\n");
        return;
    }
    if((p_z = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate p_z.\n");
        return;
    }
    if((cum_sum_p_z = calloc((nclass+1),sizeof(double))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate cum_sum_p_z.\n");
        return;
    }
    if((topics = calloc(ndocs,sizeof(int **))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate topics.\n");
        return;
    }
    if((tagtopics = calloc(ndocs,sizeof(int **))) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate tagtopics.\n");
        return;
    }
    if((temp_phi = dmatrix(nlex, nclass)) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate temp_phi.\n");
        exit(1);
    }
    if((temp_theta = dmatrix(ndocs, nclass)) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate temp_theta.\n");
        exit(1);
    }
    if((temp_psi = dmatrix(ntlex, nclass)) == NULL){
        fprintf(stderr,"corrlda_learn:: cannot allocate temp_psi.\n");
        exit(1);
    }
    
    printf("Number of documents          = %d\n",ndocs);
    printf("Number of unique words       = %d\n",nlex);
    printf("Number of unique tags        = %d\n",ntlex);
    printf("Number of latent classes     = %d\n",nclass);
    printf("Number of iteration          = %d\n",maxiter);
    
    // choose an arbitrary topic as first topic for word
    gsl_rng_set(r, random_seed);
    for(dp = data, m = 0;(dp->len) != -1;dp++, m++){
        if((topics[m] = calloc((dp->len), sizeof(int *))) == NULL){
            fprintf(stderr,"corrlda_learn:: cannot allocate topics[m].\n");
            return;
        }
        for(w = 0;w < (dp->len);w++){
            if((topics[m][w] = calloc((dp->cnt[w]), sizeof(int))) == NULL){
                fprintf(stderr,"corrlda_learn:: cannot allocate topics[m][w].\n");
                return;
            }
            word_index = dp->id[w];
            word_num = dp->cnt[w];
            for(i = 0;i < word_num;i++){
                z = (int)gsl_rng_uniform_int(r, nclass);
                n_mz[m][z] += 1;
                n_m[m] += 1;
                n_zw[z][word_index] += 1;
                n_z[z] += 1;
                topics[m][w][i] = z;
            }
        }
    }
    // setting first topics for tags using initial word-topic distributions.
    Mall = 0;
    for(tp = tag, m = 0;(tp->len) != -1;tp++, m++){
        if((tagtopics[m] = calloc((tp->len), sizeof(int *))) == NULL){
            fprintf(stderr, "corrlda_learn:: cannot allocate tagtopics[m].\n");
            return;
        }
        for(t = 0;t < (tp->len);t++){
            if((tagtopics[m][t] = calloc((tp->cnt[t]), sizeof(int))) == NULL){
                fprintf(stderr, "corrlda_learn:: cannot allocate tagtopics[m][t].\n");
                return;
            }
            tag_index = tp->id[t];
            tag_num = tp->cnt[t];
            for(i = 0;i < tag_num;i++){
                Mall += 1;
                sum_p_z = 0.0;
                for(k = 0;k < nclass;k++){
                    p_z[k] = (double)n_mz[m][k];
                    sum_p_z += p_z[k];
                }
                for(k = 0;k < nclass;k++){
                    p_z[k] /= sum_p_z;
                }
                // random sampling from p_z
                z = sampling_multinomial(r, p_z, cum_sum_p_z, nclass);
                m_mz[m][z] += 1;
                m_m[m] += 1;
                m_zt[z][tag_index] += 1;
                m_z[z] += 1;
                tagtopics[m][t][i] = z;
            }
        }
    }
    
    // learning main
    for(it = 0;it < maxiter;it++){
        printf("iteration %2d/%3d..\n", it + 1, maxiter);
        fflush(stdout);
        sum_alpha = 0.0;
        for(k = 0;k < nclass;k++)
            sum_alpha += alpha[k];
        for (dp = data, tp = tag, m = 0; (dp->len) != -1; dp++, tp++, m++){
            // for words
            for(w = 0;w < (dp->len);w++){
                word_index = dp->id[w];
                word_num = dp->cnt[w];
                for(i = 0;i < word_num;i++){
                    z = topics[m][w][i];
                    n_mz[m][z] -= 1;
                    n_m[m] -= 1;
                    n_zw[z][word_index] -= 1;
                    n_z[z] -= 1;
                    
                    // compute conditional distribution p_z
                    // p_z left
                    for(k = 0;k < nclass;k++){
                        left[k] = (double)n_mz[m][k] + alpha[k];
                        left[k] /= ((double)n_m[m] + sum_alpha);
                    }
                    // p_z center
                    for(k = 0;k < nclass;k++){
                        center[k] = (double)n_zw[k][word_index] + beta;
                        center[k] /= ((double)n_z[k] + (double)nlex * beta);
                    }
                    // p_z right
                    infinite_prob = false;
                    for(k = 0;k < nclass;k++){
                        right[k] = 1.0;
                        for(j = 0;j < m_mz[m][k];j++){
                            if(n_mz[m][k] == 0){
                                infinite_prob = true;
                                z = k;
                                break;
                            }
                            right[k] *= (((double)n_mz[m][k] + 1.0) / (double)n_mz[m][k]) * (((double)n_m[m] - 1.0) / (double)n_m[m]);
                        }
                    }
                    if(infinite_prob == false){
                        // conditional distribution p_z
                        sum_p_z = 0.0;
                        for(k = 0;k < nclass;k++){
                            p_z[k] = left[k] * center[k] * right[k];
                            sum_p_z += p_z[k];
                        }
                        for(k = 0;k < nclass;k++){
                            p_z[k] /= sum_p_z; // normalize to obtain probabilities
                        }
                        // random sampling from p_z
                        z = sampling_multinomial(r, p_z, cum_sum_p_z, nclass);
                    }
                    // update buffers
                    n_mz[m][z] += 1;
                    n_m[m] += 1;
                    n_zw[z][word_index] += 1;
                    n_z[z] += 1;
                    topics[m][w][i] = z;
                }
            }
            // for tags
            for(t = 0;t < (tp->len);t++){
                tag_index = tp->id[t];
                tag_num = tp->cnt[t];
                for(i = 0;i < tag_num;i++){
                    z = tagtopics[m][t][i];
                    m_mz[m][z] -= 1;
                    m_m[m] -= 1;
                    m_zt[z][tag_index] -= 1;
                    m_z[z] -= 1;
                    Mall -= 1;
                    // p_z left
                    for(k = 0;k < nclass;k++)
                        left[k] = ((double)m_zt[k][tag_index] + gamma) / ((double)m_z[k] + gamma * ntlex);
                    //p_z right
                    for(k = 0;k < nclass;k++)
                        right[k] = (double)n_mz[m][k] / (double)n_m[m];
                    // conditional distribution p_z
                    sum_p_z = 0.0;
                    for(k = 0;k < nclass;k++){
                        p_z[k] = left[k] * right[k];
                        sum_p_z += p_z[k];
                    }
                    for(k = 0;k < nclass;k++)
                        p_z[k] /= sum_p_z;
                    // random sampling from p_z
                    z = sampling_multinomial(r, p_z, cum_sum_p_z, nclass);
                    // update buffers
                    m_mz[m][z] += 1;
                    m_m[m] += 1;
                    m_zt[z][tag_index] += 1;
                    m_z[z] += 1;
                    tagtopics[m][t][i] = z;
                    Mall += 1;
                }
            }
        }
        
        // update hyperparameters.
        update_alpha(alpha, n_m, n_mz, ndocs, nclass);
        beta = update_beta(beta, n_z, n_zw, nclass, nlex);
        gamma = update_gamma(gamma, m_z, m_zt, nclass, ntlex);
        
        // compute likelihood.
        lik = loglikelihood(n_mz, n_zw, n_m, m_mz, m_zt, m_z, nclass, nlex, ntlex, ndocs, Mall, alpha, beta, gamma);
        printf("\tlikelihood ... %.8f\n",lik);
        printf("\talpha = \n\t");
        for(k = 0;k < nclass;k++)
            printf("%.8f ",alpha[k]);
        printf("\n\tbeta ... %.2f, gamma ... %.2f\n",beta,gamma);
        fprintf(likp,"%.8f\n",lik);
        for(k = 0;k < nclass;k++)
            fprintf(hyperp,"%.8f,",alpha[k]);
        fprintf(hyperp,"%.8f,%.8f\n",beta,gamma);
    }
    
    // compute matrix phi ([nlex, nclass] matrix)
    for(w = 0;w < nlex;w++)
        for(k = 0;k < nclass;k++)
            temp_phi[w][k] = (double)n_zw[k][w] + beta;
    normalize_matrix_col(phi, temp_phi, nlex, nclass);
    
    // compute matrix theta ([ndocs, nclass])
    for(m = 0;m < ndocs;m++)
        for(k = 0;k < nclass;k++)
            temp_theta[m][k] = (double)n_mz[m][k] + alpha[k];
    normalize_matrix_row(theta, temp_theta, ndocs, nclass);
    
    // compute matrix psi ([ntlex, nclass])
    for(t = 0;t < ntlex;t++)
        for(k = 0;k < nclass;k++)
            temp_psi[t][k] = (double)m_zt[k][t] + gamma;
    normalize_matrix_col(psi, temp_psi, ntlex, nclass);
    
    free(n_m);
    free(n_z);
    free(m_m);
    free(m_z);
    free(left);
    free(center);
    free(right);
    free(p_z);
    free(cum_sum_p_z);
    
    for(dp = data, m = 0;(dp->len) != -1;dp++, m++){
        for(w = 0;w < (dp->len);w++){
            free(topics[m][w]);
        }
        free(topics[m]);
    }
    free(topics);
    for(tp = tag, m = 0;(tp->len) != -1;tp++, m++){
        for(t = 0;t < (tp->len);t++){
            free(tagtopics[m][t]);
        }
        free(tagtopics[m]);
    }
    free(tagtopics);
    free_dmatrix(temp_phi, nlex);
    free_dmatrix(temp_theta, ndocs);
    free_dmatrix(temp_psi, ntlex);
    
    return;
}

int sampling_multinomial(gsl_rng *r, double *p, double *cum_sum_p, int len_p){
    int k, z;
    double sampling;
    
    cum_sum_p[0] = 0.0;
    for(k = 0;k < len_p;k++){
        cum_sum_p[k+1] = cum_sum_p[k] + p[k];
    }
    sampling = gsl_rng_uniform(r);
    for(k = 0;k < len_p;k++){
        if((sampling >= cum_sum_p[k]) && (sampling < cum_sum_p[k+1])){
            z = k;
            break;
        }
    }
    return z;
}
