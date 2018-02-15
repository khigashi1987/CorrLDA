/*
    writer.c
    an implementation of matrix writer
*/
#include <stdio.h>
#include "writer.h"

void corrlda_write(FILE *pp, FILE *tp, FILE *sp, FILE *n_mzp, FILE *n_wzp, FILE *m_mzp, FILE *m_tzp, double **phi, double **theta, double **psi, int **n_mz, int **n_zw, int **m_mz, int **m_zt, int nclass, int nlex, int ntlex, int ndoc){
    printf("writing model...\n"); fflush(stdout);
    write_matrix(pp, phi, nlex, nclass);
    write_matrix(tp, theta, ndoc, nclass);
    write_matrix(sp, psi, ntlex, nclass);
    write_imatrix(n_mzp, n_mz, ndoc, nclass);
    write_imatrix_transpose(n_wzp, n_zw, nclass, nlex);
    write_imatrix(m_mzp, m_mz, ndoc, nclass);
    write_imatrix_transpose(m_tzp, m_zt, nclass, ntlex);
    printf("done.\n"); fflush(stdout);
}

void write_matrix(FILE *fp, double **matrix, int rows, int cols){
    int i, j;
    for(i = 0;i < rows;i++)
        for(j = 0;j < cols;j++)
            fprintf(fp, "%.7e%s", matrix[i][j], (j == cols - 1)? "\n":" ");
}

void write_imatrix(FILE *fp, int **matrix, int rows, int cols){
    int i, j;
    for(i = 0;i < rows;i++)
        for(j = 0;j < cols;j++)
            fprintf(fp, "%d%s", matrix[i][j], (j == cols - 1)? "\n":" ");
}

void write_imatrix_transpose(FILE *fp, int **matrix, int rows, int cols){
    int i, j;
    for(j = 0;j < cols;j++)
        for(i = 0;i < rows;i++)
            fprintf(fp, "%d%s", matrix[i][j], (i == rows -1)? "\n":" ");
}
