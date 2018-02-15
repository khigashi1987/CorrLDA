/*
    writer.h
    a header file of matrix writer
*/
#ifndef WRITER_H
#define WRITER_H
#include <stdio.h>

extern void corrlda_write(FILE *pp, FILE *tp, FILE *sp,
        FILE *n_mzp, FILE *n_wzp, FILE *m_mzp, FILE *m_tzp,
        double **phi, double **theta, double **psi,
        int **n_mz, int **n_zw, int **m_mz, int **m_zt,
        int nclass, int nlex, int ntlex, int ndoc);
void write_matrix(FILE *fp, double **matrix, int rows, int cols);
void write_imatrix(FILE *fp, int **matrix, int rows, int cols);
void write_imatrix_transpose(FILE *fp, int **matrix, int rows, int cols);
#endif
