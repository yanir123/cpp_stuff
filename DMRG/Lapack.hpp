#pragma once

extern "C" {
void dgeqrf_(int *m, int *n, double *A, int *LDA, double *tau, double *WORK, int *LWORK, int *INFO);
void dorgqr_(int *m, int *n, int *k, double *A, int *LDA, double *tau, double *WORK, int *LWORK, int *INFO);
void dgelqf_(int *m, int *n, double *A, int *LDA, double *tau, double *WORK, int *LWORK, int *INFO);
void dorglq_(int *m, int *n, int *k, double *A, int *LDA, double *tau, double *WORK, int *LWORK, int *INFO);
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void dscal_(int *n, double *alpha, double *x, int *incx);
void dgemm_(char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc);
void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *X, int *incx, double *beta, double *Y, int *incy);
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda, double *W, double *work, int *lwork, int *info);
void dgesdd_(char *JOBZ, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *IWORK, int *INFO);
void dlasrt_(char *id, int *n, double *vec, int *info);
double dlansy_(char *norm, char *uplo, int *dimR, double *mx, int *lda, double *work);
double dlange_(char *norm, int *m, int *n, double *mx, int *lda, double *work);
}