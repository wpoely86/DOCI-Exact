#ifndef LAPACK_H
#define LAPACK_H

//lapack definities
extern "C" {

   void dgeqrf_(int *m,int *n,double *A,int *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dorgqr_(int *m,int *n, int *k,double *A,int  *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dgelqf_(int *m,int *n,double *A,int *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dorglq_(int *m,int *n, int *k,double *A,int  *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dcopy_(int *n,const double *x,int *incx,double *y,int *incy);
   void daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
   void dscal_(int *n,const double *alpha,double *x,int *incx);
   void dgemm_(char *transA,char *transB,const int *m,const int *n,const int *k,double *alpha,double *A,const int *lda,double *B,const int *ldb,double *beta,double *C,const int *ldc);
   void dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dgemv_(char *trans,int *m,int *n,double *alpha,double *A,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
   double ddot_(const int *n,double *x,int *incx,double *y,int *incy);
   void dsyev_(char *jobz,char *uplo,int *n,double *A,int *lda,double *W,double *work,int *lwork,int *info);
   void dpotrf_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dpotri_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dsyevr_( char* jobz, char* range, char* uplo, int* n, double* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, double* z, int* ldz, int* isuppz, double* work, int* lwork, int* iwork, int* liwork, int* info );
   void dsyevd_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* iwork, int* liwork, int* info );

   void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );

    void dsyevd_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* iwork, int* liwork, int* info);
    void dsymv_(char *uplo, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
    void dstev_( const char* jobz, const int* n, double* d, double* e, double* z, const int* ldz, double* work, int* info );
    void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );
    void dgemm_(char *transA,char *transB,const int *m,const int *n,const int *k,double *alpha,double *A,const int *lda,double *B,const int *ldb,double *beta,double *C,const int *ldc);

//    void zheevd_( const char* jobz, const char* uplo, const int* n, std::complex<double>* a, const int* lda, double* w, std::complex<double>* work, const int* lwork, double* rwork,
//            const int* lrwork, int* iwork, const int* liwork, int* info );

    void dsaupd_(int *ido, char *bmat, int *n, char *which,
            int *nev, double *tol, double *resid, int *ncv,
            double *v, int *ldv, int *iparam, int *ipntr,
            double *workd, double *workl, int *lworkl, int *info);

    void dseupd_(int *rvec, char *All, int *select, double *d,
            double *z, int *ldz, double *sigma,
            char *bmat, int *n, char *which, int *nev,
            double *tol, double *resid, int *ncv, double *v,
            int *ldv, int *iparam, int *ipntr, double *workd,
            double *workl, int *lworkl, int *info);
}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
