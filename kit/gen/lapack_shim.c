/* Sandbox only: forwards standard LAPACK/BLAS names to the scipy_-prefixed OpenBLAS wheel. */
double scipy_ddot_(int*,double*,int*,double*,int*);
void scipy_dgemv_(char*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
void scipy_dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
void scipy_dsyrk_(char*,char*,int*,int*,double*,double*,int*,double*,double*,int*);
void scipy_dsytrf_(char*,int*,double*,int*,int*,double*,int*,int*);
void scipy_dsytrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*);
double ddot_(int*n,double*x,int*ix,double*y,int*iy){return scipy_ddot_(n,x,ix,y,iy);}
void dgemv_(char*t,int*m,int*n,double*al,double*a,int*lda,double*x,int*ix,double*be,double*y,int*iy){scipy_dgemv_(t,m,n,al,a,lda,x,ix,be,y,iy);}
void dgesv_(int*n,int*nr,double*a,int*lda,int*ip,double*b,int*ldb,int*info){scipy_dgesv_(n,nr,a,lda,ip,b,ldb,info);}
void dsyrk_(char*u,char*t,int*n,int*k,double*al,double*a,int*lda,double*be,double*c,int*ldc){scipy_dsyrk_(u,t,n,k,al,a,lda,be,c,ldc);}
void dsytrf_(char*u,int*n,double*a,int*lda,int*ip,double*w,int*lw,int*info){scipy_dsytrf_(u,n,a,lda,ip,w,lw,info);}
void dsytrs_(char*u,int*n,int*nr,double*a,int*lda,int*ip,double*b,int*ldb,int*info){scipy_dsytrs_(u,n,nr,a,lda,ip,b,ldb,info);}
