#include "ocp_mex.h"

void odefun(double *f, double t, double *y, int neq, const mxArray *args[])
{
    double *t_int, *x_int, *x;
    double deltat;
    int N, i, k;
    double *f_x, *f_d;
    
    t_int = mxGetPr(args[0]);
    x_int = mxGetPr(args[1]);
    deltat = mxGetScalar(args[2]);
    neq--;
    N = mxGetN(args[1]);
    args = args+4;
    x = mxCalloc(neq,sizeof(double));
    f_x = mxCalloc(neq*neq,sizeof(double));
    f_d = mxCalloc(neq,sizeof(double));
    binterp(x,t_int,x_int,t,neq,N);
    ffun(f_d,NULL,x,args);
    fxfun(f_x,NULL,x,args);
    for(i = 0; i < neq; i++) {
        f[i] = 0;
        for(k = 0; k < neq; k++)
            f[i] = f[i]-f_x[k+i*neq]*y[k];
    }
    f[neq] = 0;
    for(k = 0; k < neq; k++)
        f[neq] = f[neq]-f_d[k]*y[k]/deltat;
    mxFree(x);
    mxFree(f_x);
    mxFree(f_d);
}