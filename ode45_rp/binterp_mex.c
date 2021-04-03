#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *t_int, *x_int, *x;
    double t;
    int M, N;
    
    t_int = mxGetPr(prhs[0]);
    x_int = mxGetPr(prhs[1]);
    t = mxGetScalar(prhs[2]);
    M = mxGetM(prhs[1]);
    N = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
    x = mxGetPr(plhs[0]);
    
    binterp(x,t_int,x_int,t,M,N);
}
