#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x, *f, *u;
    int nu, M;
    const mxArray **args;
    
    x = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    nu = mxGetScalar(prhs[1]);
    args = prhs+2;
    
    plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
    f = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(nu,1,mxREAL);
    u = mxGetPr(plhs[1]);
    
    ffun(f,u,x,args);
}