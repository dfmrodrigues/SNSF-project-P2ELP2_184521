#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x, *f_x, *u_x;
    int nu, M;
    const mxArray **args;
    
    x = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    nu = mxGetScalar(prhs[1]);
    args = prhs+2;
    
    plhs[0] = mxCreateDoubleMatrix(M,M,mxREAL);
    f_x = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(nu,M,mxREAL);
    u_x = mxGetPr(plhs[1]);
    
    fxfun(f_x,u_x,x,args);
}