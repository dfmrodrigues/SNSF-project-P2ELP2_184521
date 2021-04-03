#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double t, *x, *xi, *phi_xt1;
    int M, Mi;
    const mxArray **args;
    
    t = mxGetScalar(prhs[0]);
    x = mxGetPr(prhs[1]);
    M = mxGetM(prhs[1]);
    xi = mxGetPr(prhs[2]);
    Mi = mxGetM(prhs[2]);
    args = prhs+3;
    
    plhs[0] = mxCreateDoubleMatrix(1,M+1+Mi,mxREAL);
    phi_xt1 = mxGetPr(plhs[0]);
    
    phixt1fun(phi_xt1,t,x,xi,args);
}