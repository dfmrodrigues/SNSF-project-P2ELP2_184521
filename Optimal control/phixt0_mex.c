#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double t, *x, *xi, *phi_xt0;
    int M, Mi;
    const mxArray **args;
    
    t = mxGetScalar(prhs[0]);
    x = mxGetPr(prhs[1]);
    M = mxGetM(prhs[1]);
    xi = mxGetPr(prhs[2]);
    Mi = mxGetM(prhs[2]);
    args = prhs+3;
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    phi_xt0 = mxGetPr(plhs[0]);
    
    phixt0fun(phi_xt0,t,x,xi,args);
}