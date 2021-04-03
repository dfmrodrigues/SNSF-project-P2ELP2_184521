#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x, *f_u, *gh_u;
    int nu, ngh, M;
    const mxArray **args;
    
    x = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    nu = mxGetScalar(prhs[1]);
    ngh = mxGetScalar(prhs[2]);
    args = prhs+3;
    
    plhs[0] = mxCreateDoubleMatrix(ngh,nu,mxREAL);
    gh_u = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(M,nu,mxREAL);
    f_u = mxGetPr(plhs[1]);
    
    fufun(f_u,x,args);
    
    ghufun(gh_u,x,args);
}