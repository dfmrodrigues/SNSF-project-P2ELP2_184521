#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x, *f, *u, *gh;
    int nu, ngh, M;
    const mxArray **args;
    
    x = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    nu = mxGetScalar(prhs[1]);
    ngh = mxGetScalar(prhs[2]);
    args = prhs+3;
    
    plhs[0] = mxCreateDoubleMatrix(ngh,1,mxREAL);
    gh = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(M,1,mxREAL);
    f = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(nu,1,mxREAL);
    u = mxGetPr(plhs[2]);
    
    ffun(f,u,x,args);
    
    ghfun(gh,u,x,args);
}