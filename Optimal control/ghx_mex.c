#include "ocp_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x, *f_x, *u_x, *gh_x;
    int nu, ngh, M;
    const mxArray **args;
    
    x = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    nu = mxGetScalar(prhs[1]);
    ngh = mxGetScalar(prhs[2]);
    args = prhs+3;
    
    plhs[0] = mxCreateDoubleMatrix(ngh,M,mxREAL);
    gh_x = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(M,M,mxREAL);
    f_x = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(nu,M,mxREAL);
    u_x = mxGetPr(plhs[2]);
    
    fxfun(f_x,u_x,x,args);
    
    ghxfun(gh_x,u_x,x,args);
}