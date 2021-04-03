#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double t;
    double *y, *v, *direction;
    int nv, neq;
    const mxArray **args;
    
    t = mxGetScalar(prhs[0]);
    y = mxGetPr(prhs[1]);
    neq = mxGetScalar(prhs[2]);
    nv = mxGetScalar(prhs[3]);
    args = prhs+4;
    
    plhs[0] = mxCreateDoubleMatrix(nv,1,mxREAL);
    v = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(nv,1,mxREAL);
    direction = mxGetPr(plhs[1]);
    
    eventfun(t,y,v,direction,nv,neq,args);
}