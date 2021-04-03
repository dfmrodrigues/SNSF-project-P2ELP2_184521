#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *yinterp, tinterp, t, *y, h, *f;
    int neq;
    
    tinterp = mxGetScalar(prhs[0]);
    t = mxGetScalar(prhs[1]);
    y = mxGetPr(prhs[2]);
    h = mxGetScalar(prhs[3]);
    f = mxGetPr(prhs[4]);
    neq = mxGetM(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(neq,1,mxREAL);
    yinterp = mxGetPr(plhs[0]);
    
    ntrp45(yinterp,tinterp,t,y,h,f,neq);
}
