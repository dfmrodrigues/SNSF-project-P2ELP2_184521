#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *yinterp, tinterp, t, *y, tnew, *ynew, *yp, *ypnew;
    int neq;
    
    tinterp = mxGetScalar(prhs[0]);
    t = mxGetScalar(prhs[1]);
    y = mxGetPr(prhs[2]);
    tnew = mxGetScalar(prhs[3]);
    ynew = mxGetPr(prhs[4]);
    yp = mxGetPr(prhs[5]);
    ypnew = mxGetPr(prhs[6]);
    neq = mxGetM(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(neq,1,mxREAL);
    yinterp = mxGetPr(plhs[0]);
    
    ntrp3h(yinterp,tinterp,t,y,tnew,ynew,yp,ypnew,neq);
}
