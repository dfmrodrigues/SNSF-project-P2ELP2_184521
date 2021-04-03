#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *tint, *t, *y, *yp, *history;
    double *yout;
    int ntint, neq, nt;
    
    tint = mxGetPr(prhs[0]);
    t = mxGetPr(prhs[1]);
    y = mxGetPr(prhs[2]);
    yp = mxGetPr(prhs[3]);
    history = mxGetPr(prhs[4]);
    ntint = mxGetScalar(prhs[5]);
    neq = mxGetScalar(prhs[6]);
    nt = mxGetN(prhs[1]);
    
    yout = mxCalloc(neq*ntint,sizeof(double));
    
    rp3h(yout,tint,t,y,yp,history,ntint,neq,nt);
    
    plhs[0] = mxCreateDoubleMatrix(neq,ntint,mxREAL);
    memcpy(mxGetPr(plhs[0]),yout,neq*ntint*sizeof(double));
    
    mxFree(yout);
}
