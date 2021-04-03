#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int next, ntspan, haveeventfun, stop, tdir, nofailed;
    double *tspan, tnew, *ynew, t, *y, h, *f;
    double temp, err, rtol, absh;
    int nout_new, oldnout, nout;
    double *tout_new, *yout_new;
    int neq;
    
    next = mxGetScalar(prhs[0]);
    ntspan = mxGetScalar(prhs[1]);
    haveeventfun = mxGetScalar(prhs[2]);
    stop = mxGetScalar(prhs[3]);
    tdir = mxGetScalar(prhs[4]);
    tspan = mxGetPr(prhs[5]);
    tnew = mxGetScalar(prhs[6]);
    ynew = mxGetPr(prhs[7]);
    t = mxGetScalar(prhs[8]);
    y = mxGetPr(prhs[9]);
    h = mxGetScalar(prhs[10]);
    f = mxGetPr(prhs[11]);
    nofailed = mxGetScalar(prhs[12]);
    err = mxGetScalar(prhs[13]);
    rtol = mxGetScalar(prhs[14]);
    absh = mxGetScalar(prhs[15]);
    nout = mxGetScalar(prhs[16]);
    neq = mxGetScalar(prhs[17]);
    
    next--;
    tout_new = mxCalloc((ntspan-next),sizeof(double));
    yout_new = mxCalloc(neq*(ntspan-next),sizeof(double));
    
    rp45(&nout_new,&next,tout_new,yout_new,ntspan,haveeventfun,stop,tdir,tspan,tnew,ynew,t,y,h,f,neq);
    
    oldnout = nout;
    nout = nout+nout_new;
    next++;
    
    plhs[0] = mxCreateDoubleScalar(oldnout);
    plhs[1] = mxCreateDoubleScalar(nout);
    plhs[2] = mxCreateDoubleMatrix(1,nout_new,mxREAL);
    memcpy(mxGetPr(plhs[2]),tout_new,nout_new*sizeof(double));
    plhs[3] = mxCreateDoubleMatrix(neq,nout_new,mxREAL);
    memcpy(mxGetPr(plhs[3]),yout_new,neq*nout_new*sizeof(double));
    plhs[4] = mxCreateDoubleScalar(next);
    
    mxFree(tout_new);
    mxFree(yout_new);
}
