#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double t0, t, h, tnew, err, tfinal, htspan, threshold;
    int next, ntspan, haveeventfun;
    double rtol, absh, rh, hmax, hmin, tdir;
    double *tspan, temp;
    int done, nofailed, nfevals, nfailed, nsteps;
    int nout_new, oldnout, nout;
    double *tout_new, *yout_new, *tout, *yout;
    double *y0, *f0, *ynew, *f, *y;
    double *v0, *v, *vnew, *direction, *indzc;
    int i, j, stop, neq;
    int nv;
    const mxArray **args;
    
    t0 = mxGetScalar(prhs[0]);
    y0 = mxGetPr(prhs[1]);
    neq = mxGetScalar(prhs[2]);
    tfinal = mxGetScalar(prhs[3]);
    htspan = mxGetScalar(prhs[4]);
    threshold = mxGetScalar(prhs[5]);
    nfevals = mxGetScalar(prhs[6]);
    nfailed = mxGetScalar(prhs[7]);
    rtol = mxGetScalar(prhs[8]);
    hmax = mxGetScalar(prhs[9]);
    tdir = mxGetScalar(prhs[10]);
    nsteps = mxGetScalar(prhs[11]);
    v0 = mxGetPr(prhs[12]);
    next = mxGetScalar(prhs[13]);
    ntspan = mxGetScalar(prhs[14]);
    haveeventfun = mxGetScalar(prhs[15]);
    tspan = mxGetPr(prhs[16]);
    nv = mxGetScalar(prhs[17]);
    args = prhs+18;
    
    plhs[1] = mxCreateDoubleMatrix(neq,1,mxREAL);
    ynew = mxGetPr(plhs[1]);
    plhs[3] = mxCreateDoubleMatrix(neq,7,mxREAL);
    f = mxGetPr(plhs[3]);
    plhs[11] = mxCreateDoubleMatrix(nv,1,mxREAL);
    vnew = mxGetPr(plhs[11]);
    plhs[12] = mxCreateDoubleMatrix(nv,1,mxREAL);
    direction = mxGetPr(plhs[12]);
    plhs[13] = mxCreateDoubleMatrix(nv,1,mxREAL);
    indzc = mxGetPr(plhs[13]);
    plhs[14] = mxCreateDoubleMatrix(nv,1,mxREAL);
    v = mxGetPr(plhs[14]);
    plhs[16] = mxCreateDoubleMatrix(1,ntspan,mxREAL);
    tout = mxGetPr(plhs[16]);
    plhs[17] = mxCreateDoubleMatrix(neq,ntspan,mxREAL);
    yout = mxGetPr(plhs[17]);
    plhs[20] = mxCreateDoubleMatrix(neq,1,mxREAL);
    y = mxGetPr(plhs[20]);
    
    t = t0;
    memcpy(y,y0,neq*sizeof(double));
    odefun(f,t,y,neq,args);
    memcpy(vnew,v0,nv*sizeof(double));
    
    hmin = 16*DBL_EPSILON*pow(2,floor(log(t0)*M_LOG2E));
    absh = min(hmax,htspan);
    rh = 0;
    for(i = 0; i < neq; i++) {
        double rh_i;
        rh_i = f[i];
        rh_i = abs(rh_i)/max(abs(y[i]),threshold);
        if(rh_i > rh)
            rh = rh_i;
    }
    rh = rh/(0.8*pow(rtol,0.2));
    if(absh*rh > 1)
        absh = 1/rh;
    absh = max(absh,hmin);
    
    nout = 1;
    *tout = t;
    memcpy(yout,y,neq*sizeof(double));
    
    done = 0;
    while(!done) {
        
        hmin = 16*DBL_EPSILON*pow(2,floor(log(t)*M_LOG2E));
        absh = min(hmax, max(hmin, absh));
        h = tdir * absh;
        if(1.1*absh >= abs(tfinal - t)) {
            h = tfinal - t;
            absh = abs(h);
            done = 1;
        }
        else
            done = 0;
        nofailed = 1;
        
        while(1) {
            tnew = t+h*1.0/5.0;
            for(i = 0; i < neq; i++) {
                ynew[i] = y[i];
                ynew[i] = ynew[i]+1.0/5.0*h*f[i+0*neq];
            }
            odefun(f+neq,tnew,ynew,neq,args);
            tnew = t+h*3.0/10.0;
            for(i = 0; i < neq; i++) {
                ynew[i] = y[i];
                ynew[i] = ynew[i]+3.0/40.0*h*f[i+0*neq];
                ynew[i] = ynew[i]+9.0/40.0*h*f[i+1*neq];
            }
            odefun(f+2*neq,tnew,ynew,neq,args);
            tnew = t+h*4.0/5.0;
            for(i = 0; i < neq; i++) {
                ynew[i] = y[i];
                ynew[i] = ynew[i]+44.0/45.0*h*f[i+0*neq];
                ynew[i] = ynew[i]-56.0/15.0*h*f[i+1*neq];
                ynew[i] = ynew[i]+32.0/9.0*h*f[i+2*neq];
            }
            odefun(f+3*neq,tnew,ynew,neq,args);
            tnew = t+h*8.0/9.0;
            for(i = 0; i < neq; i++) {
                ynew[i] = y[i];
                ynew[i] = ynew[i]+19372.0/6561.0*h*f[i+0*neq];
                ynew[i] = ynew[i]-25360.0/2187.0*h*f[i+1*neq];
                ynew[i] = ynew[i]+64448.0/6561.0*h*f[i+2*neq];
                ynew[i] = ynew[i]-212.0/729.0*h*f[i+3*neq];
            }
            odefun(f+4*neq,tnew,ynew,neq,args);
            tnew = t+h;
            for(i = 0; i < neq; i++) {
                ynew[i] = y[i];
                ynew[i] = ynew[i]+9017.0/3168.0*h*f[i+0*neq];
                ynew[i] = ynew[i]-355.0/33.0*h*f[i+1*neq];
                ynew[i] = ynew[i]+46732.0/5247.0*h*f[i+2*neq];
                ynew[i] = ynew[i]+49.0/176.0*h*f[i+3*neq];
                ynew[i] = ynew[i]-5103.0/18656.0*h*f[i+4*neq];
            }
            odefun(f+5*neq,tnew,ynew,neq,args);
            
            tnew = t+h;
            if(done)
                tnew = tfinal;
            h = tnew - t;
            for(i = 0; i < neq; i++) {
                ynew[i] = y[i];
                ynew[i] = ynew[i]+35.0/384.0*h*f[i+0*neq];
                ynew[i] = ynew[i]+500.0/1113.0*h*f[i+2*neq];
                ynew[i] = ynew[i]+125.0/192.0*h*f[i+3*neq];
                ynew[i] = ynew[i]-2187.0/6784.0*h*f[i+4*neq];
                ynew[i] = ynew[i]+11.0/84.0*h*f[i+5*neq];
            }
            odefun(f+6*neq,tnew,ynew,neq,args);
            
            err = 0;
            for(i = 0; i < neq; i++) {
                double err_i;
                err_i = 0;
                err_i = err_i+71.0/57600.0*f[i+0*neq];
                err_i = err_i-71.0/16695.0*f[i+2*neq];
                err_i = err_i+71.0/1920.0*f[i+3*neq];
                err_i = err_i-17253.0/339200.0*f[i+4*neq];
                err_i = err_i+22.0/525.0*f[i+5*neq];
                err_i = err_i-1.0/40.0*f[i+6*neq];
                err_i = abs(err_i)/max(max(abs(y[i]),abs(ynew[i])),threshold);
                if(err_i > err)
                    err = err_i;
            }
            err = absh*err;
            
            nfevals = nfevals+6;
            
            if(err > rtol) {
                nfailed = nfailed + 1;
                if(nofailed) {
                    nofailed = 0;
                    absh = max(hmin,absh*max(0.1,0.8*pow(rtol/err,0.2)));
                }
                else
                    absh = max(hmin,0.5*absh);
                h = tdir * absh;
                done = 0;
            }
            else
                break;
        }
        nsteps++;
        stop = 0;
        
        if(haveeventfun) {
            memcpy(v,vnew,nv*sizeof(double));
            
            eventfun(tnew,ynew,vnew,direction,nv,neq,args);
            
            stop = 0;
            for(i = 0; i < nv; i++) {
                indzc[i] = (sign(v[i])!=sign(vnew[i]))&&(direction[i]*(vnew[i]-v[i])>=0);
                if(indzc[i])
                    stop = 1;
            }
            if(stop)
                break;
        }
        
        next--;
        tout_new = mxCalloc((ntspan-next),sizeof(double));
        yout_new = mxCalloc(neq*(ntspan-next),sizeof(double));
        
        rp45(&nout_new,&next,tout_new,yout_new,ntspan,haveeventfun,stop,tdir,tspan,tnew,ynew,t,y,h,f,neq);
        
        oldnout = nout;
        nout = nout+nout_new;
        next++;
        
        if(nofailed) {
            temp = 1.25*pow(err/rtol,0.2);
            if(temp > 0.2)
                absh = absh/temp;
            else
                absh = 5.0*absh;
        }
        t = tnew;
        memcpy(y,ynew,neq*sizeof(double));
        memcpy(f,f+6*neq,neq*sizeof(double));
        
        memcpy(tout+oldnout,tout_new,(nout-oldnout)*sizeof(double));
        memcpy(yout+oldnout*neq,yout_new,neq*(nout-oldnout)*sizeof(double));
        
        mxFree(tout_new);
        mxFree(yout_new);
    }
    
    plhs[0] = mxCreateDoubleScalar(tnew);
    plhs[2] = mxCreateDoubleScalar(h);
    plhs[4] = mxCreateDoubleScalar(err);
    plhs[5] = mxCreateDoubleScalar(nfevals);
    plhs[6] = mxCreateDoubleScalar(nfailed);
    plhs[7] = mxCreateDoubleScalar(nofailed);
    plhs[8] = mxCreateDoubleScalar(absh);
    plhs[9] = mxCreateDoubleScalar(nsteps);
    plhs[10] = mxCreateDoubleScalar(stop);
    plhs[15] = mxCreateDoubleScalar(nout);
    plhs[18] = mxCreateDoubleScalar(next);
    plhs[19] = mxCreateDoubleScalar(t);
}