#include "ocp_mex.h"

void eventfun(double t, double *y, double *v, double *direction, int nv, int neq, const mxArray **args)
{
    double *f, *u;
    int i, nu;
    
    nu = mxGetScalar(args[0]);
    args = args+1;
    
    f = mxCalloc(neq,sizeof(double));
    u = mxCalloc(nu,sizeof(double));
    
    ffun(f,u,y,args);
    
    ghfun(v,u,y,args);
    
    for(i = 0; i < nv; i++) {
        if(v[i]>-1e-9&&v[i]<=0)
            v[i] = 0;
        direction[i] = 1.0;
    }
    
    mxFree(f);
    mxFree(u);
}
