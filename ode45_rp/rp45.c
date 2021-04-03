#include "ode45_rp.h"

void rp45(int* nout_new_ptr, int* next_ptr, double *tout_new, double *yout_new, int ntspan, int haveeventfun, int stop, int tdir, double *tspan, double tnew, double *ynew, double t, double *y, double h, double *f, int neq)
{
    int nout_new, next;
    int i;
    
    nout_new = 0;
    next = *next_ptr;
    
    while(next < ntspan) {
        if(tdir * (tnew - tspan[next]) < 0) {
            if(haveeventfun && stop) {
                tout_new[nout_new] = tnew;
                for(i = 0; i < neq; i++)
                    yout_new[i+nout_new*neq] = ynew[i];
                nout_new = nout_new + 1;
            }
            break;
        }
        tout_new[nout_new] = tspan[next];
        if(tspan[next] == tnew) {
            for(i = 0; i < neq; i++)
                yout_new[i+nout_new*neq] = ynew[i];
        }
        else
            ntrp45(yout_new+nout_new*neq,tspan[next],t,y,h,f,neq);
        nout_new++;
        next++;
    }
    
    *nout_new_ptr = nout_new;
    *next_ptr = next;
}