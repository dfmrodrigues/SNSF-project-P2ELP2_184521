#include "ode45_rp.h"

void rp3h(double *yout, double *tint, double *t, double *y, double *yp, double *history, int ntint, int neq, int nt)
{
    int i, j;
    
    for(j = 0; j < ntint; j++) {
        if(tint[j] < t[0]) {
            for(i = 0; i < neq; i++) {
                yout[i+j*neq] = history[i];
            }
        }
        else {
            int n;
            for(n = 0; n+1 < nt-1 && tint[j] >= t[n+1]; n++) {
            }
            ntrp3h(yout+j*neq,tint[j],t[n],y+n*neq,t[n+1],y+(n+1)*neq,yp+n*neq,yp+(n+1)*neq,neq);
        }
    }
}