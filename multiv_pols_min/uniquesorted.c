#include <mex.h>

void uniquesorted(double *v12, int k1, int k2, double *Ik, int p, int nk)
{
    int nv12, j1, j2, j12;
    
    nv12 = 2*nk;
    j1 = 0;
    j2 = 0;
    for(j12=0;j12<nv12;j12++) {
        if(j1<nk&&(j2>=nk||Ik[k1+j1*p]<Ik[k2+j2*p])) {
            v12[j12] = Ik[k1+j1*p];
            j1++;
        }
        else {
            v12[j12] = Ik[k2+j2*p];
            j2++;
        }
    }
}
