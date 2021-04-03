#include "ocp_mex.h"

void odefun(double *f, double t, double *y, int neq, const mxArray *args[])
{
    args = args+1;
    ffun(f,NULL,y,args);
}