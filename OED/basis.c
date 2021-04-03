#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int nparam, poly_num;

double *ords, *param_basis_eval, *param_basis_deval;

double *basis_eval, *basis_deval;

int j, k, l;

nparam = mxGetScalar(prhs[0]);
poly_num = mxGetScalar(prhs[1]);
ords = mxGetData(prhs[2]);
param_basis_eval = mxGetPr(prhs[3]);
param_basis_deval = mxGetPr(prhs[4]);

plhs[0] = mxCreateDoubleMatrix(1,poly_num,mxREAL);
basis_eval = mxGetPr(plhs[0]);
plhs[1] = mxCreateDoubleMatrix(nparam,poly_num,mxREAL);
basis_deval = mxGetPr(plhs[1]);

for(j = 0; j < poly_num; j++) {
    basis_eval[j] = 1;
    for(l = 0; l < nparam; l++) {
        basis_deval[l+j*nparam] = 1;
    }
    for(k = 0; k < nparam; k++) {
        int ords_kj = (int)ords[k+j*nparam];
        basis_eval[j] = basis_eval[j]*param_basis_eval[k+ords_kj*nparam];
        for(l = 0; l < nparam; l++) {
            if(l == k)
                basis_deval[l+j*nparam] = basis_deval[l+j*nparam]*param_basis_deval[k+ords_kj*nparam];
            else
                basis_deval[l+j*nparam] = basis_deval[l+j*nparam]*param_basis_eval[k+ords_kj*nparam];
        }
    }
}
}
