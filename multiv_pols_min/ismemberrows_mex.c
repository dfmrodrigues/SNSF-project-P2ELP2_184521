#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B;
    mxLogical *lia;
    int A_s, B_s, n, i, j, k, flag;
    
    A = mxGetPr(prhs[0]);
    A_s = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    B = mxGetPr(prhs[1]);
    B_s = mxGetM(prhs[1]);
    
    plhs[0] = mxCreateLogicalMatrix(A_s,1);
    lia = mxGetLogicals(plhs[0]);
    
    j = 0;
    for(i=0;i<A_s;i++) {
        flag = 1;
        for(k=0;k<n;k++) {
            if(!flag)
                break;
            flag = (A[i+k*A_s]==B[j+k*B_s]);
        }
        if(flag) {
            lia[i] = true;
            j++;
        }
        if(j>=B_s)
            break;
    }
    
}
