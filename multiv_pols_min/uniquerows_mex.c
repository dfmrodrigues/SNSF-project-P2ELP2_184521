#include <mex.h>
#include <string.h>
void heapify(int, double *, int, double *, int, int, int *, int *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *a, *ac, *Ik, *c, *cout, *ind_a, *ind_c;
    int *ptr_a, *heap, *ptr_heap;
    int a_s, n, p, nk, ak_s, ptr_c, i, j, k, head, ptr_head, flag, Ik_j;
    
    a = mxGetPr(prhs[0]);
    a_s = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    ac = mxGetPr(prhs[1]);
    Ik = mxGetPr(prhs[2]);
    p = mxGetM(prhs[2]);
    nk = mxGetN(prhs[2]);
    ak_s = a_s/p;
    
    c = mxCalloc(a_s*n,sizeof(double));
    ind_a = mxCalloc(a_s,sizeof(double));
    plhs[2] = mxCreateDoubleMatrix(a_s,1,mxREAL);
    ind_c = mxGetPr(plhs[2]);
    
    ptr_c = 0;
    ptr_a = mxCalloc(p,sizeof(int));
    heap = mxCalloc(p,sizeof(int));
    ptr_heap = mxCalloc(p,sizeof(int));
    
    for(k=0;k<p;k++) {
        heap[k] = k*ak_s+ptr_a[k];
        ptr_heap[k] = k;
    }
    for(k=(p-1)/2;k>=0;k--) {
        heapify(k,ac,a_s,Ik,p,nk,heap,ptr_heap);
    }
    while(1) {
        head = ptr_heap[0];
        ptr_head = head*ak_s+ptr_a[head]+1;
        flag = (ptr_c<=0);
        for(j=0;j<nk;j++) {
            if(flag)
                break;
            Ik_j = (int)Ik[head+j*p];
            flag = (c[ptr_c-1+(Ik_j-1)*a_s]!=a[ptr_head-1+(Ik_j-1)*a_s]);
        }
        if(flag) {
            ptr_c++;
            for(j=0;j<nk;j++) {
                Ik_j = (int)Ik[head+j*p];
                c[ptr_c-1+(Ik_j-1)*a_s] = a[ptr_head-1+(Ik_j-1)*a_s];
            }
            ind_a[ptr_c-1] = ptr_head;
        }
        ind_c[ptr_head-1] = ptr_c;
        ptr_a[head]++;
        if(ptr_a[head]<ak_s)
            heap[0] = head*ak_s+ptr_a[head];
        else
            heap[0] = -1;
        heapify(0,ac,a_s,Ik,p,nk,heap,ptr_heap);
        flag = 0;
        for(k=0;k<p;k++) {
            if(flag)
                break;
            flag = (ptr_a[k]<ak_s);
        }
        if(!flag)
            break;
    }
    
    plhs[1] = mxCreateDoubleMatrix(ptr_c,1,mxREAL);
    memcpy(mxGetPr(plhs[1]),ind_a,ptr_c*sizeof(double));
    plhs[0] = mxCreateDoubleMatrix(ptr_c,n,mxREAL);
    cout = mxGetPr(plhs[0]);
    for(i=0;i<ptr_c;i++) {
        head = (ind_a[i]-1)/ak_s;
        ptr_head = (int)ind_a[i];
        for(j=0;j<nk;j++) {
            Ik_j = (int)Ik[head+j*p];
            cout[i+(Ik_j-1)*ptr_c] = a[ptr_head-1+(Ik_j-1)*a_s];
        }
    }
    
    mxFree(c);
    mxFree(ind_a);
    mxFree(ptr_a);
    mxFree(heap);
    mxFree(ptr_heap);
}
