#include <mex.h>
void uniquesorted(double *, int, int, double *, int, int);

void heapify(int i, double *ac, int a_s, double *Ik, int p, int nk, int *heap, int *ptr_heap)
{
    double *v12;
    int nv12, l, r, smallest, j12, v12_j12, comp, el, ptr_el;
    
    nv12 = 2*nk;
    l = 2*i+1;
    r = 2*i+2;
    smallest = i;
    if(l<p&&heap[l]>=0) {
        if(heap[i]<0)
            comp = -1;
        else {
            v12 = mxCalloc(nv12,sizeof(double));
            uniquesorted(v12,ptr_heap[l],ptr_heap[i],Ik,p,nk);
            for(j12=0;j12<nv12;j12++) {
                v12_j12 = (int)v12[j12];
                comp = ac[heap[l]+(v12_j12-1)*a_s]-ac[heap[i]+(v12_j12-1)*a_s];
                if(comp!=0)
                    break;
            }
            mxFree(v12);
        }
        if(comp<0)
            smallest = l;
    }
    if(r<p&&heap[r]>=0) {
        if(heap[smallest]<0)
            comp = -1;
        else {
            v12 = mxCalloc(nv12,sizeof(double));
            uniquesorted(v12,ptr_heap[r],ptr_heap[smallest],Ik,p,nk);
            for(j12=0;j12<nv12;j12++) {
                v12_j12 = (int)v12[j12];
                comp = ac[heap[r]+(v12_j12-1)*a_s]-ac[heap[smallest]+(v12_j12-1)*a_s];
                if(comp!=0)
                    break;
            }
            mxFree(v12);
        }
        if(comp<0)
            smallest = r;
    }
    if(smallest!=i) {
        el = heap[i];
        ptr_el = ptr_heap[i];
        heap[i] = heap[smallest];
        ptr_heap[i] = ptr_heap[smallest];
        heap[smallest] = el;
        ptr_heap[smallest] = ptr_el;
        heapify(smallest,ac,a_s,Ik,p,nk,heap,ptr_heap);
    }
}
