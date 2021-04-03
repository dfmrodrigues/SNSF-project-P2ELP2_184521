void binterp(double *yy, double *x, double *y, double xx, int M, int N)
{
    
    int idxi, idxm, idxf, i;
    double xi;
    
    idxi = 0;
    idxf = N-1;
    
    if(xx<x[idxi]) {
        for(i = 0; i < M; i++) {
            yy[i] = y[i+idxi*M];
        }
        return;
    }
    if(xx>x[idxf]) {
        for(i = 0; i < M; i++) {
            yy[i] = y[i+idxf*M];
        }
        return;
    }
    while(idxf-idxi>1) {
        idxm = (idxf+idxi)>>1;
        if(xx>=x[idxm])
            idxi = idxm;
        else
            idxf = idxm;
    }
    
    xi = x[idxi];
    for(i = 0; i < M; i++) {
        yy[i] = y[i+idxi*M];
        yy[i] = yy[i]+(xx-xi)*(y[i+idxf*M]-yy[i])/(x[idxf]-xi);
    }
    
}