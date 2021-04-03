void ntrp3h(double *yinterp, double tinterp, double t, double *y, double tnew, double *ynew, double *yp, double *ypnew, int neq)
{
    int i;
    double h, s, s2, s3, s4;
    
    h = tnew-t;
    s = (tinterp-t)/h;
    s2 = s*s;
    s3 = s2*s;
    for(i = 0; i < neq; i++) {
        double slope_i, c_i, d_i;
        slope_i = (ynew[i] - y[i])/h;
        c_i = 3*slope_i - 2*yp[i] - ypnew[i];
        d_i = yp[i] + ypnew[i] - 2*slope_i;
        yinterp[i] = y[i];
        yinterp[i] = yinterp[i]+h*d_i*s3;
        yinterp[i] = yinterp[i]+h*c_i*s2;
        yinterp[i] = yinterp[i]+h*yp[i]*s;
    }
}