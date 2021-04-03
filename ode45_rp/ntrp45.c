void ntrp45(double *yinterp, double tinterp, double t, double *y, double h, double *f, int neq)
{
    int i;
    double s, s2, s3, s4;
    double bi0, bi2, bi3, bi4, bi5, bi6;
    
    s = (tinterp-t)/h;
    s2 = s*s;
    s3 = s2*s;
    s4 = s3*s;
    bi0 = h*(1.0*s-183.0/64.0*s2+37.0/12.0*s3-145.0/128.0*s4);
    bi2 = h*(1500.0/371.0*s2-1000.0/159.0*s3+1000.0/371.0*s4);
    bi3 = h*(-125.0/32.0*s2+125.0/12.0*s3-375.0/64.0*s4);
    bi4 = h*(9477.0/3392.0*s2-729.0/106.0*s3+25515.0/6784.0*s4);
    bi5 = h*(-11.0/7.0*s2+11.0/3.0*s3-55.0/28.0*s4);
    bi6 = h*(3.0/2.0*s2-4.0*s3+5.0/2.0*s4);
    for(i = 0; i < neq; i++) {
        yinterp[i] = y[i];
        yinterp[i] = yinterp[i]+f[i+0*neq]*bi0;
        yinterp[i] = yinterp[i]+f[i+2*neq]*bi2;
        yinterp[i] = yinterp[i]+f[i+3*neq]*bi3;
        yinterp[i] = yinterp[i]+f[i+4*neq]*bi4;
        yinterp[i] = yinterp[i]+f[i+5*neq]*bi5;
        yinterp[i] = yinterp[i]+f[i+6*neq]*bi6;
    }
}