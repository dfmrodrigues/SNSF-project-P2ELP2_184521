#include <math.h>
#include <string.h>
#include <mex.h>
#define abs(x) (((x)<0) ? -(x) : (x))
#define max(x,y) (((x)>(y)) ? (x) : (y))
#define min(x,y) (((x)<(y)) ? (x) : (y))
#define sign(x) (((x)>0)-((x)<0))

void odefun(double*, double, double*, int, const mxArray**);
void eventfun(double, double*, double*, double*, int, int, const mxArray**);
void rp45(int*, int*, double*, double*, int, int, int, int, double*, double, double*, double, double*, double, double*, int);
void rp3h(double*, double*, double*, double*, double*, double*, int, int, int);
void ntrp45(double*, double, double, double*, double, double*, int);
void ntrp3h(double*, double, double, double*, double, double*, double*, double*, int);
void binterp(double*, double*, double*, double, int, int);
