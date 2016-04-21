#ifndef SPLINE_NRH
#define SPLINE_NRH

void  NRspline(double*  x,double*  y,long n,double yp1,float ypn,double* y2);
void  NRsplint(double* xa,double* ya,double* y2a,long  n,double  x,double* y);
void NRasplint(double* xa,double* ya,double* y2a,long na,double *x,double* y,long n);

#endif
