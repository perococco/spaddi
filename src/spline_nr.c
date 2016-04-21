#include <stdlib.h>
#include <stdio.h>
#include "memory_manager.h"
#include "spline_nr.h"

void NRspline(double*  x,double*  y,long n,double yp1,float ypn,double* y2)
{
  int i,k;
  double p,qn,sig,un,*u;
  
  MMV_malloc(u,n,double);
  if(!u)
    {
      printf("Not enough memory in spline.  %s at line %d\n",__FILE__,__LINE__);
      exit(1);
    }
  --x;--y;--y2;--u;
  if(yp1>0.99e30) 
    y2[1]=u[1]=0.0;
  else
    {
      y2[1]=-0.5;
      u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
  
  for(i=2;i<=n-1;i++)
    {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
  if (ypn>0.99e30)
    qn=un=0.0;
  else
    {
      qn=0.5;
      un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for(k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  ++u;++y2;++y;++x;
  MM_free(u);
}

void NRsplint(double* xa,double* ya,double* y2a,long n,double x,double* y)
{
  int klo,khi,k;
  float h,b,a;
  
  klo=1;
  khi=n;
  while(khi-klo>1)
    {
      k=(khi+klo)>>1;
      if(xa[k-1]>x) khi=k;
      else klo=k;
    }
  --khi;--klo;
  h=xa[khi]-xa[klo];
  if(h==0.0)
    {
      printf("Bad xa input to splint.  %s at line %d\n",__FILE__,__LINE__);
      exit(1);
    }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void NRasplint(double* xa,double* ya,double *y2a,long na,double* x,double* y,long n)
{
  int klo,khi,k;
  float h,b,a;
  
  klo=1;
  khi=na;
  while(khi-klo>1)
    {
      k=(khi+klo)>>1;
      if(xa[k-1]>x[0]) khi=k;
      else klo=k;
    }
  --khi;--klo;
  for(k=0;k<n;k++)
    {
      y[k]=0.;
      while((khi<na)&&((x[k]<xa[klo])||(x[k]>xa[khi]))) {++khi;++klo;}
      if(khi>=na) break;
      h=xa[khi]-xa[klo];
      if(h==0.0)
	{
	  printf("Bad xa input to splint.  %s at line %d\n",__FILE__,__LINE__);
	  exit(1);
	}
      a=(xa[khi]-x[k])/h;
      b=(x[k]-xa[klo])/h;
      y[k]=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
    }
}
