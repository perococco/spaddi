/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "outstream.h"
#include "vector.h"
#include "vector_math.h"



VC_vector _VC_span(double inf,double sup,size_t nbpt,char* file,size_t line)
{
  VC_vector rslt;
  size_t i;
  double step;

  rslt=_VC_allocate(nbpt,VCTDOUBLE,file,line);
  if(!rslt) return rslt;

  step=(sup-inf)/(nbpt-1.);
  for(i=0;i<nbpt;i++)
    rslt->ddata[i]=inf+step*i;
  return rslt;
}


VC_vector _VC_sqr(VC_vector vin,VC_vector* vout,char* file,size_t line)
{
  size_t i;
  VC_vector rslt;
  if((!vout)||((!*vout)))
    rslt=_VC_allocate(vin->size,VCTDOUBLE,file,line);
  else
    rslt=*vout;

  if(rslt)
      for(i=0;i<vin->size;i++)
	rslt->ddata[i]=(vin->ddata[i])*(vin->ddata[i]);
  if(vout) *vout=rslt;
  return rslt;
}

VC_vector _VC_sqrt(VC_vector vin,VC_vector* vout,char* file,size_t line)
{
  size_t i;
  VC_vector rslt;
  if((!vout)||((!*vout)))
    rslt=_VC_allocate(vin->size,VCTDOUBLE,file,line);
  else
    rslt=*vout;

  if(rslt)
    for(i=0;i<vin->size;i++)
      rslt->ddata[i]=sqrt((vin->ddata[i]>0)?vin->ddata[i]:0);
  if(vout) *vout=rslt;
  return rslt;
}

VC_vector _VC_binSize(VC_vector vin,VC_vector* vout,char* file,size_t line)
{
  switch (vin->kind)
    {
    case VCTNONE   :
    case VCTCHAR   : return _VC_cbinSize(vin,vout,file,line);
    case VCTSHORT  : return _VC_sbinSize(vin,vout,file,line);
    case VCTINT    : return _VC_ibinSize(vin,vout,file,line);
    case VCTLONG   : return _VC_lbinSize(vin,vout,file,line);
    case VCTFLOAT  : return _VC_fbinSize(vin,vout,file,line);
    case VCTDOUBLE : return _VC_dbinSize(vin,vout,file,line);
    }
  return NULL;
}


#define _DEFINE_TEMPLATE_FUNC_(_Fc,_Data)						\
VC_vector _Fc(VC_vector vin,VC_vector* vout,char* file,size_t line)			\
{											\
  size_t i,nbpt;									\
  VC_vector rslt;									\
											\
  rslt=NULL;										\
  nbpt=vin->size;									\
  if(nbpt>1)										\
    {											\
      rslt=_VC_allocate(VCTDOUBLE,nbpt,file,line);					\
      if(rslt)										\
	{										\
	  rslt->ddata[     0]=vin->_Data[1]-vin->_Data[0];				\
	  rslt->ddata[nbpt-1]=vin->_Data[nbpt-1]-vin->_Data[nbpt-2];			\
	  for(i=1;i<nbpt-1;i++)								\
	    rslt->ddata[i]=0.5*(vin->_Data[i-1]+vin->_Data[i+1]-2*vin->_Data[i]);	\
	}										\
    }											\
  if(vout) *vout=rslt;									\
  return rslt;										\
}

_DEFINE_TEMPLATE_FUNC_(_VC_cbinSize,cdata)
_DEFINE_TEMPLATE_FUNC_(_VC_sbinSize,sdata)
_DEFINE_TEMPLATE_FUNC_(_VC_ibinSize,idata)
_DEFINE_TEMPLATE_FUNC_(_VC_lbinSize,ldata)
_DEFINE_TEMPLATE_FUNC_(_VC_fbinSize,fdata)
_DEFINE_TEMPLATE_FUNC_(_VC_dbinSize,ddata)

#undef _DEFINE_TEMPLATE_FUNC_


VC_vector _VC_rebinl(VC_vector viny,VC_vector vinx,VC_vector voutx,VC_vector* vouty,char* file,size_t line)
{
  size_t i,j,szout,szin;
  size_t idmin,idmax;
  double *dyin,*dxin,*dxout,*dyout;
  double coef,x0,y0,x1,y1;
  double xmin,xmax,tmp;
  double cmin,cmax;
  double fmin,fmax;
  double ry;
  VC_vector rslt;

  rslt=_VC_allocate(voutx->size,VCTDOUBLE,file,line);
  if(rslt)
    {
      dxin =vinx->ddata;
      dyin =viny->ddata;
      dxout=voutx->ddata;
      dyout=rslt->ddata;
      
      szin =vinx->size;
      szout=voutx->size;

      idmin=idmax=0;
      for(i=0;i<voutx->size;i++)
	{
	  if(i==0) xmin=0.5*(3*dxout[0]-dxout[1]);
	  else     xmin=0.5*(dxout[i-1]+dxout[i]);
	  
	  if(i==(szout-1)) xmax=0.5*(3*dxout[szout-1]-dxout[szout-2]);
	  else             xmax=0.5*(dxout[i]+dxout[i+1]);

	  while((idmin<szin)&&(dxin[idmin]<xmin)) idmin++;
	  idmin=(idmin>0)?idmin-1:0;
	  while((idmax<szin)&&(dxin[idmax]<xmax)) idmax++;
	  if(idmax==szin) --idmax;
	  cmin=xmin;
	  cmax=xmax;
	  fmin=dxin[idmin];
	  fmax=dxin[idmax];
	  ry=0;
	  if(xmin<fmin) //At the beginning of vinx
	    {
	      tmp=0.5*(3*dxin[0]-dxin[1]);
	      x0=(xmin>tmp)?xmin:tmp;
	      ry+=dyin[0]*(fmin-x0)/(dxin[1]-dxin[0]);
	      cmin=fmin;
	    }
	  if(xmax>fmax) //At the end of vinx
	    {
	      tmp=0.5*(3*dxin[szin-1]-dxin[szin-2]);
	      x0=(xmax<tmp)?xmax:tmp;
	      ry+=dyin[szin-1]*(x0-fmax)/(dxin[szin-1]-dxin[szin-2]);
	      cmax=fmax;
	    }


	  for(j=idmin;j<idmax;j++)
	    {
	      x0=(cmin>dxin[  j])?cmin:dxin[  j];
	      x1=(cmax<dxin[1+j])?cmax:dxin[1+j];
	      
	      coef=(dyin[j+1]-dyin[j])/(dxin[j+1]-dxin[j]);
	      y0=coef*(x0-dxin[j])+dyin[j];
	      y1=coef*(x1-dxin[j])+dyin[j];

	      ry+=0.5*(y0+y1)*(x1-x0)/(dxin[j+1]-dxin[j]);
	    }
	  rslt->ddata[i]=ry;
	}
    }
  if(vouty) *vouty=rslt;
  return rslt;
}

VC_vector _VC_rebins(VC_vector viny,VC_vector vinx,VC_vector voutx,VC_vector* vouty,char* file,size_t line)
{
  size_t i,nbout;
  double *dxout,*dyout;
  double xmin,xmax;
  VC_vector rslt,vinyd2;

  rslt=_VC_allocate(voutx->size,VCTDOUBLE,file,line);
  if(rslt)
    {
      vinyd2=VC_spline(viny,vinx,1e30,1e30,NULL);
      if(vinyd2)
	{
	  dxout=voutx->ddata;
	  dyout=rslt->ddata;
	  nbout=voutx->size;
	  for(i=0;i<nbout;i++)
	    {
	      if(i==0) xmin=0.5*(3*dxout[  0]-dxout[1]);
	      else     xmin=0.5*(  dxout[i-1]+dxout[i]);
	      
	      if(i==(nbout-1)) xmax=0.5*(3*dxout[nbout-1]-dxout[nbout-2]);
	      else             xmax=0.5*(  dxout[      i]+dxout[    i+1]);
	      
	      dyout[i]=_qgaus(splint,xmin,xmax,viny,vinx,vinyd2)/(xmax-xmin);
	    }
	  VC_free(vinyd2);
	  
	}
      else
	{
	  VC_free(rslt);
	  rslt=NULL;
	}
    }
  if(vouty) *vouty=rslt;
  return rslt;
}


VC_vector _VC_spline(VC_vector y,VC_vector x,double yp1,double ypn,
		     VC_vector *vout,char* file,size_t line)
{
  size_t i,k,n;
  double p,qn,sig,un;
  VC_vector u,y2;
  double *dx,*dy,*dy2,*du;

  n=y->size;
  y2=_VC_allocate(n,VCTDOUBLE,file,line);
  if(y2)
    {
      u=VC_allocate(n-1,VCTDOUBLE);
      if(u)
	{
	  dx=x->ddata-1;
	  dy=y->ddata-1;
	  dy2=y2->ddata-1;
	  du=u->ddata-1;
	  
	  if (yp1 > 0.99e30)
	    dy2[1]=du[1]=0.0;
	  else {
	    dy2[1] = -0.5;
	    du[1]=(3.0/(dx[2]-dx[1]))*((dy[2]-dy[1])/(dx[2]-dx[1])-yp1);
	  }
	  for (i=2;i<=n-1;i++) 
	    {
	      sig=(dx[i]-dx[i-1])/(dx[i+1]-dx[i-1]);
	      p=sig*dy2[i-1]+2.0;
	      dy2[i]=(sig-1.0)/p;
	      du[i]=(dy[i+1]-dy[i])/(dx[i+1]-dx[i]) - (dy[i]-dy[i-1])/(dx[i]-dx[i-1]);
	      du[i]=(6.0*du[i]/(dx[i+1]-dx[i-1])-sig*du[i-1])/p;
	    }
	  if (ypn > 0.99e30)
	    qn=un=0.0;
	  else {
	    qn=0.5;
	    un=(3.0/(dx[n]-dx[n-1]))*(ypn-(dy[n]-dy[n-1])/(dx[n]-dx[n-1]));
	  }
	  dy2[n]=(un-qn*du[n-1])/(qn*dy2[n-1]+1.0);
	  for (k=n-1;k>=1;k--)
	    dy2[k]=dy2[k]*dy2[k+1]+du[k];
	  VC_free(u);
	}
      else
	{
	  VC_free(y2);
	  y2=NULL;
	}
    }
  if(vout) *vout=y2;
  return y2;
}


VC_vector _VC_splint(VC_vector ya, VC_vector xa,VC_vector x,VC_vector y2a,
		     VC_vector* vout,char* file,size_t line)
{
  size_t i;
  VC_vector rslt;

  rslt=_VC_allocate(x->size,VCTDOUBLE,file,line);
  if(rslt)
    {
      for(i=0;i<rslt->size;i++)
	rslt->ddata[i]=splint(ya,xa,x->ddata[i],y2a);
    }
  if(vout) *vout=rslt;
  return rslt;
}


double splint(VC_vector ya,VC_vector xa,double x,VC_vector y2a)
{
  size_t klo,khi,k,n;
  double h,b,a;
  double *dxa,*dya,*dy2a;

  dxa =xa->ddata-1;
  dya =ya->ddata-1;
  dy2a=y2a->ddata-1;
  n=xa->size;

  klo=1;
  khi=n;
  while(khi-klo>1) 
    {
      k=(khi+klo)>>1;
      if(dxa[k]>x) khi=k;
      else klo=k;
    }
  h=dxa[khi]-dxa[klo];
  if(h==0.0)
    {
      OS_message(0,OS_ERR,"[SPLINT] Bad xa input !!\n");
      exit(1);
    }
    
  a=(dxa[khi]-x)/h;
  b=(x-dxa[klo])/h;
  return a*dya[klo]+b*dya[khi]+((a*a*a-a)*dy2a[klo]+(b*b*b-b)*dy2a[khi])*(h*h)/6.0;
}


double _qgaus(double (*func)(VC_vector ya, VC_vector xa,double x, VC_vector y2a),
	     double a, double b,VC_vector ya, VC_vector xa, VC_vector y2a)
{
	int j;
	double xr,xm,dx,s;
	static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(ya,xa,xm+dx,y2a)+(*func)(ya,xa,xm-dx,y2a));
	}
	return s *= xr;
}
/* (C) Copr. 1986-92 Numerical Recipes Software )!0. */
