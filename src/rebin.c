/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdlib.h>
#include <math.h>
#include "outstream.h"
#include "vector.h"
#include "rebin.h"


#define _SQ(a) ((a)*(a))

int _rebinl(VC_vector signal,VC_vector noise,VC_vector xinput,
	    VC_vector xoutput,VC_vector* sout,VC_vector* nout,
	    char* file,size_t line)
{
  VC_vector sgout,nsout;
  size_t    i,szout,szin;
  long int  kmin,kmax;
  double    xmin,xmax;
  double    *xout,*xin;

  szout=xoutput->size;
  szin =xinput->size;
  
  xout=xoutput->ddata;
  xin=xinput->ddata;

  sgout=_VC_allocate(szout,VCTDOUBLE,file,line);
  if(!sgout) 
    return 1;
  
  nsout=_VC_allocate(szout,VCTDOUBLE,file,line);
  if(!nsout)
    {
      VC_free(sgout);
      return 1;
    }
  /***********************************/
  /*Initialisation of sgout and nsout*/
  /***********************************/
  kmin=-1;kmax=-1;
  for(i=0;i<szout;i++) 
    {
      sgout->ddata[i]=0;
      nsout->ddata[i]=0;

      if(0==i)         xmin=0.5*(3*xout[0]-xout[1]);
      else             xmin=0.5*(xout[i-1]+xout[i]);
  
      if(i==(szout-1)) xmax=0.5*(3*xout[szout-1]-xout[szout-2]);
      else             xmax=0.5*(  xout[i+1]    +xout[i]);

      if((xmax<=xin[0])||(xmin>=xin[szin-1])) continue;
      
      rebinLocal(signal->ddata,noise->ddata,xin,szin,
		 xmin,xmax,&kmin,&kmax,
		 sgout->ddata+i,nsout->ddata+i);
    }
  *sout=sgout;
  *nout=nsout;
   
  return 0;
}


int rebinLocal(double *signal,double  *noise,double* xin,size_t szin,
	       double    xmin,double    xmax,
	       long int *kmin,long int *kmax,
	       double   *sout,double   *nout)
{
  long   k;
  double fact1,fact2,tot;
  double ximin,ximax;
  double xcmin,xcmax;
  double delt;
  double sig,noi;
  long   imin,imax;

  if((xin[0]>xmax)||(xin[szin-1]<xmin)) return 1;


  if((*kmin<0)||(*kmax<0)) 
    {
      imin=0;
      imax=0;
    }
  else 
    {
      imin=*kmin;
      imax=*kmax;
    }
  
  imin=(imin)?imin-1:0;
  while((imin<szin)&&(xmin>=xin[imin])) imin++;

  if(imin==0)
    xmin=(xmin>xin[0])?xmin:xin[0];
  else
    imin--;
      
  imax=(imax)?imax-1:0;
  while((imax<szin)&&(xmax>=xin[imax])) imax++;
  imax--;
  if(imax==(szin-1))
    {
      xmax=(xmax<xin[szin-1])?xmax:xin[szin-1];
      imax--;
    }
  *kmin=imin;
  *kmax=imax;

  *sout=*nout=0;
  if(imax<imin) return 1;

  fact1=fact2=tot=0.;
  noi=sig=0.;
  for(k=imin;k<=imax;k++)
    {
      ximin=xin[  k];
      ximax=xin[1+k];
      delt=0.5/(ximax-ximin);

      xcmin=(ximin>xmin)?ximin:xmin;
      xcmax=(ximax<xmax)?ximax:xmax;

      fact1+= delt* (-_SQ(ximax - xcmax) + _SQ(ximax - xcmin));
      fact2 = delt*( _SQ(xcmax - ximin) - _SQ(xcmin - ximin));

      tot+=fact1;
      sig+=fact1*signal[k];
      noi+=_SQ(fact1*noise[k]);

      fact1=fact2;
    }
  tot+=fact1;
  sig+=fact1*signal[imax+1];
  noi+=_SQ(fact1*noise[imax+1]);

  if(tot)
    delt=1./tot;//delt=1./(xmax-xmin);
  else
    delt=0.;
  
  sig=delt*sig;
  noi=delt*sqrt(noi);

  *sout=sig;
  *nout=noi;
  
  return 0;
}




#undef _SQ
