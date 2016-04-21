/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdlib.h>
#include <math.h>
#include "outstream.h"
#include "memory_manager.h"
#include "vector.h"
#include "vector_math.h"
#include "vector_sort.h"
#include "table.h"
#include "spectrum.h"
#include "astro.h"

#define _SQ(x) ((x)*(x))


VC_vector _cosmicMask(VC_vector flux,VC_vector sigm,double clip,size_t* nbc,VC_vector* mask,char* file,size_t line)
{
  VC_vector rslt;
  double    kernel[]={0.0625,0.25,0.375};
  size_t i,nf;
  double   *dr,*df,*ds;
  
  *nbc=0;
  nf=flux->size;
  if(mask&&(*mask)&&((*mask)->size==nf)&&((*mask)->kind==VCTDOUBLE))
    rslt=*mask;
  else
    rslt=_VC_allocate(nf,VCTDOUBLE,file,line);

  if(rslt)
    {
      dr=rslt->ddata;
      df=flux->ddata;
      ds=sigm->ddata;
      for(i=0;i<rslt->size;i++)
	dr[i]=df[i]*(1-kernel[2]);
      
      dr[0]-=kernel[0]*(df[0]+df[2])+kernel[1]*(df[0]+df[1]);
      dr[1]-=kernel[0]*(df[0]+df[3])+kernel[1]*(df[0]+df[2]);
      
      dr[nf-2]-=kernel[0]*(df[nf-4]+dr[nf-1])+kernel[1]*(df[nf-3]+df[nf-1]);
      dr[nf-1]-=kernel[0]*(df[nf-3]+dr[nf-1])+kernel[1]*(df[nf-2]+df[nf-1]);

      for(i=2;i<nf-2;i++)
	dr[i]-=kernel[0]*(df[i-2]+df[i+2])+kernel[1]*(df[i-1]+df[i+1]);

      for(i=0;i<nf;i++)
	if(ds[i]<=0) {dr[i]=0.;continue;}
	else
	  {
	    if(dr[i]>clip*ds[i])
	      {
		dr[i]=1.;
		++(*nbc);
	      }
	    else
	      dr[i]=0.;
	  }
    }
  if(mask) *mask=rslt;
  return rslt;

}

VC_vector correctFlux(VC_vector* flux,VC_vector mask)
{
  long i,inf,sup,nbi,nbs;
  double *dm;
  double yinf,ysup;
  long width=3;
  long *idxinf,*idxsup;

  idxinf=(long *)MM_malloc(width*sizeof(long));
  idxsup=(long *)MM_malloc(width*sizeof(long));

  dm=mask->ddata;
  for(i=0;i<(*flux)->size;i++)
    {
      if(!dm[i]) continue;
      inf=i;nbi=0;
      while((--inf>=0)&&(nbi<width)) if(!dm[inf]) {idxinf[nbi++]=inf;}

      sup=i;nbs=0;
      while((++sup<(*flux)->size)&&(nbs<width)) if(!dm[sup]) {idxsup[nbs++]=sup;}

      if(0==nbi)
	{
	  if(nbs>0) (*flux)->ddata[i]=(*flux)->ddata[idxsup[nbs-1]];
	  continue;
	}

      if(0==nbs)
	{
	  if(nbi>0) (*flux)->ddata[i]=(*flux)->ddata[idxinf[nbi-1]];
	  continue;
	}
      nbi=idxinf[nbi-1];
      nbs=idxsup[nbs-1];
      yinf=(*flux)->ddata[nbi];
      ysup=(*flux)->ddata[nbs];
      (*flux)->ddata[i]=(ysup-yinf)/(double)(nbs-nbi)*(i-nbi)+yinf;
    }
  MM_free(idxsup);
  MM_free(idxinf);

  return *flux;
}


VC_vector _addiScaleFactors(SP_spectrum* spectra,size_t nbsp,double width,VC_vector* vout,char* file,size_t line)
{
  int done;
  size_t i,j,k,l,sz,msz,swork;
  double *df;
  double *dw;
  double *work;
  VC_vector rslt,idmin,idmax,lwave,fctrs;
  size_t spref,nbchunk,cmin,cmax,rmin,rmax;
  double medref,aux;
  double wmin,wmax;
  double chkmin,chkmax;
  
  
  rslt  =_VC_allocate(nbsp,VCTDOUBLE,file,line);
  VCV_allocate(idmin,nbsp,VCTLONG);
  VCV_allocate(idmax,nbsp,VCTLONG);
/*   printf("%d\n",(long)idmin); */
  if(rslt)
    {
      /*************************/
      /* Create working space */
      /*************************/
      msz=0;
      for(i=0;i<nbsp;i++)
	{
	  sz=spectra[i]->npix;
	  if(sz>msz) msz=sz;
	}
      MMV_malloc(work,msz,double);
      swork=0;
      /****************************************/
      /*First find the limits of each spectrum*/
      /****************************************/
      for(i=0;i<nbsp;i++)
	{
	  VC_getIdxEdge(spectra[i]->data[spectra[i]->colF],idmin->ldata+i,idmax->ldata+i);

	  if(idmin->ldata[i]>=idmax->ldata[i])
	    {
	      OS_message(0,OS_ERR,"[addiScaleFactors] No data in spectrum #%ld\n",i);
	      exit(1);
	    }
	}
      
      /*********************************************/
      /*Find the spectrum with the greatest median */
      /*********************************************/
      OS_message(0,OS_STD,"[addiScaleFactors] Search for spectra with the greatest median\n");
      spref=0;
      medref=1;
      for(i=0;i<nbsp;i++)
	{
	  sz=spectra[i]->npix;
	  df=spectra[i]->data[spectra[i]->colF]->ddata;
	  for(j=idmin->ldata[i],k=0;j<idmax->ldata[i];j++,k++)
	    work[k]=df[j];
	  swork=idmax->ldata[i]-idmin->ldata[i];
	  
	  aux=median(work,swork);
	  OS_message(0,OS_STD,"[addiScaleFactors] #%2d : %f\n",i+1,aux);
	  if((0==i)||(aux>medref))
	    {
	      medref=aux;
	      spref=i;
	    }
	}
      OS_message(0,OS_STD,"[addiScaleFactors] ------------------------------------------\n");
      OS_message(0,OS_STD,"[addiScaleFactors] Use spectrum #%d as reference\n",spref+1);

      wmin=spectra[spref]->data[spectra[spref]->colW]->ddata[idmin->ldata[spref]];
      wmax=spectra[spref]->data[spectra[spref]->colW]->ddata[idmax->ldata[spref]-1];

      nbchunk=(size_t)((wmax-wmin)/width);
      nbchunk=(nbchunk<1)?1:nbchunk;
      VCV_allocate(fctrs,nbchunk,VCTDOUBLE);
      lwave=VC_span(wmin,wmax,nbchunk+1);

      OS_message(0,OS_STD,"[addiScaleFactors]\n");
      OS_message(0,OS_STD,"[addiScaleFactors] Compute scale factors for each spectrum \n");
      for(i=0;i<nbsp;i++)
	{
	  dw=spectra[i]->data[spectra[i]->colW]->ddata;
	  df=spectra[i]->data[spectra[i]->colF]->ddata;
	  wmin=dw[idmin->ldata[i]];
	  wmax=dw[idmax->ldata[i]-1];
	  if(i!=spref)
	    {
	      done=0;
	      for(j=0;j<nbchunk;j++)
		{
		  chkmin=(wmin>lwave->ddata[  j])?wmin:lwave->ddata[  j];
		  chkmax=(wmax<lwave->ddata[1+j])?wmax:lwave->ddata[1+j];
		  
		  if(chkmin>=chkmax) continue;

		  rmin=findNear(spectra[spref]->data[spectra[spref]->colW],chkmin);
		  rmax=findNear(spectra[spref]->data[spectra[spref]->colW],chkmax)+1;
		  cmin=findNear(spectra[    i]->data[spectra[    i]->colW],chkmin);
		  cmax=findNear(spectra[    i]->data[spectra[    i]->colW],chkmax)+1;
		  for(k=rmin,l=0;k<rmax;k++,l++)
		    work[l]=spectra[spref]->data[spectra[spref]->colF]->ddata[k];
		  medref=median(work,rmax-rmin);

		  for(k=cmin,l=0;k<cmax;k++,l++)
		    work[l]=df[k];
		  aux=median(work,cmax-cmin);

		  if((aux>0)&&(medref>0))
		    fctrs->ddata[done++]=medref/aux;
		  if((medref<=0)&&(aux<=0))
		    fctrs->ddata[done++]=1.;
		}
	      if(done<=0)
		{
		  OS_message(0,OS_STD,"[addiScaleFactors] no common region. Set factor to 1\n");
		  rslt->ddata[i]=1.;
		}
	      else
		rslt->ddata[i]=median(fctrs->ddata,done);
	    }
	  else
	    rslt->ddata[i]=1.;

	  OS_message(0,OS_STD,"[addiScaleFactors] Factor for spectrum #%2d : %f\n",i+1,rslt->ddata[i]);
	}
      VC_free(lwave);
      VC_free(fctrs);
      MM_free(work);
    }
  VC_free(idmax);
  VC_free(idmin);
  if(vout) *vout=rslt;
  return rslt;
}

VC_vector _mergeScaleFactors(SP_spectrum* spectra,size_t nbsp,long width,VC_vector* vout,char* file,size_t line)
{
  size_t    i,j,k,sz,msz,swork,i1,i2;
  double    *df1,*dw1,*df2,*dw2,*work;
  VC_vector rslt,idmin,idmax,sidx;
  size_t    spref,r1min,r1max,r2min,r2max;
  double    med1,med2;
  double    wmin1,wmax1,wmin2,wmax2,wmin,wmax;
  
  
  rslt  =_VC_allocate(nbsp,VCTDOUBLE,file,line);
  for(i=0;i<rslt->size;i++) rslt->ddata[i]=1;

  if(nbsp>1)
    {
      idmin =VC_allocate(nbsp,VCTLONG);
      idmax =VC_allocate(nbsp,VCTLONG);
      sidx  =VC_allocate(nbsp,VCTLONG);
      if(rslt)
	{
	  /************************/
	  /* Create working space */
	  /************************/
	  msz=0;
	  for(i=0;i<nbsp;i++)
	    {
	      sz=spectra[i]->npix;
	      if(sz>msz) msz=sz;
	    }
	  work=(double*)MM_malloc(msz*sizeof(double));
	  swork=0;
	  /****************************************/
	  /*First find the limits of each spectrum*/
	  /****************************************/
	  for(i=0;i<nbsp;i++)
	    {
	      VC_getIdxEdge(spectra[i]->data[spectra[i]->colF],idmin->ldata+i,idmax->ldata+i);
	      
	      if(idmin->ldata[i]>=idmax->ldata[i])
		{
		  OS_message(0,OS_ERR,"[mergeScaleFactors] No data in spectrum #%ld\n",i);
		  exit(1);
		}
	      sidx->ldata[i]=i; //initialise sort array
	    }
	  
	  /****************************************/
	  /*Sort spectra in increasing wavelength */
	  /****************************************/
	  OS_message(0,OS_STD,"[mergeScaleFactors] Sort spectra in increasing wavelength\n");
	  for(i=0;i<nbsp-1;i++)
	    for(j=i+1;j<nbsp;j++)
	      {
		i1=sidx->ldata[i];
		i2=sidx->ldata[j];
		if((spectra[i1]->data[spectra[i1]->colW]->ddata[idmin->ldata[i1]]) <=
		   (spectra[i2]->data[spectra[i2]->colW]->ddata[idmin->ldata[i2]]))
		  continue;
		sidx->ldata[i]=i2;
		sidx->ldata[j]=i1;
	      }

	  spref=sidx->ldata[0];
	  OS_message(0,OS_STD,"[mergeScaleFactors] ------------------------------------------\n");
	  OS_message(0,OS_STD,"[mergeScaleFactors] Use spectrum #%d as reference\n",spref+1);
	  
	  rslt->ddata[spref]=1.;
	  OS_message(0,OS_STD,"[mergeScaleFactors]\n");
	  OS_message(0,OS_STD,"[mergeScaleFactors] Compute scale factors for each spectrum \n");
	  for(i=1;i<nbsp;i++)
	    {
	      i1=sidx->ldata[i-1];
	      i2=sidx->ldata[  i];

	      dw1=spectra[i1]->data[spectra[i1]->colW]->ddata;
	      df1=spectra[i1]->data[spectra[i1]->colF]->ddata;

	      dw2=spectra[i2]->data[spectra[i2]->colW]->ddata;
	      df2=spectra[i2]->data[spectra[i2]->colF]->ddata;

	      /************************/
	      /* Find overlap regions */
	      /************************/
	      wmin1=dw1[idmin->ldata[i1]];
	      wmax1=dw1[idmax->ldata[i1]-1];

	      wmin2=dw2[idmin->ldata[i2]];
	      wmax2=dw2[idmax->ldata[i2]-1];

	      wmin=(wmin1>wmin2)?wmin1:wmin2;
	      wmax=(wmax1<wmax2)?wmax1:wmax2;

	      /***************************/
	      /* Get index of the limits */
	      /* of the overlap          */
	      /***************************/

	      if(wmin>=wmax)
		{
		  r1max=idmax->ldata[i1];
		  r1min=r1max-10*width;
		  r2min=idmin->ldata[i2];
		  r2max=r2min+10*width;
		}
	      else
		{
		  r1min=findNear(spectra[i1]->data[spectra[i1]->colW],wmin);
		  r1max=findNear(spectra[i1]->data[spectra[i1]->colW],wmax);
		  r2min=findNear(spectra[i2]->data[spectra[i2]->colW],wmin);
		  r2max=findNear(spectra[i2]->data[spectra[i2]->colW],wmax);
		  r1max++;
		  r2max++;
		}
	      
	      /**********************/
	      /* Compute median for */
	      /* the overlap        */
	      /**********************/
	      for(j=0,k=r1min;k<r1max;j++,k++) work[j]=df1[k];
	      med1=median(work,r1max-r1min);

	      for(j=0,k=r2min;k<r2max;j++,k++) work[j]=df2[k];
	      med2=median(work,r2max-r2min);
	      
	      if((med1>0)&&(med2>0)) rslt->ddata[i2]=med1/med2*rslt->ddata[i1];

	      OS_message(0,OS_STD,"[mergeScaleFactors] Factor for spectrum #%2d : %f\n",i+1,rslt->ddata[i]);
	    }
	  MM_free(work);
	}
      VC_free(sidx);
      VC_free(idmax);
      VC_free(idmin);
    }
  if(vout) *vout=rslt;
  return rslt;
}




size_t findNear(VC_vector vin,double val)
{
  size_t kmin,kmax,kmid;
  kmin=0;kmax=vin->size;
  while((kmax-kmin)>1)
    {
      kmid=(kmax+kmin)>>1;
      if(vin->ddata[kmid]>val) kmax=kmid;
      else kmin=kmid;
    }
  return kmin;
}


double median(double* tab,size_t size)
{
  size_t half;
  if(size==0) return 0;
  if(size==1) return tab[0];
  if(size==2) return 0.5*(tab[0]+tab[1]);

  qsort(tab,size,sizeof(double),_super);
  half=size>>1;
  if(size%2)
    return tab[half];
  else
    return 0.5*(tab[half]+tab[half-1]);
}

int fitLine(double* y,double* x,size_t size,double* a,double* b)
{
  double mx,mx2,my,mxy,n;
  size_t i;
  
  mx=mx2=my=mxy=0;
  n=(double)size;
  for(i=0;i<size;i++)
    {
      mx+=x[i];
      my+=y[i];
      mx2+=x[i]*x[i];
      mxy+=x[i]*y[i];
    }
  mx/=n;my/=n;mx2/=n;mxy/=n;
  *a=(mxy-mx*my)/(mx2-mx*mx);
  *b=my-*a*mx;
  return 1;
}

int fitLineIter(double* y,double* x,size_t size,double* a,double* b)
{
  size_t i,nbok,oldnbok,nbiter,nbitermax;
  double rms,mean,yy,norm;

  nbok=size;
  norm=1./(double)size;
  nbiter=0;nbitermax=100;
  while(++nbiter<=nbitermax)
    {
      oldnbok=nbok;
      fitLine(y,x,size,a,b);
      rms=mean=0.;
      for(i=0;i<size;i++)
	{
	  yy=y[i]-x[i]**a-*b;
	  mean+=yy;
	  rms+=yy*yy;
	}
      mean*=norm;
      rms*=norm;
      rms=sqrt(rms-mean*mean);
      nbok=0;
      for(i=0;i<size;i++)
	{
	  yy=y[i]-x[i]**a-*b;
	  yy=(yy<0)?-yy:yy;
	  if(yy<=3*rms) nbok++;
	}
      if(nbok==oldnbok) break;
    }
  return 1;
}

int _super(const void* a,const void* b)
{
  double ra,rb;
  ra=*(double*)a;
  rb=*(double*)b;
  if(ra>rb) return 1;
  if(ra<rb) return -1;
  return 0;
}

VC_vector _sigmaClipping(VC_vector flux,VC_vector* mask,double clip,int verb,char *file,size_t line)
{
  size_t i;
  VC_vector rslt;
  double mean,rms,aux;
  double dev,devmax;
  long idevmax,nsel;

  if(mask&&*mask)
    rslt=*mask;
  else
    {
      rslt=VC_allocate(flux->size,VCTDOUBLE);
      if(rslt) for(i=0;i<rslt->size;i++) rslt->ddata[i]=1.;
    }
  if(!rslt) 
    {
      return NULL;
    }
  nsel=0;
  for(i=0;i<rslt->size;i++) if(rslt->ddata[i]) nsel++;  
  
  while(nsel>2)
    {
      devmax=0.;
      idevmax=0;
      for(i=0;i<rslt->size;i++)
	{
	  aux=rslt->ddata[i];
	  rslt->ddata[i]=0.;
	  statistique(flux->ddata,flux->size,rslt->ddata,&mean,&rms);
	  rslt->ddata[i]=aux;
	  if(rms<=0) continue;
	  dev=fabs((flux->ddata[i]-mean)/rms);
	  if(dev>devmax)
	    {
	      devmax=dev;
	      idevmax=i;
	    }
	}
      if(devmax<clip) break;
      rslt->ddata[idevmax]=0.;
      nsel--;
    }
  return rslt;
}


VC_vector _dataClipping(VC_vector flux,VC_vector sigm,VC_vector* mask,double clip,char* file,size_t line)
{
  size_t i;
  VC_vector rslt;
  size_t idmax;
  double fsum,wsum,devmax,dev,weight,err2,flx;
  size_t nsel;
  int pass,mean;

  mean=0;
  if(mask&&*mask)
    rslt=*mask;
  else
    {
      rslt=VC_allocate(flux->size,VCTDOUBLE);
      if(rslt) for(i=0;i<rslt->size;i++) rslt->ddata[i]=1;
    }
  if(!rslt) return NULL;


  nsel=0;for(i=0;i<rslt->size;i++) if(rslt->ddata[i]) nsel++;
  
  if(!nsel)
    {
      for(i=0;i<rslt->size;i++) rslt->ddata[i]=1;
      return rslt;
    }
  
  pass=0;
  while(nsel>2)
    {
      fsum=wsum=0;
      nsel=0;
      for(i=0;i<rslt->size;i++)
	{
	  if(pass==0) rslt->ddata[i]=1;
	  if(rslt->ddata[i]) 	    //	    lidx[nsel++]=i;
	    if(sigm->ddata[i])
	      {
		if(mean)
		  {
		    weight=flux->ddata[i];
		    fsum+=weight;
		    wsum+=weight*weight;
		  }
		else
		  {
		    weight=1./(sigm->ddata[i]*sigm->ddata[i]);
		    wsum+=weight;
		    fsum+=flux->ddata[i]*weight;
		  }
		nsel++;
	      }
	}
      if(pass==0) {nsel=rslt->size;pass=1;}
      

      if(wsum<=0)
	{
	  for(i=0;i<rslt->size;i++) rslt->ddata[i]=0;
	  return rslt;
	}

      if(mean)
	{
	  flx=fsum/(double)nsel;
	  err2=sqrt((wsum/(double)nsel)-flx*flx);
	  if(err2<=0) err2=1e-10;
	}
      else
	{
	  flx=fsum/wsum;
	  err2=1./wsum;
	}
      
      devmax=0.;
      idmax=0;
      for(i=0;i<rslt->size;i++)
	if(rslt->ddata[i])
	  {
	    if(mean)
	      dev=fabs(flux->ddata[i]-flx)/err2;
	    else
	      dev=fabs(flux->ddata[i]-flx)/sqrt(err2+sigm->ddata[i]*sigm->ddata[i]);
	    if(dev>devmax)
	      {
		devmax=dev;
		idmax=i;
	      }
	  }
      if(devmax<clip) break;
      rslt->ddata[idmax]=0;
      nsel--;
    }
  return rslt;
}

int VC_getIdxEdge(VC_vector data,long* idmin,long* idmax)
{
  long   imin,imax;
  double *df;

  df=data->ddata;
  imin=0;
  while((imin<data->size)&&(df[imin]<=0)) imin++;
  
  imax=data->size-1;
  while((imax>=0)&&(df[imax]<=0)) imax--;
  imax++;

  *idmin=imin;
  *idmax=imax;
  return 0;
}


int statistique(double* data,long n,double* mask,double* mean,double* rms)
{
  long i,nsel;
  double dsum,d2sum;
  
  nsel=0;
  dsum=d2sum=0.;
  for(i=0;i<n;i++)
    {
      if(mask[i]<=0.) continue;
      dsum +=data[i];
      d2sum+=data[i]*data[i];
      nsel++;
    }
  if(nsel)
    {
      dsum/=(double)nsel;
      d2sum=sqrt(d2sum/(double)nsel-dsum*dsum);
    }
  *mean=dsum;
  *rms=d2sum;
  return nsel;
}
