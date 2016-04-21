/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "memory_manager.h"
#include "outstream.h"
#include "vector.h"
#include "vector_math.h"
#include "vector_sort.h"
#include "table.h"
#include "table_tfits.h"
#include "spectrum.h"
#include "astro.h"
#include "rebin.h"
#include "spline_nr.h"

#define IS_NAN(x) (x!=x)
#define _SQ(a) ((a)*(a))




SP_spectrum _SP_addSpectra(SP_spectrum* spectra,size_t nbsp,double width,double clip,double prc,int error_type,double rebin_pixfactor,double rebin_pixshift,int spline,int tclip,int keepnbsp,char* file,size_t line)
{

  VC_vector   group,factors;
  double*     work,fct;
  double*     tmp;
  double      *dw,*df,*ds;
  double      fact1,fact2;
  double      *wmin,*wmax;
  VC_vector    wlim,wave;
  SP_spectrum rslt;
  VC_vector   eflux,esigm,emask;
  VC_vector   flx1,flx2,wav1,pds1,pds2;
  long        *npix;
  double      y1,y2;
  double       tmin, tmax,tcen,wpas,old_wpas;
  double       w1min,w1max,w2min,w2max;
  long         k1min,k1max,k2min,k2max;
  long         imin,imax,jmin,jmax;
  long        *idmin,*idmax; 
  long        *kmin,*kmax; 
  long        nbgroup,igroup,nbcol;
  long        i,j;
  size_t      si,sj,k,ngood;
  size_t      npix_max,npix_tot;
  int         pourc,old_pourc;
  double      fsum,wsum,weight,sigm_srianand,mean_srianand;
  long        nbg_srianand;
  long        nbg_bastien;
  double      a1,b1;
  SP_spectrum* spline_spectra;

  if(nbsp<2) return NULL;

  /****************************/
  /*  Sort and group by setup */
  /****************************/
  SP_sortSpectrum(spectra,nbsp); 
  OS_message(0,OS_STD,"[SP_addSpectra] new order for spectra :\n");
  for(i=0;i<nbsp;i++)
    {
      OS_message(0,OS_STD,"[SP_addSpectra] %4d : %s\n",i+1,spectra[i]->name);
      
    }
  group=SP_groupSpectrum(spectra,nbsp,prc);
  if(!group)
    {
      OS_message(0,OS_ERR,"not enough memory !!\n");
      exit(1);
    }
  nbgroup=1+group->ldata[nbsp-1];
  OS_message(0,OS_STD,"[SP_addSpectra] found %d group%s\n",nbgroup,(nbgroup>1)?"s":" ");
  
  
  npix_max=0; /*will contain the number of pixel of the biggest spectrum*/
  for(si=0;si<nbsp;si++) 
    npix_max=(spectra[si]->npix>npix_max)?spectra[si]->npix:npix_max;
   
  MMV_malloc(tmp,npix_max,double); /*working buffer*/
  MMV_malloc(idmin,nbsp,long);
  MMV_malloc(idmax,nbsp,long);
  MMV_malloc( kmin,nbsp,long);
  MMV_malloc( kmax,nbsp,long);
  MMV_malloc( wmin,nbsp,double);
  MMV_malloc( wmax,nbsp,double);
  for(i=0;i<nbsp;i++) 
    {
      /**********************************************************/
      /*Find first and last pixel that not have null flux values*/
      /*sort the corresponding wavelength in wmin and wmax and  */
      /*the pixel index in idmin and idmax.                     */
      /**********************************************************/
      VC_getIdxEdge(spectra[i]->data[spectra[i]->colF],idmin+i,idmax+i);
      wmin[i]=spectra[i]->wave[idmin[i]];
      wmax[i]=spectra[i]->wave[idmax[i]-1];
      kmin[i]=kmax[i]=-1;
    }
  for(i=0;i<nbgroup;i++)
    {
      k=0;
      for(j=0;j<nbsp;j++) if(group->ldata[j]==i) k++;
      OS_message(0,OS_STD,"[SP_addSpectra] group #%-2d : %d spectr%s\n",i+1,k,(k>1)?"a":"um");
    }

  /************************************************/
  /*Compute scale factors by group and apply them */
  /************************************************/
  
  imax=imin=igroup=0;
  while(imax<nbsp)
    {
      imax=imin;
      while((imax<nbsp)&&(group->ldata[imax]==igroup)) ++imax;
      //      printf("igroup: %2ld imin: %2ld  imax: %2ld\n",igroup,imin,imax);
      factors=addiScaleFactors(spectra+imin,imax-imin,width,NULL);
      for(i=0,si=imin;si<imax;i++,si++)
	for(sj=0;sj<spectra[si]->npix;sj++)
	  {
	    spectra[si]->flux[sj]*=factors->ddata[i];
	    spectra[si]->sigm[sj]*=factors->ddata[i];
	  }
      VC_free(factors);
      imin=imax;
      igroup++;
    }

  /*************/
  /* Find edge */
  /*************/


  /****************************************/
  /* Compute merge factors and apply them */
  /****************************************/
  /* Two steps for the computation.             */
  /* First find scale factor for the edges      */
  /* Then rebin each edges and compute their    */
  /* ratio. Fit the ratio by a line and take    */
  /* max of 1 and the line value for correction */
  
  if(nbgroup>1)
    {    
      OS_message(0,OS_STD,"[SP_addSpectra] Compute factor for the merging\n");
      MMV_malloc(work,nbsp,double);
      imin=igroup=0;
      while(igroup<nbgroup-1)
	{
	  imax=imin;
	  while((imax<nbsp)&&(group->ldata[imax]==igroup)) imax++;
	  //	  printf("igroup: %2ld imin: %2ld  imax: %2ld\n",igroup,imin,imax);
	  ++igroup;
	  if(imax==nbsp)
	    {
	      OS_message(0,OS_ERR,"[SP_addSpectra] BUG at line %d of file %s\n",__LINE__,__FILE__);
	      exit(1);
	    }
	  jmin=imax;
	  jmax=jmin;
	  while((jmax<nbsp)&&(group->ldata[jmax]==igroup)) jmax++;
	  //printf("igroup: %2ld jmin: %2ld  jmax: %2ld\n",igroup,jmin,jmax);

	  w1min=wmin[imin];
	  w1max=wmax[imin];
	  for(si=imin+1;si<imax;si++)
	    {
	      w1min=(w1min<wmin[si])?w1min:wmin[si];
	      w1max=(w1max>wmax[si])?w1max:wmax[si];
	    }

	  w2min=wmin[jmin];
	  w2max=wmax[jmin];
	  for(sj=jmin+1;sj<jmax;sj++)
	    {
	      w2min=(w2min<wmin[sj])?w2min:wmin[sj];
	      w2max=(w2max>wmax[sj])?w2max:wmax[sj];
	    }

	  tmin=(w1min>w2min)?w1min:w2min;
	  tmax=(w1max<w2max)?w1max:w2max;
	  if(tmin>=tmax) /* no overlap*/
	    {
	      w1min=w1max-width/10.;
	      w2max=w2min+width/10.;
	    }
	  else
	    {
	      w1min=w2min=tmin;
	      w1max=w2max=tmax;
	    }

	  for(si=imin;si<imax;si++)
	    {
	      k1min=findNear(spectra[si]->data[spectra[si]->colW],w1min);
	      k1max=findNear(spectra[si]->data[spectra[si]->colW],w1max)+1;
	      for(i=0,j=k1min;j<k1max;j++,i++) 
		tmp[i]=spectra[si]->flux[j];
	      work[si]=median(tmp,k1max-k1min);
	    }

	  for(sj=jmin;sj<jmax;sj++)
	    {
	      k2min=findNear(spectra[sj]->data[spectra[sj]->colW],w2min);
	      k2max=findNear(spectra[sj]->data[spectra[sj]->colW],w2max)+1;
	      for(i=0,j=k2min;j<k2max;j++,i++) 
		tmp[i]=spectra[sj]->flux[j];
	      work[sj]=median(tmp,k2max-k2min);
	    }

	  fact1=median(work+imin,imax-imin);
	  fact2=median(work+jmin,jmax-jmin);
	  if(fact1<=0) fact1=1;
	  if(fact2<=0) fact2=1;


	  /* Factor are computed. Now rebin edges*/
	  /* imin->imax index of spectra in group 1 */
	  /* jmin->jmax index of spectra in group 2 */
	  /* Compute flx1 and flx2 */
	  
	  if(tmin>=tmax)
	    {
	      fact1/=fact2;
	      OS_message(0,OS_STD,"[SP_addSpectra] Factor for merging for group #%d : %f\n",1+igroup,fact1);
	      for(sj=jmin;sj<jmax;sj++)
		for(j=0;j<spectra[sj]->npix;j++)
		  {
		    spectra[sj]->flux[j]*=fact1;
		    spectra[sj]->sigm[j]*=fact1;
		  }
	  
	    }
	  else
	    {
	      /*First crude addition of the edge 2*/
	      k1min=findNear(spectra[imin]->data[spectra[imin]->colW],w1min);
	      k1max=findNear(spectra[imin]->data[spectra[imin]->colW],w1max)+1;
	      //	      printf("%d  %d\n",k1min,k1max);
	      VCV_allocate(flx1,k1max-k1min,VCTDOUBLE);
	      VCV_allocate(wav1,k1max-k1min,VCTDOUBLE);
	      VCV_allocate(flx2,k1max-k1min,VCTDOUBLE);
	      VCV_allocate(pds1,k1max-k1min,VCTDOUBLE);
	      VCV_allocate(pds2,k1max-k1min,VCTDOUBLE);

	      //	      printf("flx1 %x \n",flx1);
	      //	      printf("flx2 %x \n",flx2);
	      //	      printf("wav1 %x \n",wav1);

	      for(i=k1min,j=0;i<k1max;i++,j++)
		{
		  wav1->ddata[j]=spectra[imin]->wave[i];
		  flx1->ddata[j]=0.;
		  pds1->ddata[j]=0.;
		  flx2->ddata[j]=0.;
		  pds2->ddata[j]=0.;
		}

	      for(si=imin;si<imax;si++)
		{
		  k2min=k2max=-1;
		  dw=spectra[si]->wave+idmin[si];
		  df=spectra[si]->flux+idmin[si];
		  ds=spectra[si]->sigm+idmin[si];

		  for(i=k1min,j=0;i<k1max;i++,j++)
		    {
		      if(j==0) tmin=0.5*(3*wav1->ddata[0]-wav1->ddata[1]);
		      else     tmin=0.5*(wav1->ddata[j-1]+wav1->ddata[j]);
		      
		      if(j==(wav1->size-1)) tmax=0.5*(3*wav1->ddata[j]-wav1->ddata[j-1]);
		      else                  tmax=0.5*(wav1->ddata[j+1]+wav1->ddata[j]);
		      
		      rebinLocal(df,ds,dw,idmax[si]-idmin[si],
				 tmin,tmax,&k2min,&k2max,
				 &fsum,&wsum);
		      if(wsum<=0) continue;
		      weight=fsum/wsum;
		      weight*=weight;
		      pds1->ddata[j]+=weight;
		      flx1->ddata[j]+=weight*fsum;
		    }
		}
	      
	      for(sj=jmin;sj<jmax;sj++)
		{
		  k2min=k2max=-1;
		  dw=spectra[sj]->wave+idmin[sj];
		  df=spectra[sj]->flux+idmin[sj];
		  ds=spectra[sj]->sigm+idmin[sj];

		  for(i=k1min,j=0;i<k1max;i++,j++)
		    {
		      if(j==0) tmin=0.5*(3*wav1->ddata[0]-wav1->ddata[1]);
		      else     tmin=0.5*(wav1->ddata[j-1]+wav1->ddata[j]);
		      
		      if(j==(wav1->size-1)) tmax=0.5*(3*wav1->ddata[j]-wav1->ddata[j-1]);
		      else                  tmax=0.5*(wav1->ddata[j+1]+wav1->ddata[j]);
		      
		      rebinLocal(df,ds,dw,idmax[sj]-idmin[sj],
				 tmin,tmax,&k2min,&k2max,
				 &fsum,&wsum);
		      if(wsum<=0) continue;
		      weight=fsum/wsum;
		      weight*=weight;
		      pds2->ddata[j]+=weight;
		      flx2->ddata[j]+=weight*fsum;
		    }
		}

	      ngood=0;
	      for(j=0;j<wav1->size;j++)
		{
		  if((!pds1->ddata[j])||(!pds2->ddata[j])) continue;
		  y1=flx1->ddata[j]/(fact1*pds1->ddata[j]);
		  y2=flx2->ddata[j]/(fact2*pds2->ddata[j]);
		  
		  if((y1<0.2)||(y2<0.2)) continue;
		  flx1->ddata[ngood]=y1/y2;
		  wav1->ddata[ngood]=wav1->ddata[j];
		  ngood++;
		}
	      if(ngood)
		fitLineIter(flx1->ddata,wav1->ddata,ngood,&a1,&b1);
	      else
		{
		  a1=0;
		  b1=1;
		}
	      wpas=0.1/(wav1->ddata[1]-wav1->ddata[0]);
	      fact1/=fact2;
	      /* apply correction */
	      for(sj=jmin;sj<jmax;sj++)
		for(j=0;j<spectra[sj]->npix;j++)
		  {
		    wsum=spectra[sj]->wave[j];
		    weight=(fact1*(a1*wsum+b1)*expFilter(wsum,w1max,-wpas))+
		      (fact1*(a1*w1max+b1)*expFilter(wsum,w1max,wpas));
		      
		    spectra[sj]->flux[j]*=weight;
		    spectra[sj]->sigm[j]*=weight;
		  }
	      VC_free(pds2);
	      VC_free(pds1);
	      VC_free(flx1);
	      VC_free(flx2);
	      VC_free(wav1);

	    }
	  imin=jmin;
	}
      MM_free(work);
    }


  /**************************/
  /* Compute the wave range */
  /**************************/

  OS_message(0,OS_STD,"[SP_addSpectra] Compute wavelength range\n");

  /*****************************************/
  /* merge all the wave limits and compute */
  /* average pixel size for each group     */
  /*****************************************/
  VCV_allocate(wlim,2*nbgroup,VCTDOUBLE);
  imin=igroup=0;
  while(imin<nbsp)
    {
      imax=imin;
      tmin=wmin[imin];
      tmax=wmax[imin];
      while((imax<nbsp)&&(group->ldata[imax]==igroup))
	{
	  tmin=(tmin<wmin[imax])?tmin:wmin[imax];
	  tmax=(tmax>wmax[imax])?tmax:wmax[imax];
	  imax++;
	}
      wlim->ddata[  2*igroup]=tmin;
      wlim->ddata[1+2*igroup]=tmax;
      imin=imax;
      igroup++;
    }
  VC_hpsort(wlim);
  /* Compute pixel size */
  
  MMV_malloc(npix,(2*nbgroup-1),long);
  
  wpas=0;
  for(i=0;i<2*nbgroup-1;i++)
    {
      k=0;
      npix[i]=0;
      tcen=0.5*(wlim->ddata[i]+wlim->ddata[i+1]);
      old_wpas=wpas;
      wpas=0;
      for(j=0;j<nbsp;j++)
	{
	  if((tcen>=wmax[j])||(tcen<=wmin[j])) continue;
	  k++;
	  wpas+=(wmax[j]-wmin[j])/(float)(idmax[j]-idmin[j]);
	}
      //      printf("chunck %2ld  k: %2ld   swpas: %f   rap:%f\n",i,k,wpas,(k>0)?(wpas/(double)k):0);
      if(k)
	wpas/=(double)k;
      else
	if(!i)
	  {
	    OS_message(0,OS_ERR,"Bug at %s line %d !!!\n",__FILE__,__LINE__);
	    exit(1);
	  }
	else
	  wpas=old_wpas;
      wpas*=rebin_pixfactor;
      npix[i]=(long)((wlim->ddata[i+1]-wlim->ddata[i])/wpas+0.5);
      if(!npix[i])
	{
	  OS_message(0,OS_ERR,"Bug at %s line %d !!!\n",__FILE__,__LINE__);
	  exit(1);
	}
    }
  npix_tot=0;
  for(i=0;i<2*nbgroup-1;i++) npix_tot+=npix[i]-1;
  npix_tot++;

  /* Fill the final wave array */
  VCV_allocate(wave,npix_tot,VCTDOUBLE);
  k=0;
  for(i=0;i<2*nbgroup-1;i++)
    {
      wpas=(wlim->ddata[i+1]-wlim->ddata[i])/(double)(npix[i]);
      for(j=0;j<npix[i]-1;j++)
	wave->ddata[k++]=wlim->ddata[i]+((float)j+rebin_pixshift)*wpas;
    }
  wave->ddata[k++]=wlim->ddata[2*nbgroup-1]+rebin_pixshift*wpas;

  /* display some information */
  OS_message(0,OS_STD,"\n[SP_mergeSpectra] New wavelength range (%d chunks) \n",2*nbgroup-1);
  for(i=0;i<2*nbgroup-1;i++)
    OS_message(0,OS_STD,"[SP_mergeSpectra]  chunks #%-2d : %10.2f %10.2f %9.4f\n",
	       i+1,wlim->ddata[i],wlim->ddata[i+1],
	       (wlim->ddata[i+1]-wlim->ddata[i])/(double)npix[i]);


  //  for(i=0;i<nbsp;i++) printf("i:%2d   wmin:%f   wmax:%f\n",(int)i,wmin[i],wmax[i]);
  MM_free( npix);
  VC_free( wlim);
  MM_free(wmin);
  MM_free(wmax);


  /**********************************/
  /* if spline is set rebin spectra */
  /* with cubic spline              */
  /**********************************/
  if(spline)
    {
      OS_message(0,OS_STD,"[SP_mergeSpectra]  Doing spline rebinning : \n");
      MMV_malloc(spline_spectra,nbsp,SP_spectrum);
      if(!spline_spectra)
	{
	  OS_message(0,OS_ERR,"not enough memory !!\n");
	  MM_free( kmax);
	  MM_free( kmin);
	  MM_free(idmax);
	  MM_free(idmin);
	  MM_free(  tmp);
	  return NULL;
	}
      for(i=0;i<nbsp;i++)
	{
	  tmin=spectra[i]->wave[idmin[i]];
	  tmax=spectra[i]->wave[idmax[i]-1];
	  if(wave->ddata[0]>tmin)
	    imin=0;
	  else
	    while((imin<wave->size)&&(wave->ddata[imin]<tmin)) ++imin;

	  imax=imin;
	  while((imax<wave->size)&&(wave->ddata[imax]<tmax)) ++imax;

	  spline_spectra[i]=SP_allocate(imax-imin,3);
	  if(!spline_spectra[i])
	    {
	      OS_message(0,OS_ERR,"not enough memory !!\n");
	      MM_free(spline_spectra);
	      MM_free( kmax);
	      MM_free( kmin);
	      MM_free(idmax);
	      MM_free(idmin);
	      MM_free(  tmp);
	      return NULL;
	    }
	  /*Copy Wavelength array*/
	  for(k=0,j=imin;j<imax;j++,k++) 
	    spline_spectra[i]->wave[k]=wave->ddata[j];
	  /*********************************/
	  /* do the cubic spline estimation*/
	  /* for the flux                  */
	  NRspline(spectra[i]->wave,spectra[i]->flux,spectra[i]->npix,1e30,1e30,tmp);
	  NRasplint(spectra[i]->wave,spectra[i]->flux,tmp,spectra[i]->npix,
		    spline_spectra[i]->wave,spline_spectra[i]->flux,spline_spectra[i]->npix);
	  /* for the error                 */
	  /* (strange, maybe to fix latter)*/
	  NRspline(spectra[i]->wave,spectra[i]->sigm,spectra[i]->npix,1e30,1e30,tmp);
	  NRasplint(spectra[i]->wave,spectra[i]->sigm,tmp,spectra[i]->npix,
		    spline_spectra[i]->wave,spline_spectra[i]->sigm,spline_spectra[i]->npix);
	  /* Rebin Done                    */
	  /* Put index instead of wavelenght */
	  for(k=0,j=imin;j<imax;j++,k++) 
	    spline_spectra[i]->wave[k]=j;
	  
	  /*********************************/
	  /*Old spectra[i] useless, set it to spline_spectra[i] */
	  SP_free(spectra[i]);
	  spectra[i]=spline_spectra[i];
	}
    }

  /***************************************/
  /* Allocate memory for output spectrum */
  /***************************************/
  nbcol=3;
  if(keepnbsp) ++nbcol;
  if(error_type==_BOTH_ERRORS_) ++nbcol;

  rslt=_SP_allocate(wave->size,nbcol,file,line);

  if(!rslt) 
    {
      OS_message(0,OS_ERR,"not enough memory !!\n");
      MM_free( kmax);
      MM_free( kmin);
      MM_free(idmax);
      MM_free(idmin);
      MM_free(  tmp);
      return NULL;
    }
  
  VC_free(rslt->data[rslt->colW]);
  rslt->data[rslt->colW]=wave;

  /*******************/
  /* Do the addition */
  /*******************/
  VCV_allocate(eflux,nbsp,VCTDOUBLE);
  VCV_allocate(esigm,nbsp,VCTDOUBLE);
  VCV_allocate(emask,nbsp,VCTDOUBLE);

  OS_message(0,OS_STD,"\n[SP_addSpectra] Do the addition :   0%%");
  old_pourc=pourc=0;


  for(j=0;j<rslt->npix;j++)
    {
      pourc=(int)(j*100/(rslt->npix-1));
      if(pourc>old_pourc)
	{
	  old_pourc=pourc;
	  OS_message(0,OS_STD,"\b\b\b\b%3ld%%",pourc);
	  OS_flush();
	}
      
      if(spline)
	{
	  ngood=0;
	  for(i=0;i<nbsp;i++)
	    {
	      eflux->ddata[ngood]=esigm->ddata[ngood]=0.;
	      dw=spectra[i]->wave;
	      df=spectra[i]->flux;
	      ds=spectra[i]->sigm;
	      if(getSplineInter(df,ds,dw,spectra[i]->npix,j,kmin+i,
				eflux->ddata+ngood,esigm->ddata+ngood)) continue;
	      emask->ddata[ngood]=1.;
	      ++ngood;
	    }
	}
      else
	{
	  if(j==0)  tmin=0.5*(3*wave->ddata[0]-wave->ddata[1]);
	  else      tmin=0.5*(wave->ddata[j-1]+wave->ddata[j]);
	  
	  if(j==(rslt->npix-1)) tmax=0.5*(3*wave->ddata[rslt->npix-1]-wave->ddata[rslt->npix-2]);
	  else                  tmax=0.5*(  wave->ddata[j+1]+         wave->ddata[j]);
	  
	  ngood=0;
	  for(i=0;i<nbsp;i++)
	    {
	      eflux->ddata[ngood]=esigm->ddata[ngood]=0.;
	      dw=spectra[i]->wave+idmin[i];
	      df=spectra[i]->flux+idmin[i];
	      ds=spectra[i]->sigm+idmin[i];
	      
	      if(rebinLocal(df,ds,dw,idmax[i]-idmin[i],
			    tmin,tmax,kmin+i,kmax+i,
			    eflux->ddata+ngood,esigm->ddata+ngood)) continue;
	      emask->ddata[ngood]=1.;
	      ++ngood;
	    }
	}
/*       if((tmin>5156)&&(tmax<5158)) */
/* 	printf("tmin:%f  tmax:%f   ngood: %d\n",tmin,tmax,ngood); */

      mean_srianand=sigm_srianand=fsum=wsum=0.0;
      nbg_srianand=nbg_bastien=0;
      if(ngood)
	{
	  eflux->size=ngood;
	  esigm->size=ngood;
	  emask->size=ngood;
	  
	  if(tclip==_DATACLIPPING_)
	    dataClipping(eflux,esigm,&emask,clip);
	  else if (tclip!=_NOCLIPPING_)
	    sigmaClipping(eflux,&emask,clip,0);

	  for(i=0;i<ngood;i++)
	    {
	      if(emask->ddata[i])
		{
		  weight=esigm->ddata[i];
		  mean_srianand+=eflux->ddata[i];
		  sigm_srianand+=eflux->ddata[i]*eflux->ddata[i];
		  ++nbg_srianand;
		  if(weight<=0) continue;
		  ++nbg_bastien;
		  weight=1./(weight*weight);
		  wsum+=weight;
		  fsum+=eflux->ddata[i]*weight;
		}
	    }
	}
      if(nbg_srianand>0)
	{
/* 	  mean_srianand*=mean_srianand; */
/* 	  mean_srianand/=(double)(nbg_srianand); /\* n.<X>^2 *\/ */
/* 	  sigm_srianand-=mean_srianand;  // n<X^2>-n<X>^2 */
/* 	  sigm_srianand/=(double)(nbg_srianand);  //(n/n-1)(<X^2>-<X>^2) */
	  mean_srianand/=(double)(nbg_srianand);
	  mean_srianand*=mean_srianand;
	  
	  sigm_srianand/=(double)(nbg_srianand);
	  sigm_srianand=(sigm_srianand-mean_srianand)/(double)(nbg_srianand);
	}
      else
	sigm_srianand=0.;
      if(wsum>0)
	{
	  wsum=1./wsum;
	  fsum*=wsum;
	  switch (error_type)
	    {
	    case _BOTH_ERRORS_ :
	    case _BASTIEN_     : 
	      {
		wsum=sqrt(wsum);
		sigm_srianand=sqrt(sigm_srianand);
	      }
	      break;
	    case _SRIANAND_                 : wsum=sqrt(sigm_srianand);
	      break;
	    case _BASTIEN_AND_SRIANAND_     : wsum=sqrt(wsum+sigm_srianand);
	      break;
	    case _MAX_BASTIEN_AND_SRIANAND_ : wsum=sqrt((sigm_srianand>wsum)?sigm_srianand:wsum);
	      break;
	    case _MAX_BASTIEN_AND_MAX_SRIANAND_ : 
	      if(nbg_srianand>1)
		fct=(1+sqrt(2./(nbg_srianand-1)));
	      else
		fct=1.;
	      wsum=sqrt((sigm_srianand>wsum)?sigm_srianand*fct:wsum);
	      break;
	    default : wsum=sqrt((sigm_srianand>wsum)?sigm_srianand:wsum);
	    }
	}
      else
	fsum=wsum=0;

      if((!IS_NAN(fsum))&&(!IS_NAN(wsum)))
	{
	  rslt->flux[j]=fsum;
	  rslt->sigm[j]=wsum;
	}
      if((error_type==_BOTH_ERRORS_)&&(!IS_NAN(sigm_srianand)))
	rslt->cont[j]=sigm_srianand;
      if(keepnbsp)
	rslt->data[rslt->ndata-1]->ddata[j]=nbg_bastien;

      eflux->size=nbsp;
      esigm->size=nbsp;
      emask->size=nbsp;
    }
  OS_message(0,OS_STD,"\b\b\b\b%3ld%%\n",100);
  VC_free(emask);
  VC_free(esigm);
  VC_free(eflux);
  MM_free( kmin);
  MM_free( kmax);
  MM_free(idmin);
  MM_free(idmax);
  MM_free(tmp);
  return rslt;
}


SP_spectrum _SP_allocate(size_t npix,size_t nd,char *file,size_t line)
{
  SP_spectrum rslt;
  size_t i,j;
  
  nd=(nd<2)?2:nd;
  


  rslt=(SP_spectrum)_MM_malloc(sizeof(SP_spectrum_ref),file,line);
  if(!rslt) return NULL;

  rslt->wave=NULL;
  rslt->flux=NULL;
  rslt->sigm=NULL;
  rslt->cont=NULL;

  rslt->colW=-1;
  rslt->colF=-1;
  rslt->colS=-1;
  rslt->colC=-1;

  rslt->data=(VC_vector*)MM_malloc(nd*sizeof(VC_vector));
  if(!rslt->data)
    {
      MM_free(rslt);
      return NULL;
    }
  
  for(i=0;i<nd;i++)
    {
      rslt->data[i]=VC_allocate(npix,VCTDOUBLE);
      if(!rslt->data[i])
	{
	  for(j=0;j<i;j++) VC_free(rslt->data[j]);
	  MM_free(rslt->data);
	  MM_free(rslt);
	  return NULL;
	}
    }

  rslt->npix=npix;
  rslt->ndata=nd;

  rslt->colW=(nd>0)?0:-1;
  rslt->colF=(nd>1)?1:-1;
  rslt->colS=(nd>2)?2:-1;
  rslt->colC=(nd>3)?3:-1;
  
  rslt->wave=(rslt->colW>=0)?rslt->data[rslt->colW]->ddata:NULL;
  rslt->flux=(rslt->colF>=0)?rslt->data[rslt->colF]->ddata:NULL;
  rslt->sigm=(rslt->colS>=0)?rslt->data[rslt->colS]->ddata:NULL;
  rslt->cont=(rslt->colC>=0)?rslt->data[rslt->colC]->ddata:NULL;
  rslt->name=NULL;

  return rslt;
}


void _SP_free(SP_spectrum sp,char* file,size_t line)
{
  size_t i;
  for(i=0;i<sp->ndata;i++) VC_free(sp->data[i]);
  MM_free(sp->data);
  if(sp->name) MM_free(sp->name);
  MM_free(sp);
}

SP_spectrum _SP_readTfits(char* filename,VC_vector which,char* file,size_t line)
{
  size_t i;
  TB_table tmp;
  SP_spectrum rslt;

  tmp=readTfits(filename,which,NULL);

  if(!tmp) return NULL;
  
  rslt=(SP_spectrum)_MM_malloc(sizeof(SP_spectrum_ref),file,line);
  if(!rslt)
    {
      TB_free(tmp);
      return NULL;
    }
  rslt->name=tmp->name;
  rslt->npix=tmp->nbRow;
  rslt->ndata=tmp->nbCol;
  
  rslt->data=tmp->data;
  for(i=0;i<tmp->nbCol;i++)
    if(tmp->labels[i]) MM_free(tmp->labels[i]);
  MM_free(tmp->labels);
  MM_free(tmp);

  rslt->colW=(rslt->ndata>0)?0:-1;
  rslt->colF=(rslt->ndata>1)?1:-1;
  rslt->colS=(rslt->ndata>2)?2:-1;
  rslt->colC=(rslt->ndata>3)?3:-1;

  rslt->wave=(rslt->colW>=0)?rslt->data[rslt->colW]->ddata:NULL;
  rslt->flux=(rslt->colF>=0)?rslt->data[rslt->colF]->ddata:NULL;
  rslt->sigm=(rslt->colS>=0)?rslt->data[rslt->colS]->ddata:NULL;
  rslt->cont=(rslt->colC>=0)?rslt->data[rslt->colC]->ddata:NULL;

  return rslt;

}

SP_spectrum _SP_readAscii(char* filename,VC_vector which,char* file,size_t line)
{
  FILE *fstrm;
  size_t i,j;
  int ok;
  long nbl,nbw;
  long nbc;
  SP_spectrum rslt;

  if(wordCount(filename,&nbl,&nbw,NULL)) return NULL;
  if((fstrm=fopen(filename,"r"))==NULL) return NULL;

  if(nbw%nbl) return NULL;
  nbc=nbw/nbl;

  rslt=_SP_allocate(nbl,nbc,file,line);
  if(!rslt) return NULL;

  MMV_malloc(rslt->name,strlen(filename)+1,char);
  strcpy(rslt->name,filename);
  ok=1;

  for(i=0;(i<nbl)&&ok;i++)
    for(j=0;(j<nbc)&&ok;j++)
      {
	if(gotoNextWord(fstrm)) 
	  {
	    ok=0;
	    break;
	  }
	if(readNumberD(fstrm,rslt->data[j]->ddata+i))
	  {
	    fprintf(stderr,"Error when reading file '%s' (line: %ld column: %ld)\n",filename,i,j);
	    exit(1);
	  }
      }
  
  if(ok==0) fprintf(stderr,"Reach end of file !!\n");
  return rslt;
}



int _SP_writeTfits(SP_spectrum spin,char* filename,VC_vector which,char* file,size_t line)
{
  size_t i;
  TB_table tmp;
  int rtn;

  tmp=(TB_table)MM_malloc(sizeof(TB_table_ref));
  if(!tmp) return 1;
  tmp->name=spin->name;
  tmp->nbCol=spin->ndata;
  tmp->nbRow=spin->npix;
  tmp->labels=(char**)MM_malloc(tmp->nbCol*sizeof(char*));
  if(!tmp->labels)
    {
      MM_free(tmp);
      return 1;
    }
  if(spin->ndata>0) tmp->labels[0]="Wave";
  if(spin->ndata>1) tmp->labels[1]="Flux";
  if(spin->ndata>2) tmp->labels[2]="Sigm";
  if(spin->ndata>3) tmp->labels[3]="Cont";

  for(i=4;i<spin->ndata;i++)
    tmp->labels[i]="Unknown";
  
  tmp->data=spin->data;

  rtn=writeTfits(tmp,filename,which);
  MM_free(tmp->labels);
  MM_free(tmp);
  return rtn;
}

int _SP_writeAscii(SP_spectrum spin,char* filename,VC_vector which,int prec,char* file,size_t line)
{
  FILE* stream;
  size_t i,j,jj,nbc,nbr;

  if(!(stream=fopen(filename,"w")))
    {
      OS_message(0,OS_ERR,"[SP_writeAscii] Cannot open file '%s' !!\n",filename);
      exit(1);
    }
  nbc=(which)?which->size:spin->ndata;
  nbr=spin->npix;
  for(i=0;i<nbr;i++)
    for(j=0;j<nbc;j++)
      {
	jj=(which)?which->ldata[j]:j;
	fprintf(stream,"%.*f%1s",prec,spin->data[jj]->ddata[i],((j==(nbc-1))?"\n":" "));
      }
  fclose(stream);
  return 0;
}

void SP_sortSpectrum(SP_spectrum* spectra,size_t nbsp)
{
  size_t i,j;
  double wmini,wminj;
  long  idmin,idmax;
  SP_spectrum aux;

  for(i=0;i<nbsp-1;i++)
    {
      VC_getIdxEdge(spectra[i]->data[spectra[i]->colF],&idmin,&idmax);
      wmini=spectra[i]->data[spectra[i]->colW]->ddata[idmin];
      for(j=i+1;j<nbsp;j++)
	{
	  VC_getIdxEdge(spectra[j]->data[spectra[j]->colF],&idmin,&idmax);
	  wminj=spectra[j]->data[spectra[j]->colW]->ddata[idmin];
	  if(wminj>=wmini) continue;

	  aux=spectra[i];
	  spectra[i]=spectra[j];
	  spectra[j]=aux;
	  wmini=wminj;
	}
    }
  return ;
}

VC_vector _SP_groupSpectrum(SP_spectrum* spectra,size_t nbsp,double pourc,char* file,size_t line)
{
  VC_vector group;
  long i,j,idx;
  long idmin,idmax;
  double wmin,wmax,wtmin,wtmax;
  double prc;

  VCV_allocate(group,nbsp,VCTLONG);
  for(i=0;i<nbsp;i++) group->ldata[i]=-1;
  idx=0;
  for(i=0;i<nbsp;i++)
    {
      if(group->ldata[i]>=0) continue;
      group->ldata[i]=idx;
      VC_getIdxEdge(spectra[i]->data[spectra[i]->colF],&idmin,&idmax);
      wmin=spectra[i]->data[spectra[i]->colW]->ddata[idmin];
      wmax=spectra[i]->data[spectra[i]->colW]->ddata[idmax-1];
      for(j=0;j<nbsp;j++)
	{
	  if(group->ldata[j]>=0) continue;
	  VC_getIdxEdge(spectra[j]->data[spectra[j]->colF],&idmin,&idmax);	  
	  wtmin=spectra[j]->data[spectra[j]->colW]->ddata[idmin];
	  wtmax=spectra[j]->data[spectra[j]->colW]->ddata[idmax-1];
	  prc=pourcOver(wmin,wmax,wtmin,wtmax);
	  if(prc>pourc)
	    group->ldata[j]=idx;
	}
      idx++;
    }
  return group;
}



double pourcOver(double a1,double b1,double a2,double b2)
{
  double a,b;
  
  a=(a1>a2)?a1:a2;
  b=(b1<b2)?b1:b2;
  
  if(a>=b) return 0;
  
  return (b-a)/(b2-a2);
}

double expFilter(double x,double x0,double a)
{
  double xx;

  xx=-a*(x-x0);
  if(xx>100) return 0.;
  return 1./(1+exp(xx));
}

int getSplineInter(double* flux,double* sigm,double* index,long npix,
		   long idx,long* kpos,double* fout,double* sout)
{
  long ipos;

  if((index[     0]-idx)>0.5) return 1;
  if((idx-index[npix-1])>0.5) return 1;
  
  if(*kpos<0) ipos=0;
  else ipos=*kpos;
  while(fabs(index[ipos]-idx)>0.5) ipos++;
  *kpos=ipos;
  *fout=flux[ipos];
  *sout=sigm[ipos];
  return 0;
}

