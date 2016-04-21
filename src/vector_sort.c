/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdlib.h>
#include "memory_manager.h"
#include "vector_sort.h"




int VC_hpsort(VC_vector vcin)
{
  switch (vcin->kind)
    {
    case VCTNONE   :
    case VCTCHAR   : return VC_chpsort(vcin);
    case VCTSHORT  : return VC_shpsort(vcin);
    case VCTINT    : return VC_ihpsort(vcin);
    case VCTLONG   : return VC_lhpsort(vcin);
    case VCTFLOAT  : return VC_fhpsort(vcin);
    case VCTDOUBLE : return VC_dhpsort(vcin);
    }
  return 1;
}

VC_vector _VC_hpsort_idx(VC_vector vcin,char* file,size_t line)
{
  switch (vcin->kind)
    {
    case VCTNONE   :
    case VCTCHAR   : return _VC_chpsort_idx(vcin,file,line);
    case VCTSHORT  : return _VC_shpsort_idx(vcin,file,line);
    case VCTINT    : return _VC_ihpsort_idx(vcin,file,line);
    case VCTLONG   : return _VC_lhpsort_idx(vcin,file,line);
    case VCTFLOAT  : return _VC_fhpsort_idx(vcin,file,line);
    case VCTDOUBLE : return _VC_dhpsort_idx(vcin,file,line);
    }
  return NULL;
}





#define _DEFINE_FUNC_TEMPLATE_(_Fc,_Data,_Tp)			\
int _Fc(VC_vector vcin)          				\
{								\
  size_t i,ir,j,l,n;						\
  _Tp rra;							\
  _Tp *ra;							\
								\
  n=vcin->size;							\
								\
  if(n<2) return 0;						\
  ra=vcin->_Data;				                \
  ra--; /*To have the Num-Rec convention*/                      \
  l=(n >> 1)+1;							\
  ir=n;								\
  while(1)							\
    {								\
      if (l > 1)						\
	rra=ra[--l];						\
      else							\
	{							\
	  rra=ra[ir];						\
	  ra[ir]=ra[1];						\
	  if (--ir == 1) {					\
	    ra[1]=rra;						\
	    break;						\
	  }							\
	}							\
      i=l;							\
      j=l+l;							\
      while (j <= ir)						\
	{							\
	  if (j < ir && ra[j] < ra[j+1]) j++;			\
	  if (rra < ra[j])					\
	    {							\
	      ra[i]=ra[j];					\
	      i=j;						\
	      j <<= 1;						\
	    }							\
	  else							\
	    j=ir+1;						\
	}							\
      ra[i]=rra;						\
    }								\
  return 0;							\
}

_DEFINE_FUNC_TEMPLATE_(VC_chpsort,cdata,char);
_DEFINE_FUNC_TEMPLATE_(VC_shpsort,sdata,short);
_DEFINE_FUNC_TEMPLATE_(VC_ihpsort,idata,int);
_DEFINE_FUNC_TEMPLATE_(VC_lhpsort,ldata,long);
_DEFINE_FUNC_TEMPLATE_(VC_fhpsort,fdata,float);
_DEFINE_FUNC_TEMPLATE_(VC_dhpsort,ddata,double);

#undef _DEFINE_FUNC_TEMPLATE_

#define _DEFINE_FUNC_TEMPLATE_(_Fc,_Data)				\
VC_vector _Fc(VC_vector vcin,char* file,size_t line)			\
{									\
  size_t i,ir,j,l,n;							\
  VC_vector rslt;							\
  long rra;								\
  long *ra;								\
									\
  n=vcin->size;								\
  rslt=_VC_allocate(n,VCTLONG,file,line);						\
									\
  if(!rslt) return NULL;						\
  for(i=0;i<=n;i++)							\
    rslt->ldata[i]=i;							\
									\
  if(n<2) return rslt;							\
									\
									\
  ra=rslt->ldata;							\
  ra--; /*To have the Num-Rec convention*/				\
  l=(n >> 1)+1;								\
  ir=n;									\
  while(1)								\
    {									\
      if (l > 1)							\
	rra=ra[--l];							\
      else								\
	{								\
	  rra=ra[ir];							\
	  ra[ir]=ra[1];							\
	  if (--ir == 1) {						\
	    ra[1]=rra;							\
	    break;							\
	  }								\
	}								\
      i=l;								\
      j=l+l;								\
      while (j <= ir)							\
	{								\
	  if (j < ir && vcin->_Data[ra[j]] < vcin->_Data[ra[j+1]]) j++;	\
	  if (vcin->_Data[rra] < vcin->_Data[ra[j]])			\
	    {								\
	      ra[i]=ra[j];						\
	      i=j;							\
	      j <<= 1;							\
	    }								\
	  else								\
	    j=ir+1;							\
	}								\
      ra[i]=rra;							\
    }									\
  return rslt;								\
}

_DEFINE_FUNC_TEMPLATE_(_VC_chpsort_idx,cdata);
_DEFINE_FUNC_TEMPLATE_(_VC_shpsort_idx,sdata);
_DEFINE_FUNC_TEMPLATE_(_VC_ihpsort_idx,idata);
_DEFINE_FUNC_TEMPLATE_(_VC_lhpsort_idx,ldata);
_DEFINE_FUNC_TEMPLATE_(_VC_fhpsort_idx,fdata);
_DEFINE_FUNC_TEMPLATE_(_VC_dhpsort_idx,ddata);
