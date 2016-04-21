/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdlib.h>
#include <string.h>
#include "outstream.h"
#include "memory_manager.h"
#include "vector.h"

#define  Balise(a) printf("## [%s at %d] %s\n",__FILE__,__LINE__,a)

size_t VC_sizeof(VC_TYPE kind)
{
  switch (kind)
    {
    case VCTNONE   : return sizeof(char);
    case VCTCHAR   : return sizeof(char);
    case VCTSHORT  : return sizeof(short);
    case VCTINT    : return sizeof(int);
    case VCTLONG   : return sizeof(long);
    case VCTFLOAT  : return sizeof(float);
    case VCTDOUBLE : return sizeof(double);
  }
  OS_message(1,OS_ERR,"[VC_sizeof] Unknown type !!\n");
  return sizeof(char);
}

void VC_update(VC_vector vector)
{
  vector->cdata=(  char*)vector->ndata;
  vector->cdata=(  char*)vector->ndata;
  vector->sdata=( short*)vector->ndata;
  vector->idata=(   int*)vector->ndata;
  vector->ldata=(  long*)vector->ndata;
  vector->fdata=( float*)vector->ndata;
  vector->ddata=(double*)vector->ndata;
  return;
}


VC_vector _VC_allocate(size_t size,VC_TYPE kind,char* file ,size_t line)
{
  VC_vector rslt;

  rslt=(VC_vector)_MM_malloc(sizeof(VC_vector_ref),file,line);
  if(!rslt) return NULL;

  rslt->kind=kind;
  rslt->size=size;

  rslt->ndata=MM_malloc(size*VC_sizeof(kind));
  if(!rslt->ndata) {_MM_free(rslt,file,line);return NULL;}

  VC_update(rslt);

  return rslt;
}

VC_vector _VCV_allocate(size_t size,VC_TYPE kind,char* var,char* file ,size_t line)
{
  VC_vector rslt;

  rslt=(VC_vector)_MMV_malloc(sizeof(VC_vector_ref),var,file,line);
  if(!rslt) return NULL;

  rslt->kind=kind;
  rslt->size=size;

  rslt->ndata=_MMV_malloc(size*VC_sizeof(kind),"rslt->ndata",__FILE__,__LINE__);
  if(!rslt->ndata) {MM_free(rslt);return NULL;}

  VC_update(rslt);
  return rslt;
}

size_t _VC_resize(VC_vector vector,size_t nsize,char* file,size_t line)
{
  void* tmp;
  if(vector==NULL) return 0;

  tmp=MM_realloc((void *)vector->ndata,nsize*VC_sizeof(vector->kind));
  if(!tmp) return vector->size;

  vector->size=nsize;
  vector->ndata=tmp;

  VC_update(vector);

  return vector->size;
}

void _VC_free(VC_vector vector,char* file,size_t line)
{
  MM_free(vector->ndata);
  _MM_free(vector,file,line);
}

int VC_get(VC_vector vct,void* values,VC_TYPE kind,long start,long end,long stride)
{
  switch (kind)
    {
    case VCTNONE   :
    case VCTCHAR   : return VC_getc(vct,values,start,end,stride);
    case VCTSHORT  : return VC_gets(vct,values,start,end,stride);
    case VCTINT    : return VC_geti(vct,values,start,end,stride);
    case VCTLONG   : return VC_getl(vct,values,start,end,stride);
    case VCTFLOAT  : return VC_getf(vct,values,start,end,stride);
    case VCTDOUBLE : return VC_getd(vct,values,start,end,stride);
    }
  OS_message(1,OS_ERR,"[VC_get] Unknown type !!\n");
  return 1;
}

int VC_set(VC_vector vct,void* values,VC_TYPE kind,long start,long end,long stride)
{
  switch (kind)
    {
    case VCTNONE   :
    case VCTCHAR   : return VC_setc(vct,values,start,end,stride);
    case VCTSHORT  : return VC_sets(vct,values,start,end,stride);
    case VCTINT    : return VC_seti(vct,values,start,end,stride);
    case VCTLONG   : return VC_setl(vct,values,start,end,stride);
    case VCTFLOAT  : return VC_setf(vct,values,start,end,stride);
    case VCTDOUBLE : return VC_setd(vct,values,start,end,stride);
    }
  OS_message(1,OS_ERR,"[VC_get] Unknown type !!\n");
  return 1;
}


#define xstr(x) str(x)
#define str(x) #x
#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c
#define GLUE4(a,b,c,d) a##b##c##d

#define _DEFINE_TEMPLATE_FUNC_FCTX(_Fc,_Car,_Tp)						\
int GLUE3(VC_,_Fc,_Car)(VC_vector vct,_Tp* values,long start,long end,long stride)	\
{												\
  switch (vct->kind)										\
    {												\
    case VCTNONE   :										\
    case VCTCHAR   : return GLUE4(VC_,c,_Fc,_Car)(vct,values,start,end,stride);			\
    case VCTSHORT  : return GLUE4(VC_,s,_Fc,_Car)(vct,values,start,end,stride);			\
    case VCTINT    : return GLUE4(VC_,i,_Fc,_Car)(vct,values,start,end,stride);			\
    case VCTLONG   : return GLUE4(VC_,l,_Fc,_Car)(vct,values,start,end,stride);			\
    case VCTFLOAT  : return GLUE4(VC_,f,_Fc,_Car)(vct,values,start,end,stride);			\
    case VCTDOUBLE : return GLUE4(VC_,d,_Fc,_Car)(vct,values,start,end,stride);			\
    }												\
  OS_message(1,OS_ERR,"[" xstr(GLUE3(VC_,_Fc,_Car)) "] Unknown type !!\n");			\
  return 1;											\
}

#define _DEFINE_TEMPLATE_FUNC_YGETX(_Car1,_Data,_Car2,_Tp)					\
int GLUE4(VC_,_Car1,get,_Car2)(VC_vector vct,_Tp* values,long start,long end,long stride)	\
{												\
  size_t i,j;											\
  if(stride>0)											\
    {												\
      for(i=start,j=0;i<end;i+=stride,j++)							\
	values[j]=(_Tp)vct->_Data[i];								\
      return 0;											\
    }												\
  if(stride<0)											\
    {												\
      for(i=start-1,j=0;i>=end;i+=stride,j++)							\
	values[j]=(_Tp)vct->_Data[i];								\
      return 0;											\
    }												\
  return 1;											\
}

#define _DEFINE_TEMPLATE_FUNC_YSETX(_Car1,_Data,_Car2,_Tp)					\
int GLUE4(VC_,_Car1,set,_Car2)(VC_vector vct,_Tp* values,long start,long end,long stride)	\
{												\
  long i,j;											\
  if(stride>0)											\
    {												\
      for(i=start,j=0;i<end;i+=stride,j++)							\
	vct->_Data[i]=values[j];								\
      return 0;											\
    }												\
  if(stride<0)											\
    {												\
      for(i=start-1,j=0;i>=end;i+=stride,j++)							\
	vct->_Data[i]=values[j];								\
      return 0;											\
    }												\
  return 1;											\
}

#define _DEFINE_TEMPLATE_FUNC_YFCTX(_Car1,_Data1,_Car2,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YGETX(_Car1,_Data1,_Car2,_Tp)		\
_DEFINE_TEMPLATE_FUNC_YSETX(_Car1,_Data1,_Car2,_Tp)

#define _DEFINE_TEMPLATE_FUNC_(_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_FCTX(get,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_FCTX(set,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YFCTX(c,cdata,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YFCTX(s,sdata,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YFCTX(i,idata,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YFCTX(l,ldata,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YFCTX(f,fdata,_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_YFCTX(d,ddata,_Car,_Tp)

_DEFINE_TEMPLATE_FUNC_(c,  char);
_DEFINE_TEMPLATE_FUNC_(s, short);
_DEFINE_TEMPLATE_FUNC_(i,   int);
_DEFINE_TEMPLATE_FUNC_(l,  long);
_DEFINE_TEMPLATE_FUNC_(f, float);
_DEFINE_TEMPLATE_FUNC_(d,double);


#undef _DEFINE_TEMPLATE_FUNC_
#undef _DEFINE_TEMPLATE_FUNC_YFCTX
#undef _DEFINE_TEMPLATE_FUNC_YSETX
#undef _DEFINE_TEMPLATE_FUNC_YGETX
#undef _DEFINE_TEMPLATE_FUNC_FCTX

#undef GLUE2
#undef GLUE3
#undef GLUE4

VC_vector _VC_init(VC_TYPE kind,void* values,size_t size,char* file,size_t line)
{
  VC_vector rslt;

  rslt=_VC_allocate(size,kind,file,line);
  if(!rslt) return NULL;
  VC_set(rslt,values,kind,1,size,1);
  return rslt;
}


VC_vector  _VC_copy(VC_vector vin,VC_vector* vout,char* file,size_t line)
{
  void* tmp;
  VC_vector rslt;
  if(vout==NULL||*vout==NULL)
    {
      rslt=_VC_allocate(vin->size,vin->kind,file,line);
      if(vout) *vout=rslt;
      if(!rslt) return NULL;
    }
  else
    {
      tmp=MM_realloc((*vout)->ndata,vin->size*VC_sizeof(vin->kind));
      if(!tmp) return NULL;
      (*vout)->kind=vin->kind;
      (*vout)->size=vin->size;
      (*vout)->ndata=tmp;
      VC_update(*vout);
      rslt=*vout;
    }
  bcopy(vin->ndata,rslt->ndata,vin->size*VC_sizeof(vin->kind));
  return rslt;
}

int VC_print(VC_vector vin,int level)
{

  switch (vin->kind)
    {
    case VCTNONE   : return VC_cprint(vin,level);
    case VCTCHAR   : return VC_cprint(vin,level);
    case VCTSHORT  : return VC_sprint(vin,level);
    case VCTINT    : return VC_iprint(vin,level);
    case VCTLONG   : return VC_lprint(vin,level);
    case VCTFLOAT  : return VC_fprint(vin,level);
    case VCTDOUBLE : return VC_dprint(vin,level);
    }
  return 1;
}

#define _DEFINE_FUNC_TEMPLATE_(_Fc,_Data,_Fmt)			\
int _Fc(VC_vector vin,int level)				\
{								\
  size_t i;							\
  OS_message(level,OS_STD,"[");					\
  if(vin->size)						\
    {								\
      for(i=0;i<(vin->size-1);i++)				\
	OS_message(level,OS_STD,#_Fmt ",",vin->_Data[i]);	\
      OS_message(level,OS_STD,#_Fmt,vin->_Data[vin->size-1]);	\
    }								\
  OS_message(level,OS_STD,"]");					\
 return 0;							\
}

_DEFINE_FUNC_TEMPLATE_(VC_cprint,cdata,%c)
_DEFINE_FUNC_TEMPLATE_(VC_sprint,sdata,%d)
_DEFINE_FUNC_TEMPLATE_(VC_iprint,idata,%d)
_DEFINE_FUNC_TEMPLATE_(VC_lprint,ldata,%d)
_DEFINE_FUNC_TEMPLATE_(VC_fprint,fdata,%f)
_DEFINE_FUNC_TEMPLATE_(VC_dprint,ddata,%f)

#undef _DEFINE_FUNC_TEMPLATE_
