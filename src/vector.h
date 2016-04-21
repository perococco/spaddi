/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef VECTORH
#define VECTORH

#define VC_allocate(_size,_kind)  _VC_allocate(_size,_kind,__FILE__,__LINE__)
#define VC_resize(_vct,_nsize)    _VC_resize(_vct,_nsize,__FILE__,__LINE__)
#define VC_free(_vct)             _VC_free(_vct,__FILE__,__LINE__)
#define VC_init(_kind,_val,_size) _VC_init(_kind,_val,_size,__FILE__,__LINE__)
#define VC_copy(_vin,_vout)       _VC_copy(_vin,_vout,__FILE__,__LINE__)
#define VCV_allocate(var,size,kind) var=_VCV_allocate(size,kind,#var,__FILE__,__LINE__)

typedef enum {VCTNONE,VCTCHAR,VCTSHORT,VCTINT,VCTLONG,VCTFLOAT,VCTDOUBLE} VC_TYPE;
  
typedef struct {
  VC_TYPE kind;
  size_t  size;

  void*    ndata;
  char*    cdata;
  short*   sdata;
  int*     idata;
  long*    ldata;
  float*   fdata;
  double*  ddata;
} VC_vector_ref;

typedef VC_vector_ref* VC_vector;



size_t    VC_sizeof(VC_TYPE kind);
void      VC_update(VC_vector vector);
VC_vector _VC_allocate(size_t size,VC_TYPE kind,char* file,size_t line);
size_t    _VC_resize(VC_vector vector,size_t nsize,char* file,size_t line);
void      _VC_free(VC_vector vector,char* file,size_t line);
int       VC_get(VC_vector vct,void* value,VC_TYPE kind,long start,long end,long stride);
int       VC_set(VC_vector vct,void* value,VC_TYPE kind,long start,long end,long stride);
VC_vector _VCV_allocate(size_t size,VC_TYPE kind,char* var,char* file ,size_t line);

#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c
#define GLUE4(a,b,c,d) a##b##c##d

#define _DEFINE_TEMPLATE_FUNC_(_Fc,_Car,_Tp)						       \
int GLUE3(VC_,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride);  \
int GLUE3(VC_c,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride); \
int GLUE3(VC_s,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride); \
int GLUE3(VC_i,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride); \
int GLUE3(VC_l,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride); \
int GLUE3(VC_f,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride); \
int GLUE3(VC_d,_Fc,_Car)(VC_vector vct,_Tp* value,long start,long end,long stride);


#define _DEFINE_TEMPLATE_FUNC2_(_Car,_Tp)	\
_DEFINE_TEMPLATE_FUNC_(get,_Car,_Tp);		\
_DEFINE_TEMPLATE_FUNC_(set,_Car,_Tp);

_DEFINE_TEMPLATE_FUNC2_(c, char  );
_DEFINE_TEMPLATE_FUNC2_(s, short );
_DEFINE_TEMPLATE_FUNC2_(i,   int);
_DEFINE_TEMPLATE_FUNC2_(l,  long);
_DEFINE_TEMPLATE_FUNC2_(f, float);
_DEFINE_TEMPLATE_FUNC2_(d,double);


#undef _DEFINE_TEMPLATE_FUNC2_
#undef _DEFINE_TEMPLATE_FUNC_

#undef GLUE2
#undef GLUE3
#undef GLUE4

VC_vector _VC_init(VC_TYPE kind,void* values,size_t size,char* file,size_t line);
VC_vector _VC_copy(VC_vector vin,VC_vector *vout,char* file,size_t line);
int       VC_print(VC_vector vin,int level);
int       VC_cprint(VC_vector vin,int level);
int       VC_sprint(VC_vector vin,int level);
int       VC_iprint(VC_vector vin,int level);
int       VC_lprint(VC_vector vin,int level);
int       VC_fprint(VC_vector vin,int level);
int       VC_dprint(VC_vector vin,int level);

#endif
