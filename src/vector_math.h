/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef VECTORMATHH
#define VECTORMATHH

#include "vector.h"

#define VC_span(_inf,_sup,_nbpt) _VC_span(_inf,_sup,_nbpt,__FILE__,__LINE__)
#define VC_binSize(_vin,_vout)   _VC_binSize(_vin,_vout,__FILE__,__LINE__)
#define VC_rebinl(_vy,_vx,_nvx,_nvy)  _VC_rebinl(_vy,_vx,_nvx,_nvy,__FILE__,__LINE__)
#define VC_rebins(_vy,_vx,_nvx,_nvy)  _VC_rebins(_vy,_vx,_nvx,_nvy,__FILE__,__LINE__)
#define VC_sqr(_vin,_vout)       _VC_sqr(_vin,_vout,__FILE__,__LINE__)
#define VC_sqrt(_vin,_vout)       _VC_sqrt(_vin,_vout,__FILE__,__LINE__)
#define VC_spline(_vy,_vx,_yp1,_ypn,_vout) _VC_spline(_vy,_vx,_yp1,_ypn,_vout,__FILE__,__LINE__)
#define VC_splint(_vya,_vxa,_vx,_vyda,_vout) _VC_splint(_vya,_vxa,_vx,_vyda,_vout,__FILE__,__LINE__)

VC_vector _VC_span(double inf,double sup,size_t nbpt,char* file,size_t line);
VC_vector _VC_binSize(VC_vector vin,VC_vector* vout,char* file,size_t line);

VC_vector _VC_sqr(VC_vector vin,VC_vector* vout,char* file,size_t line);
VC_vector _VC_sqrt(VC_vector vin,VC_vector* vout,char* file,size_t line);

VC_vector _VC_cbinSize(VC_vector vin,VC_vector* vout,char* file,size_t line);
VC_vector _VC_sbinSize(VC_vector vin,VC_vector* vout,char* file,size_t line);
VC_vector _VC_ibinSize(VC_vector vin,VC_vector* vout,char* file,size_t line);
VC_vector _VC_lbinSize(VC_vector vin,VC_vector* vout,char* file,size_t line);
VC_vector _VC_fbinSize(VC_vector vin,VC_vector* vout,char* file,size_t line);
VC_vector _VC_dbinSize(VC_vector vin,VC_vector* vout,char* file,size_t line);

VC_vector _VC_rebinl(VC_vector viny,VC_vector vinx,VC_vector voutx,VC_vector* vouty,char* file,size_t line);

VC_vector _VC_spline(VC_vector y,VC_vector x,double yp1,double ypn,VC_vector* vout,char* file,size_t line);
VC_vector _VC_splint(VC_vector ya, VC_vector xa,VC_vector x,VC_vector y2a,VC_vector* vout,char* file,size_t line);
double splint(VC_vector ya, VC_vector xa,double x, VC_vector y2a);

double _qgaus(double (*func)(VC_vector ya, VC_vector xa,double x, VC_vector y2a),
	      double a, double b,VC_vector ya, VC_vector xa, VC_vector y2a);

VC_vector _VC_rebins(VC_vector viny,VC_vector vinx,VC_vector voutx,VC_vector* vouty,char* file,size_t line);

#endif
