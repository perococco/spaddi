/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef SPECTRUMH
#define SPECTRUMH
#include "vector.h"

#define _BASTIEN_ 1
#define _SRIANAND_ 2
#define _BASTIEN_AND_SRIANAND_ 3
#define _MAX_BASTIEN_AND_SRIANAND_ 4
#define _MAX_BASTIEN_AND_MAX_SRIANAND_ 5
#define _BOTH_ERRORS_ 6

#define _DATACLIPPING_ 1
#define _SIGMCLIPPING_ 2
#define _NOCLIPPING_ 3

#define SP_allocate(_a,_b) _SP_allocate(_a,_b,__FILE__,__LINE__)
#define SP_free(_a) _SP_free(_a,__FILE__,__LINE__)
#define SP_readTfits(_a,_b) _SP_readTfits(_a,_b,__FILE__,__LINE__)
#define SP_readAscii(_a,_b) _SP_readAscii(_a,_b,__FILE__,__LINE__)
#define SP_writeTfits(_a,_b,_c) _SP_writeTfits(_a,_b,_c,__FILE__,__LINE__)
#define SP_writeAscii(_a,_b,_c,_d) _SP_writeAscii(_a,_b,_c,_d,__FILE__,__LINE__)
#define SP_addSpectra(_a,_b,_c,_d,_e,_f,_g,_h,_i,_j,_k) _SP_addSpectra(_a,_b,_c,_d,_e,_f,_g,_h,_i,_j,_k,__FILE__,__LINE__)
#define SP_groupSpectrum(_a,_b,_c) _SP_groupSpectrum(_a,_b,_c,__FILE__,__LINE__)

typedef struct 
{

  char* name;
  size_t npix;
  
  long colW;
  long colF;
  long colS;
  long colC;

  size_t    ndata;
  VC_vector* data;

  double* wave;
  double* flux;
  double* sigm;
  double* cont;
} SP_spectrum_ref;

typedef SP_spectrum_ref* SP_spectrum;

SP_spectrum _SP_addSpectra(SP_spectrum* spectra,size_t nbsp,double width,double clip,double pourc,int error_type,double rebin_pixfactor,double rebin_pixshift,int spline,int tclip,int keepnbsp,char* file,size_t line);

SP_spectrum _SP_allocate(size_t npix,size_t nd,char *file,size_t line);
void        _SP_free(SP_spectrum sp,char* file,size_t line);
SP_spectrum _SP_readTfits(char* filename,VC_vector which,char* file,size_t line);
SP_spectrum _SP_readAscii(char* filename,VC_vector which,char* file,size_t line);
int         _SP_writeTfits(SP_spectrum spin,char* filename,VC_vector which,char* file,size_t line);
int         _SP_writeAscii(SP_spectrum spin,char* filename,VC_vector which,int prec,char* file,size_t line);
VC_vector   _SP_groupSpectrum(SP_spectrum* spectra,size_t nbsp,double pourc,char* file,size_t line);
double      pourcOver(double a1,double b1,double a2,double b2);
double expFilter(double x,double x0,double a);
void        SP_sortSpectrum(SP_spectrum* spectra,size_t nbsp);
int getSplineInter(double* flux,double* sigm,double* index,long npix,long idx,long* kpos,double* fout,double* sout);

#endif
