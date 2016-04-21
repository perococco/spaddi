/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef VECTORSORTH
#define VECTORSORTH
#include "vector.h"

int VC_hpsort(VC_vector vcin);
int VC_chpsort(VC_vector vcin);
int VC_shpsort(VC_vector vcin);
int VC_ihpsort(VC_vector vcin);
int VC_lhpsort(VC_vector vcin);
int VC_fhpsort(VC_vector vcin);
int VC_dhpsort(VC_vector vcin);

#define VC_hpsort_idx(_vcin)  _VC_hpsort_idx(_vcin,__FILE__,__LINE__)
#define VC_chpsort_idx(_vcin) _VC_chpsort_idx(_vcin,__FILE__,__LINE__)
#define VC_shpsort_idx(_vcin) _VC_shpsort_idx(_vcin,__FILE__,__LINE__)
#define VC_ihpsort_idx(_vcin) _VC_ihpsort_idx(_vcin,__FILE__,__LINE__)
#define VC_lhpsort_idx(_vcin) _VC_lhpsort_idx(_vcin,__FILE__,__LINE__)
#define VC_fhpsort_idx(_vcin) _VC_fhpsort_idx(_vcin,__FILE__,__LINE__)
#define VC_dhpsort_idx(_vcin) _VC_dhpsort_idx(_vcin,__FILE__,__LINE__)


VC_vector _VC_hpsort_idx(VC_vector vcin,char* file,size_t line);
VC_vector _VC_chpsort_idx(VC_vector vcin,char* file,size_t line);
VC_vector _VC_shpsort_idx(VC_vector vcin,char* file,size_t line);
VC_vector _VC_ihpsort_idx(VC_vector vcin,char* file,size_t line);
VC_vector _VC_lhpsort_idx(VC_vector vcin,char* file,size_t line);
VC_vector _VC_fhpsort_idx(VC_vector vcin,char* file,size_t line);
VC_vector _VC_dhpsort_idx(VC_vector vcin,char* file,size_t line);


#endif
