/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef ASTROH
#define ASTROH

#define cosmicMask(_flux,_sigm,_clip,_nbc,_mask) _cosmicMask(_flux,_sigm,_clip,_nbc,_mask,__FILE__,__LINE__);
#define addiScaleFactors(_spc,_nbsp,_width,_vout) _addiScaleFactors(_spc,_nbsp,_width,_vout,__FILE__,__LINE__);
#define mergeScaleFactors(_spc,_nbsp,_width,_vout) _mergeScaleFactors(_spc,_nbsp,_width,_vout,__FILE__,__LINE__);
#define sigmaClipping(_a,_b,_c,_d) _sigmaClipping(_a,_b,_c,_d,__FILE__,__LINE__)
#define dataClipping(_a,_b,_c,_d) _dataClipping(_a,_b,_c,_d,__FILE__,__LINE__)

VC_vector _cosmicMask(VC_vector flux,VC_vector sigm,double clip,size_t* nbc,VC_vector* mask,char* file,size_t line);
VC_vector correctFlux(VC_vector* flux,VC_vector mask);
VC_vector _addiScaleFactors(SP_spectrum* spectra,size_t nbsp,double width,VC_vector* vout,char* file,size_t line);
VC_vector _mergeScaleFactors(SP_spectrum* spectra,size_t nbsp,long width,VC_vector* vout,char* file,size_t line);
double    median(double* tab,size_t size);
int fitLine(double* y,double* x,size_t size,double* a,double* b);

int       _super(const void* a,const void* b);
size_t    findNear(VC_vector vin,double val);
VC_vector _sigmaClipping(VC_vector flux,VC_vector* mask,double clip,int verb,char *file,size_t line);
VC_vector _dataClipping(VC_vector flux,VC_vector sigm,VC_vector* mask,double clip,char* file,size_t line);
int VC_getIdxEdge(VC_vector data,long* idmin,long* idmax);
int fitLineIter(double* y,double* x,size_t size,double* a,double* b);
int statistique(double* data,long n,double* mask,double* mean,double* rms);


#endif
