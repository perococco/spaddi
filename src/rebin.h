/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef REBINH
#define REBINH

#define rebinl(a,b,c,d,e,f) _rebinl(a,b,c,d,e,f,__FILE__,__LINE__)

int _rebinl(VC_vector signal,VC_vector noise,VC_vector xinput,
	    VC_vector xoutput,VC_vector* sout,VC_vector* nout,
	    char* file,size_t line);

int  rebinLocal(double *signal,double  *noise,double* xin,size_t szin,
		double    xmin,double    xmax,
		long int *kmin,long int *kmax,
		double   *sout,double   *nout);

#endif
