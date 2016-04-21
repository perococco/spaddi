/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef TABLETFITSH
#define TABLETFITSH

#define readTfits(_file,_which,_nullval) _readTfits(_file,_which,_nullval,__FILE__,__LINE__)

TB_table _readTfits(char* filename,VC_vector which,double* nullval,char* file,size_t line);
int      writeTfits(TB_table tbin,char* filename,VC_vector which);

#endif
