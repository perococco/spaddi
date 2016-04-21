/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef UTILSH
#define UTILSH
#include <stdio.h>

int wordCount(char* filename,long* nbl,long* nbw,long* nbc);
int readNumberL(FILE* fstrm,long*   n,long* nd);
int readNumberD(FILE* fstrm,double* n);
//int readNumberF(FILE* fstrm,float*  n);
int gotoNextWord(FILE* fstrm);

#endif
