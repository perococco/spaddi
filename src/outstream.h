/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef OUTSTREAMH
#define OUTSTREAMH

#include <stdio.h>
#include <stdarg.h>

typedef enum  {OS_TOUT_NONE,OS_TOUT_STD,OS_TOUT_FILE} OS_TOUT;
typedef enum  {OS_STD,OS_WARN,OS_ERR} OS_TMES; 

extern int   outputVerbosityLevel;
extern int   outputType;
extern FILE* outputStream;


void OS_init(OS_TOUT ostype,int level,char* filename);
void OS_setLevel(int level);
void OS_message(int level,OS_TMES type,char* fmt, ...);
void OS_finish();
void OS_flush();

#endif
