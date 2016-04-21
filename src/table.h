/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef TABLEH
#define TABLEH

#include "vector.h"

#define TB_allocate(_name,_nbCol,_nbRow) _TB_allocate(_name,_nbCol,_nbRow,__FILE__,__LINE__)
#define TB_free(_tabl)                   _TB_free(_tabl,__FILE__,__LINE__)


typedef struct {
  char*      name;
  size_t     nbCol;
  size_t     nbRow;
  char**     labels;
  VC_vector* data;
} TB_table_ref;
  
typedef TB_table_ref* TB_table;


TB_table _TB_allocate(char* name,size_t nbCol,size_t nbRow,char* file,size_t line);
void     _TB_free(TB_table table,char* file,size_t line);



#endif
