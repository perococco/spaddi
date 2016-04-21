/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdlib.h>
#include <string.h>
#include "memory_manager.h"
#include "table.h"


TB_table _TB_allocate(char* name,size_t nbCol,size_t nbRow,char* file,size_t line)
{
  TB_table rslt;
  size_t i,j;
  
  rslt=(TB_table)_MM_malloc(sizeof(TB_table_ref),file,line);
  if(!rslt) return NULL;

  rslt->nbCol=nbCol;
  rslt->nbRow=nbRow;
  rslt->name=(char*)MM_malloc((strlen(name)+1)*sizeof(char));
  if(!rslt->name)
    {
      _MM_free(rslt,file,line);
      return NULL;
    }
  strcpy(rslt->name,name);

  rslt->labels=(char**)MM_malloc(nbCol*sizeof(char*));
  if(!rslt->labels)
    {
      _MM_free(rslt->name,file,line);
      _MM_free(rslt,file,line);
      return NULL;
    }

  rslt->data=(VC_vector*)MM_malloc(nbCol*sizeof(VC_vector));
  if(!rslt->data)
    {
      _MM_free(rslt->labels,file,line);
      _MM_free(rslt->name,file,line);
      _MM_free(rslt,file,line);
      return NULL;
    }

  for(i=0;i<nbCol;i++)
    {
      rslt->labels[i]=NULL;
      rslt->data[i]=VC_allocate(nbRow,VCTDOUBLE);
      if(!rslt->data[i])
	{
	  for(j=0;j<i;j++)
	      VC_free(rslt->data[j]);
	  _MM_free(rslt->data,file,line);
	  _MM_free(rslt->labels,file,line);
	  _MM_free(rslt->name,file,line);
	  _MM_free(rslt,file,line);
	  return NULL;
	}
    }
  return rslt;
}

void _TB_free(TB_table table,char* file,size_t line)
{
  size_t i;
  
  for(i=0;i<table->nbCol;i++)
    {
      _VC_free(table->data[i],file,line);
      if(table->labels[i]) _MM_free(table->labels[i],file,line);
    }
  _MM_free(table->data,file,line);
  _MM_free(table->labels,file,line);
  _MM_free(table->name,file,line);
  _MM_free(table,file,line);
}
