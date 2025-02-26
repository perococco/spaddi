/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "memory_manager.h"
#include "outstream.h"
#include "vector.h"
#include "table.h"
#include "table_ascii.h"

void writeTAscii(TB_table tbin,char* filename,VC_vector which,int prec)
{
  FILE* stream;
  size_t i,j,nbc,nbr;


  if(strcmp(filename,"-"))
    {
      if(!(stream=fopen(filename,"w")))
	{
	  OS_message(0,OS_ERR,"[writeTascii] Cannot open file '%f' !!\n",filename);
	  exit(1);
	}
    }
  else
    stream=stdout;

  nbc=(which)?which->size:tbin->nbCol;
  nbr=tbin->nbRow;
  for(i=0;i<nbr;i++)
    for(j=0;j<nbc;j++)
      {
	fprintf(stream,"%.*f%1s",prec,tbin->data[j]->ddata[i],((j==(nbc-1))?"\n":" "));
      }
  fclose(stream);
  return;
}
