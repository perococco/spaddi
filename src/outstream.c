/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include "outstream.h"

int   outputVerbosityLevel=-1;
int   outputType=OS_TOUT_STD;
FILE* outputStream=NULL;

void OS_init(OS_TOUT ostype,int level,char* filename)
{
  outputVerbosityLevel=level;

  if((ostype==OS_TOUT_NONE)||(ostype==OS_TOUT_STD))
    {
      outputType=ostype;
      return;
    }
  if(ostype==OS_TOUT_FILE)
    {
      outputType=ostype;
      if(filename)
	{
	  outputStream=fopen(filename,"a");
	  if(!outputStream)
	    {
	      printf("Cannot open file '%s' ...\nUse Standard Output ...\n",filename);
	      outputType=OS_TOUT_STD;
	    }
	}
      else
	{
	  outputStream=fopen("output.log","a");
	  if(!outputStream)
	    {
	      printf("Cannot open file 'output.log' ...\nUse Standard Output ...\n");
	      outputType=OS_TOUT_STD;
	    }
	}
      return;
    }
  outputType=OS_TOUT_STD;
}

void OS_setLevel(int level)
{
  outputVerbosityLevel=level;
}

void OS_message(int level,OS_TMES type,char *fmt, ...)
{
  va_list argp;

  level=(level<0)?0:level;
  if((outputVerbosityLevel>=0)&&(outputVerbosityLevel<=level)) return;
  if(OS_TOUT_NONE) return;
  va_start(argp,fmt);
  switch (outputType) 
    {
    case OS_TOUT_STD  : vprintf(fmt,argp);break;
    case OS_TOUT_FILE : vfprintf(outputStream,fmt,argp);break;
    }
  va_end(argp);
  return;
}

void OS_finish()
{
  if(outputType==OS_TOUT_FILE)
    fclose(outputStream);
}

void OS_flush()
{
  switch (outputType) 
    {
    case OS_TOUT_STD  : fflush(stdout);break;
    case OS_TOUT_FILE : fflush(outputStream);break;
    }
}
