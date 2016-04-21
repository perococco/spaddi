/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "outstream.h"
#include "table.h"
#include "table_tfits.h"
#include "table_ascii.h"


int main(int argc,char** argv)
{
  TB_table table;
  int prec,init_nulval;
  double nulval;

  prec=10;
  init_nulval=0;
  OS_init(OS_TOUT_STD,-1,NULL);

  if((argc<3)||(argc==4)||(argc==6)||(argc>8))
    {
      OS_message(0,OS_ERR,"Usage : %s  tfits_filein  ascii_fileout [-p prec] [-n nulval] \n",argv[0]);
      return 1;
    }

  if(argc==5)
    {
      if(strcmp(argv[3],"-p")==0)
	prec=atoi(argv[4]);
      if(strcmp(argv[3],"-n")==0)
	{
	  init_nulval=1;
	  nulval=atof(argv[4]);
	}
    }

  if(argc==7)
    {
      if(strcmp(argv[5],"-p")==0)
	prec=atoi(argv[6]);
      if(strcmp(argv[5],"-n")==0)
	{
	  init_nulval=1;
	  nulval=atof(argv[4]);
	}
    }

  table=readTfits(argv[1],NULL,((init_nulval)?&nulval:NULL));
  writeTAscii(table,argv[2],NULL,prec);
  
  TB_free(table);
  return 0;
}
