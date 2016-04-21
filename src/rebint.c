/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "vector_yor.h"
#include "rebin.h"


int main(int argc,char** argv)
{
  VC_vector signal,noise,xinput;
  VC_vector sout,nout,xoutput;

  if (argc<7)
    {
      fprintf(stderr,"Usage %s signal noise xinput xouput sout nout\n",argv[0]);
      return 1;
    }

  signal =binRead(argv[1],NULL);
  noise  =binRead(argv[2],NULL);
  xinput =binRead(argv[3],NULL);
  xoutput=binRead(argv[4],NULL);

  rebinl(signal,noise,xinput,xoutput,&sout,&nout);
  
  binWrite(argv[5],sout);
  binWrite(argv[6],nout);
  
  return 0;
  

}
