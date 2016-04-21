/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "version.h"
#include "outstream.h"
#include "memory_manager.h"
#include "vector.h"
#include "table.h"
#include "spectrum.h"

int usageMessage(char* progname,int help)
{
  if(help)
    {
      OS_message(0,OS_ERR,"Spaddi version %4.2f\n",VERSION);
      OS_message(0,OS_ERR,"Usage : %s tfits_or_dat_spectra\n",progname);
      OS_message(0,OS_ERR,"Options:\n");
      OS_message(0,OS_ERR,"\t-h,-help     : print this help message\n");
      OS_message(0,OS_ERR,"\t-o file      : use 'file' as output file [added_spectrum.tfits]\n");
      OS_message(0,OS_ERR,"\t-ascii       : make a copy of the output in .dat format\n");
      OS_message(0,OS_ERR,"\t-nbsp        : put in a fourth column the number of\n");
      OS_message(0,OS_ERR,"\t               spectra used for the addition.\n");
      OS_message(0,OS_ERR,"\t-tclip nb    : type of clipping [2]\n");
      OS_message(0,OS_ERR,"\t               1 -> (flx-mean)/rms\n");
      OS_message(0,OS_ERR,"\t               2 -> (flx-mean)/rms without checked point\n");
      OS_message(0,OS_ERR,"\t               3 -> no clipping\n");
      OS_message(0,OS_ERR,"\t-clip clip   : clip value for the cosmic ray detection [3]\n"); 
      OS_message(0,OS_ERR,"\t-width width : size of the chunk for the scale factors computation [500]\n");
      OS_message(0,OS_ERR,"\t-terr nb     : Variance type\n");
      OS_message(0,OS_ERR,"\t               1 -> weigthed noise\n");
      OS_message(0,OS_ERR,"\t               2 -> rms of input / numberof of input\n");
      OS_message(0,OS_ERR,"\t               3 -> (#1^2+#2^2)\n");
      OS_message(0,OS_ERR,"\t               4 -> max(#1,#2)     [*]\n");
      OS_message(0,OS_ERR,"\t               5 -> max(#1,max #2)\n");
      OS_message(0,OS_ERR,"\t               6 -> #1 and #2\n");
      OS_message(0,OS_ERR,"\t-spline      : spline rebinning    [*]\n");
      OS_message(0,OS_ERR,"\t-linear      : linear rebinning\n");
      exit(1);
    }
  else
    {
      OS_message(0,OS_ERR,"Usage : %s tfits_or_dat_spectra [-o outfile] [-ascii] [-help|-h] [-clip clip] [-width width] [-terr error_type] [-spline | -linear] [-nbsp]\n",progname);
      exit(1);

    }
}

int main(int argc,char** argv)
{
  double width,clip,pourc;
  double rfct,rshift;
  size_t nbsp,i,j;
  char *tfitsout,*asciiout;
  size_t* nameidx;
  char  terminator[]=".tfits";
  long k,good;
  int ascii,prec,terr,spline,help,tclip,keepnbsp;
  SP_spectrum result;
  SP_spectrum* spectra;
  VC_vector which;


  MM_init();
  nbsp=0;
  width=500;
  pourc=0.90;
  tclip=_SIGMCLIPPING_;
  clip=3;
  ascii=0;
  tfitsout=NULL;
  i=1;
  prec=10;
  terr=_MAX_BASTIEN_AND_SRIANAND_;
  rfct=1.;
  rshift=0.;
  spline=1;
  help=0;
  keepnbsp=0;
  while(i<argc)
    {
      if(strcmp(argv[i],"-o")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  tfitsout=argv[i++];continue;
	}
      if(strcmp(argv[i],"-h")==0)
	{
	  usageMessage(argv[0],1);
	}
      if(strcmp(argv[i],"-help")==0)
	{
	  usageMessage(argv[0],1);
	}
      if(strcmp(argv[i],"-spline")==0)
	{
	  i++;
	  spline=1;continue;
	}
      if(strcmp(argv[i],"-linear")==0)
	{
	  i++;
	  spline=0;continue;
	}
      if(strcmp(argv[i],"-ascii")==0)
	{
	  i++;
	  ascii=1;continue;
	}      
      if(strcmp(argv[i],"-nbsp")==0)
	{
	  i++;
	  keepnbsp=1;continue;
	}      
      if(strcmp(argv[i],"-clip")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  clip=atof(argv[i++]);continue;
	}
      if(strcmp(argv[i],"-pourc")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  pourc=atof(argv[i++]);continue;
	}
      if(strcmp(argv[i],"-width")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  width=atof(argv[i++]);continue;
	}

      if(strcmp(argv[i],"-prec")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  prec=atoi(argv[i++]);continue;
	}


      if(strcmp(argv[i],"-terr")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  terr=atoi(argv[i++]);continue;
	}

      if(strcmp(argv[i],"-tclip")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  tclip=atoi(argv[i++]);continue;
	}

      if(strcmp(argv[i],"-rfct")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  rfct=atof(argv[i++]);continue;
	}

      if(strcmp(argv[i],"-rshift")==0)
	{
	  i++;if(i>=argc) usageMessage(argv[0],help);
	  rshift=atof(argv[i++]);continue;
	}
      nbsp++;
      i++;
    }
  
  if(nbsp<2) usageMessage(argv[0],help);

  nameidx=(size_t*)MM_malloc(nbsp*sizeof(size_t));
  i=1;j=0;
  while(i<argc)
    {
      if(strcmp(argv[i],     "-o")==0) {i++;i++;continue;}
      if(strcmp(argv[i], "-ascii")==0) {i++;continue;}
      if(strcmp(argv[i],  "-nbsp")==0) {i++;continue;}
      if(strcmp(argv[i],"-spline")==0) {i++;continue;}
      if(strcmp(argv[i],"-linear")==0) {i++;continue;}
      if(strcmp(argv[i],  "-clip")==0) {i++;i++;continue;}
      if(strcmp(argv[i], "-width")==0) {i++;i++;continue;}
      if(strcmp(argv[i],  "-prec")==0) {i++;i++;continue;}
      if(strcmp(argv[i],  "-terr")==0) {i++;i++;continue;}
      if(strcmp(argv[i], "-tclip")==0) {i++;i++;continue;}
      if(strcmp(argv[i],  "-rfct")==0) {i++;i++;continue;}
      if(strcmp(argv[i],"-rshift")==0) {i++;i++;continue;}
      nameidx[j++]=i;
      i++;
    }

  
  if(j!=nbsp)
    {
      OS_message(0,OS_ERR,"Euh bug ... at %s line %d\n",__FILE__,__LINE__);
      return 1;
    }

  
  if(!tfitsout) tfitsout="added_spectrum.tfits";
  good=5;
  for(k=strlen(tfitsout)-1;k>=0;k--)
    {
      if(tfitsout[k]==terminator[good])
	if(good==0) break;
	else --good;
      else break;
    }
  if(good==0)
    {
      k=strlen(tfitsout);
      asciiout=(char*)MM_malloc((k-1)*sizeof(char));
      for(i=0;i<k-5;i++) asciiout[i]=tfitsout[i];
      strcpy(asciiout+k-5,"dat");
    }
  else
    {
      k=strlen(tfitsout);
      asciiout=(char*)MM_malloc((k+5)*sizeof(char));
      for(i=0;i<k;i++) asciiout[i]=tfitsout[i];
      strcpy(asciiout+k,".dat");
    }

  which=VC_allocate(3,VCTLONG);
  which->ldata[0]=0;
  which->ldata[1]=1;
  which->ldata[2]=2;

  spectra=(SP_spectrum*)MM_malloc(nbsp*sizeof(SP_spectrum));
  for(i=0;i<nbsp;i++)
    {
      k=strlen(argv[nameidx[i]]);
      if((k<6)||strcmp((argv[nameidx[i]]+k-6),terminator))
	spectra[i]=SP_readAscii(argv[nameidx[i]],NULL);
      else
	spectra[i]=SP_readTfits(argv[nameidx[i]],which);
    }
  
  printf("\n");
  printf("---------- Spaddi options ----------\n");
  printf("tclip    : %d\n",tclip);
  printf("clip     : %f\n",clip);
  printf("terr     : %d\n",terr);
  printf("rebinnig : %s\n",(spline?"spline":"linear"));
  printf("width    : %f\n",width);
  printf("ascii    : %d\n",ascii);
  printf("tfitsout : %s\n",tfitsout);
  printf("asciiout : %s\n",asciiout);
  printf("------------------------------------\n");
  printf("\n");

  result=SP_addSpectra(spectra,nbsp,width,clip,pourc,terr,rfct,rshift,spline,tclip,keepnbsp);
  
/*   for(i=0;i<result->npix;i++) */
/*     printf("%06d %f %f %f\n",i,result->data[0]->ddata[i], */
/* 	   result->data[1]->ddata[i], */
/* 	   result->data[2]->ddata[i]); */

  SP_writeTfits(result,tfitsout,NULL);
  if(ascii) SP_writeAscii(result,asciiout,NULL,prec);
  SP_free(result);
  
  for(i=0;i<nbsp;i++)
    SP_free(spectra[i]);
  MM_free(spectra);
  VC_free(which);
  MM_free(nameidx);
  MM_free(asciiout);

  MM_finish();
  return 0;
}
