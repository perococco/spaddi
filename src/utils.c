/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <errno.h>
#include "utils.h"

int wordCount(char* filename,long* nbl,long* nbw,long* nbc)
{
  long  nlgn,nwrd,nchr;
  int   spc,tab,nwl,strl,strw,rchr;
  FILE *fstrm;

  if((fstrm=fopen(filename,"r"))==NULL) return 1;
  nlgn=nwrd=nchr=0;
  errno=0;
  spc=(int)' ';
  tab=(int)'\t';
  nwl=(int)'\n';
  strl=0;
  strw=0;
  while(1)
    {
      rchr=fgetc(fstrm);
      if((rchr==spc)||(rchr==tab))
	{
	  if(strw) {++nwrd;strw=0;}
	  continue;
	}
      if(rchr==nwl)
	{
	  ++nlgn;
	  if(strw) {++nwrd;strw=0;}
	  continue;
	}
      if(rchr==EOF) break;
      strw=1;
      ++nchr;
    }
  fclose(fstrm);
  if(errno) return 1;
  nlgn+=strw;
  if(nbl) *nbl=nlgn;
  if(nbw) *nbw=nwrd;
  if(nbc) *nbc=nchr;
  return 0;
}


int readNumberL(FILE* fstrm,long* n,long *nd)
{
  long v;
  int sgn,rchr,dgt,ndgt;

  sgn=0;
  ndgt=0;
  v=0;
  while(1)
    {
      rchr=fgetc(fstrm);
      //      printf("   rchr: %c\n",rchr);
      if(rchr==EOF) break;
      if((rchr=='-')||(rchr=='+'))
	{
	  if(sgn) return 1;
	  sgn=(rchr=='-')?-1:1;
	  continue;
	}
      dgt=rchr-'0';
      if((dgt>=0)&&(dgt<=9))
	{
	  ++ndgt;
	  v=v*10+dgt;
	  continue;
	}
      fseek(fstrm,-1,SEEK_CUR);
      break;
    }
  if(sgn<0) v=-v;
  if( n) *n=v;
  if(nd) *nd=ndgt;
  return 0;
}


int readNumberD(FILE* fstrm,double* n)
{
  long  ipL,fpL,xp,ndgL,ndgF;
  double v,ipD,fpD;

  int rchr,sgn;

  ipD=fpD=0.;
  v=0;sgn=1;
  if((rchr=fgetc(fstrm))==EOF) return 1;
  if((rchr=='-')||(rchr=='+'))
    {
      sgn=(rchr=='-')?-1:1;
      if((rchr=fgetc(fstrm))==EOF) return 1;
      if((rchr<'0')||(rchr>'9')) return 1;
    }
  fseek(fstrm,-1,SEEK_CUR);

  if(readNumberL(fstrm,&ipL,&ndgL)) return 1;
  //  printf("ipL : %d\n",ipL);
  ipD=(double)ipL;
  rchr=fgetc(fstrm);
  
  if(rchr==EOF)
    {
      if(n) *n=ipD;
      return 0;
    }
  
  if(rchr=='.')
    {
      if(readNumberL(fstrm,&fpL,&ndgF)) return 1;
      //            printf("fpL : %d\n",fpL);
      if(fpL<0) return 1;
      fpD=(double)fpL;
      while(ndgF--) fpD*=0.1;
      rchr=fgetc(fstrm);
    }
  v=fpD;
  v+=ipD;
  if(sgn<0) v=-v;

  if((rchr=='E')||(rchr=='e'))
    {
      if(readNumberL(fstrm,&xp,NULL)) return 1;
      if (xp<0)
	{
	  v*=0.1;
	  while(++xp<0) v*=0.1;
	}
      else if(xp>0)
	{
	  v*=10;
	  while(--xp>0) v*=10;
	}
      rchr=fgetc(fstrm);
    }

  if(n) *n=v;
  fseek(fstrm,-1,SEEK_CUR);
  if((rchr!=' ')&&(rchr!='\t')&&(rchr!='\n')&&(rchr!=EOF)) return 1;
  return 0;
}

int gotoNextWord(FILE* fstrm)
{
  int rchr;

  while(1)
    {
      rchr=fgetc(fstrm);
      if(rchr==EOF) return 1;
      if((rchr!=' ')&&(rchr!='\t')&&(rchr!='\n')) break;
      continue;
    }
  fseek(fstrm,-1,SEEK_CUR);
  return 0;
}
