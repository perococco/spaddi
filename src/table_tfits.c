/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdlib.h>
#include <string.h>

#include <fitsio.h>

#include "memory_manager.h"
#include "outstream.h"
#include "vector.h"
#include "table.h"
#include "table_tfits.h"






TB_table _readTfits(char* filename,VC_vector which,double* nullval,char* file,size_t line)
{
  fitsfile *fptr;
  int status,anynul,i;
  int naxis,hdunb,hdutype;
  int nbcol,tot_col,curCol;
  VC_vector naxes;
  char keyw[81];
  TB_table rslt;
  char value[FLEN_VALUE];

  /* First Open the file */
  status=0;
  status=fits_open_file(&fptr,filename,READONLY,&status);
  if(status) 
    {
      OS_message(0,OS_ERR,"[readTfits:%d] Cannot open file '%s' !!\n",status,filename);
      exit(1);
    }
  OS_message(1,OS_STD,"[ReadTfits] Open file '%s'\n",filename);

  /* Find the good header */
  hdunb=0;
  do
    {
      hdunb++;
      status=0;
      status=fits_movabs_hdu(fptr,hdunb,&hdutype,&status);
      if(status)
	{
	  OS_message(0,OS_ERR,"[readTfits:%d] Error when searchind hdu number \n",status);
	  exit(1);
	}
    } while((hdutype!=1)&&(hdutype!=2));
  OS_message(1,OS_STD,"[readTfits] Select HDU number #%d\n",hdunb);
  /*Read NAXIS KEYWORD*/
  status=0;
  status=fits_read_key(fptr,TINT,"NAXIS",&naxis,NULL,&status);
  if(status)
    {
      OS_message(0,OS_ERR,"[readTfits:%d] Cannot get NAXIS keyword \n",status);
      exit(1);
    }
  
  naxes=VC_allocate(naxis,VCTINT);
  for(i=0;i<naxis;i++)
    {
      sprintf(keyw,"NAXIS%d",i+1);
      status=0;
      status=fits_read_key(fptr,TINT,keyw,&(naxes->idata[i]),NULL,&status);
      if(status)
	{
	  VC_free(naxes);
	  OS_message(0,OS_ERR,"[readTfits:%d] Cannot get %s keyword \n",status,keyw);
	  exit(1);
	}
    }
  status=0;
  status=fits_read_key(fptr,TINT,"TFIELDS",&nbcol,NULL,&status);
  if(status)
    {
      VC_free(naxes);
      OS_message(0,OS_ERR,"[readTfits:%d] Cannot get TFIELDS keyword \n",status);
      exit(1);
    }

  OS_message(1,OS_STD,"[readTfits] Number of columns : %6d\n",nbcol);
  OS_message(1,OS_STD,"[readTfits] Number of rows    : %6d\n",naxes->idata[1]);

  tot_col=(which)?which->size:nbcol;
  rslt=_TB_allocate(filename,tot_col,naxes->idata[1],file,line);
  
  for(i=0;i<tot_col;i++)
    {
      curCol=(which)?which->ldata[i]:i;
      ++curCol; //starting index for tfits is 1 !!
      sprintf(keyw,"TTYPE%d",curCol);
      status=0;
      status=fits_read_key_str(fptr,keyw,value,NULL,&status);
      if(status)
	{
	  VC_free(naxes);
	  OS_message(0,OS_ERR,"[readTFits] Cannot get keyword %s\n",keyw);
	  exit(1);
	}
      rslt->labels[i]=(char*)MM_malloc((strlen(value)+1)*sizeof(char));
      strcpy(rslt->labels[i],value);
      
      OS_message(1,OS_STD,"[readTfits] Read column %3d   :  %s\n",curCol,value);
      status=0;
      status=fits_read_col(fptr,TDOUBLE,curCol,1,1,naxes->idata[1],nullval,
			   rslt->data[i]->ddata,&anynul,&status);
      if(anynul)
	{
	  OS_message(0,OS_WARN,"[readTFits] Some undefined value !\n");
	}
      if(status)
	{
	  VC_free(naxes);
	  OS_message(0,OS_ERR,"[readTfits] Cannot read column #%d\n",i);
	  exit(1);
	}
    }
  VC_free(naxes);
  status=0;
  status=fits_close_file(fptr,&status);
  if(status)
    {
      OS_message(0,OS_ERR,"[readTfits] Cannot close file '%s'!!\n",filename);
      exit(1);
    }
  return rslt;
}


int writeTfits(TB_table tbin,char* filename,VC_vector which)
{
  fitsfile *fptr;
  int      status;
  size_t   i,nb_columns,curCol;
  int      naxis;
  long     naxes,tmp;

  char**            ctform;
  char**            cttype;
  char**            ctunit;

  nb_columns=(which)?which->size:tbin->nbCol;
  status=0;
  status=fits_create_file(&fptr,filename,&status);
  if(status)
    {
      OS_message(0,OS_ERR,"[writeTfits] Cannot create file '%s'!!\n",filename);
      exit(1);
    }

  
  naxis=0;
  naxes=0;
  status=0;
  status=fits_create_img(fptr,8,naxis,&naxes,&status);
  if(status)
    {
      OS_message(0,OS_ERR,"[writeTfits] Cannot create image for file '%s'!!\n",filename);
      exit(1);
    }
  
  ctform=(char **)MM_malloc(nb_columns*sizeof(char*));
  cttype=(char **)MM_malloc(nb_columns*sizeof(char*));
  ctunit=(char **)MM_malloc(nb_columns*sizeof(char*));
  
  tmp=sizeof(char);
  for(i=0;i<nb_columns;i++)
    {
      curCol=(which)?which->ldata[i]:i;
      ctform[i]=(char*)MM_malloc(tmp*3);
      ctunit[i]=(char*)MM_malloc(tmp*2);
      cttype[i]=(char*)MM_malloc(tmp*(strlen(tbin->labels[curCol])+1));
      strcpy(ctform[i],"1D");
      strcpy(ctunit[i]," ");
      strcpy(cttype[i],tbin->labels[curCol]);
    }

  status=0;
  status=fits_create_tbl(fptr,BINARY_TBL,tbin->nbRow,nb_columns,cttype,
			 ctform,ctunit,NULL,&status);

  for(i=0;i<nb_columns;i++)
    {
      MM_free(cttype[i]);
      MM_free(ctunit[i]);
      MM_free(ctform[i]);
    }
  MM_free(ctunit);
  MM_free(cttype);
  MM_free(ctform);

  if(status)
    {
      OS_message(0,OS_ERR,"[writeTfits] Cannot create table for file '%s'!!\n",filename);
      exit(1);
    }

  for(i=0;i<nb_columns;i++)
    {
      curCol=(which)?which->ldata[i]:i;
      status=0;
      status=fits_write_col(fptr,TDOUBLE,(int)(i+1),1,1,tbin->nbRow,
			    tbin->data[curCol]->ddata,&status);
      if(status)
	{
	  OS_message(0,OS_ERR,"[writeTfits] Cannot write column %s !!\n",tbin->labels[curCol]);
	  exit(1);
	}
    }
  status=0;
  status=fits_close_file(fptr,&status);
  if(status)
    {
      OS_message(0,OS_ERR,"[writeTfits] Cannot close file '%s'!!\n",filename);
      exit(1);
    }
  return 0;
}
