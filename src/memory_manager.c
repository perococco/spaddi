/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#include <stdio.h>
#include <stdlib.h>
#include "memory_manager.h"

MM_memItem* MM_list=NULL;
size_t      MM_list_size=0;
size_t      MM_list_last=0;

#define  Balise(a) printf("## [%s at %d] %s\n",__FILE__,__LINE__,a);fflush(stdout);

void MM_info()
{
  size_t i;
  printf("-- MEMORY MANAGER INFORMATIONS --\n");
  printf(" MM_list : adr %lx , size %lu , length %lu\n",(size_t)MM_list,MM_list_size,MM_list_last);
  for(i=0;i<MM_list_last;i++)
    printf(" MM_list[%4lu] : adr=%lx %s at %s line %lu .\n",i,(size_t)MM_list[i].adr,(MM_list[i].varname)?(MM_list[i].varname):"unknown",MM_list[i].file,MM_list[i].line);
  printf("---------------------------------\n");
}

int MM_deleteItems()
{
  void* tmp;
  size_t newl;

  if(MM_list_last==0) newl=1;
  else newl=(MM_list_last-1)/MM_list_stride+1;

  if((MM_list_last+MM_list_stride+MM_list_stride/2)>=MM_list_size) return 0;
  tmp=realloc(MM_list,(newl*MM_list_stride)*sizeof(MM_memItem));
  if(tmp)
    {
      MM_list=tmp;
      MM_list_size=newl;
      return 0;
    }
  fprintf(stderr,"[MEMORY MANAGER] !! Cannot reallocate MM_list !!\n");
  return 1;
}

int MM_addItem()
{
  void* tmp;
  size_t i;
  if(MM_list_last<MM_list_size) return 0;

  
  tmp=(MM_memItem*)realloc(MM_list,(MM_list_size+MM_list_stride)*sizeof(MM_memItem));
  if(tmp)
    {
      MM_list=tmp;
      for(i=MM_list_size;i< MM_list_size+MM_list_stride;i++)
	MM_list[i].varname="unknown";
      MM_list_size+=MM_list_stride;
      return 0;
    }
  fprintf(stderr,"[MEMORY MANAGER] !! No more place left for MM_list !!\n");
  return 1;
}


MM_memItem* MM_init()
{
  MM_list_size=0;
  MM_list_last=0;
  MM_list=NULL;

  MM_addItem();
  return MM_list;
}

void *_MM_calloc(size_t num_of_elts,size_t elt_size,char* file,size_t line)
{
  void* tmp;
  tmp=calloc(num_of_elts,elt_size);
  if(!tmp) return NULL;
  if(!MM_addItem())
    {
      MM_list[MM_list_last].adr=tmp;
      MM_list[MM_list_last].line=line;
      MM_list[MM_list_last].file=file;
      ++MM_list_last;
    }
  return tmp;
}

void  _MM_free(void *pointer,char* file,size_t line)
{
  size_t i,j;
  char ok;
  ok=0;
  for(i=0;i<MM_list_last;i++)
    if((ok=(pointer==MM_list[i].adr))) break;

  if(ok)
    {
      for(j=i+1;j<MM_list_last;j++)
	MM_list[j-1]=MM_list[j];
      MM_list_last--;
      MM_deleteItems();
    }
  else
    fprintf(stderr,"[MEMORY MANAGER] !! Warning !! Pointer not in MM_list at '%s' line %ld \n",file,line);

  free(pointer);
}

void *_MM_malloc(size_t size,char* file,size_t line)
{
  void* tmp;
  tmp=malloc(size);
  if(!tmp) return NULL;
  if(!MM_addItem())
    {
      MM_list[MM_list_last].adr=tmp;
      MM_list[MM_list_last].line=line;
      MM_list[MM_list_last].file=file;
      ++MM_list_last;
    }
  return tmp;  
}

void *_MMV_malloc(size_t size,char* var,char* file,size_t line)
{
  void* tmp;
  tmp=malloc(size);
  if(!tmp) return NULL;
  if(!MM_addItem())
    {
      MM_list[MM_list_last].adr=tmp;
      MM_list[MM_list_last].line=line;
      MM_list[MM_list_last].file=file;
      MM_list[MM_list_last].varname=var;
      ++MM_list_last;
    }
  return tmp;  
}


void *_MM_realloc(void *pointer,size_t size,char* file,size_t line)
{
  void* tmp;
  size_t i;
  char ok;

  tmp=realloc(pointer,size);
  if(!tmp) return NULL;

  ok=0;
  for(i=0;i<MM_list_last;i++)
    if((ok=(pointer==MM_list[i].adr))) break;

  if(ok) 
    {
      MM_list[i].adr=tmp;
      MM_list[i].line=line;
      MM_list[i].file=file;
    }
  else
    fprintf(stderr,"[MEMORY MANAGER] !! Warning !! Pointer not in MM_list at %s line %ld\n",file,line);
  return tmp;
}

void *_MM_valloc(size_t size,char* file,size_t line)
{
  void* tmp;
  tmp=valloc(size);
  if(!tmp) return NULL;
  if(!MM_addItem())
    {
      MM_list[MM_list_last].adr=tmp;
      MM_list[MM_list_last].line=line;
      MM_list[MM_list_last].file=file;
      ++MM_list_last;
    }
  return tmp;  
}

void MM_clean()
{
  size_t i;
  for(i=MM_list_last;i>=1;i--)
    free(MM_list[i-1].adr);
  MM_list_last=0;
  MM_deleteItems();
}
     
void MM_finish()
{
  MM_clean();
  free(MM_list);
  MM_list_size=0;
  MM_list_last=0;
  MM_list=NULL;
  return;
}
