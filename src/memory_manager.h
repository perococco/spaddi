/*******************************/
/*  Spectra Addition Procedure */
/*     Bastien Aracil 2002     */
/*       aracil@iap.fr         */
/*******************************/   
#ifndef MEMORYMANAGERH
#define MEMORYMANAGERH

#define MM_list_stride 200

#define MM_calloc(nelt,elts) _MM_calloc(nelt,elts,__FILE__,__LINE__)
#define MM_free(ptr)         _MM_free(ptr,__FILE__,__LINE__)
#define MM_malloc(size)      _MM_malloc(size,__FILE__,__LINE__)
#define MM_realloc(ptr,size) _MM_realloc(ptr,size,__FILE__,__LINE__)
#define MM_valloc(size)      _MM_valloc(size,__FILE__,__LINE__)
#define MMV_malloc(var,size,type)    var=(type*)_MMV_malloc((size)*sizeof(type),#var,__FILE__,__LINE__)


typedef struct {
  void* adr;
  size_t line;
  char* file;
  char* varname;
} MM_memItem;

extern MM_memItem* MM_list;
extern size_t      MM_list_size;
extern size_t      MM_list_last;

int MM_deleteItems();
int MM_AddItem();
void MM_info();

MM_memItem* MM_init();

void *_MM_calloc(size_t num_of_elts,size_t elt_size,char* file,size_t line);
void  _MM_free(void *pointer,char* file,size_t line);
void *_MM_malloc(size_t size,char* file,size_t line);
void *_MMV_malloc(size_t size,char* var,char* file,size_t line);
void *_MM_realloc(void *pointer,size_t size,char* file,size_t line);
void *_MM_valloc(size_t size,char* file,size_t line);

void  MM_clean();
void  MM_finish();

#endif
