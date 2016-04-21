
#define DEBUG_BALISE

#ifdef DEBUG_BALISE

#define Balise(a) _Balise(a,__FILE__,__LINE__)

int _Balise(char* a,char* f,size_t l)
{
  printf("##[%s at %ld] %s##\n",f,l,a);fflush(stdout);
  return 0;
}
#else

#define Balise(a) /* a */

#endif
