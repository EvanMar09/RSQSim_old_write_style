#include <stdlib.h>
#include <stdio.h>
#include "xmalloc.h"

void *xmalloc(size_t size)
{
   void *ptr;

   if (size == 0) return((void *)NULL);

   ptr = malloc(size);
   if (ptr == (void *)NULL)
   {
      fprintf(stderr, "error:  NULL returned from malloc (size = %zd b)\n",
              size);
      exit(-1);
   }
   return(ptr);
}   

