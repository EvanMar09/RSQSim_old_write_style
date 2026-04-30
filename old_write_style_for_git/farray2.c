/* allocate contiguous memory for 2D arrays of floats */
#include "xmalloc.h"

float **farray2(int m, int n)
{
   float **A,*tmp;
   int i;

   tmp = (float *) xmalloc( (unsigned)(m*n*sizeof(float)) );

   A = (float **) xmalloc( (unsigned) (m*sizeof(float *)) );
   for (i = 0; i <= m-1; i++)
      A[i] = tmp + i*n;
   
   return(A);
}

void ffree2(float **A, int m, int n)
{
  int i;

  if (A != (float **)NULL)
  {
    if (m!=0) free(A[0]);
    free(A);
  }

  return;
}
