/* allocate contiguous memory for 2D arrays of doubles */
#include "xmalloc.h"

double **darray2(int m, int n)
{
   double **A,*tmp;
   int i;

   tmp = (double *) xmalloc( (unsigned)(m*n*sizeof(double)) );

   A = (double **) xmalloc( (unsigned) (m*sizeof(double *)) );
   for (i = 0; i <= m-1; i++)
      A[i] = tmp + i*n;
   
   return(A);
}

void dfree2(double **A, int m, int n)
{
  if (A != (double **)NULL)
  {
    free(A[0]);
    free(A);
  }

  return;
}
