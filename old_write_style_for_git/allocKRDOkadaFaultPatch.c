#include "KRDOkada.h"

/*
   Allocates space for the various vectors in a KRDOkadaFaultPatch
*/
int allocKRDOkadaFaultPatch(KRDOkadaFaultPatch *f, int N, int allK, int anyK)
{

  if (anyK)
  {
    if (allK != 0) f->K = darray2(N, 6);
    else f->K = NULL;
    f->Ksigma = xmalloc(N*sizeof(double));
    f->Ktau = xmalloc(N*sizeof(double));
  }
  else
  {
    f->K = NULL;
    f->Ksigma = NULL;
    f->Ktau = NULL;
  }
  f->N = N;
  
  return(0);
}



/*
   Frees space for the various vectors in a KRDOkadaFaultPatch
*/
int freeKRDOkadaFaultPatch(KRDOkadaFaultPatch *f)
{
  if (f->K != NULL) dfree2(f->K, f->N, 6);
  if (f->Ksigma != NULL) free(f->Ksigma);
  if (f->Ktau != NULL) free(f->Ktau);
  f->N = 0;
  
  return(0);
}
