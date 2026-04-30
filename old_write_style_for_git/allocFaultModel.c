#include "RSQSim.h"

static char svnId[] = {"$Id: allocFaultModel.c 755 2014-05-15 22:31:26Z dinger $"};

void allocFaultModel(FaultModel *m)
{
  int i;
  
  m->p = xmalloc(m->N*sizeof(FaultPatch));

  m->memorySize+=m->N*sizeof(FaultPatch);

  for (i=0; i<m->N; i++) allocFaultPatch(m, &(m->p[i]), i);
  
  return;
}

void freeFaultModel(FaultModel *m)
{
  int i;
  
  if (m->p != NULL)
  {
    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      freeKRDOkadaFaultPatch(&(m->p[i].op));
    }
    free(m->p);
  }
  m->N = 0;
  
  return;
}
  
/*
   Allocates space for the various vectors in a FaultPatch
*/
void allocFaultPatch(FaultModel *m, FaultPatch *f, int iPatch)
{
  allocKRDOkadaFaultPatch(&(f->op), m->N, ALLK, 
                          (iPatch >= m->startPatch && iPatch <= m->stopPatch ? TRUE : FALSE));

  if (iPatch >= m->startPatch && iPatch <= m->stopPatch)
    m->memorySize+=2*m->N*sizeof(double);
 
  return;
}



/*
   Frees space for the various vectors in a FaultPatch
*/
int freeFaultPatch(FaultPatch *f)
{
  freeKRDOkadaFaultPatch(&(f->op));
  
  return(0);
}
  
