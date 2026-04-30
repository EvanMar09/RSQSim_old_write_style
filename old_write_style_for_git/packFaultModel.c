/* functions in this file pack up most of the entries of a FaultModel struct
   to be sent around to all the subprocesses 
   
   this is called from within initFaultModel() after the master process has 
   read in all the parameters, fault geometry, and initial conditions, but 
   before the stiffness matrices are calculated (each process will calculate
   the entries of the matrix for its patches as receiving patches)
*/

#include "RSQSim.h"

static char svnId[] = {"$Id$"};

void packFaultModel(FaultModel *m, void **mbuffer, int *bufsize)
{
#ifdef USE_MPI
  int position = 0;
  int ip;  
   
  /* the below expression should be reasonable tight upper bound (or maybe need to use MPI_Pack_size()?? */
  *bufsize = sizeof(FaultModel) + m->N*(sizeof(FaultPatch) + m->maxNeighbors*sizeof(int) + sizeof(KRDOkadaFaultPatch));
  *mbuffer = xmalloc(*bufsize);

  /* pack FaultPatches */
  for (ip = 0; ip < m->N; ip++)
    packFaultPatch(&(m->p[ip]), *mbuffer, &position, *bufsize); 

  MPI_Pack(&(m->ddotEQ),                  1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->maxThetaPin),             1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->stressOvershootFactor),   1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->stressRateSpecification), 1, MPI_INT,  *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->nRupture),                1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->fA),                      1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->nTransitions),            1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->minDtWriteCoseismic),     1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->minDtWriteInterseismic),  1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->tLastWrite),              1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->minMagWrite),             1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->dMu3),                    1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->DEBUG),                   1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->ZBrentUpperBracket),      1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->ZBrentTol),               1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->state3patches),           1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->lowSigmaAction),          1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->slipInState1),            1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->resetStressingRates),     1, MPI_INT,    *mbuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(m->maxTimeStep),             1, MPI_DOUBLE, *mbuffer, *bufsize, &position, MPI_COMM_WORLD);

  if (m->DEBUG) fprintf(stderr, "bufsize = %d;  position = %d\n", *bufsize, position);   
#endif

  return;
}

void unpackFaultModel(FaultModel *m, void *mbuffer, int bufsize)
{
#ifdef USE_MPI
  int ip, position = 0;

  for (ip = 0; ip < m->N; ip++)
    unpackFaultPatch(&(m->p[ip]), mbuffer, &position, bufsize);

  MPI_Unpack(mbuffer, bufsize, &position, &(m->ddotEQ),                1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->maxThetaPin),           1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->stressOvershootFactor), 1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->stressRateSpecification), 1, MPI_INT,  MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->nRupture),              1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->fA),                    1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->nTransitions),          1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->minDtWriteCoseismic),   1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->minDtWriteInterseismic),1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->tLastWrite),            1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->minMagWrite),           1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->dMu3),                  1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->DEBUG),                 1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->ZBrentUpperBracket),    1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->ZBrentTol),             1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->state3patches),         1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->lowSigmaAction),        1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->slipInState1),          1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->resetStressingRates),   1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(mbuffer, bufsize, &position, &(m->maxTimeStep),           1, MPI_DOUBLE, MPI_COMM_WORLD);

  if (m->DEBUG) fprintf(stderr, "bufsize = %d;  position = %d\n", bufsize, position);   
#endif
  
  return;
}

void packFaultPatch(FaultPatch *p, void *buffer, int *position, int bufsize)
{
#ifdef USE_MPI
  packKRDOkadaFaultPatch(&(p->op), buffer, position, bufsize);
  MPI_Pack(&(p->tau),           1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->sigma),         1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->theta),         1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->d),             1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->ddot),          1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->mu0),           1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->ddotStar),      1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->ddotAB),        1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->ddot2A),        1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->alpha),         1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->ddotEQ),        1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->A),             1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->B),             1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->Dc),            1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->H),             1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->slowSlip),      1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->writePatch),    1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->tauDot),        1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->sigmaDot),      1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->sigmaPin),      1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->muPin),         1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->dDotDrive),     1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->tauDotDrive),   1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->sigmaDotDrive), 1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->tFail),         1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->tauFail),       1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->muTrans2),      1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->state),         1, MPI_CHAR,   buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->APFail),        1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->pinned),        1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->nNeighbors),    1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
  if (p->nNeighbors > 0) 
    MPI_Pack(&(p->neighbors[0]),  p->nNeighbors, MPI_INT, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->neighborInRuptureState), 1, MPI_INT, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->tauDotDrive),   1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(p->sigmaDotDrive), 1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
#endif
  
  return;
}

void unpackFaultPatch(FaultPatch *p, void *buffer, int *position, int bufsize)
{
#ifdef USE_MPI
  unpackKRDOkadaFaultPatch(&(p->op), buffer, position, bufsize);
  MPI_Unpack(buffer, bufsize, position, &(p->tau),           1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->sigma),         1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->theta),         1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->d),             1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->ddot),          1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->mu0),           1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->ddotStar),      1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->ddotAB),        1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->ddot2A),        1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->alpha),         1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->ddotEQ),        1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->A),             1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->B),             1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->Dc),            1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->H),             1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->slowSlip),      1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->writePatch),    1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->tauDot),        1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->sigmaDot),      1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->sigmaPin),      1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->muPin),         1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->dDotDrive),     1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->tauDotDrive),   1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->sigmaDotDrive), 1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->tFail),         1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->tauFail),       1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->muTrans2),      1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->state),         1, MPI_CHAR,   MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->APFail),        1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->pinned),        1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->nNeighbors),    1, MPI_INT,    MPI_COMM_WORLD);
  if (p->nNeighbors > 0)
  {
    p->neighbors = xmalloc(p->nNeighbors * sizeof(int));
    MPI_Unpack(buffer, bufsize, position, &(p->neighbors[0]),  p->nNeighbors, MPI_INT, MPI_COMM_WORLD);
  }
  else
    p->neighbors = NULL;
  MPI_Unpack(buffer, bufsize, position, &(p->neighborInRuptureState), 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->tauDotDrive),      1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(p->sigmaDotDrive),    1, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
  
  return;
}


void packKRDOkadaFaultPatch(KRDOkadaFaultPatch *op, void *buffer, int *position, int bufsize)
{
#ifdef USE_MPI
  MPI_Pack(&(op->corner), 12, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->center),  3, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->L),       1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->W),       1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->area),    1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->strike),  1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->dip),     1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->rake),    1, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->nu),      3, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->u),       3, MPI_DOUBLE, buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->N),       1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
  MPI_Pack(&(op->type),    1, MPI_INT,    buffer, bufsize, position, MPI_COMM_WORLD);
#endif
  
  return;
}


void unpackKRDOkadaFaultPatch(KRDOkadaFaultPatch *op, void *buffer, int *position, int bufsize)
{
#ifdef USE_MPI
  MPI_Unpack(buffer, bufsize, position, &(op->corner), 12, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->center),  3, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->L),       1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->W),       1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->area),    1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->strike),  1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->dip),     1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->rake),    1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->nu),      3, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->u),       3, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->N),       1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufsize, position, &(op->type),    1, MPI_INT,    MPI_COMM_WORLD);
#endif
  
  return;
}
