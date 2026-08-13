#include "RSQSim.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif


static char svnId[] = {"$Id: transitions.c 1065 2017-04-07 21:57:01Z dinger $"};

/*

 In this particular file are the functions that update the system
 when a patch has gone through a state transition

 makeTransition() is a wrapper function which calls all three of the
 individual transitionSNtoSN() functions.


*/

void makeTransition(FaultModel *m, Transition tr)
{
  if (tr.patch < 0 && tr.nextState != MAX_TIME_STEP)
  {
    if (m->master)
      fprintf(stderr, "makeTransition(): Warning: tr.patch < 0; doing nothing\n");

    return;
  }


  if (tr.nextState == LOW_SIGMA || tr.nextState == LOW_TAU || tr.nextState == HIGH_THETA ||
      tr.nextState == HIGH_TAU)
  {
    transitionLowSigma(m, tr.patch, tr.nextState);
  }
  else
  {
    if (tr.patch >= m->startPatch && tr.patch <= m->stopPatch && m->p[tr.patch].APFail == NOW)
    {
      aPrioriFail(m, tr.patch);
    }

    if (tr.nextState == LOW_SIGMA || tr.nextState == LOW_TAU || tr.nextState == HIGH_THETA ||
        tr.nextState == HIGH_TAU)
    {
      transitionLowSigma(m, tr.patch, tr.nextState);
    }
    else if (tr.nextState == MAX_TIME_STEP) transitionMaxTimeStep(m);
    else if (tr.nextState == STRESS_RATE_STEP) transitionStressRateStep(m);
    else
    {
      switch(m->p[tr.patch].state)
      {
        case 0:
          if (tr.nextState == 1) transition0to1(m, tr.patch);
          else if (tr.nextState == 2 && (tr.patch < m->startPatch || tr.patch > m->stopPatch))
          {   /* then this must be another process's patch which did 0 -> 1 in timeStep() */
            m->p[tr.patch].state = 1;
            transition1to2(m, tr.patch);
          }
          else if (tr.nextState == SLOWSLIP_2A && (tr.patch < m->startPatch || tr.patch > m->stopPatch))
          {   /* then this must be another process's patch which did 0 -> 1 in timeStep() */
            m->p[tr.patch].state = 1;
            transition1to2a(m, tr.patch);
          }
          else if (tr.nextState == 0 && (tr.patch < m->startPatch || tr.patch > m->stopPatch))
          {   /* then this is likely another process's patch which did a 0 -> 1 in timeStep() and now is going back to 0 */
            m->p[tr.patch].state = 0;  /* don't really need to do anything here as this is already 0 */
          } else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, 0, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        case 1:
          if (tr.nextState == 0) transition1to0(m, tr.patch);
          else if (tr.nextState == 2) transition1to2(m, tr.patch);
          else if (tr.nextState == SLOWSLIP_2A) transition1to2a(m, tr.patch);
          else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, 1, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        case 2:
          if (tr.nextState == 0) transition2to0(m, tr.patch);
          else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, 2, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        case 3:
          if (tr.nextState == 3) transition3to3(m, tr.patch);
          else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, 3, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        case SLOWSLIP_2A:
          if (tr.nextState == 0)
            transition2ato0(m, tr.patch);
          else if (tr.nextState == SLOWSLIP_2B)
            transition2ato2b(m, tr.patch);
          else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, SLOWSLIP_2A, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        case SLOWSLIP_2B:
          if (tr.nextState == SLOWSLIP_2A)
            transition2bto2a(m, tr.patch);
          else if (tr.nextState == SLOWSLIP_2B)
            transition2bto2b(m, tr.patch);
          else if (tr.nextState == SLOWSLIP_2C)
            transition2bto2c(m, tr.patch);
          else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, SLOWSLIP_2B, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        case SLOWSLIP_2C:
          if (tr.nextState == LOCK)
            transition2cto0(m, tr.patch);
          else
          {
            fprintf(stderr, "makeTransition(): error: on patch %d, state %d -> state %d\n", tr.patch, SLOWSLIP_2C, tr.nextState);
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            exit(-1);
          }
          break;

        default:
          fprintf(stderr, "makeTransition():  error: m.p[%d].state = %d\n",
                  tr.patch, m->p[tr.patch].state);
#ifdef USE_MPI
          MPI_Abort(MPI_COMM_WORLD, -1);
#endif
          exit(-1);
      }
    }
  }

  m->nTransitions++;

  return;
}


void aPrioriFail(FaultModel *m, int patch)
{
  FaultPatch *p;

  p = &(m->p[patch]); /* for brevity below */

  if (p->tauFail > 0) /* set shear stress to tauFail before setting theta for immediate failure */
    p->tau = p->tauFail;

  p->theta = (p->Dc/p->ddotStar) * pow(p->ddotEQ/p->ddotStar, -p->A/p->B)*exp((p->tau/p->sigma - p->mu0)/p->B);

  p->ddot = p->ddotEQ; /* I think this will get overwritten if p is in state 0,
                          but that it is needed if p is in state 1 */


  p->APFail = ALREADY;  /* indicate that this patch has already had its shear
                           stress set at its a priori failure time */
  return;
}

void transition0to1(FaultModel *m, int patch)
{
  FaultPatch *p;

  p = &(m->p[patch]); /* for brevity below */

  /* if not one of my patches just update the state */
  if (patch < m->startPatch || patch > m->stopPatch)
  {
    p->state = 1;
    return;
  }

  /* only the transitioning element has any changes */
  m->p[patch].state = 1;
/*  m->p[patch].ddot = (p->Dc/p->theta)*(1 - (p->alpha/p->B)*(p->sigmaDot/p->sigma));
   now using the expression below which should give the same answer for a patch which
   actually transitioned 0 -> 1 and should also be right if it immediately transitions
   upon model startup. Note that these expressions for ddot have to be consistent with
   what criteria are used for the 0->1 transition time, for example with the currently
   used approximation of tau == sigma*(mu0 + (B-A)*ln(theta)), the above should really
   just be (Dc/theta) */
  p->ddot = p->ddotStar*exp((p->tau - p->mu0*p->sigma)/(p->A*p->sigma))*
                     pow(p->theta*p->ddotStar/p->Dc, -p->B/p->A);
  p->ddot = MIN(m->p[patch].ddot, p->ddotEQ); /* may need to transition straight to state 2 */

  if (m->slipInState1 != 0)
    p->H = p->B/p->Dc + (p->Kself - (p->tau/p->sigma)*p->op.Ksigma[patch])/p->sigma;
  else
    p->H = p->B/p->Dc;

  checkA(m, patch);  /* will reduce A if appropriate */

  p->t01 = T01_NA_VAL;

  return;
}

void transition1to0(FaultModel *m, int patch)
{
  FaultPatch *p;

  p = &(m->p[patch]); /* for brevity below */

  /* if not one of our patches, just update the state */
  if (patch < m->startPatch || patch > m->stopPatch)
  {
    p->state = 0;
    return;
  }

  p->state = 0;
/*  if (p->neighborInRuptureState > 0) p->A /= m->fA; this is not allowed to happen anymore */
  checkA(m, patch); /* may allow this to happen again */

  if (m->DEBUG > 2) fprintf(stdout, "patch %d transitioning from 1 to 0\n", patch);

}

void transition1to2(FaultModel *m, int patch)
{
  int i, iNeighbor,ipatch;
  FaultPatch *p, *ptrans;

  m->nRupture++;

  p = &(m->p[patch]);
  p->state = 2;
  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    /* first make changes to transitioning element */
    p->tauDot += p->ddotEQ*p->Kself; /* now could just let the general tauDot loop take care of this */
    if (!p->sigmaClipped) p->sigmaDot += p->ddotEQ*p->op.Ksigma[patch];
    p->ddot = p->ddotEQ;
    p->theta = p->Dc/p->ddotEQ; /* actually, will quickly evolve to about this */

    checkA(m, patch); /* will unreduce A if necessary */

    p->muTrans2 = (p->mu0 + (p->A - p->B)*log(p->ddotEQ/p->ddotStar))*(1+m->stressOvershootFactor) -
                   m->stressOvershootFactor*p->tau/p->sigma;
  } /* if  (patch >= m->startPatch && patch <= m->stopPatch) */
  else /* if not one of our patches, just make note of the state change for the transitioning patch */
  {
    p->state = 2;
  }

  ptrans = p; /* the 1->2 transitioning patch */
  /* only tauDot (and possibly A) changes on the other elements */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,ptrans) private(i,p,iNeighbor)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    if (i == patch) continue;

    p->tauDot += ptrans->ddotEQ * p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += ptrans->ddotEQ * p->op.Ksigma[patch];

    if (p->neighborInRuptureState == 0) /* then check to see if the 1->2 patch            */
    {                             /* is one of this patch's neighbors, for A reduction purposes; now done below instead */
      if (p->nNeighbors < 0) /* then all patches are neighbors of this one */
      {
        p->neighborInRuptureState = 1;
      }
    }

    checkA(m, i);  /* will reduce A if appropriate */

  } /* for (i = m->startPatch; i <= m->stopPatch; i++) */

  for(iNeighbor=0; iNeighbor < ptrans->nNeighbors; iNeighbor++)
  {
    ipatch = ptrans->neighbors[iNeighbor];
    if(ipatch >= m->startPatch && ipatch <= m->stopPatch)
    {
      p = &(m->p[ipatch]);
      p->neighborInRuptureState = 1;
      checkA(m, ipatch);
    }
  }

  return;
}

void transition1to2a(FaultModel *m, int patch)
{
  double DeltaDdot;
  int i, iNeighbor;
  FaultPatch *p; /* for brevity */
  FaultPatch *ptrans;

  m->nRupture++;

  p = &(m->p[patch]);
  p->state = SLOWSLIP_2A;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    checkA(m, patch); /* will unreduce A if necessary */
    p->ddot = p->ddot2A;
    p->theta = p->Dc/p->ddot;
    DeltaDdot = p->ddot;
    p->muTrans2aLow = p->mu0*(1+m->stressOvershootFactor) - m->stressOvershootFactor * p->tau/p->sigma;
    p->muTrans2High = (p->mu0 + p->A*log(p->ddotAB/p->ddotStar))*(1.0 + MU_EPS);

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
    p->state = SLOWSLIP_2A;
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  ptrans = p;

#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,iNeighbor,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);
    p->tauDot += DeltaDdot*p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot*p->op.Ksigma[patch];

    p->t01 = T01_NA_VAL;

    if (i == patch) continue;

    if (p->neighborInRuptureState == 0) /* then check to see if the 1->2 patch            */
    {                                  /* is one of this patch's neighbors, for A reduction purposes */
      if (p->nNeighbors < 0) /* then all patches are neighbors of this one */
      {
        p->neighborInRuptureState = 1;
      }
      else /* need to look through list of neighbors: note possibly more efficient way */
      {    /* of doing this above in transition1to2() */
        for (iNeighbor=0; iNeighbor < p->nNeighbors; iNeighbor++)
          if (p->neighbors[iNeighbor] == patch) break;

        if (iNeighbor < p->nNeighbors) /* then the transitioning patch is a neighbor */
        {
          p->neighborInRuptureState = 1;
        }
      }
    }

    checkA(m, i);  /* will reduce A if appropriate */


  }

  return;
}


void transition2to0(FaultModel *m, int patch)
{
  int i;
  FaultPatch *p, *ptrans;

  m->nRupture--;

  p = &(m->p[patch]);
  p->state = 0;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    /* first make changes to transitioning element */
    p->tauDot -= p->ddotEQ * p->Kself;
    if (!p->sigmaClipped) p->sigmaDot -= p->ddotEQ * p->op.Ksigma[patch];
    p->ddot = 0;

    p->theta = 2.0/(p->ddotEQ/p->Dc + (p->alpha/p->B)*(p->sigmaDot/p->sigma));
    /* probably could just ignore (alpha/B)*(sigmaDot/sigma) relative to ddotEQ/Dc.
       in sigmaDot == 0 case, needs to be slightly bigger than Dc/ddotEQ,
       now with the new parameter stressOvershootFactor, probably don't need that 2 in the numerator
       Turns out do need the 2 in the numerator for using large D_c, small patches, and slipInState1 = 0 */
  } /* if (patch >= m->startPatch && patch <= m->stopPatch) */
  else
  {
    p->state = 0;
  }
  ptrans = p;

  /* only tauDot (and possibly A) changes on the other elements; unless an event just stopped,
     then see next section below */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,ptrans) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    if (i == patch) continue;

    p->tauDot -= ptrans->ddotEQ * p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot -= ptrans->ddotEQ * p->op.Ksigma[patch];
    /* if an eq just ended, reset A for state 1 patches back to original and adjust theta  */
    if (m->nRupture == 0)
      p->neighborInRuptureState = 0; /* this will get set to zero again in updateEqs() */

    checkA(m, i);
  }


  return;
}

void transition2ato0(FaultModel *m, int patch)
{
  int i;
  FaultPatch *p;
  double DeltaDdot;

  m->nRupture--;

  p = &(m->p[patch]);
  p->state = 0;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    /* first make changes to transitioning element */
    p->tauDot -= p->ddot * p->Kself;
    if (!p->sigmaClipped) p->sigmaDot -= p->ddot * p->op.Ksigma[patch];
    DeltaDdot = -p->ddot;
    p->ddot = 0;

    p->theta = p->Dc/p->ddotStar;

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  } /* if (patch >= m->startPatch && patch <= m->stopPatch) */
  else
  {
    p->state = 0;
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  /* only tauDot (and possibly A) changes on the other elements; unless an event just stopped,
     then see next section below */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    if (i == patch) continue;

    p->tauDot += DeltaDdot * p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot * p->op.Ksigma[patch];
    /* if an eq just ended, reset A for state 1 patches back to original and adjust theta  */
    if (m->nRupture == 0)
      p->neighborInRuptureState = 0; /* this will get set to zero again in updateEqs() */

    checkA(m, i);
  }


  return;
}



void transition2ato2b(FaultModel *m, int patch)
{
  double dDotOld, DeltaDdot;
  int i;
  FaultPatch *p; /* for brevity */

  p = &(m->p[patch]);
  p->state = SLOWSLIP_2B;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    dDotOld = p->ddot;
    setDdot2b(&(m->p[patch]));
    p->theta = p->Dc/p->ddot;
    DeltaDdot = p->ddot - dDotOld;
    p->muTrans2Low = MAX(p->tau/p->sigma - m->dMu3, (p->mu0 + p->A*log(p->ddotAB/p->ddotStar))*(1.0 - MU_EPS));
    p->muTrans2High = MIN(p->tau/p->sigma + m->dMu3, (p->mu0 + p->A*log(p->ddotEQ/p->ddotStar))*(1.0 + MU_EPS));
    p->muTrans2aLow = MIN(p->muTrans2aLow, p->mu0 - m->stressOvershootFactor * (p->tau/p->sigma - p->mu0));

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    p->tauDot += DeltaDdot*p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot*p->op.Ksigma[patch];

    checkA(m, i);
  }

  return;
}

void transition2bto2b(FaultModel *m, int patch)
{
  double dDotOld, DeltaDdot;
  int i;
  FaultPatch *p; /* for brevity */

  p = &(m->p[patch]);

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    dDotOld = p->ddot;
    setDdot2b(&(m->p[patch]));
    p->theta = p->Dc/p->ddot;
    DeltaDdot = p->ddot - dDotOld;
    p->muTrans2Low = MAX(p->tau/p->sigma - m->dMu3, (p->mu0 + p->A*log(p->ddotAB/p->ddotStar))*(1.0 - MU_EPS));
    p->muTrans2High = MIN(p->tau/p->sigma + m->dMu3, (p->mu0 + p->A*log(p->ddotEQ/p->ddotStar))*(1.0 + MU_EPS));
    p->muTrans2aLow = MIN(p->muTrans2aLow, p->mu0 - m->stressOvershootFactor * (p->tau/p->sigma - p->mu0));

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    p->tauDot += DeltaDdot*p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot*p->op.Ksigma[patch];

    checkA(m, i);
  }

  return;
}




void transition2bto2a(FaultModel *m, int patch)
{
  double dDotOld, DeltaDdot;
  int i;
  FaultPatch *p; /* for brevity */

  p = &(m->p[patch]);
  p->state = SLOWSLIP_2A;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    dDotOld = p->ddot;
    p->ddot = p->ddotStar;
    p->theta = p->Dc/p->ddot;
    DeltaDdot = p->ddot - dDotOld;
    p->muTrans2aLow = MIN(p->muTrans2aLow, p->mu0 - m->stressOvershootFactor * (p->tau/p->sigma - p->mu0));
    p->muTrans2High = (p->mu0 + p->A*log(p->ddotAB/p->ddotStar))*(1.0 + MU_EPS);

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    p->tauDot += DeltaDdot*p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot*p->op.Ksigma[patch];

    checkA(m, i);
  }

  return;
}



void transition2bto2c(FaultModel *m, int patch)
{
  double dDotOld, DeltaDdot;
  int i;
  FaultPatch *p; /* for brevity */

  p = &(m->p[patch]);
  p->state = SLOWSLIP_2C;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    dDotOld = p->ddot;
    p->ddot = p->ddotEQ;
    p->theta = p->Dc/p->ddot;
    DeltaDdot = p->ddot - dDotOld;
    p->muTrans2aLow = MIN(p->muTrans2aLow, p->mu0 - m->stressOvershootFactor * (p->tau/p->sigma - p->mu0));
    p->muTrans2Low = (p->mu0 + p->A*log(p->ddotEQ/p->ddotStar))*(1.0 - MU_EPS);
    p->muTrans2High = INFINITY;

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    p->tauDot += DeltaDdot*p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot*p->op.Ksigma[patch];

    checkA(m, i);
  }

  return;
}

void transition2cto0(FaultModel *m, int patch)
{
  double dDotOld, DeltaDdot;
  int i;
  FaultPatch *p; /* for brevity */

  m->nRupture--;

  p = &(m->p[patch]);
  p->state = LOCK;

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    /* first make changes to transitioning element */
    p->tauDot -= p->ddot * p->Kself;
    if (!p->sigmaClipped) p->sigmaDot -= p->ddot * p->op.Ksigma[patch];
    p->theta = p->Dc/p->ddot;
    DeltaDdot = -p->ddot;
    p->ddot = 0;

    p->muTrans2Low = NAN;
    p->muTrans2High = NAN;
    p->muTrans2aLow = NAN;

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  /* only tauDot (and possibly A) changes on the other elements; unless an event just stopped,
     then see next section below */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif
  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    if (i == patch) continue;

    p->tauDot += DeltaDdot * p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot * p->op.Ksigma[patch];
    /* if an eq just ended, reset A for state 1 patches back to original and adjust theta  */
    if (m->nRupture == 0)
      p->neighborInRuptureState = 0; /* this will get set to zero again in updateEqs() */

    checkA(m, i);
  }

  return;
}



void transition3to3(FaultModel *m, int patch)
{
  double dDotOld, DeltaDdot;
  int i;
  FaultPatch *p; /* for brevity */

  p = &(m->p[patch]);

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    dDotOld = p->ddot;
    setDdot3(&(m->p[patch]));
    p->theta = p->Dc/p->ddot;
    DeltaDdot = p->ddot - dDotOld;
    p->muTrans3Low = p->tau/p->sigma - m->dMu3;
    p->muTrans3High = p->tau/p->sigma + m->dMu3;

#ifdef USE_MPI
    if (!m->master) MPI_Send(&DeltaDdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
#endif
  }
  else
  {
#ifdef USE_MPI
    if (m->master) MPI_Recv(&DeltaDdot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
  }

#ifdef USE_MPI
  MPI_Bcast(&DeltaDdot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch,DeltaDdot) private(i,p)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    p = &(m->p[i]);

    p->t01 = T01_NA_VAL;

    p->tauDot += DeltaDdot*p->op.Ktau[patch];
    if (!p->sigmaClipped) p->sigmaDot += DeltaDdot*p->op.Ksigma[patch];

    checkA(m, i);
  }

  return;
}

/* do nothing for now, later will update stresses and/or stressing rates
   from external sources */
void transitionMaxTimeStep(FaultModel *m)
{
  return;
}

void transitionStressRateStep(FaultModel *m)
{
  int i;

  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    m->p[i].tauDot += m->p[i].tauDotExt[m->iExtStressRateStep];
    m->p[i].sigmaDot += m->p[i].sigmaDotExt[m->iExtStressRateStep];
    if (m->iExtStressRateStep > 0)
    {
      m->p[i].tauDot -= m->p[i].tauDotExt[m->iExtStressRateStep - 1];
      m->p[i].sigmaDot -= m->p[i].sigmaDotExt[m->iExtStressRateStep - 1];
    }

    checkA(m, i);
  }

  m->iExtStressRateStep++;

  return;
}


/* for now, we are just locking such patches permanently */
void transitionLowSigma(FaultModel *m, int patch, int nextState)
{
  int i;

  if (m->lowSigmaAction == CLIP && nextState == LOW_SIGMA)
  {
    m->p[patch].sigmaClipped = 1;
    m->p[patch].sigmaDot = 0.0;
    return;
  }

  if (m->p[patch].state == 2)
    transition2to0(m, patch);
  else
    m->p[patch].state = 0;

  m->p[patch].pinned = 1;

  /* adjust backslip stress rates for the fact that this patch is now locked */
  if (m->stressRateSpecification == BACKSLIP)
  {
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,patch) private(i)
#endif
    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      m->p[i].tauDotDrive += m->p[i].op.Ktau[patch]*m->p[patch].dDotDrive;
      if (!m->p[i].sigmaClipped) m->p[i].sigmaDotDrive += m->p[i].op.Ksigma[patch]*m->p[patch].dDotDrive;
    }
  }

  if (m->master)
  {
    fprintf(stderr, "transitionLowSigma(): Eliminating element %d due to ", patch+1);
    switch(nextState)
    {
      case LOW_SIGMA:
        fprintf(stderr, "low sigma\n");
        break;
      case LOW_TAU:
        fprintf(stderr, "low tau\n");
        break;
      case HIGH_THETA:
        fprintf(stderr, "high theta\n");
        break;
      case HIGH_TAU:
        fprintf(stderr, "high tau\n");
        break;
    }
  }

  if (patch >= m->startPatch && patch <= m->stopPatch)
  {
    fprintf(stderr, "tau = %lg;  sigma = %lg;  theta = %lg\n",
            m->p[patch].tau, m->p[patch].sigma, m->p[patch].theta);
  }

  return;
}

void checkA(FaultModel *m, int patch)
{
  double SDotExt; /* modified Coulomb stress */

  FaultPatch *p; /* for brevity below */

  p = &(m->p[patch]);

  SDotExt = p->tauDot - ((p->tau/p->sigma) - p->alpha)*p->sigmaDot;
  if (p->state == 1 && SDotExt > 0 && p->neighborInRuptureState > 0)
  {
    if (!(p->aReduced)) reduceA(m, patch);
  }
  else if (p->aReduced)
    unreduceA(m, patch);

  return;
}

void reduceA(FaultModel *m, int patch)
{
  if (m->aReduceN2 > 0) /* don't use aReduceN1, aReduceN2, or fAmin if aReduceN2 <= 0 */
  {
    double fA;

    if (m->nRupture < m->aReduceN1)
      fA = m->fAmin;
    else if (m->nRupture <= m->aReduceN2)
      fA = m->fAmin + fAfunc((m->nRupture - m->aReduceN1)/
                             (m->aReduceN2 - m->aReduceN1)) * (1 - m->fAmin);
    else
      fA = 1;

    m->p[patch].A = m->p[patch].Asv * fA;
  }
/*
  if (m->aReduceNmax > 0)
  {
    if (m->nRupture < m->aReduceNmax)
      m->p[patch].A = m->p[patch].Asv * sqrt(m->nRupture/m->aReduceNmax);
    else
      m->p[patch].A = m->p[patch].Asv;
  }
*/
  else  m->p[patch].A = m->p[patch].Asv * m->fA;

  m->p[patch].aReduced = 1;

  if (m->DEBUG > 7) fprintf(stderr, "reducing A on patch %d\n", patch);

  return;
}

void unreduceA(FaultModel *m, int patch)
{

  FaultPatch *p;

  p = &(m->p[patch]); /* for brevity below */

  p->A = p->Asv;
  p->aReduced = 0;
  if (p->state == 1)
  {
    if (m->unreduceAFix == SLIPSPEED)
      p->theta = theta_1(p->tau, p->sigma, p->ddot, p->mu0, p->ddotStar, p->A, p->B, p->Dc);
    else
      p->ddot = ddotRS(p->ddotStar, p->tau, p->sigma, p->mu0, p->A, p->B, p->theta, p->Dc);
  }

  if (m->DEBUG > 7) fprintf(stderr, "unreducing A on patch %d\n", patch);

  return;
}
