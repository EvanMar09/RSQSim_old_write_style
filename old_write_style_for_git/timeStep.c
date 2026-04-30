#include "RSQSim.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif


static char svnId[] = {"$Id: timeStep.c 809 2014-08-20 23:23:31Z dinger $"};

/*

   In this particular file are the functions that update the system
   through a time period in which no transitions take place

   timeStep() is a wrapper function which calls all three of the
   individual updateStateN() functions.

   The input argument p for all these functions is a list that
   needs to have at least the following elements:
    tau, tauDot, theta, state, A, B, ddot, H, Kself, d
   each of which should be a vector that has the corresponding
   quantity for each patch


*/

void timeStep(FaultModel *m, double dt)
{
  int i;
  
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,dt) private(i)
#endif

  for (i = m->startPatch; i <= m->stopPatch; i++)
  { 
    switch(m->p[i].state) 
    {
      case 0: 
        if (m->p[i].t01 != T01_NA_VAL && dt >= m->p[i].t01)
        {
          double t01sv;
          t01sv = m->p[i].t01;
          updateState0(m, i, m->p[i].t01);
          transition0to1(m, i);
          updateState1(m, i, dt - t01sv );
        } else
	  updateState0(m, i, dt);
	break;
      case 1:
	updateState1(m, i, dt);
	break;
      case 2:
	updateState2(m, i, dt);
	break;
      case 3:
        updateState3(m, i, dt);
        break; 
      case SLOWSLIP_2A:
      case SLOWSLIP_2B:
      case SLOWSLIP_2C:
        updateState2x(m, i, dt);
        break;
      default:
	fprintf(stderr, "timeStep():  error: m.p[%d].state = %d\n",
        	i, m->p[i].state);
    }
  }

  return;
}

void updateState0(FaultModel *m, int patch, double dt)
{
  FaultPatch *p;
  
  p = &(m->p[patch]);
  
  p->theta = theta_0(dt, p->sigmaDot, p->sigma, p->alpha, p->B, p->theta);
  p->tau += p->tauDot*dt;
  p->sigma += p->sigmaDot*dt;
  
  if (p->t01 != T01_NA_VAL) p->t01 -= dt; 
    
  return;
}

void updateState1(FaultModel *m, int patch, double dt)
{
  double dd; /* an increment of displacement */
  double SDotExt; /* modified Coulomb stressing rate */
  double ddotsv; /* for debugging purposes */
  double A; /* may need to use p->A or original A for updating theta depending whether a
               rupture is in progress */
  FaultPatch *p;
  
  p = &(m->p[patch]);
  
  
  SDotExt = p->tauDot - ((p->tau/p->sigma) - p->alpha)*p->sigmaDot;

  if (m->slipInState1 != 0)
  {
    p->H = p->B/p->Dc + (p->Kself - (p->tau/p->sigma - p->alpha)*p->op.Ksigma[patch])/p->sigma;
    dd = d1(dt, p->A, p->H, p->ddot, p->sigma, SDotExt);
  }
  else
  {
    p->H = p->B/p->Dc;
    dd = 0;
  }
  ddotsv = p->ddot;
  p->ddot = ddot1(dt, p->ddot, SDotExt, p->A, p->sigma, p->H, p->ddotEQ);
  if (ddotsv != 0 && p->ddot == 0.0) 
  {
    fprintf(stderr, "Warning: updateState1(); ddot is zero on patch %d\n", patch);
  }
  p->tau += p->tauDot*dt + p->Kself*dd;
  p->sigma += p->sigmaDot*dt + p->op.Ksigma[patch]*dd;
  p->d += dd;
  if (p->aReduced > 0) A = p->A/m->fA; /* use original A for updating theta */
  else A = p->A;
/*  p->theta = theta_1(p->tau, p->sigma, p->ddot, p->mu0, p->ddotStar, A, p->B, p->Dc); */
  p->theta = theta_1a(p->theta, dd, p->Dc, p->alpha, p->B, 
                      p->sigma - p->sigmaDot*dt - p->op.Ksigma[patch]*dd, p->sigma);
  
  return;
}

void updateState2(FaultModel *m, int patch, double dt)
{
  FaultPatch *p;
  
  p = &(m->p[patch]);
  
  p->tau += p->tauDot*dt;
  p->sigma += p->sigmaDot*dt;
  p->d += p->ddot*dt;
  
  return;
}
  
  
/* note that, for now, we are only updating the stresses, not the
   slip speeds, so the stresses and slip speeds will be inconsistent
   with each other until we update the slip speeds */  
void updateState3(FaultModel *m, int patch, double dt)
{
  FaultPatch *p;
  
  p = &(m->p[patch]);
  
  p->tau += p->tauDot*dt;
  p->sigma += p->sigmaDot*dt;
  p->d += p->ddot*dt;
    
  return;
}
    
void updateState2x(FaultModel *m, int patch, double dt)
{
  FaultPatch *p;
  
  p = &(m->p[patch]);
  
  p->tau += p->tauDot*dt;
  p->sigma += p->sigmaDot*dt;
  p->d += p->ddot*dt;
  
  
  return;
}
    
