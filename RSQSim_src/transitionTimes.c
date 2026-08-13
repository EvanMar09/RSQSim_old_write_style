#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "RSQSim.h"


static char svnId[] = {"$Id: transitionTimes.c 1083 2017-04-11 23:29:47Z dinger $"};
static void checkIfMinT(Transition *T, Transition *minT, FaultModel *m, int i);

/*

   In this particular file are the functions that calculate the
   transition times.  The SNtoSN() take all the patches but then
   just work on the ones in their particular starting state.
   They all return a list with two elements: dt is the minimum
   transition time, and patch is the index of the patch for
   which it occurred.

   The input argument p for all these functions is a list that
   needs to have at least the following elements:
    tau, tauDot, theta, state, A, B, ddot, H
   each of which should be a vector that has the corresponding
   quantity for each patch

   findMinTransitionTime() is just a wrapper function that calls
   all three of the individual SNtoSN() functions and returns
   the overall minimum dt and the patch for which it occurred.

   What happens with ties?  Hopefully the first one found will
   get returned and the next time in, another one of the ties
   will be found

   Will any of these ever give negative dt's?
*/

Transition findMinTransitionTime(FaultModel *m)
{
/*
 *
 */
#ifdef USE_OPENMP
  int omp_num_procs=omp_get_num_procs();
#else
  int omp_num_procs=1;
#endif
/*
 *
 */

  int i,tid;
  double thetasv, tausv, sigmasv, t01sv;

#ifdef USE_OPENMP
/* ---llnl openMP
  Fork a team of threads with each thread having a private tid variable */
#pragma omp parallel private(tid)
  {

  /* Obtain and print thread id */
  tid = omp_get_thread_num();

  /* Only master thread does this */
  if (tid == 0)
    {
    omp_num_procs = omp_get_num_threads();
//         fprintf(stdout,"Number of threads per MPI (OPENMP) %d\n",omp_num_procs);
    }

  }  /* All threads join master thread and terminate */

#endif

  Transition minT[omp_num_procs], T, T1;

/*
 *
 */  
 for (i=0;i<omp_num_procs;i++) {
          minT[i].dt = INFINITY;/* smallest dt found so far */
          minT[i].patch = -1;   /* patch for which dt == minT.dt */
          minT[i].nextState = 0;
  }
#ifdef USE_OPENMP
#pragma omp parallel for shared(minT,m,omp_num_procs) private(T,T1,i,thetasv,tausv,sigmasv,t01sv)
#endif
/*
 *
 */

//  minT.dt = INFINITY;/* smallest dt found so far */
//  minT.patch = -1;   /* patch for which dt == minT.dt */
  
  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
/*
 *
 */
#ifdef USE_OPENMP
        int omp_rank=omp_get_thread_num();
  //int mpi_rank;
  //  //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  //    //printf("mpi rank=%d, openmp rank=%d\n", mpi_rank, omp_rank);
#else
  int omp_rank=0;
#endif
/*
 *
 */


    if (m->p[i].pinned) continue; /* if patch is pinned, never transitions out of its initial state of 0 */
    
    T = checkForLowStress(&(m->p[i]));
    checkIfMinT(&T, minT+omp_rank, m, i);
    
    T = checkForHighTheta(&(m->p[i]), m->t, m->maxThetaPin);
    checkIfMinT(&T, minT+omp_rank, m, i);
    
    switch(m->p[i].state)
    {
      case 0:
        T = S0TransTime(m, i, minT[omp_rank].dt);
        if (isfinite(T.dt) && m->p[i].APFail != NOW)
        {
          thetasv = m->p[i].theta;
          tausv = m->p[i].tau;
          sigmasv = m->p[i].sigma;
          t01sv = m->p[i].t01;
          updateState0(m, i, T.dt);
          transition0to1(m, i);
          T1 = S1TransTime(m, i);
          m->p[i].theta = thetasv;
          m->p[i].tau = tausv;
          m->p[i].sigma = sigmasv;
          m->p[i].t01 = t01sv;
          m->p[i].state = 0;
          m->p[i].ddot = 0;
          checkA(m, i);
          if (T1.dt <= 0.0)
	  {
	    static int count = 0; /* hopefully only initialized once per run */
	    
	    count++;
	    
	    fprintf(stderr, "Warning: T1.dt <= 0 for patch %d\n", i);
	    
	    if (count >= m->maxDt1Zero) 
	    {
	      FILE *tmpfp;
	      
	      fprintf(stderr, "Error on process %d: T1.dt <= 0 too many times, stopping\n", 
	              m->iProc);
              tmpfp = fopen("stopnow", "a");
	      fclose(tmpfp);
	    }
	    
	  }
          T.dt += T1.dt;
          T.nextState = T1.nextState;
          if (T1.nextState == 0)
          {
            fprintf(stderr, "Warning: 1 -> 0 -> 1 on patch %d during time of constant stressing rates\n", i);
	    if (m->DEBUG > 2)
	      fprintf(stderr, "   tau = %lg;  tauDot = %lg; sigma = %lg;  sigmaDot = %lg; theta = %lg; slip = %lg; slipSpeed = %lg\n",
	              m->p[i].tau, m->p[i].tauDot, m->p[i].sigma, m->p[i].sigmaDot, m->p[i].theta, m->p[i].d, m->p[i].ddot);
            T1.dt = INFINITY;
          }
        }
	break;
      case 1:
        T = S1TransTime(m, i);
	break;
      case 2:
        T = S2TransTime(m, i);
	break;
      case 3:
        T = S3TransTime(&(m->p[i]));
        break;
      case SLOWSLIP_2A:
      case SLOWSLIP_2B:
      case SLOWSLIP_2C:
        T = S2xTransTime(&(m->p[i]));
        break;
      default:
	fprintf(stderr, "makeTransition():  error: m.p[%d].state = %d\n",
        	i, m->p[i].state);
        T.dt = INFINITY;
    }
    checkIfMinT(&T, minT+omp_rank, m, i);
    
    
  } /* for (i=m->startPatch; i<=m->stopPatch; i++) */
  
/*
 *
 */

  int omp_min_index=0;
#ifdef USE_OPENMP
  for (i=1;i<omp_num_procs;i++) {
          if (minT[i].dt < minT[omp_min_index].dt)
                  omp_min_index=i;
  }
#endif
  
  /* check if max time step is smaller than time to next actual state transition */
  if ( (m->maxTimeStep > 0) && (m->maxTimeStep < minT[omp_min_index].dt) )
  {
    minT[omp_min_index].dt = m->maxTimeStep;
    minT[omp_min_index].patch = 0;
    minT[omp_min_index].nextState = MAX_TIME_STEP;
  }
  
  /* check if next step change in external stressing rate is sooner */
  if (m->nExtStressRateSteps > 0)
  {
    if (m->iExtStressRateStep < m->nExtStressRateSteps) 
    {
      if (m->tExtStressRateSteps[m->iExtStressRateStep] - m->t < minT[omp_min_index].dt)
      {
        minT[omp_min_index].dt = m->tExtStressRateSteps[m->iExtStressRateStep] - m->t;
        minT[omp_min_index].patch = 0;
        minT[omp_min_index].nextState = STRESS_RATE_STEP;
      }
    }
  }
  
  return(minT[omp_min_index]);
}

static void checkIfMinT(Transition *T, Transition *minT, FaultModel *m, int i)
{
  if (T->dt < minT->dt)
  {
    if (minT->patch >= 0 && minT->patch != i)
      if (m->p[minT->patch].APFail == NOW) m->p[minT->patch].APFail = NOTYET;
    minT->dt = T->dt;
    minT->nextState = T->nextState;
    minT->patch = i;
  }
  else
  {
    if (m->p[i].APFail == NOW) m->p[i].APFail = NOTYET;
  }
  
  return;
}


Transition S2TransTime(FaultModel *m, int patch)
{
  Transition T;
  
  FaultPatch *p;
  
  p = &(m->p[patch]);
  
  if (p->tauDot - p->muTrans2*p->sigmaDot >= 0.0)
  { T.dt = INFINITY; T.nextState = 2; return(T);} /* patch will stay in state 2 */

  /* tau should never be less than muTrans2*sigma, if it is transition out immediately */ 
  if (p->tau < p->muTrans2*p->sigma) 
  {
    fprintf(stderr, "tau < muTrans2*sigma in S2toS0Time() on patch %d\n", patch+1);
    T.dt = 0.0; T.nextState = 0;
    return(T);
  }
  
  T.dt = (p->sigma*p->muTrans2 - p->tau)/(p->tauDot - p->muTrans2*p->sigmaDot);
  T.nextState = 0;
  return(T);
} 

Transition S1TransTime(FaultModel *m, int patch) /* also looks for S1toS0 transitions */
{
  Transition T;
  double SDotExt;
  double tol, tolFrac;
  double f1;
  double tTrans;
  double ddotTrans;
  
  FaultPatch *p = &(m->p[patch]);

  tolFrac = 2*(fabs(p->B/(p->A*p->theta)) + fabs((p->tauDot - p->mu0*p->sigmaDot)/(p->A*p->sigma)) +
               fabs(((p->tau - p->mu0*p->sigma)/(p->A*p->sigma))*(p->sigmaDot/p->sigma)));
  if (tolFrac < 10) tolFrac = 10;
  tol = tolFrac*m->ZBrentTol; /* to avoid spurious 0 -> 1 -> 0 -> 1, etc. */

  ddotTrans = (p->slowSlip ? p->ddot2A : p->ddotEQ);
  
  if (p->ddot < (1-tol)*p->Dc/p->theta && p->ddot < ddotTrans && !p->neighborInRuptureState) /* then transition immediately back to state 0 */
  {
    T.dt = 0;
    T.nextState = 0;
    return(T);  
  }

  SDotExt = p->tauDot - ((p->tau/p->sigma) - p->alpha)*p->sigmaDot;
  
  if (m->slipInState1)
    p->H = p->B*p->iDc + (p->Kself - (p->tau/p->sigma - p->alpha)*p->op.Ksigma[patch])/p->sigma;
  else
    p->H = p->B*p->iDc;
    
  if (SDotExt <= -p->H*p->sigma*p->ddot) tTrans = INFINITY; /* really should look for transition back to state 0 */
  else if (fabs(SDotExt/(p->H*p->sigma*ddotTrans)) < tol) 
  {
    f1 = 1 + (SDotExt/(p->H*p->sigma))*(1/p->ddot - 1/ddotTrans);
    if (fabs(1 - f1) < tol)
    {
      tTrans = (p->A/p->H)*(1/p->ddot - 1/ddotTrans);
    }
    else
    {
      tTrans = (p->A*p->sigma/SDotExt)*log(f1);
    }
  }
  else
    tTrans = (p->A*p->sigma/SDotExt)*log( (1 + SDotExt/(p->H*p->sigma*p->ddot)) /
                                      (1 + SDotExt/(p->H*p->sigma*ddotTrans)) );

  /* now check if a priori failure time is sooner that tTrans */
  if (p->tFail - m->t >= 0 && p->tFail - m->t < tTrans && p->APFail == NOTYET) 
  {
    tTrans = p->tFail - m->t;
    p->APFail = NOW; 
  }
  
  T.dt = tTrans;
  T.nextState = (p->slowSlip ? SLOWSLIP_2A : 2);
  return(T);
}

/**** macros, function and variables used by S0toS1Time()        ****/
/**** need to declare some variables out here so that can pass   ****/
/**** a function of one variable to fixedPoint()/zbrent()        ****/
/*static double tau0, tauDot, sigma0, sigmaDot, theta0, A, B, 
              mu0, alpha, ddotStar, Dc;
static int slowSlip;
*/

/* f is the function of which we want to find a zero.
   for patches in state 0, it is guaranteed to be less
   than zero at t = 0.  */
static double f(double t,FaultPatch *p)
{                                                        
  double sigma, tau, theta, thetaDot;
	//setVariables;
   double tau0, tauDot, sigma0, sigmaDot, theta0, A, B, mu0, alpha, ddotStar, iDc;
	int slowSlip;
        tau0 = p->tau;\
        tauDot = p->tauDot;\
        sigma0 = p->sigma;\
        sigmaDot = p->sigmaDot;\
        theta0 = p->theta;\
        A = p->A;\
        B = p->B;\
        mu0 = p->mu0;\
        alpha = p->alpha;\
        ddotStar = p->ddotStar;\
        iDc = p->iDc;\
	//        slowSlip=p->slowSlip;

  sigma = sigma0 + sigmaDot*t;
  tau = tau0 + tauDot*t;
  theta = theta_0(t, sigmaDot, sigma0, alpha, B, theta0);
  /*  thetaDot = 1 - (alpha/B)*(sigmaDot/sigma)*theta; */
  
/*  return( tau -  sigma*(mu0 + (B - A)*log(theta) + A*log(fabs(thetaDot))) ); 
    below seems to work as well on average and simplifies many things */

/*  return(exp(((tau/sigma)-mu0)/(B-A))-theta*ddotStar/Dc); */

  return( tau -  sigma*(mu0 + (B - A)*log(theta*ddotStar*iDc)) );
}

/* g is a function whose fixed point should be the zero
   of f, and because of the nature of f and the fact
   that f(0) < 0, starting the FixedPoint() of g at t = 0
   should always converge <- not so sure of this now with sigmaDot != 0 
   
   No longer used, would need to be updated for dimensional form of 
   constitutive relation */
static double g(double t,FaultPatch *p)
{
  double sigma, tau, theta, thetaDot;
	
//	setVariables;
   double tau0, tauDot, sigma0, sigmaDot, theta0, A, B, mu0, alpha, ddotStar, Dc;
	int slowSlip;
        tau0 = p->tau;\
        tauDot = p->tauDot;\
        sigma0 = p->sigma;\
        sigmaDot = p->sigmaDot;\
        theta0 = p->theta;\
        A = p->A;\
        B = p->B;\
        mu0 = p->mu0;\
        alpha = p->alpha;\
        ddotStar = p->ddotStar;\
	//        Dc = p->Dc;  \
        //slowSlip=p->slowSlip;

  sigma = sigma0 + sigmaDot*t;
  tau = tau0 + tauDot*t;
  theta = theta_0(t, sigmaDot, sigma0, alpha, B, theta0);
  thetaDot = 1 - (alpha/B)*(sigmaDot/sigma)*theta; 

  return( (1/tauDot)*( sigma*(mu0 + (B-A)*log(theta) + A*log(fabs(thetaDot))) - tau0 ) );
}


Transition S0TransTime(FaultModel *m, int i, double minDt)
{
  double f0;
  double tTrans=-1; /* -1 is just a flag so we can tell if it has been reset
                       by any of the if blocks below */
  double zBrentUpperBracket;
  
  Transition T;
// 	setVariables;
  double tau0, tauDot, sigma0, sigmaDot, theta0, A, B, mu0, alpha, ddotStar, Dc;
  int slowSlip;
  double t0,t1;
  
  FaultPatch *p;
  
  p = &(m->p[i]);

  tau0 = p->tau;
  tauDot = p->tauDot;
  sigma0 = p->sigma;
  sigmaDot = p->sigmaDot;
  theta0 = p->theta;
  A = p->A;
  B = p->B;
  mu0 = p->mu0;
  alpha = p->alpha;
  ddotStar = p->ddotStar;
  //  Dc = p->Dc;
  //slowSlip = p->slowSlip;


  /* check sign of f() at dt = 0, if positive, then this patch should
     already have transitioned out of state 0 (this can happen at model
     startup, depending on the initial values of tau and theta)  */
  if (f(0,p) > 0.0) {T.dt = 0.0;  T.nextState = 1; p->t01 = 0.0; return(T);}

  /* if sigmaDot == 0 and tauDot is negative, then this patch
     will stay in state 0  */
  if (p->sigmaDot == 0 && p->tauDot <= 0) 
    tTrans = INFINITY;
  
  /* if sigmaDot < 0 and tau at -sigma0/sigmaDot is less than zero
     then the patch will stay in state 0 */
  else if (p->sigmaDot < 0 && p->tau + p->tauDot*(-p->sigma/p->sigmaDot) < 0.0) 
    tTrans = INFINITY;
  

  /* next check the sign of f() what we will try to use as the upper bracket for zbrent(),
     does this still work if sigmaDot != 0 ??? */
  else 
  {
      if (m->ZBrentUpperBracket == 0.0)
        zBrentUpperBracket = MIN(MIN(minDt, fabs(0.99*sigma0/sigmaDot)), 1e30);
      else
        zBrentUpperBracket = MIN(fabs(0.99*sigma0/sigmaDot), m->ZBrentUpperBracket);

      if (f(zBrentUpperBracket,p) < 0.0) tTrans = INFINITY;
  }
  
  if (tTrans == -1) /* that is, if we haven't reset it in any of the if blocks above, then root find */ 
  {
    if (p->t01 != T01_NA_VAL) 
    {
/*      fprintf(stderr, "Saved a rootfinding episode\n"); */
      tTrans = p->t01;
    } else
    {
/*      fprintf(stderr, "Didn't save a rootfinding episode\n"); */

      t0 = MPI_Wtime();

      tTrans = zbrent(f, 0.0, zBrentUpperBracket, m->ZBrentTol, p);
      
      t1 = MPI_Wtime();
      m->zbrentTime += t1-t0;

      if (!isfinite(tTrans) && m->DEBUG > 1)
      {  
        fprintf(stderr, "zbrent() failed on patch %d\n", i);
        fprintf(stderr, "tau0 = %lg; tauDot = %lg; sigma0 = %lg; sigmaDot = %lg; theta0 = %lg\n",
                tau0, tauDot, sigma0, sigmaDot, theta0);
      }
    }
  }
  
  /* now check if a priori failure time is sooner that tTrans */
  if (p->tFail - m->t >= 0 && p->tFail - m->t < tTrans && p->APFail == NOTYET)
  {
    tTrans = p->tFail - m->t;
    p->APFail = NOW;
  }
  
  T.dt = tTrans;
  T.nextState = 1;
  
  p->t01 = (isfinite(tTrans) ? tTrans : T01_NA_VAL);
  
  return(T);
}


Transition S3TransTime(FaultPatch *p)
{
  double tLow, tHigh;
  Transition T;
  
  /* check if need to immediately transition; sometimes need to due to roundoff */
  if (p->tau/p->sigma <= p->muTrans3Low || p->tau/p->sigma >= p->muTrans3High)
  {
    T.dt = 0;
    T.nextState = 3;
    return(T);
  }
  
  tLow = (p->muTrans3Low*p->sigma - p->tau)/(p->tauDot - p->muTrans3Low*p->sigmaDot);
  tHigh = (p->muTrans3High*p->sigma - p->tau)/(p->tauDot - p->muTrans3High*p->sigmaDot);
  T.dt = INFINITY;
  if (tLow >= 0) T.dt = tLow;
  if (tHigh >= 0 && tHigh < T.dt) T.dt = tHigh;
  
  T.nextState = 3;
  
  
   /* turns out that this is not a problem, just means that tauDot/sigmaDot is between
      muTrans3Low and muTrans3High, so that at the current stressing rates, this
      patch will never need to have its slip speed adjusted 
  if (tLow < 0 && tHigh < 0)
  {
    fprintf(stderr, "Warning in S3toS3Time(): tLow and tHigh both < 0\n");
  } */
  
  
  return(T);
}


Transition S2xTransTime(FaultPatch *p)
{
  double tLow, tHigh;
  Transition T;
  
  /* check if need to transition immediately */
  if (p->state == SLOWSLIP_2A)
  {
    if (p->tau/p->sigma > p->muTrans2High)
      {T.dt = 0; T.nextState = SLOWSLIP_2B; return(T);}
    if (p->tau/p->sigma < p->muTrans2aLow)
      {T.dt = 0; T.nextState = LOCK; return(T);}
  }
  if (p->state == SLOWSLIP_2B)
  {
    if (p->tau/p->sigma >= (p->mu0 + p->A*log(p->ddotEQ/p->ddotStar))*(1.0 /* + MU_EPS */))
      {T.dt = 0; T.nextState = SLOWSLIP_2C; return(T);}
    if (p->tau/p->sigma < (p->mu0 + p->A*log(p->ddotAB/p->ddotStar))*(1.0 - MU_EPS))
      {T.dt = 0; T.nextState = SLOWSLIP_2A; return(T);}
  }
  if (p->state == SLOWSLIP_2C)
  {
    if (p->tau/p->sigma < p->muTrans2aLow)
      {T.dt = 0; T.nextState = LOCK; return(T);}
  }
  
  
  if (p->state == SLOWSLIP_2A || p->state == SLOWSLIP_2C)
    tLow = (p->muTrans2aLow*p->sigma - p->tau)/(p->tauDot - p->muTrans2aLow*p->sigmaDot);
  else   /* p->state == SLOWSLIP_2B */
    tLow = (p->muTrans2Low*p->sigma - p->tau)/(p->tauDot - p->muTrans2Low*p->sigmaDot);
  tHigh = (p->muTrans2High*p->sigma - p->tau)/(p->tauDot - p->muTrans2High*p->sigmaDot);
  T.dt = INFINITY;
  if (tLow > 0) T.dt = tLow;
  if (tHigh > 0 && tHigh < T.dt) T.dt = tHigh;
  
  if (p->state == SLOWSLIP_2A)
  {
    if (T.dt == tLow) {T.nextState = 0;} 
    else if (T.dt == tHigh) {T.nextState = SLOWSLIP_2B;}
    else {T.nextState = SLOWSLIP_2A;} /* shouldn't happen */
  }
  else if (p->state == SLOWSLIP_2B)
  {
    if (T.dt == tLow && p->muTrans2Low == (p->mu0 + p->A*log(p->ddotAB/p->ddotStar))*(1.0 - MU_EPS)) {T.nextState = SLOWSLIP_2A;}
    else if (T.dt == tHigh && p->muTrans2High == (p->mu0 + p->A*log(p->ddotEQ/p->ddotStar))*(1.0 + MU_EPS)) {T.nextState = SLOWSLIP_2C;}
    else {T.nextState = SLOWSLIP_2B;}  /* will happen a lot */
  }
  else /* p->state == SLOWSLIP_2C */
  {
    if (T.dt == tLow) {T.nextState = LOCK;}
    else {T.nextState = SLOWSLIP_2C;} /* shouldn't happen */
  }
  
  return(T);
}


Transition checkForLowStress(FaultPatch *p)
{
  Transition T;
  
  T.dt = INFINITY;
  
  if (p->sigma < 0.0)
  {
    T.dt = 0;
    T.nextState = LOW_SIGMA;
    return(T);
  }
  
  if (p->tau < 0)
  {
    T.dt = 0;
    T.nextState = LOW_TAU;
    return(T);
  }

  if (p->sigmaDot < 0  && !p->sigmaClipped) 
  {  
    T.dt = (p->sigmaPin - p->sigma)/p->sigmaDot;
    T.nextState = LOW_SIGMA;
  }
  
  if (p->tauDot < 0 && -p->tau/p->tauDot < T.dt)
  {
    double dt;
    dt = -p->tau/p->tauDot;
    if (dt < T.dt)
    {
      T.dt = dt;
      T.nextState = LOW_TAU;
    }
  }
  
  if (p->muPin > 0) /* check for mu exceeding muPin */
  {
    double dt;
    dt = (p->muPin * p->sigma - p->tau)/(p->tauDot - p->muPin * p->sigmaDot);
    if (dt >= 0 && dt < T.dt)
    {
      T.dt = dt;
      T.nextState = HIGH_TAU;
    }
  }
  
  if (T.dt < 0)
  {
     fprintf(stderr, "checkForLowStress(): warning: T.dt < 0 (%lf): setting to zero\n", T.dt);
     T.dt = 0;
  }
  
  return(T);
}

Transition checkForHighTheta(FaultPatch *p, double t, double maxThetaPin)
{
  Transition T;
  
  T.dt = INFINITY;
  
  if (p->state == CREEP) return(T);  /* don't apply this to creeping elements */
  
  if (p->theta > 1e3*t + maxThetaPin) 
  {  
    T.dt = 0;
    T.nextState = HIGH_THETA;
  }
    
  return(T);
}
