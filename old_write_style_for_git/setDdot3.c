#include "RSQSim.h"

static char svnId[] = {"$Id$"};

void setDdot3(FaultPatch *p)
{
  p->ddot = p->ddotStar * exp( (p->tau/p->sigma - p->mu0)/(p->A - p->B) );
  p->ddot = MIN(p->ddot, p->ddotEQ);

  return;
}

void setDdot2b(FaultPatch *p)
{
  p->ddot = p->ddotStar * exp( (p->tau/p->sigma - p->mu0)/p->A );
  p->ddot = MIN(p->ddot, p->ddotEQ);
  
  return;
}  

  
