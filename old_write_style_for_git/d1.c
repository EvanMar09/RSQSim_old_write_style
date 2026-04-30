#include "RSQSim.h"

/* displacement in state 1 */

double d1(double t, double A, double H, double ddot0, double sigma0,
          double SDotExt)
{
  double tol=1e-7;
  double f1;
  
  if (fabs((SDotExt/(A*sigma0))*t) < tol) f1 = -(ddot0*H/A)*t;
  else f1 = (ddot0*H*sigma0/SDotExt)*( 1 - exp((SDotExt/(A*sigma0))*t) );
  
  if (fabs(f1) < tol) return(-(A/H)*f1);
  else if (f1 <= -1) return(-(A/H)*log(DBL_EPSILON)); /* COMPLETE HACK  FIXME */
  else return(-(A/H)*log(1 + f1));
}

