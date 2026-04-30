#include "RSQSim.h"

static char svnId[] = {"$Id: fixedPoint.c 161 2009-01-06 00:07:10Z dinger $"};

/* 
   repeatedly applies the function f to the result of f(x0) until the
   change from one iteration to another is less than acc + abs(x)*prec
   
   need to add: max number of iterations
                checks for if x goes to +-Inf
*/

double fixedPoint(double (*f)(double x), double x0, double acc, double prec,
                  double xmin, double xmax)
{
  double dx = INFINITY;
  double xnew;
  
  while (fabs(dx) > acc + fabs(x0)*prec)
  {  
    if (x0 < xmin || x0 > xmax) return(NAN);
    xnew = f(x0);
    dx = xnew - x0;
    x0 = xnew;
    if (!isfinite(dx)) return(NAN);
  }
  return(xnew);
}
