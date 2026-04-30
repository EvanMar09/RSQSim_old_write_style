#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifdef USE_MPI
#include "mpi.h"
#endif
#include "KRDOkada.h"
#include "RSQSim_dataStructures.h"

#define ITMAX 100 /* Maximum allowed number of iterations. */
#define EPS DBL_EPSILON /* Machine doubleing-point precision. */
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double zbrent(double (*func)(double,FaultPatch*), double x1, double x2, double tol,FaultPatch *pFault)
/* Using Brent's method, find the root of a function func known to lie between x1 and x2. The
root, returned as zbrent, will be refined until its accuracy is tol. */
{
  int iter;
  double a=x1,b=x2,c=x2,d,e,min1,min2;
/* double fa=(*func)(a);
  double fb=(*func)(b);
*/
 double fa=(*func)(a,pFault);
 double fb=(*func)(b,pFault);
  
double fc,p,q,r,s,tol1,xm;
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
  {
/*    fprintf(stderr, "Root must be bracketed in zbrent\n"); */
    return(INFINITY);
  }
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
    c=a; /* Rename a, b, c and adjust bounding interval */
    fc=fa; 
    e=d=b-a;
  }
  if (fabs(fc) < fabs(fb)) {
    a=b;
    b=c;
    c=a;
    fa=fb;
    fb=fc;
    fc=fa;
  }
  tol1=2.0*EPS*fabs(b)+0.5*tol; /* Convergence check. */
  xm=0.5*(c-b);
  if (fabs(xm) <= tol1 || fb == 0.0) return b;
  if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
    s=fb/fa; /* Attempt inverse quadratic interpolation. */
    if (a == c) {
      p=2.0*xm*s;
      q=1.0-s;
    } else {
      q=fa/fc;
      r=fb/fc;
      p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
      q=(q-1.0)*(r-1.0)*(s-1.0);
    }
    if (p > 0.0) q = -q; /* Check whether in bounds. */
    p=fabs(p);
    min1=3.0*xm*q-fabs(tol1*q);
    min2=fabs(e*q);
    if (2.0*p < (min1 < min2 ? min1 : min2)) {
      e=d; /* Accept interpolation. */
      d=p/q;
    } else {
      d=xm; /* Interpolation failed, use bisection. */
      e=d;
    }
    } else { /* Bounds decreasing too slowly, use bisection. */
      d=xm;
      e=d;
    }
    a=b; /* Move last best guess to a. */
    fa=fb;
    if (fabs(d) > tol1) /* Evaluate new trial root. */
      b += d;
    else
      b += SIGN(tol1,xm);
//    fb=(*func)(b);
  fb=(*func)(b,pFault);

  }
  fprintf(stdout, "Maximum number of iterations exceeded in zbrent\n");
  return(INFINITY);
}
