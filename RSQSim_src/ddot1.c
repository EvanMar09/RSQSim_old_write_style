#include "RSQSim.h"

/* slip speed evolution in state 1 */

double ddot1(double t, double ddot0, double SDotExt, double A,
             double sigma0, double H, double ddotEQ)
{
  double tol = 1e-7;
  
  double ddot1r;  /* value to be returned */
  
  if (fabs((SDotExt/(A*sigma0))*t) < tol)  
    ddot1r = 1/( (1/ddot0)*exp(-(SDotExt/(A*sigma0))*t) - (H/A)*t);
  else
    ddot1r = 1/( (1/ddot0 + H*sigma0/SDotExt)*exp(-(SDotExt/(A*sigma0))*t) -
                H*sigma0/SDotExt);

  /* assume that if ddot1r is too big or negative or Inf or NaN that it was supposed 
     to go to ddotEQ but missed due to roundoff error.  FIXME */
  if (!isfinite(ddot1r) || ddot1r < 0 || ddot1r > ddotEQ) ddot1r = ddotEQ;
  
  return(ddot1r);
}         

/* returns the slip speed as a function of stresses, theta, etc.
   from solving the constitutive relation for V */
double ddotRS(double ddotStar, double tau, double sigma, double mu0, 
              double a, double b, double theta, double Dc)
{
  return(ddotStar * exp( (tau - mu0*sigma)/(a*sigma) - (b/a)*log(theta * ddotStar / Dc)));
}
