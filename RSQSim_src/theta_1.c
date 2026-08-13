#include "RSQSim.h"

/* theta as a function of tau, sigma, ddot, etc. - from solving  the
   R&S constitutive relation for theta; used in state 1
   
   used for resetting theta on state 1 patches upon termination of
   rupture
 */
 
double theta_1(double tau, double sigma, double ddot, double mu0, double ddotStar,
               double A, double B, double Dc)
{
  return((Dc/ddotStar) * exp( (tau - mu0*sigma)/(B*sigma))*pow(ddot/ddotStar, -A/B));
}

/* alternate expression to use for theta in state 1;  from integrating
   the state 1 approx for thetaDot 
   
   
   used for evolving theta (which is for informational purposes
   only, actually also for checking if should transition back to 
   state 0) in state 1 
*/
double theta_1a(double theta0, double d, double Dc, double alpha, 
                double B, double sigma0, double sigma)
{
  return( theta0*exp(-d/Dc - (alpha/B)*log(sigma/sigma0)) );
}
