#include "RSQSim.h"

/* theta evolution in state 0 */
double theta_0_no_use(double t, double sigmaDot, double sigma0, double alpha,
               double B, double theta0)
{
  double tol = 1e-10;
  double temp,aB,pW;
	temp=sigmaDot*t/sigma0;
	aB= alpha/B;

  if (fabs(temp) < tol)
  {
    return(theta0*(1 - aB*temp) + t);
  }
  else if ((1 + temp ) <= 0)
  {
    return(INFINITY);
  }
  else
  {
   return( pow(1+temp, -aB)*(theta0 +
             (1/(1+aB))*(sigma0/sigmaDot)*(pow(1+temp, 1+aB) - 1)) );

  }
}

/* theta evolution in state 0 */
double theta_0(double t, double sigmaDot, double sigma0, double alpha,
               double B, double theta0)
{
  double tol = 1e-10;
  double tmp1,tmp2;

  tmp1 = sigmaDot*t/sigma0;
  tmp2 = alpha/B;

  if (fabs(tmp1) < tol)
  {
    return(theta0*(1 - tmp2*tmp1) + t);
  }
  else if (1 + tmp1 <= 0)
  {
    return(INFINITY);
  }
  else
  {
    return( pow(1+tmp1, -tmp2)*(theta0 + 
		(sigma0/(sigmaDot*(1+tmp2)))*(pow(1+tmp1, 1+tmp2) - 1)) );
  }
}
