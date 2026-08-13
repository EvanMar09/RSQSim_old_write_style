#include "KRDOkada.h"

/* 
   Given a stress tensor Sigma (as a 6-component vector like: (0 1 2)
                                                              (  3 4)
							      (    5)),
   and a unit normal vector nu to a fault plane and a unit slip vector u (of the
   material into which nu points), calculates the normal stress on the
   plane (tension positive) and the shear stress in the direction of u
*/

int projectStress(double *Sigma, /* stress tensor as a 6-component vector */
                  double *nu, /* unit normal */
		  double *u, /* slip direction (of material into which nu points) */
		  double *sigma, /* scalar normal stress (output) */
		  double *tau /* scalar shear stress (output) */
		  )
{
  double ss[3][3];  /* stress tensor as a 3x3 array */
  int i, j, k;
  
  /* expand Sigma to a 3x3 array */
  k = 0;
  for (i=0; i<3; i++) /* upper triangle */
  {
    for (j=i; j<3; j++) ss[i][j] = Sigma[k++];
  }
  for (i=0; i<3; i++) /* lower triangle */
  {
    for (j=0; j<i; j++) ss[i][j] = ss[j][i];
  }
  
  *sigma = 0.0;
  *tau = 0.0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      *sigma += nu[i]*ss[i][j]*nu[j];
      *tau += u[i]*ss[i][j]*nu[j];
    }
  }
    
  return(0);
}
  
  
  
  


