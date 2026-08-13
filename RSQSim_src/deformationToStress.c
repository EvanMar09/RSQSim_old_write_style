#include "KRDOkada.h"

/* 
   given a deformation tensor (i.e. du[i][j] = du_i/dx_j) calculates the
   stress tensor.  The stress tensor is stored as a vector with the components
   ordered thusly: (0 1 2)
		   (  3 4)
		   (    5) 
		   
   lambda and mu are Lame's parameters
   
   space needs to already have been allocated for sigma
*/

int deformationToStress(double du[3][3], /* deformation tensor, see above */
                        double lambda, double mu, /* Lame's constants */
			double *sigma  /* stress tensor as a vector, see above */
			)
{
  double ss[3][3]; /* the stress tensor as a 3x3 array */
  double tre; /* the trace of the strain tensor (which is equal to the trace of du */
  int i, j, k;
  
  /* first just calc the mu part */
  for (i=0; i<3; i++)
  {
    for (j=i; j<3; j++) /* only need to calc upper triangle */
      ss[i][j] = mu*(du[i][j] + du[j][i]);
  }
  
  /* find the trace of the strain tensor (= trace of du), needed for 
     lambda part of stress tensor */
  tre = 0.0;
  for (i=0; i<3; i++) tre += du[i][i];
  
  /* add on lambda part (only affects diagonal) */
  for (i=0; i<3; i++) ss[i][i] += lambda*tre;
    
  /* store stress as 6-component vector */
  k = 0;
  for (i=0; i<3; i++)
  {
    for (j=i; j<3; j++) /* only need to do upper triangle */
      sigma[k++] = ss[i][j];
  }
  
  return(0);
}

      
