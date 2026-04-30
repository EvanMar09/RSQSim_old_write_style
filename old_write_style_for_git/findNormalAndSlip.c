#include "KRDOkada.h"

/*
   for type == OKADA, does what the name of the function implies:
   populates the unit normal and unit slip vector elements of the KRDOkadaFaultPatch
   using the strike, dip, and rake
   
   for type == TRIANGULAR, does the normal and slip elements as well as the
   strike, dip, and area from the corners and the rake
*/

int findNormalAndSlip(KRDOkadaFaultPatch *f)
{
  double st[3]; /* unit vector in strike direction */
  double ud[3]; /* unit vector in up-dip direction */
  double cosr, sinr;  /* cosine and sine of rake angle */
  int ic;
  
  if (f->type == OKADA)
  {
    st[0] = sin(f->strike);  st[1] = cos(f->strike); st[2] = 0.0;
    ud[0] = -cos(f->dip)*cos(f->strike);  
    ud[1] = cos(f->dip)*sin(f->strike);
    ud[2] = sin(f->dip);

    cross(st, ud, f->nu);

    cosr = cos(f->rake);  sinr = sin(f->rake);
    for (ic=0; ic<3; ic++) f->u[ic] = cosr*st[ic] + sinr*ud[ic];
  }
  else if (f->type == TRIANGULAR)
  {
    double dx1[3], dx2[3];
    
    /* get normal from cross product of two legs of triangle */
    for (ic=0; ic<3; ic++) 
    {
      dx1[ic] = f->corner[1][ic] - f->corner[0][ic];
      dx2[ic] = f->corner[2][ic] - f->corner[0][ic];
    }
    cross(dx1, dx2, f->nu);
    f->area = 0.5*norm(f->nu);
    normalize(f->nu);
    
    f->strike = -atan2(f->nu[1], f->nu[0]); /* see man page for atan2() for if nu[0] == nu[1] == 0 */
    
    f->dip = M_PI/2 - atan2(f->nu[2], sqrt(f->nu[0]*f->nu[0] + f->nu[1]*f->nu[1]));
    
    /* with strike and dip now calculated, we can use the same method
       as above to calculate the slip vector */
    st[0] = sin(f->strike);  st[1] = cos(f->strike); st[2] = 0.0;
    ud[0] = -cos(f->dip)*cos(f->strike);  
    ud[1] = cos(f->dip)*sin(f->strike);
    ud[2] = sin(f->dip);

    cosr = cos(f->rake);  sinr = sin(f->rake);
    for (ic=0; ic<3; ic++) f->u[ic] = cosr*st[ic] + sinr*ud[ic];
    
  }  
  
  return(0);
}

