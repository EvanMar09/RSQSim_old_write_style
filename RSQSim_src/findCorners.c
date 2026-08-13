#include "KRDOkada.h"

/* 
   finds the x, y, and z coords of the corners of the KRDOkadaFaultPatch
   from the coords of the center, length, width, strike, and dip.
   With the corners numbered in the order described in fault.h
   
   If the latter entries of the KRDOkadaFaultPatch are not valid, you will
   get garbage and/or NaN results - should add code to check
*/

int findCorners(KRDOkadaFaultPatch *f)
{
  double st[3]; /* unit vector in strike direction */
  double dd[3]; /* unit vector in down-dip direction */
  int i;
  
  st[0] = sin(f->strike);  st[1] = cos(f->strike); st[2] = 0.0;
  dd[0] = cos(f->dip)*cos(f->strike);  
  dd[1] = -cos(f->dip)*sin(f->strike);
  dd[2] = -sin(f->dip);
  
  for (i=0; i<3; i++) 
  {
    f->corner[0][i] = f->center[i] - (f->L/2)*st[i] + (f->W/2)*dd[i];
    f->corner[1][i] = f->center[i] + (f->L/2)*st[i] + (f->W/2)*dd[i];
    f->corner[2][i] = f->center[i] + (f->L/2)*st[i] - (f->W/2)*dd[i];
    f->corner[3][i] = f->center[i] - (f->L/2)*st[i] - (f->W/2)*dd[i];
  }
  
  return(0);
}
  

/* findCenter() finds the center of the TRIANGULAR KRDOKadaFaultPatch
   from the coords of the corners (it is just the average of the 
   corners) */

int findCenter(KRDOkadaFaultPatch *f)
{
  int ii, ic;
  
  if (f->type != TRIANGULAR)
  {
    fprintf(stderr, "findCenter(): Error, only for TRIANGULAR patches\n");
    return(-1);
  }
  
  for (ii=0; ii<3; ii++)
  {
    f->center[ii] = 0;
    for (ic=0; ic<3; ic++) f->center[ii] += f->corner[ic][ii]/3;
  }
  
  return(0);
}
