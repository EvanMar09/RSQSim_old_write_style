#ifndef KRDOKADA_H
#define KRDOKADA_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dc3d.h"
#include "array.h"
#include "xmalloc.h"

#define OKADA      0
#define TRIANGULAR 1

/* 
   for type == OKADA, patch must be rectangular with two edges parallel to xy-plane 
   
   coord system is East, North, Up
*/
typedef struct
{
  /* several redundant ways of specifying the geometry of the
     patch and the slip direction on that patch */
  double corner[4][3]; /* corner[i][j] is the x_j coord. of the i^th corner 
                         for type == OKADA:
                         corner[0] should be the lower corner in the
			 away-from-strike direction, corner[1] the lower
			 corner in the along-strike direction, corner[2]
			 the upper corner in the along-strike direction and
			 corner[3] the upper corner in the away-from-strike
			 direction 
                         for type == TRIANGULAR:
                         the corners need to be in an order such that as you curl the
                         fingers of your right hand around the perimeter of the triangle
                         in order, the thumb of your right hand points from the footwall
                         to the hanging wall */
  double center[3];  
  double L;            /* length in the along-strike direction (only for type == OKADA) */
  double W;            /* width in the along-dip direction (only for type == OKADA) */ 
  double area;  
  double strike;       /* as in Aki and Richards; in radians */
  double dip;          /* as in Aki and Richards; in radians */
  double rake;         /* as in Aki and Richards; in radians */
  double nu[3];         /* upward pointing unit normal vector, i.e. unit
                          vector pointing into hanging wall */
  double u[3];         /* unit slip vector (of hanging wall relative to footwall).  u.nu should probably be zero */ 
  

  double **K;    /* the rows of the "stiffness" matrices for this patch as a
                   receiving patch.  That is, the change in the alpha^th component
		   of the stress tensor at the center of this patch due to
		   slips d[j] on the other patches is given by
		   K[j][alpha]*d[j] (with a sum over j, obviously).  The components
		   of the stress tensor (and therefore the stiffness matrix)
		   are stored like so: (0 1 2)
		                       (  3 4)
				       (    5) 
	           Note:  I ordered the j and alpha indices this way for ease
		          in passing to, e.g., deformationToStress() */
  double *Ksigma;  /* the row of the stiffness matrix for the normal stress on this 
                     patch as a receiving patch.  That is, the change in normal 
		     stress at the center of this patch due to slips d[j] on the 
		     other patches is given by Ksigma[j]*d[j] */
  double *Ktau;    /* the row of the stiffness matrix for the shear stress (IN THE
                     PRESCRIBED SLIP DIRECTION FOR THIS PATCH) for this patch
		     as a receiving patch */
  
  
  int N; /* number of fault patches for which space in the K vectors has
            been allocated */

  int type;  /* is this an Okada-type patch (OKADA) or a triangular one (TRIANGULAR)? */
} KRDOkadaFaultPatch;

/* cross() stores the cross product of x and y in z */
static void cross(double *x, double *y, double *z)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  
  return;
}

/* makes x a unit vector */
static void normalize(double *x)
{
  double norm;
  
  norm = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  x[0] /= norm;
  x[1] /= norm;
  x[2] /= norm;
  
  return;
}

/* returns norm of a 3-vector */
static double norm(double *x)
{
  return(sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
}  
  

int initKRDOkadaFaultPatch(KRDOkadaFaultPatch *f);
int calcStiffnessMatrix(void *f, size_t offset, size_t size, int N, double lambda, double mu, int allK,
                        int startPatch, int stopPatch, int pointSource, double dist1, double dist2);
int deformationToStress(double du[3][3], double lambda, double mu, double *sigma);
int projectStress(double *Sigma, double *nu, double *u, double *sigma, double *tau);
int findNormalAndSlip(KRDOkadaFaultPatch *f);
int findCorners(KRDOkadaFaultPatch *f);
int findCenter(KRDOkadaFaultPatch *f);
int allocKRDOkadaFaultPatch(KRDOkadaFaultPatch *f, int N, int allK, int anyK);
int freeKRDOkadaFaultPatch(KRDOkadaFaultPatch *f);
int approxTgf(double lambda, double mu, double triangle[3][3], double target[][3], int ntarg, double u1, double u2, int npt,
              double nu[3], double area, double strike, double dip,double u[][3], double du[][3][3]);
void approxTgf_passByRef(double *lambda_p, double *mu_p,                          
                         double *trianglev1, double *trianglev2, double *trianglev3,
                         double *targetx, double *targety, double *targetz,
                         int *ntarg_p, double *u1_p, double *u2_p, int *npts_p, double *nu, double *area_p,
                         double *strike_p, double *dip_p, 
                         double *ux, double *uy, double *uz,
                         double *uxx, double *uxy, double *uxz,
                         double *uyx, double *uyy, double *uyz,
                         double *uzx, double *uzy, double *uzz);
int getDunavantPtsWeights(int n, double *w, double *alpha, double *beta, 
                          double *gamma);

extern void 
  tg_elh3dtriadirecttarg_(double *rlam,
                          double *rmu,
                          double triangle[][3][3],
                          double trinorm[][3],
                          int *nsource,
                          double source[][3],
                          int *ifsingle,
                          double sigma_sl[][3],
                          int *ifdouble,
                          double sigma_dl[][3],
                          int *ifptfrc,
                          double ptfrc[][3],
                          int *ifstrain,
                          double strain[][3][3],
                          int *ntarget,
                          double target[][3],
                          int *ifptfrctarg,
                          double ptfrctarg[][3],
                          int *ifstraintarg,
                          double straintarg[][3][3]);


#endif /* ndef KRDOKADA_H */
