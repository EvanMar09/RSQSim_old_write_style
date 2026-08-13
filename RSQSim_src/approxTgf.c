#include "RSQSim.h"

/*
###############################################################
##### approximate tgf with n Dunavant point sources
#####
##### u1 and u2 are the strike-slip and dip-slip components 
##### of the motion of the "+" face (rh rule) of the triangle 
##### relative to the "-" face
###############################################################
*/

int approxTgf(double lambda, double mu, double triangle[3][3], double target[][3], int ntarg, double u1, double u2, int npts,
              double nu[3], double area, double strike, double dip, double u[][3], double du[][3][3])
{
  double *w, *alpha, *beta, *gamma;
  double theta, cost, sint, xsv, ysv; /* xsv and ysv for use in rotating */
  double u0[3], du0[3][3];
  double alpha1; /* Okada's alpha */
  int icorner, i, j, ipt, itarg;
  double xpt, ypt, zpt;
  double ltriangle[3][3], (*ltarget)[3]; /* local copies of triangle and targets */
      
  w = xmalloc(npts*sizeof(double));
  alpha = xmalloc(npts*sizeof(double));
  beta = xmalloc(npts*sizeof(double));
  gamma = xmalloc(npts*sizeof(double));
  ltarget = xmalloc(ntarg*3*sizeof(double));
  
  if (getDunavantPtsWeights(npts, w, alpha, beta, gamma) < 0)  
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(-1);
#endif
   
  alpha1 = (lambda + mu)/(lambda + 2*mu); /* this is the alpha used by Okada */
  

  theta = M_PI/2 - strike;
  cost = cos(theta);
  sint = sin(theta);
  
/* targets in coord. system rotated so that strike direction of triangle is along +x axis */
  for (itarg=0; itarg<ntarg; itarg++)
  {
    xsv = target[itarg][0]; ysv = target[itarg][1];
    ltarget[itarg][0] = cost*target[itarg][0] + sint*target[itarg][1]; 
    ltarget[itarg][1] = -sint*target[itarg][0] + cost*target[itarg][1];
    ltarget[itarg][2] = target[itarg][2];
  }

/* triangle in coord. system rotated so that strike direction of triangle is along +x axis */
  for (icorner=0; icorner<3; icorner++)
  {
    ltriangle[icorner][0] = cost*triangle[icorner][0] + sint*triangle[icorner][1];
    ltriangle[icorner][1] = -sint*triangle[icorner][0] + cost*triangle[icorner][1];
    ltriangle[icorner][2] = triangle[icorner][2];
  }
  
  for (itarg=0; itarg<ntarg; itarg++) 
  {
    for (i=0; i<3; i++) 
    {
      u[itarg][i] = 0.0; 
      for (j=0; j<3; j++) du[itarg][i][j] = 0.0;
    }
  }
  for (ipt=0; ipt<npts; ipt++)
  {
    xpt = alpha[ipt]*ltriangle[0][0] + beta[ipt]*ltriangle[1][0] + gamma[ipt]*ltriangle[2][0];
    ypt = alpha[ipt]*ltriangle[0][1] + beta[ipt]*ltriangle[1][1] + gamma[ipt]*ltriangle[2][1];
    zpt = alpha[ipt]*ltriangle[0][2] + beta[ipt]*ltriangle[1][2] + gamma[ipt]*ltriangle[2][2];
    
    for (itarg=0; itarg<ntarg; itarg++)
    {
      dc3d0(alpha1, ltarget[itarg][0] - xpt, ltarget[itarg][1] - ypt, ltarget[itarg][2], -zpt, 
            (180.0/M_PI)*dip, area*u1, area*u2, 0, 0,
            &(u0[0]), &(u0[1]), &(u0[2]),
            &(du0[0][0]), &(du0[1][0]), &(du0[2][0]),
            &(du0[0][1]), &(du0[1][1]), &(du0[2][1]),
            &(du0[0][2]), &(du0[1][2]), &(du0[2][2]));

      for (i=0; i<3; i++)
      {
        u[itarg][i] += w[ipt]*u0[i];
        for (j=0; j<3; j++) du[itarg][i][j] += w[ipt]*du0[i][j];
      }
    }
  }

/* rotate results back to orig coord system */
  for (itarg=0; itarg<ntarg; itarg++)
  {
    xsv = u[itarg][0]; ysv = u[itarg][1];
    u[itarg][0] = cost*xsv - sint*ysv;
    u[itarg][1] = sint*xsv + cost*ysv;
    
    for (i=0; i<3; i++)
    {
      xsv = du[itarg][i][0]; ysv = du[itarg][i][1];
      du[itarg][i][0] = cost*xsv - sint*ysv;
      du[itarg][i][1] = sint*xsv + cost*ysv;
    }
    for (i=0; i<3; i++)
    {
      xsv = du[itarg][0][i]; ysv = du[itarg][1][i];
      du[itarg][0][i] = cost*xsv - sint*ysv;
      du[itarg][1][i] = sint*xsv + cost*ysv;
    }
  }
  
  free(w); free(alpha); free(beta); free(gamma);
  free(ltarget); 

  return(0);
}

/* this version is so R can call it */
void approxTgf_passByRef(double *lambda_p, double *mu_p, 
                         double *trianglev1, double *trianglev2, double *trianglev3,
                         double *targetx, double *targety, double *targetz,
                         int *ntarg_p, double *u1_p, double *u2_p, int *npts_p, double *nu, double *area_p,
                         double *strike_p, double *dip_p, 
                         double *ux, double *uy, double *uz,
                         double *uxx, double *uxy, double *uxz,
                         double *uyx, double *uyy, double *uyz,
                         double *uzx, double *uzy, double *uzz)
{
  double lambda, mu, area, strike, dip, u1, u2;
  int ntarg, npts, err;
  int itarg, ixyz;
  
  double u[1][3], du[1][3][3];
  double triangle[3][3], target[1][3];
  
  
  for (ixyz=0; ixyz<3; ixyz++)
  {
    triangle[0][ixyz] = trianglev1[ixyz];
    triangle[1][ixyz] = trianglev2[ixyz];
    triangle[2][ixyz] = trianglev3[ixyz];
  }
  
  lambda = *lambda_p;
  mu = *mu_p;
  area = *area_p;
  strike = *strike_p;
  dip = *dip_p;
  u1 = *u1_p;
  u2 = *u2_p;
  ntarg = *ntarg_p;
  npts = *npts_p;
  
  
  for (itarg=0; itarg < ntarg; itarg++)
  {
    target[0][0] = targetx[itarg];
    target[0][1] = targety[itarg];
    target[0][2] = targetz[itarg];
    
    err = approxTgf(lambda, mu, triangle, target, 1, u1, u2, npts, nu, area,
                    strike, dip, u, du);
    
    ux[itarg] = u[0][0];
    uy[itarg] = u[0][1];
    uz[itarg] = u[0][2];
    uxx[itarg] = du[0][0][0];
    uxy[itarg] = du[0][0][1];
    uxz[itarg] = du[0][0][2];
    uyx[itarg] = du[0][1][0];
    uyy[itarg] = du[0][1][1];
    uyz[itarg] = du[0][1][2];
    uzx[itarg] = du[0][2][0];
    uzy[itarg] = du[0][2][1];
    uzz[itarg] = du[0][2][2];
        
  }
  
                  
  return;
}
