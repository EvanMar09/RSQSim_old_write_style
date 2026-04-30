#include "RSQSim.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif

/* for now, am assuming that the corner[0], center, L, W, strike, dip, and
   rake entries of all the KRDOkadaFaultPatchs are valid (and that the lower
   away-from-strike direction corner is corner[0]); and also nu and u 
   
   assumes that there are N KRDOkadaFaultPatch's at f + offset + i*size with
   i running from 0 to N-1 
   
   if allK is non-zero, will fill the K[6] vectors 
   
   only calcs entries for receiving patches startPatch through stopPatch
   
   patches whose centers are within dist1 (as a fraction of the source patch's
   size) will use the full analytic gf, those between dist1 and dist2 will
   use some number of point sources (12 for triangular source patches, for
   now), and those greater than dist2 will use some smaller number of 
   point sources (6 for triangular patches for now) - only implemented for
   triangular source patches so far.
*/

int calcStiffnessMatrix(void *f, size_t offset, size_t size, int N, double lambda, double mu,
                        int allK, int startPatch, int stopPatch, int pointSource, 
                        double dist1, double dist2)
{
  KRDOkadaFaultPatch *fi, *fj;
  
  int i, j, iK, icorner, icomp, jcomp;
  double x, y, z;  /* the center of the receiving patch, in
                     the translated coord. system */
  double xx, yy, zz;  /* the center of the receiving patch, in
                        the translated and rotated coord. system */
  double theta, cost, sint;  /* angle by which the coord system needs to be rotated
                               ccl about the +z axis to put the strike direction along
			       the +x axis.  And its cosine and sine */
  double dip;  /* need in degrees to pass to dc3d() */
  
  double alpha;  /* (lambda+mu)/(lambda+2*mu) */
  double c; /* depth of lower away-from-strike direction corner of source patch,
              need separate variable because Fortran passes by reference */
  double depth; /* depth of center of source patch */
  double u1, u2; /* strike-slip and dip-slip components of slip vector, again
                   need separate variable because Fortran passes by ref */
  double pot1, pot2; /* strike-slip and dip-slip potencies (for when using point sources), again
                   need separate variable because Fortran passes by ref */

  /* to hold the output from dc3d_() */
  double u[3], du[3][3];
    
  double Sigma[6]; /* stress tensor as a 6-vector in Okada coord. system */
  
  double K[6]; /* stiffness coeff for current pair of patches, ordered to match Sigma */
  
  double zero = 0.0;  /* not sure why this seems to be necessary */
  int izero = 0;
  int one = 1;
  
  /* for passing to tg_elh3dtriadirecttarg_() */
  double dummy2[1][3], dummy3[1][3][3];
  double corner[1][3][3], center[1][3], nu[1][3], disp[1][3]; 
/*
  double (*target)[3], (*strain)[3][3];
  double dist, (*uu)[3]; 
  */
	double target[1][3],strain[1][3][3],uu[1][3],dist;
  int np;
  
  np = stopPatch - startPatch + 1;
  /*target = xmalloc(1*3*sizeof(double));
  strain = xmalloc(1*3*3*sizeof(double));
  uu = xmalloc(1*3*sizeof(double));
  */
  
  alpha = (lambda + mu)/(lambda + 2*mu);
		   
#ifdef USE_OPENMP
#pragma omp parallel for  private(i,j,fi,fj,theta,cost,sint,dip,depth,c,pot1,pot2,u1,u2,x,y,z,xx,yy,zz,u,du,Sigma,K,iK, icomp,center,nu,disp,icorner,corner,dist,jcomp,strain,target,uu)
#endif

  for (i=0; i<N; i++) /* loop over source patches */
  {
    
    fi = f + offset + i*size;
    
    if (pointSource)
    {
      if (fi->type == TRIANGULAR)
      {
        fprintf(stderr, "calcStiffnessMatrix(): pointSource not yet supported for triangular elements\n");
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#else
        exit(-1);
#endif
      }
      
      theta = M_PI/2.0 - fi->strike;
      cost = cos(theta);
      sint = sin(theta);
      dip = fi->dip*180/M_PI;
      depth = -fi->center[2]; /* for point source, dc3d0() needs the depth of center of patch */
      c = -fi->corner[0][2];
      pot1 = fi->area * cos(fi->rake);
      pot2 = fi->area * sin(fi->rake);
      u1 = cos(fi->rake);
      u2 = sin(fi->rake);

      for (j=startPatch; j<=stopPatch; j++) /* loop over receiver patches */
      {
        fj = f + offset + j*size;
        
        if (j == i) /* need to use finite Greens functions for diagonal elements */
        {
          /* center of j^th patch in coord. system translated so that
             f[i].corner[0] is at x=0, y=0 (depth stays same) */
          x = fj->center[0] - fi->corner[0][0];
          y = fj->center[1] - fi->corner[0][1];
          z = fj->center[2];

          /* center of j^th patch in coord. system rotated so that 
             strike direction (of fault i) is along +x axis */
          xx = x*cost + y*sint;
          yy = -x*sint + y*cost;
          zz = z;

          dc3d(alpha, xx, yy, zz, c, dip, 0.0, fi->L, 0.0, fi->W,
                u1, u2, 0.0, 
	        &u[0], &u[1], &u[2],
	        &du[0][0], &du[1][0], &du[2][0], 
	        &du[0][1], &du[1][1], &du[2][1], 
	        &du[0][2], &du[1][2], &du[2][2]);
        } else   
        {     
          /* center of j^th patch in coord. system translated so that
             f[i].center is at x=0, y=0 (depth stays same) */
          x = fj->center[0] - fi->center[0];
          y = fj->center[1] - fi->center[1];
          z = fj->center[2];

          /* center of j^th patch in coord. system rotated so that 
             strike direction (of fault i) is along +x axis */
          xx = x*cost + y*sint;
          yy = -x*sint + y*cost;
          zz = z;
          dc3d0(alpha, xx, yy, zz, depth, dip, pot1, pot2, 0.0, 0.0,
 	       &u[0], &u[1], &u[2],
	       &du[0][0], &du[1][0], &du[2][0], 
	       &du[0][1], &du[1][1], &du[2][1], 
	       &du[0][2], &du[1][2], &du[2][2]);

        }
        
        deformationToStress(du, lambda, mu, Sigma);

        /* now need to rotate Sigma back to global coords. and store in f[j].K[i], this
           is just R^T * Sigma * R, which I've multiplied out by hand */
        K[0] = Sigma[0]*cost*cost - 2*Sigma[1]*sint*cost + Sigma[3]*sint*sint;
        K[1] = Sigma[0]*sint*cost + Sigma[1]*(cost*cost - sint*sint) - Sigma[3]*sint*cost;
        K[2] = Sigma[2]*cost - Sigma[4]*sint;
        K[3] = Sigma[0]*sint*sint + 2*Sigma[1]*sint*cost + Sigma[3]*cost*cost;
        K[4] = Sigma[2]*sint + Sigma[4]*cost;
        K[5] = Sigma[5];

        projectStress(K, fj->nu, fj->u, &fj->Ksigma[i], &fj->Ktau[i]);

        /* up until here, extension has been reckoned positive, but in the simulator
           code, compression is reckoned positive, so need to switch sign of Ksigma */
        fj->Ksigma[i] = -fj->Ksigma[i];

        if (allK != 0)
        {
          for (iK=0; iK<6; iK++) fj->K[i][iK] = K[iK];
        }

      } /*  for (j=0; j<N; j++)  */
    } else  /* e.g. not point source */
    {
      if (fi->type == TRIANGULAR)
      {

        for (icomp=0; icomp<3; icomp++)
        {
          center[0][icomp] = fi->center[icomp];
          nu[0][icomp] = fi->nu[icomp];
          disp[0][icomp] = fi->u[icomp];
          for (icorner=0; icorner<3; icorner++)
            corner[0][icorner][icomp] = fi->corner[icorner][icomp];
        }

//	fprintf(stderr," pass one %d %d \n",startPatch,stopPatch);

        for (j=startPatch; j<=stopPatch; j++) /* loop over receiver patches */
        {
          fj = f + offset + j*size;
          
          dist =  sqrt( (fi->center[0] - fj->center[0])*(fi->center[0] - fj->center[0]) +
                        (fi->center[1] - fj->center[1])*(fi->center[1] - fj->center[1]) +
                        (fi->center[2] - fj->center[2])*(fi->center[2] - fj->center[2]) )/sqrt(fi->area);

//		 fprintf(stderr," j= %d id = %d\n", j, omp_get_thread_num());

          if (dist <= dist1)
          {
            
            for (icomp=0; icomp<3; icomp++) target[0][icomp] = fj->center[icomp];
            /* use full triangular greens functions */
            tg_elh3dtriadirecttarg_(&lambda, &mu, corner, nu, &one, center,          /* elh3dtriadirecttarg_() actually returns strain instead of du */
                       &izero, dummy2, &one, disp, &izero, dummy2,       /* but this will still work when passed to deformationToStress */
                       &izero, dummy3, &one, target, &izero, dummy2,   /* because of symmetry of c_{ijkl}: */
                       &one, strain);

            for (icomp=0; icomp<3; icomp++)
              for (jcomp=0; jcomp<3; jcomp++)
                strain[0][icomp][jcomp] /= 4*M_PI;  /* actually tg_elh3dtriadirecttarg_() returns 4*pi*strain!!, argh! */

          } 
          else if (dist > dist1 && dist <= dist2)
          {
            for (icomp=0; icomp<3; icomp++) target[0][icomp] = fj->center[icomp];
            /* use 12 point sources for patches farther than dist2 */
            approxTgf(lambda, mu, corner[0], target, 1, cos(fi->rake), sin(fi->rake), 12,
                      fi->nu, fi->area, fi->strike, fi->dip, uu, strain);                                                

          } 
          else
          {
            for (icomp=0; icomp<3; icomp++) target[0][icomp] = fj->center[icomp];
            /* use 6 point sources for patches farther than dist2 */
            approxTgf(lambda, mu, corner[0], target, 1, cos(fi->rake), sin(fi->rake), 6,
                      fi->nu, fi->area, fi->strike, fi->dip, uu, strain);                                                

          }
          
          deformationToStress(strain[0], lambda, mu, K);
          projectStress(K, fj->nu, fj->u, &fj->Ksigma[i], &fj->Ktau[i]);

          /* up until here, extension has been reckoned positive, but in the simulator
             code, compression is reckoned positive, so need to switch sign of Ksigma */
          fj->Ksigma[i] = -fj->Ksigma[i];

          if (allK != 0)
          {
            for (iK=0; iK<6; iK++) fj->K[i][iK] = K[iK];
          }
          
        } /* loop over receiver patches */
 //               fprintf(stderr," pass two \n");
		
        
      } /* if fi->type == TRIANGULAR */
      else   /* else use Okada, 1992 */
      {
        theta = M_PI/2.0 - fi->strike;
        cost = cos(theta);
        sint = sin(theta);
        dip = fi->dip*180/M_PI;
        c = -fi->corner[0][2];
        u1 = cos(fi->rake);
        u2 = sin(fi->rake);


        for (j=startPatch; j<=stopPatch; j++) /* loop over receiver patches */
        {
          fj = f + offset + j*size;
          /* center of j^th patch in coord. system translated so that
             f[i].corner[0] is at x=0, y=0 (depth stays same) */
          x = fj->center[0] - fi->corner[0][0];
          y = fj->center[1] - fi->corner[0][1];
          z = fj->center[2];

          /* center of j^th patch in coord. system rotated so that 
             strike direction (of fault i) is along +x axis */
          xx = x*cost + y*sint;
          yy = -x*sint + y*cost;
          zz = z;

          dc3d(alpha, xx, yy, zz, c, dip, 0.0, fi->L, 0.0, fi->W,
                u1, u2, 0.0, 
	        &u[0], &u[1], &u[2],
	        &du[0][0], &du[1][0], &du[2][0], 
	        &du[0][1], &du[1][1], &du[2][1], 
	        &du[0][2], &du[1][2], &du[2][2]);

          deformationToStress(du, lambda, mu, Sigma);

          /* now need to rotate Sigma back to global coords. and store in f[j].K[i], this
             is just R^T * Sigma * R, which I've multiplied out by hand */
          K[0] = Sigma[0]*cost*cost - 2*Sigma[1]*sint*cost + Sigma[3]*sint*sint;
          K[1] = Sigma[0]*sint*cost + Sigma[1]*(cost*cost - sint*sint) - Sigma[3]*sint*cost;
          K[2] = Sigma[2]*cost - Sigma[4]*sint;
          K[3] = Sigma[0]*sint*sint + 2*Sigma[1]*sint*cost + Sigma[3]*cost*cost;
          K[4] = Sigma[2]*sint + Sigma[4]*cost;
          K[5] = Sigma[5];

          projectStress(K, fj->nu, fj->u, &fj->Ksigma[i], &fj->Ktau[i]);

          /* up until here, extension has been reckoned positive, but in the simulator
             code, compression is reckoned positive, so need to switch sign of Ksigma */
          fj->Ksigma[i] = -fj->Ksigma[i];

          if (allK != 0)
          {
            for (iK=0; iK<6; iK++) fj->K[i][iK] = K[iK];
          }

        } /*  for (j=0; j<N; j++)  */

      } /* fi->type == OKADA */
    } /* for (i=startPatch; i <= stopPatch; i++) */
  }
  
/*  free(target);
  free(strain);
  free(uu);
  */

  return(0);
}
