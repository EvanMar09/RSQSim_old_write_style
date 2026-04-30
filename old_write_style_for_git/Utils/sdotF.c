#include "../RSQSim.h"
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

int main(int argc, char **argv)
{
  FaultModel m;
  Params params;

  char comment[] = {"#"}; /* comment char in faultFname */
  int i, icorner, icomp, jcomp;
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
  double u1, u2; /* strike-slip and dip-slip components of slip vector, again
                   need separate variable because Fortran passes by ref */

  /* to hold the output from dc3d_() */
  double u[3], du[3][3];
    
  double Sigma[6]; /* stress tensor as a 6-vector in Okada coord. system */
  
  double K[6]; /* stiffness coeff for current pair of patches, ordered to match Sigma */
  
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
  
  double pt[3], sdotF[6];
  double lambda, mu, dist1, dist2;
  void *f;
  size_t offset, size;
  
  FILE *fpfp;
  
  static char *usage[] =
    {"sdotF faultFname fieldPtsFname > sdotF.txt",
     " ",
     "faultFname is a usual RSQSim fault file",
     "fieldPtsFname is a 3-column ascii file giving the x, y, z coords of the field points",
     " ",
     "outputs a 6-column ascii file giving the 6 components of the long-term",
     "stressing rate tensor (in the usual order) at each of the field points",
     "(in MPa/s)",
     ""};
  
  if (argc == 1)
  {
    i = 0;
    while(usage[i][0] != '\0') fprintf(stderr, "%s\n", usage[i++]);
    
    return(0);
  }
  
  params.faultFname = argv[1];
  m.N = nlines(params.faultFname, comment);
  m.startPatch = -1;
  m.stopPatch = -1;
  allocFaultModel(&m);
  readFault(&m, &params);
  
  lambda = 30e3;
  mu = 30e3;
  dist1 = 20;
  dist2 = 100;
  f = &(m.p);
  offset = offsetof(FaultPatch, op);
  size = sizeof(FaultPatch);
  
  KRDOkadaFaultPatch *fi;
  
  
  
  alpha = (lambda + mu)/(lambda + 2*mu);
		   
#ifdef USE_OPENMP
#pragma omp parallel for  private(i,j,fi,fj,theta,cost,sint,dip,depth,c,pot1,pot2,u1,u2,x,y,z,xx,yy,zz,u,du,Sigma,K,iK, icomp,center,nu,disp,icorner,corner,dist,jcomp,strain,target,uu)
#endif

  fpfp = xfopen(argv[2], "r");
  while (fscanf(fpfp, "%lf %lf %lf", &(pt[0]), &(pt[1]), &(pt[2])) == 3)
  {
    for (icomp=0; icomp<6; icomp++) sdotF[icomp] = 0.0;
    
    for (i=0; i<m.N; i++) /* loop over source patches */
    {

/*      fi = f + offset + i*size; */
      fi = &(m.p[i].op);

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

        dist =  sqrt( (fi->center[0] - pt[0])*(fi->center[0] - pt[0]) +
                      (fi->center[1] - pt[1])*(fi->center[1] - pt[1]) +
                      (fi->center[2] - pt[2])*(fi->center[2] - pt[2]) )/sqrt(fi->area);


        if (dist <= dist1)
        {

          for (icomp=0; icomp<3; icomp++) target[0][icomp] = pt[icomp];
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
          for (icomp=0; icomp<3; icomp++) target[0][icomp] = pt[icomp];
          /* use 12 point sources for patches farther than dist2 */
          approxTgf(lambda, mu, corner[0], target, 1, cos(fi->rake), sin(fi->rake), 12,
                    fi->nu, fi->area, fi->strike, fi->dip, uu, strain);                                                

        } 
        else
        {
          for (icomp=0; icomp<3; icomp++) target[0][icomp] = pt[icomp];
          /* use 6 point sources for patches farther than dist2 */
          approxTgf(lambda, mu, corner[0], target, 1, cos(fi->rake), sin(fi->rake), 6,
                    fi->nu, fi->area, fi->strike, fi->dip, uu, strain);                                                

        }

        deformationToStress(strain[0], lambda, mu, K);

        for (icomp=0; icomp<6; icomp++) sdotF[icomp] += K[icomp]*m.p[i].dDotDrive;

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


         /* center of j^th patch in coord. system translated so that
            f[i].corner[0] is at x=0, y=0 (depth stays same) */
         x = pt[0] - fi->corner[0][0];
         y = pt[1] - fi->corner[0][1];
         z = pt[2];

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

         for (icomp=0; icomp<6; icomp++) sdotF[icomp] += K[icomp]*m.p[i].dDotDrive;


      } /* fi->type == OKADA */
    } /* for (i=startPatch; i <= stopPatch; i++) */
    fprintf(stdout, "%lg %lg %lg %lg %lg %lg\n", sdotF[0], sdotF[1], sdotF[2], sdotF[3], sdotF[4], sdotF[5]);
   
  } /* while fscanf(fpfp, "%lf %lf %lf", &pt[0], &pt[1], &pt[2]) == 3) */
  

  return(0);
}
