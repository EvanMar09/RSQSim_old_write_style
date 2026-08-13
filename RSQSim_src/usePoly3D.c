#include "KRDOkada.h"

#ifdef USEPOLY3D
void usePoly3D(int i, void *f, size_t offset, size_t size, int N, double lambda, double mu)
{

  KRDOkadaFaultPatch *fi, *fj;
  char *tmpfname *tmpPoly3doutfname;
  FILE *fp;
  int j;
  char command[RSQSIM_MAXLINE];

  fi = f + offset + i*size;
  
  tmpfname = tempnam(".", NULL);
  fp = fopen(tmpfname, "w");
  
  /* write out poly3d input file for the i^th element as a source */
  /* Section 1 of poly3d input file */
  fprintf(fp, "shear_mod   =  %.10lg\n", mu*1e6);      /* for Poly3D, everything will be in Pa and m */
  fprintf(fp, "lame_lambda =  %.10lg\n", lambda*1e6);
  fprintf(fp, "end\n");
  /* Section 2 of poly3d input file */
  fprintf(tmpfp, "end\n");
  /* Section 3 of poly3d input file */
  for (j=0; j<N; j++)
  {
    fj = f + offset + j*size;
    fprintf(fp, "obs 0 s global global global %.10lg %.10lg %.10lg\n", 
            fj->center[0], fj->center[1], fj->center[2]);
  }
  fprintf(fp, "end\n");
  /* Section 4 of poly3d input file */
  for (k=0; k<3; k++) 
    fprintf(fp, "v %d global %.10lg %.10lg %.10lg\n", 
            fi->corner[k][0], fi->corner[k][1], fi->corner[k][2]);
  fprintf(fp, "e 3 global bbb %.10lg %.10lg %.10lg\n", 
          fi->u[0], fi->u[1], fi->u[2]);
  fprintf(fp, "end\n");
  fclose(fp);
  
  /* run poly3d */
  tmpPoly3doutfname = tempnam(".", NULL);
  sprintf(command, "poly3d -i %s -o %s", tmpfname, tmpPoly3doutfname);
  
  
  


}
#endif
