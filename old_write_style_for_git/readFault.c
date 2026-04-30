#include "RSQSim.h"

static void writeFaultToOutputFile(FaultModel *m, Params *params);

void readFault(FaultModel *m, Params *params)
{
  FILE *faultFp;
  int nOKADA = 0, nTRIANGULAR = 0;
  int i;
  char tmpLine[RSQSIM_MAXLINE];
  double z;

  
if (m->master)  writeFaultToOutputFile(m, params);

  faultFp = xfopen(params->faultFname, "r");

  for (i=0; i<m->N; i++)
  {
    fgets(tmpLine, RSQSIM_MAXLINE, faultFp);
    if (tmpLine[0] == '#') {i--; continue;}  /* comment line in fault file */

    if (sscanf(tmpLine,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &z, &z, &z, &z, &z, &z, &z, &z, &z, &z, &z)  == 9 ||   /* 9 or 10 because it might or might not have the optional segment number */
        sscanf(tmpLine,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &z, &z, &z, &z, &z, &z, &z, &z, &z, &z, &z)  == 10) /* this is an OKADA patch */
    {
      m->p[i].op.type = OKADA;
      sscanf(tmpLine, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
             &(m->p[i].op.center[0]), &(m->p[i].op.center[1]), &(m->p[i].op.center[2]),
             &(m->p[i].op.L), &(m->p[i].op.W),
             &(m->p[i].op.strike), &(m->p[i].op.dip), &(m->p[i].op.rake),
             &(m->p[i].dDotDrive)); /* note that dDotDrive might actually be tauDot depending on 
                                       params->stressRateSpecification */
      /* convert strike, dip, and rake to radians */
      m->p[i].op.strike *= M_PI/180.0;
      m->p[i].op.dip *= M_PI/180.0;
      m->p[i].op.rake *= M_PI/180.0;
      m->p[i].op.area = m->p[i].op.L*m->p[i].op.W;
      findCorners(&(m->p[i].op));
      findNormalAndSlip(&(m->p[i].op));
      nOKADA++;
//	if (m->master) fprintf(stdout," readFault i=%d dDot =%e %lf \n", i,m->p[i].dDotDrive,m->p[i].op.strike);
    }
    else if (sscanf(tmpLine,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &z, &z, &z, &z, &z, &z, &z, &z, &z, &z, &z)  ==  11) /* this is a TRIANGULAR patch */
    {
      m->p[i].op.type = TRIANGULAR;
      sscanf(tmpLine, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
             &(m->p[i].op.corner[0][0]), &(m->p[i].op.corner[0][1]), &(m->p[i].op.corner[0][2]),
             &(m->p[i].op.corner[1][0]), &(m->p[i].op.corner[1][1]), &(m->p[i].op.corner[1][2]),
             &(m->p[i].op.corner[2][0]), &(m->p[i].op.corner[2][1]), &(m->p[i].op.corner[2][2]),
             &(m->p[i].op.rake), &(m->p[i].dDotDrive)); /* see note above about dDotDrive */
      m->p[i].op.rake *= M_PI/180.0;
      findCenter(&(m->p[i].op));
      findNormalAndSlip(&(m->p[i].op));
      nTRIANGULAR++;
    }
    else
    {
      fprintf(stderr, "\nError: initFaultModel(): unrecognized fault patch line: \n%s\n", tmpLine);
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
      exit(-1);
    }
  }

  return;
}

static void writeFaultToOutputFile(FaultModel *m, Params *params)
{
  FILE *faultfp;
  char fline[RSQSIM_MAXLINE];
  
  faultfp = xfopen(params->faultFname, "r");
   
  fprintf(m->eqsfp, "%%%% faultFname = %s\n", params->faultFname); 
  while (fgets(fline, RSQSIM_MAXLINE, faultfp) != NULL)
  { 
    fprintf(m->eqsfp, "%% ");
    fputs(fline, m->eqsfp);
  }
  
  fprintf(m->eqsfp, "%%%% end faultFile\n");
  
  
  fclose(faultfp);
  
  return;
}

