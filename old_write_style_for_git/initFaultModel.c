#include "RSQSim.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif

static char svnId[] = {"$Id: initFaultModel.c 1136 2017-06-14 23:29:05Z dinger $"};

void initFaultModel(FaultModel *m, Params *params)
{  
 
  int i, j, iwhen, jwhen, iproc;
  FILE *tmpfp, *stressRateFp;
  FILE *Fps[NFPS], *kzfp;
  char tmpFname[RSQSIM_MAXLINE], tmpLine[RSQSIM_MAXLINE];
  MPI_Status status;
  
#ifdef USE_MPI  
  void *mbuffer;
  int bufsize, root;
#endif

  double stressRate[6];
  
  double KCoul, KCoulSelf;

#ifndef USE_MPI
  int ip;
#endif 
  
  allocFaultModel(m);
//  if (m->master)
  if(1) 
  {
   if (m->master)   fprintf(m->eqsfp, "%%%%%% begin input files\n");
    
    readFault(m, params);
    xferParamsToM(m, params);
    allocSnapBuffers(m);

    openFiles(m, params, Fps);
    
    if (params->stressRateSpecification == TENSOR_RATE)
    {
      stressRateFp = xfopen(params->stressRateFname, "r");
      if (m->master)
      {   
        fprintf(m->eqsfp, "%%%% stressRateFname = %s\n", params->stressRateFname);
        while (fgets(tmpLine, RSQSIM_MAXLINE, stressRateFp) != NULL)
        {
          fprintf(m->eqsfp, "%% ");
          fputs(tmpLine, m->eqsfp);
        }
        if (m->master) fprintf(m->eqsfp, "%%%% end stressRateFile\n");
        fseek(stressRateFp, 0L, SEEK_SET);
      }
      
      fscanf(stressRateFp, "%lf %lf %lf %lf %lf %lf", &(stressRate[0]), &(stressRate[1]),
             &(stressRate[2]), &(stressRate[3]), &(stressRate[4]), &(stressRate[5]));
      fclose(stressRateFp);
    }

    /* initialize the patches in the FaultModel */
    m->state3patches = FALSE;
    m->maxNeighbors = 0;

    for (i=0; i<m->N; i++)
      initializeFaultPatch(m, params, Fps, i, stressRate);
    
    closeFiles(Fps, NFPS);
     if (m->master) fprintf(m->eqsfp, "%%%%%% end input files\n");
    
#ifdef USE_MPI_XX
    packFaultModel(m, &mbuffer, &bufsize); /* pack the relevant parts of the fault model up for MPI_Bcast() */
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(mbuffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);
    free(mbuffer);
#endif
  } /* if (m->master) */
  else
  {
#ifdef USE_MPI_XX   
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    mbuffer = xmalloc(bufsize);
    MPI_Bcast(mbuffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);
    unpackFaultModel(m, mbuffer, bufsize); /* unpack fault model */
    free(mbuffer);
#endif
  }
 
  
  /* calculate or read stiffness matrix */
  if (strlen(params->KTauFname) == 0 || strlen(params->KSigmaFname) == 0)
  {
    if (params->nEq == 0 && params->writeStiffness == 0)
    {
      fprintf(stderr, "Warning: nEq = 0 and writeStiffness = 0, so not calculating stiffness matrix\n");
    }
    else
    {
      fprintf(stderr, "Calculating stiffness matrix ...");
      calcStiffnessMatrix(m->p, offsetof(FaultPatch, op), sizeof(FaultPatch),
                          m->N, params->lameLambda, params->lameMu, ALLK, 
                          m->startPatch, m->stopPatch, params->pointSource, 
                          params->tgfDist1, params->tgfDist2);
      fprintf(stderr, "done\n");
    }
  } 
  else  
  {
    
    tmpfp = xfopen(params->KTauFname, "r");
    fseek(tmpfp, m->startPatch*m->N*sizeof(double), SEEK_SET);
    for (i = m->startPatch; i <= m->stopPatch; i++)
      fread(m->p[i].op.Ktau, sizeof(double), m->N, tmpfp);

    fclose(tmpfp);
    
    tmpfp = xfopen(params->KSigmaFname, "r");
    fseek(tmpfp, m->startPatch*m->N*sizeof(double), SEEK_SET);
    for (i = m->startPatch; i <= m->stopPatch; i++)
      fread(m->p[i].op.Ksigma, sizeof(double), m->N, tmpfp);

    fclose(tmpfp);
  } 
  
  
  /* write out stiffness matrices, if desired.  */
  if (params->writeStiffness == KTAU || params->writeStiffness == BOTH)
  {
    if (m->master)
    {
      fprintf(stderr, "writing out shear stress stiffness matrix ... ");
      sprintf(tmpFname, "Ktau.%s.out", params->outFnameInfix);
      
      /* master needs to write its part first */
      tmpfp = xfopen(tmpFname, "w");
      for (i = m->startPatch; i <= m->stopPatch; i++)
        fwrite(m->p[i].op.Ktau, sizeof(double), m->N, tmpfp);
      fclose(tmpfp);

      for (iproc=1; iproc < m->nProc; iproc++) /* let one process at a time write out its part */
      {
        if (m->DEBUG > 5) {fprintf(stderr, "iproc = %d\n", iproc); fflush(stderr);}
        MPI_Send(NULL, 0, MPI_INT, iproc, 0, MPI_COMM_WORLD);
	MPI_Recv(NULL, 0, MPI_INT, iproc, 0, MPI_COMM_WORLD, &status);
      }
      fprintf(stderr, "done\n");
      
      
    } else /* slaves wait to receive a message and then write out their part */
    {
      MPI_Recv(NULL, 0, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      sprintf(tmpFname, "Ktau.%s.out", params->outFnameInfix);
      tmpfp = xfopen(tmpFname, "a");
      for (i = m->startPatch; i <= m->stopPatch; i++)
        fwrite(m->p[i].op.Ktau, sizeof(double), m->N, tmpfp);
      fclose(tmpfp);
      MPI_Send(NULL, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

  }

  if (params->writeStiffness == KSIGMA || params->writeStiffness == BOTH)
  {
    if (m->master)
    {
      fprintf(stderr, "writing out normal stress stiffness matrix ... ");
      sprintf(tmpFname, "Ksigma.%s.out", params->outFnameInfix);
     
      /* master needs to write its part out first */ 
      tmpfp = xfopen(tmpFname, "w");
      for (i = m->startPatch; i <= m->stopPatch; i++)
        fwrite(m->p[i].op.Ksigma, sizeof(double), m->N, tmpfp);
      fclose(tmpfp);

      for (iproc=1; iproc < m->nProc; iproc++) /* let one process at a time write out its part */
      {
        if (m->DEBUG > 5) {fprintf(stderr, "iproc = %d\n", iproc); fflush(stderr);}
        MPI_Send(NULL, 0, MPI_INT, iproc, 0, MPI_COMM_WORLD);
	MPI_Recv(NULL, 0, MPI_INT, iproc, 0, MPI_COMM_WORLD, &status);
      }
      fprintf(stderr, "done\n");
      
      
    } else /* slaves wait to receive a message and then write out their part */
    {
      MPI_Recv(NULL, 0, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      sprintf(tmpFname, "Ksigma.%s.out", params->outFnameInfix);
      tmpfp = xfopen(tmpFname, "a");
      for (i = m->startPatch; i <= m->stopPatch; i++)
        fwrite(m->p[i].op.Ksigma, sizeof(double), m->N, tmpfp);
      fclose(tmpfp);
      MPI_Send(NULL, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

  } 

  /* set Kself */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i)
#endif
  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    m->p[i].Kself = m->p[i].op.Ktau[i];
  } 
  
 /* possibly zero out elements of the stiffness matrices that seem pathologically large */
  if (params->KZeroFrac != 0.0)
  {
    for (i=m->startPatch; i<= m->stopPatch; i++)
    {
      KCoulSelf = fabs(m->p[i].op.Ktau[i]) + fabs(m->p[i].mu0*m->p[i].op.Ksigma[i]);
      for (j=0; j<m->N; j++)
      {
        if (j == i) continue;

        KCoul = fabs(m->p[i].op.Ktau[j]) + fabs(m->p[i].mu0*m->p[i].op.Ksigma[j]);
        if (fabs(KCoul) > fabs(params->KZeroFrac*KCoulSelf))
        {
          m->p[i].op.Ktau[j] = 0;
          m->p[i].op.Ksigma[j] = 0;
          if (params->DEBUG) fprintf(stderr, "initFaultModel(): setting K[ %d ][ %d ] to zero\n", i, j);
        }
      }
    }
  }
  
  /* zero out elements of the stiffness matrices specified in KZeroFile */
  if (strlen(params->KZeroFname) > 0)
  {
    FILE *kzfp;
    int jbad, ibad;
    
    kzfp = xfopen(params->KZeroFname, "r");
    while(fscanf(kzfp, "%d %d", &jbad, &ibad) == 2)
    {
      jbad--;  ibad--;  /* patch #s in input and output files are 1-based */
      if (ibad >= m->startPatch && ibad <= m->stopPatch)
      {
        m->p[ibad].op.Ktau[jbad] = 0;
        m->p[ibad].op.Ksigma[jbad] = 0;
      }
    }
    fclose(kzfp);
  }

  /* calc driving stressing rates according to whether tauDotFname and/or sigmaDotFname were supplied */
  if (params->stressRateSpecification == BACKSLIP) 
  {
    fprintf(stderr, "calculating backslip stressing rates... ");
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i,j)
#endif
    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      m->p[i].tauDotDrive = 0.0;
      m->p[i].sigmaDotDrive = 0.0;
      for (j=0; j<m->N; j++)
      {
        m->p[i].tauDotDrive += -m->p[i].op.Ktau[j]*m->p[j].dDotDrive;
        m->p[i].sigmaDotDrive += -m->p[i].op.Ksigma[j]*m->p[j].dDotDrive;
      }
      m->p[i].tauDot = m->p[i].tauDotDrive;
      m->p[i].sigmaDot = m->p[i].sigmaDotDrive;


    }
    
            
    fprintf(stderr, "done\n");
  }
  else /* stressing rates specified directly or from tensor stress rate; just initialize tauDot and sigmaDot */
  {
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i)
#endif
    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      m->p[i].tauDot = m->p[i].tauDotDrive;
      m->p[i].sigmaDot = m->p[i].sigmaDotDrive;
      m->p[i].dDotDrive = NAN; /* dDotDrive not applicable when specifying stressing rates directly */
    }
  }
  
  /* set initial tau on state 3 patches to that value which will make each
     state 3 patch slip at its long-term average slip rate.  Only do this
     if stressRateSpecification is BACKSLIP and if any of tau, sigma, 
     theta, or ddot were not read in from a file */
  if (params->stressRateSpecification == BACKSLIP && 
      (Fps[INITTAU] == NULL || Fps[INITSIGMA] == NULL || Fps[INITTHETA] == NULL || Fps[INITSLIPSPEED] == NULL))
  {
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i)
#endif
    for (i=m->startPatch; i <= m->stopPatch; i++)
    {
      if (m->p[i].state == 3)
      {
        m->p[i].tau = MAX(0.0, m->p[i].sigma*(m->p[i].mu0 + (m->p[i].A - m->p[i].B)*log(m->p[i].dDotDrive/m->p[i].ddotStar)));
      } /* if (m->p[i].state == 3) */
    } /* for (i=m->startPatch; i <= m->stopPatch; i++) */
  } /* if (params->stressRateSpecification == BACKSLIP) */

  /* now set slip speeds for any state 3 patches (if initial slip speeds weren't
     read in from a file), and then update stressing rates */
  if (Fps[INITSLIPSPEED] == NULL)
  {
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i)
#endif
    for (i = m->startPatch; i <= m->stopPatch; i++)
      if (m->p[i].state == 3) setDdot3(&(m->p[i]));
  }
  
  /* need to send all slip speeds to all processes 
     FIXME:  should only send slip speeds of state 3 patches
             and then also use only those elements of the
             stiffness matrix that correspond to state 3 patches 
             More importantly, should also do this in updateEqs() */
#ifdef USE_MPI
  /* can't get MPI_Allgatherv() to work so doing a series of MPI_Gatherv()s for now 
  MPI_Allgatherv(&(m->p[m->startPatch].ddot), 1, m->SendFaultDouble,
                 m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
                 MPI_COMM_WORLD); */
  for (root = 0; root < m->nProc; root++)
    MPI_Gatherv(&(m->p[m->startPatch].ddot), 1, m->SendFaultDouble,
                   m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
                   root, MPI_COMM_WORLD);
#else
  for (ip = 0;  ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].ddot;
#endif

  /* adjust initial tauDot and sigmaDot for state 3 creeping patches */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i,j)
#endif 
 for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    for (j = 0; j < m->N; j++)
    {
      if (m->p[j].state == 3)
      {
        m->p[i].tauDot += m->p[i].op.Ktau[j]*m->allFaultPropDouble[j];
        m->p[i].sigmaDot += m->p[i].op.Ksigma[j]*m->allFaultPropDouble[j];
      }
    }
  }
  
  /* adjust initial tauDot and sigmaDot for external stressing rates */
  if (m->iExtStressRateStep > 0)
  {
    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      m->p[i].tauDot += m->p[i].tauDotExt[m->iExtStressRateStep - 1];
      m->p[i].sigmaDot += m->p[i].sigmaDotExt[m->iExtStressRateStep - 1];
    }
  }
  
  /* initialize muTrans3* on state 3 patches */
#ifdef USE_OPENMP
#pragma omp parallel for shared(m) private(i)
#endif
  for (i = m->startPatch; i <= m->stopPatch; i++)
  {
    if (m->p[i].state == 3)
    {
      m->p[i].muTrans3Low = m->p[i].tau/m->p[i].sigma - m->dMu3;
      m->p[i].muTrans3High = m->p[i].tau/m->p[i].sigma + m->dMu3;
    }
  }

  
  
  /* write out tectonic tauDot and sigmaDot */
#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].tauDotDrive), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (i=0; i<m->N; i++) m->allFaultPropDouble[i] = m->p[i].tauDotDrive;
#endif

  if (m->master)
  {
    sprintf(tmpFname, "tauDot.%s.out", params->outFnameInfix);
    tmpfp = xfopen(tmpFname, "w");
    for (i=0; i<m->N; i++) fprintf(tmpfp, "%lg %lg %lg %lg\n", 
                                  m->p[i].op.center[0], m->p[i].op.center[1], 
                                  m->p[i].op.center[2], m->allFaultPropDouble[i]);
    fclose(tmpfp);
  }
  

#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].sigmaDotDrive), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (i=0; i<m->N; i++) m->allFaultPropDouble[i] = m->p[i].sigmaDot;
#endif

  if (m->master)
  {
    sprintf(tmpFname, "sigmaDot.%s.out", params->outFnameInfix);
    tmpfp = xfopen(tmpFname, "w");
    for (i=0; i<m->N; i++) fprintf(tmpfp, "%lg %lg %lg %lg\n", 
                                  m->p[i].op.center[0], m->p[i].op.center[1], 
                                  m->p[i].op.center[2], m->allFaultPropDouble[i]);
    fclose(tmpfp);
  }
  
#ifdef USE_SERIAL_IO
 if (m->master)
#else
 if(1)
#endif
  {

    if (params->writePED)
    {
#ifdef USE_SERIAL_IO
      /* open eList, pList, and dList files */
      
      sprintf(tmpFname, "%s.eList", params->outFnameInfix);
      m->eListfp = xfopen(tmpFname, "w");  
      sprintf(tmpFname, "%s.pList", params->outFnameInfix);
      m->pListfp = xfopen(tmpFname, "w");  
      sprintf(tmpFname, "%s.dList", params->outFnameInfix);
      m->dListfp = xfopen(tmpFname, "w");  
      sprintf(tmpFname, "%s.dtauList", params->outFnameInfix);
      m->dtauListfp = xfopen(tmpFname, "w");  
      sprintf(tmpFname, "%s.dsigmaList", params->outFnameInfix);
      m->dsigmaListfp = xfopen(tmpFname, "w");
      sprintf(tmpFname, "%s.dmuList", params->outFnameInfix);
      m->dmuListfp = xfopen(tmpFname, "w");

      sprintf(tmpFname, "%s.tList", params->outFnameInfix);
      m->tListfp = xfopen(tmpFname, "w");  
      sprintf(tmpFname, "%s.taupList", params->outFnameInfix);
      m->taupListfp = xfopen(tmpFname, "w");  
      sprintf(tmpFname, "%s.sigmapList", params->outFnameInfix);
      m->sigmapListfp = xfopen(tmpFname, "w");  
      
#else
      MPI_File eFile,dFile,pFile,tauFile,sigmaFile,muFile,tFile,taupFile;
      /* open eList, pList, and dList files */
      sprintf(tmpFname, "%s.eList", params->outFnameInfix);
      //m->eListfp = xfopen(tmpFname, "w");

        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&eFile);
        m->eListfp = eFile;


      sprintf(tmpFname, "%s.pList", params->outFnameInfix);
      //m->pListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&pFile);
        m->pListfp = pFile;

      sprintf(tmpFname, "%s.dList", params->outFnameInfix);
      //m->dListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&dFile);
        m->dListfp = dFile;

      sprintf(tmpFname, "%s.dtauList", params->outFnameInfix);
      //m->dtauListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&tauFile);
        m->dtauListfp = tauFile;

      sprintf(tmpFname, "%s.dsigmaList", params->outFnameInfix);
      //m->dsigmaListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&sigmaFile);
        m->dsigmaListfp = sigmaFile;

      sprintf(tmpFname, "%s.dmuList", params->outFnameInfix);
      //m->dmuListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&muFile);
        m->dmuListfp = muFile;

      sprintf(tmpFname, "%s.tList", params->outFnameInfix);
      //m->tListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&tFile);
        m->tListfp = tFile;

      sprintf(tmpFname, "%s.taupList", params->outFnameInfix);
      //m->taupListfp = xfopen(tmpFname, "w");
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&taupFile);
        m->taupListfp = taupFile;
#endif
    }
    
#ifdef USE_SERIAL_IO
    /* open desired tau, sigma, slip, state output files */
    for (iwhen=0; iwhen<3; iwhen++)
    {
      jwhen = rint(pow(2, iwhen));
      if ((params->writeTau & jwhen) > 0)
      {
        sprintf(tmpFname, "tau.%s.out.%d", params->outFnameInfix, jwhen);
        m->taufp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&(m->t), sizeof(m->t), 1, m->taufp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].tau), sizeof(m->p[i].tau), 1, m->taufp[iwhen]);
      }
      else
      {
        m->taufp[iwhen] = NULL;
      }
      /* always write out tau at end of run */
      sprintf(tmpFname, "tau.%s.out.final", params->outFnameInfix);
      m->taufp[3] = xfopen(tmpFname, "w");
      
      /* if writing out sigma, open file here */
      if ((params->writeSigma & jwhen) > 0)
      {
        sprintf(tmpFname, "sigma.%s.out.%d", params->outFnameInfix, jwhen);
        m->sigmafp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&(m->t), sizeof(m->t), 1, m->sigmafp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].sigma), sizeof(m->p[i].sigma), 1, m->sigmafp[iwhen]);
      }
      else
      {
        m->sigmafp[iwhen] = NULL;
      }
      /* always write out sigma at end of run */
      sprintf(tmpFname, "sigma.%s.out.final", params->outFnameInfix);
      m->sigmafp[3] = xfopen(tmpFname, "w");

      /* if writing out slip, open file here */
      if ((params->writeSlip & jwhen) > 0)
      {
        sprintf(tmpFname, "slip.%s.out.%d", params->outFnameInfix, jwhen);
        m->slipfp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&(m->t), sizeof(m->t), 1, m->slipfp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].d), sizeof(m->p[i].d), 1, m->slipfp[iwhen]);
      }
      else
      {
        m->slipfp[iwhen] = NULL;
      }
      /* always write out slip at end of run */
      sprintf(tmpFname, "slip.%s.out.final", params->outFnameInfix);
      m->slipfp[3] = xfopen(tmpFname, "w");

      /* if writing out slip speed, open file here */
      if ((params->writeSlipSpeed & jwhen) > 0)
      {
        sprintf(tmpFname, "slipSpeed.%s.out.%d", params->outFnameInfix, jwhen);
        m->slipSpeedfp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&(m->t), sizeof(m->t), 1, m->slipSpeedfp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].ddot), sizeof(m->p[i].ddot), 1, m->slipSpeedfp[iwhen]);
      }
      else
      {
        m->slipSpeedfp[iwhen] = NULL;
      }
      /* always write out slipSpeed at end of run */
      sprintf(tmpFname, "slipSpeed.%s.out.final", params->outFnameInfix);
      m->slipSpeedfp[3] = xfopen(tmpFname, "w");

      /* if writing out state, open file here */
      if ((params->writeState & jwhen) > 0)
      {
        sprintf(tmpFname, "state.%s.out.%d", params->outFnameInfix, jwhen);
        m->statefp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&(m->t), sizeof(m->t), 1, m->statefp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].state), sizeof(m->p[i].state), 1, m->statefp[iwhen]);
      }
      else
      {
        m->statefp[iwhen] = NULL;
      }
      /* always write out state at end of run */
      sprintf(tmpFname, "state.%s.out.final", params->outFnameInfix);
      m->statefp[3] = xfopen(tmpFname, "w");

      /* if writing out theta, open file here */
      if ((params->writeTheta & jwhen) > 0)
      {
        sprintf(tmpFname, "theta.%s.out.%d", params->outFnameInfix, jwhen);
        m->thetafp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&(m->t), sizeof(m->t), 1, m->thetafp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].theta), sizeof(m->p[i].theta), 1, m->thetafp[iwhen]);
      }
      else
      {
        m->thetafp[iwhen] = NULL;
      }
      /* always write out theta at end of run */
      sprintf(tmpFname, "theta.%s.out.final", params->outFnameInfix);
      m->thetafp[3] = xfopen(tmpFname, "w");
    }
#else
MPI_File File_tau_trans[4];
MPI_File File_sigma_trans[4];
MPI_File File_slip_trans[4];
MPI_File File_ddot_trans[4];
MPI_File File_theta_trans[4];
MPI_File File_state_trans[4];

    /* open desired tau, sigma, slip, state output files */
    for (iwhen=0; iwhen<3; iwhen++)
    {
      jwhen = rint(pow(2, iwhen));
      if ((params->writeTau & jwhen) > 0)
      {
        sprintf(tmpFname, "tau.%s.out.%d", params->outFnameInfix, jwhen);
      /*  m->taufp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&t, sizeof(t), 1, m->taufp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].tau), sizeof(m->p[i].tau), 1, m->taufp[iwhen]);
        */
//      MPI_File File_tau_trans[3];
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_tau_trans[iwhen]);
        m->taufp[iwhen]=File_tau_trans[iwhen];

      }
      else
      {
        m->taufp[iwhen] = NULL;
      }
      /* if writing out sigma, open file here */
      if ((params->writeSigma & jwhen) > 0)
      {
        sprintf(tmpFname, "sigma.%s.out.%d", params->outFnameInfix, jwhen);
       /* m->sigmafp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&t, sizeof(t), 1, m->sigmafp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].sigma), sizeof(m->p[i].sigma), 1, m->sigmafp[iwhen]);
        */

//      MPI_File File_sigma_trans[3];
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_sigma_trans[iwhen]);
        m->sigmafp[iwhen]=File_sigma_trans[iwhen];

      }
      else
      {
        m->sigmafp[iwhen] = NULL;
      }
      /* if writing out slip, open file here */
      if ((params->writeSlip & jwhen) > 0)
      {
        sprintf(tmpFname, "slip.%s.out.%d", params->outFnameInfix, jwhen);
        /*m->slipfp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&t, sizeof(t), 1, m->slipfp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].d), sizeof(m->p[i].d), 1, m->slipfp[iwhen]);
        */

        //MPI_File File_slip_trans[3];
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_slip_trans[iwhen]);
        m->slipfp[iwhen]=File_slip_trans[iwhen];

      }
      else
      {
        m->slipfp[iwhen] = NULL;
      }

     /* if writing out slip speed, open file here */
      if ((params->writeSlipSpeed & jwhen) > 0)
      {
        sprintf(tmpFname, "slipSpeed.%s.out.%d", params->outFnameInfix, jwhen);
        /*m->slipSpeedfp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&t, sizeof(t), 1, m->slipSpeedfp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].ddot), sizeof(m->p[i].ddot), 1, m->slipSpeedfp[iwhen]);
        */

        //MPI_File File_ddot_trans[3];
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_ddot_trans[iwhen]);
        m->slipSpeedfp[iwhen]=File_ddot_trans[iwhen];

      }
      else
      {
        m->slipSpeedfp[iwhen] = NULL;
      }
      /* if writing out state, open file here */
      if ((params->writeState & jwhen) > 0)
      {
        sprintf(tmpFname, "state.%s.out.%d", params->outFnameInfix, jwhen);
        /*m->statefp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&t, sizeof(t), 1, m->statefp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].state), sizeof(m->p[i].state), 1, m->statefp[iwhen]);
        */

//      MPI_File File_state_trans[3];
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_state_trans[iwhen]);
        m->statefp[iwhen]=File_state_trans[iwhen];

      }
      else
      {
        m->statefp[iwhen] = NULL;
      }
      /* if writing out theta, open file here */
      if ((params->writeTheta & jwhen) > 0)
      {
        sprintf(tmpFname, "theta.%s.out.%d", params->outFnameInfix, jwhen);
      /*  m->thetafp[iwhen] = xfopen(tmpFname, "w");
        fwrite(&t, sizeof(t), 1, m->thetafp[iwhen]);
        for (i=0; i<m->N; i++) if (m->p[i].writePatch !=0) fwrite(&(m->p[i].theta), sizeof(m->p[i].theta), 1, m->thetafp[iwhen]);
        */

//      MPI_File File_theta_trans[3];
        MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_theta_trans[iwhen]);
        m->thetafp[iwhen]=File_theta_trans[iwhen];

      }
      else
      {
        m->thetafp[iwhen] = NULL;
      }
    }

	iwhen=3;
     /* always write out theta at end of run */
      sprintf(tmpFname, "tau.%s.out.final", params->outFnameInfix);
      MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_tau_trans[iwhen]);
      m->taufp[iwhen] = File_tau_trans[iwhen]; ;
      
      sprintf(tmpFname, "sigma.%s.out.final", params->outFnameInfix);
      MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_sigma_trans[iwhen]);
      m->sigmafp[iwhen] = File_sigma_trans[iwhen]; ;
      
      sprintf(tmpFname, "slip.%s.out.final", params->outFnameInfix);
      MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_slip_trans[iwhen]);
      m->slipfp[iwhen] = File_slip_trans[iwhen]; ;
      
      sprintf(tmpFname, "slipSpeed.%s.out.final", params->outFnameInfix);
      MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_ddot_trans[iwhen]);
      m->slipSpeedfp[iwhen] = File_ddot_trans[iwhen]; ;
      
      sprintf(tmpFname, "theta.%s.out.final", params->outFnameInfix);
      MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_theta_trans[iwhen]);
      m->thetafp[iwhen] = File_theta_trans[iwhen]; ;
      
      sprintf(tmpFname, "state.%s.out.final", params->outFnameInfix);
      MPI_File_open(MPI_COMM_WORLD,tmpFname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&File_state_trans[iwhen]);
      m->statefp[iwhen] = File_state_trans[iwhen]; ;
    

#endif
    if (m->master)
    if (params->writeTransitions)
    {
      sprintf(tmpFname, "trans.%s.out", params->outFnameInfix);
      m->transfp = xfopen(tmpFname, "w");
    }
    else
    {
      m->transfp = NULL;
    }
  } /* if (m->master) */  
  
  return;
}



void checkDdot(FaultPatch *p, int ip)
{
  if (p->slowSlip && (p->ddotAB < p->ddotStar || p->ddotEQ < p->ddotAB))
  {
    fprintf(stderr, "Error: on patch %d, ddotStar = %lf, ddotAB = %lf, and ddotEQ = %lf\n",
            ip, p->ddotStar, p->ddotAB, p->ddotEQ);
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
  }
  
  return;
}

void openFiles(FaultModel *m, Params *params, FILE **Fps)
{
  char fline[RSQSIM_MAXLINE];
  
  
  /* check if there was a ddotEQ file supplied; if so, open it and write to eqs.*.out */
  if (strlen(params->ddotEQFname) > 0)
  {
    Fps[DDOTEQ] = xfopen(params->ddotEQFname, "r");
	if (m->master)    fprintf(m->eqsfp, "%%%% ddotEQFname = %s\n", params->ddotEQFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[DDOTEQ]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
	
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end ddotEQFile\n");
    fseek(Fps[DDOTEQ], 0L, SEEK_SET);
  }
  else
    Fps[DDOTEQ] = NULL;

  /* check if there was a initial shear stress file supplied */
  if (strlen(params->initTauFname) > 0)
  {
    Fps[INITTAU] = xfopen(params->initTauFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% initTauFname = %s\n", params->initTauFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[INITTAU]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
	
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end initTauFile\n");
    fseek(Fps[INITTAU], 0L, SEEK_SET);
  }
  else
    Fps[INITTAU] = NULL;

  /* check if there was a initial normal stress file supplied */
  if (strlen(params->initSigmaFname) > 0)
  {
    Fps[INITSIGMA] = xfopen(params->initSigmaFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% initSigmaFname = %s\n", params->initSigmaFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[INITSIGMA]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end initSigmaFile\n");
    fseek(Fps[INITSIGMA], 0L, SEEK_SET);
  }
  else
    Fps[INITSIGMA] = NULL;

  /* check if there was a initial state file supplied */
  if (strlen(params->initThetaFname) > 0)
  {
    Fps[INITTHETA] = xfopen(params->initThetaFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% initThetaFname = %s\n", params->initThetaFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[INITTHETA]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end initThetaFile\n");
    fseek(Fps[INITTHETA], 0L, SEEK_SET);
  }
  else
    Fps[INITTHETA] = NULL;

  /* check if there was a initial slip speed file supplied */
  if (strlen(params->initSlipSpeedFname) > 0)
  {
    Fps[INITSLIPSPEED] = xfopen(params->initSlipSpeedFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% initSlipSpeedFname = %s\n", params->initSlipSpeedFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[INITSLIPSPEED]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end initSlipSpeedFile\n");
    fseek(Fps[INITSLIPSPEED], 0L, SEEK_SET);
  }
  else
    Fps[INITSLIPSPEED] = NULL;

  /* check if there was a rate coeff file supplied */
  if (strlen(params->AFname) > 0)
  {
    Fps[FP_A] = xfopen(params->AFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% AFname = %s\n", params->AFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[FP_A]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end AFile\n");
    fseek(Fps[FP_A], 0L, SEEK_SET);
  }
  else
    Fps[FP_A] = NULL;

  /* check if there was a state coeff file supplied */
  if (strlen(params->BFname) > 0)
  {
    Fps[FP_B] = xfopen(params->BFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% BFname = %s\n", params->BFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[FP_B]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end BFile\n");
    fseek(Fps[FP_B], 0L, SEEK_SET);
  }
  else
    Fps[FP_B] = NULL;

  /* check if there was a Dc file supplied */
  if (strlen(params->DcFname) > 0)
  {
    Fps[DC] = xfopen(params->DcFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% DcFname = %s\n", params->DcFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[DC]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end DcFile\n");
    fseek(Fps[DC], 0L, SEEK_SET);
  }
  else
    Fps[DC] = NULL;

  /* check if there was a mu0 file supplied */
  if (strlen(params->mu0Fname) > 0)
  {
    Fps[MU0] = xfopen(params->mu0Fname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% mu0Fname = %s\n", params->mu0Fname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[MU0]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end mu0File\n");
    fseek(Fps[MU0], 0L, SEEK_SET);
  }
  else
    Fps[MU0] = NULL;

  /* check if there was a ddotStar file supplied */
  if (strlen(params->ddotStarFname) > 0)
  {
    Fps[DDOTSTAR] = xfopen(params->ddotStarFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% ddotStarFname = %s\n", params->ddotStarFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[DDOTSTAR]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end ddotStarFile\n");
    fseek(Fps[DDOTSTAR], 0L, SEEK_SET);
  }
  else
    Fps[DDOTSTAR] = NULL;

  /* check if there was a ddotAB file supplied */
  if (strlen(params->ddotABFname) > 0)
  {
    Fps[DDOTAB] = xfopen(params->ddotABFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% ddotABFname = %s\n", params->ddotABFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[DDOTAB]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end ddotABFile\n");
    fseek(Fps[DDOTAB], 0L, SEEK_SET);
  }
  else
    Fps[DDOTAB] = NULL;

  /* check if there was a alpha file supplied */
  if (strlen(params->alphaFname) > 0)
  {
    Fps[ALPHA] = xfopen(params->alphaFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% alphaFname = %s\n", params->alphaFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[ALPHA]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end alphaFile\n");
    fseek(Fps[ALPHA], 0L, SEEK_SET);
  }
  else
    Fps[ALPHA] = NULL;

  /* check if there was a tFail file supplied */
  if (strlen(params->tFailFname) > 0)
  {
    Fps[TFAIL] = xfopen(params->tFailFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% tFailFname = %s\n", params->tFailFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[TFAIL]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end tFailFile\n");
    fseek(Fps[TFAIL], 0L, SEEK_SET);
  }
  else
    Fps[TFAIL] = NULL;

  /* check if there was a tauFail file supplied */
  if (strlen(params->tauFailFname) > 0)
  {
    Fps[TAUFAIL] = xfopen(params->tauFailFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% tauFailFname = %s\n", params->tauFailFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[TAUFAIL]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end tauFailFile\n");
    fseek(Fps[TAUFAIL], 0L, SEEK_SET);
  }
  else
    Fps[TAUFAIL] = NULL;

  /* check if there was a pinned file supplied */
  if (strlen(params->pinnedFname) > 0)
  {
    Fps[PINNED] = xfopen(params->pinnedFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% pinnedFname = %s\n", params->pinnedFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[PINNED]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end pinnedFile\n");
    fseek(Fps[PINNED], 0L, SEEK_SET);
  }
  else
    Fps[PINNED] = NULL;

  /* check if there was a neighbor file supplied */
  if (strlen(params->neighborFname) > 0)
  {
    Fps[NEIGHBOR] = xfopen(params->neighborFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% neighborFname = %s\n", params->neighborFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[NEIGHBOR]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end neighborFile\n");
    fseek(Fps[NEIGHBOR], 0L, SEEK_SET);
  }
  else
    Fps[NEIGHBOR] = NULL;

  /* check if there was a slowSlip file supplied */
  if (strlen(params->slowSlipFname) > 0)
  {
    Fps[SLOWSLIP] = xfopen(params->slowSlipFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% slowSlipFname = %s\n", params->slowSlipFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[SLOWSLIP]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end slowSlipFile\n");
    fseek(Fps[SLOWSLIP], 0L, SEEK_SET);
  }
  else
    Fps[SLOWSLIP] = NULL;

  /* check if there was a writePatch file supplied */
  if (strlen(params->writePatchFname) > 0)
  {
    Fps[WRITEPATCH] = xfopen(params->writePatchFname, "r");
    if (m->master)  fprintf(m->eqsfp, "%%%% writePatchFname = %s\n", params->writePatchFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[WRITEPATCH]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end writePatchFile\n");
    fseek(Fps[WRITEPATCH], 0L, SEEK_SET);
  }
  else
  {
    Fps[WRITEPATCH] = NULL;
    m->npatchWrite = m->N; /* otherwise gets incremented in initializeFaultPatch() */
  }

  Fps[TAUDOT] = NULL;
  Fps[SIGMADOT] = NULL;
  if (params->stressRateSpecification == DIRECT)
  {
    /* check if there was a tauDot file supplied */
    if (strlen(params->tauDotFname) > 0)
    {
      Fps[TAUDOT] = xfopen(params->tauDotFname, "r");
      if (m->master)  fprintf(m->eqsfp, "%%%% tauDotFname = %s\n", params->tauDotFname);
      while (fgets(fline, RSQSIM_MAXLINE, Fps[TAUDOT]) != NULL)
      {
        if (m->master)  fprintf(m->eqsfp, "%% ");
        if (m->master) fputs(fline, m->eqsfp);
      }
      if (m->master)  fprintf(m->eqsfp, "%%%% end tauDotFile\n");
      fseek(Fps[TAUDOT], 0L, SEEK_SET);
    }
      
    /* check if there was a sigmaDot file supplied */
    if (strlen(params->sigmaDotFname) > 0)
    {
      Fps[SIGMADOT] = xfopen(params->sigmaDotFname, "r");
      if (m->master)  fprintf(m->eqsfp, "%%%% sigmaDotFname = %s\n", params->sigmaDotFname);
      while (fgets(fline, RSQSIM_MAXLINE, Fps[SIGMADOT]) != NULL)
      {
        if (m->master)  fprintf(m->eqsfp, "%% ");
        if (m->master) fputs(fline, m->eqsfp);
      }
      if (m->master)  fprintf(m->eqsfp, "%%%% end sigmaDotFile\n");
      fseek(Fps[SIGMADOT], 0L, SEEK_SET);
    }
  }

  /* check if there was a external stress rate file supplied */
  if (strlen(params->extStressRateFname) > 0)
  {
    Fps[EXTSTRESSRATE] = xfopen(params->extStressRateFname, "r");
    /**  hmm, can't really write out this to eqs.*.out because it is a binary file....
    if (m->master)  fprintf(m->eqsfp, "%%%% slowSlipFname = %s\n", params->slowSlipFname);
    while (fgets(fline, RSQSIM_MAXLINE, Fps[SLOWSLIP]) != NULL)
    {
      if (m->master)  fprintf(m->eqsfp, "%% ");
      if (m->master) fputs(fline, m->eqsfp);
    }
    if (m->master)  fprintf(m->eqsfp, "%%%% end slowSlipFile\n");
    fseek(Fps[SLOWSLIP], 0L, SEEK_SET); */
    
    fread(&(m->nExtStressRateSteps), sizeof(m->nExtStressRateSteps), 1, Fps[EXTSTRESSRATE]);
    m->tExtStressRateSteps = xmalloc(m->nExtStressRateSteps * sizeof(double));
    fread(m->tExtStressRateSteps, sizeof(double), m->nExtStressRateSteps, Fps[EXTSTRESSRATE]);
    fseek(Fps[EXTSTRESSRATE], m->startPatch*m->nExtStressRateSteps*2*sizeof(double), SEEK_CUR);
    m->iExtStressRateStep = 0;
    while(m->tExtStressRateSteps[m->iExtStressRateStep] < params->tStart) m->iExtStressRateStep++;
  }
  else
  {
    Fps[EXTSTRESSRATE] = NULL;
    m->nExtStressRateSteps = 0;
    m->iExtStressRateStep = 0;
    m->tExtStressRateSteps = NULL;
  }

  return;
}


void initializeFaultPatch(FaultModel *m, Params *params, FILE **Fps, int i, double *stressRate)
{
  int iNeighbor;
  char *neighborPtr, neighborString[RSQSIM_MAXLINE];

  m->p[i].d = 0.0;
  m->p[i].ddot = 0.0;
  m->p[i].APFail = 0;
  m->p[i].t01 = T01_NA_VAL;

  /* ddot EQ */
  if (Fps[DDOTEQ] == NULL)
    m->p[i].ddotEQ = m->ddotEQ;
  else
    fscanf(Fps[DDOTEQ], "%lf", &(m->p[i].ddotEQ));

  /* initTau */
  if (Fps[INITTAU] == NULL)
    m->p[i].tau = params->tau0_1;
  else
    fscanf(Fps[INITTAU], "%lf", &(m->p[i].tau));

  /* initSigma */
  if (Fps[INITSIGMA] == NULL)
    m->p[i].sigma = params->sigma0_1;
  else
    fscanf(Fps[INITSIGMA], "%lf", &(m->p[i].sigma));

  m->p[i].sigmaPin = params->sigmaFracPin*m->p[i].sigma;
  m->p[i].muPin = params->muPin;

  /* initTheta */
  if (Fps[INITTHETA] == NULL)
    m->p[i].theta = params->theta0_1;
  else
    fscanf(Fps[INITTHETA], "%lf", &(m->p[i].theta));

  /* initSlipSpeed */
  if (Fps[INITSLIPSPEED] == NULL)
    m->p[i].ddot = 0.0;
  else
    fscanf(Fps[INITSLIPSPEED], "%lf", &(m->p[i].ddot));

  /* A */
  if (Fps[FP_A] == NULL)  
    m->p[i].A = params->A_1;
  else
    fscanf(Fps[FP_A], "%lf", &(m->p[i].A));
  m->p[i].Asv = m->p[i].A;

  /* B */
  if (Fps[FP_B] == NULL)  
    m->p[i].B = params->B_1;
  else
    fscanf(Fps[FP_B], "%lf", &(m->p[i].B));

  if (m->p[i].B > m->p[i].A) 
    m->p[i].state = 0;
  else
  {
    m->p[i].state = CREEP;  /* if rate-strengthening, treat separately as state 3 */
    m->state3patches = TRUE;
  }

  /* Dc */
  if (Fps[DC] == NULL)  {
    m->p[i].Dc = params->Dc_1;
    m->p[i].iDc = 1.0/params->Dc_1;
  }
  else
    fscanf(Fps[DC], "%lf", &(m->p[i].Dc));

  /* mu0 */
  if (Fps[MU0] == NULL)  
    m->p[i].mu0 = params->mu0_1;
  else
    fscanf(Fps[MU0], "%lf", &(m->p[i].mu0));


  /* ddotStar */
  if (Fps[DDOTSTAR] == NULL)  
    m->p[i].ddotStar = params->ddotStar_1;
  else
    fscanf(Fps[DDOTSTAR], "%lf", &(m->p[i].ddotStar));
  
  /* ddot2A */
  m->p[i].ddot2A = m->p[i].ddotStar;
  
  /* ddotAB */
  if (Fps[DDOTAB] == NULL)  
    m->p[i].ddotAB = params->ddotAB_1;
  else
    fscanf(Fps[DDOTAB], "%lf", &(m->p[i].ddotAB));
  if (m->p[i].ddotAB < m->p[i].ddot2A)
  {
    m->p[i].ddotAB = m->p[i].ddot2A;
    if (m->p[i].ddotAB > m->p[i].ddotEQ) m->p[i].ddotAB = m->p[i].ddotEQ;
  }

  /* alpha */
  if (Fps[ALPHA] == NULL)  
    m->p[i].alpha = params->alpha_1;
  else
    fscanf(Fps[ALPHA], "%lf", &(m->p[i].alpha));

  /* tFail */
  if (Fps[TFAIL] == NULL)  
    m->p[i].tFail = -1;
  else
    fscanf(Fps[TFAIL], "%lf", &(m->p[i].tFail));

  /* tauFail */
  if (Fps[TAUFAIL] == NULL)  
    m->p[i].tauFail = -1;
  else
    fscanf(Fps[TAUFAIL], "%lf", &(m->p[i].tauFail));

  /* pinned */
  if (Fps[PINNED] == NULL)  
    m->p[i].pinned = 0;
  else
  {
    fscanf(Fps[PINNED], "%d", &(m->p[i].pinned));
    if (m->p[i].pinned == 1) m->p[i].dDotDrive = 0; /* if pinned then shouldn't contribute to backslipping stressing rates */
  }

  /* slowSlip */
  if (Fps[SLOWSLIP] == NULL)
    m->p[i].slowSlip = params->slowSlip_1;
  else
    fscanf(Fps[SLOWSLIP], "%d", &(m->p[i].slowSlip));
    
  /* writePatch */
  if (Fps[WRITEPATCH] == NULL)
    m->p[i].writePatch = 1;
  else
  {
    fscanf(Fps[WRITEPATCH], "%d", &(m->p[i].writePatch));
    if (m->p[i].writePatch != 0) m->npatchWrite += 1;
  }
    
  if (params->stressRateSpecification == DIRECT)
  {
     m->p[i].dDotDrive = NAN; /* dDotDrive isn't applicable when stressing rates are specified directly */

     /* tauDot */
     if (Fps[TAUDOT] == NULL)  
       m->p[i].tauDotDrive = 0;
     else
       fscanf(Fps[TAUDOT], "%lf", &(m->p[i].tauDotDrive));

     /* sigmaDot */
     if (Fps[SIGMADOT] == NULL)  
       m->p[i].sigmaDotDrive = 0;
     else
       fscanf(Fps[SIGMADOT], "%lf", &(m->p[i].sigmaDotDrive));
   }
   
      /* neighbors */
  if (Fps[NEIGHBOR] == NULL)  
  {
    m->p[i].nNeighbors = -1;
    m->p[i].neighbors = NULL;
  }
  else
  {

    fgets(neighborString, RSQSIM_MAXLINE, Fps[NEIGHBOR]);
    m->p[i].nNeighbors = strflds(neighborString);
    m->maxNeighbors = MAX(m->maxNeighbors, m->p[i].nNeighbors);
    m->p[i].neighbors = xmalloc(m->p[i].nNeighbors * sizeof(int));
    neighborPtr = strtok(neighborString, " \t\n");
    for (iNeighbor=0; iNeighbor < m->p[i].nNeighbors; iNeighbor++)
    {
      sscanf(neighborPtr, "%d", &(m->p[i].neighbors[iNeighbor]));
      neighborPtr = strtok(NULL, " \t\n");
    }        
  }
  m->p[i].neighborInRuptureState = 0;

  m->p[i].aReduced = 0;

  m->p[i].sigmaClipped = 0;

  if (params->stressRateSpecification == TENSOR_RATE)
  {
    projectStress(stressRate, m->p[i].op.nu, m->p[i].op.u, &(m->p[i].sigmaDotDrive), &(m->p[i].tauDotDrive));
    m->p[i].sigmaDotDrive = -m->p[i].sigmaDotDrive; /* change to compressive stress positive */
  }
  
  /* extStressRate */
  m->p[i].tauDotExt = NULL;
  m->p[i].sigmaDotExt = NULL;
  if (Fps[EXTSTRESSRATE] != NULL && i >= m->startPatch && i <= m->stopPatch)
  {
    m->p[i].tauDotExt = xmalloc(m->nExtStressRateSteps*sizeof(double));
    m->p[i].sigmaDotExt = xmalloc(m->nExtStressRateSteps*sizeof(double));
    fread(m->p[i].tauDotExt, sizeof(double), m->nExtStressRateSteps, Fps[EXTSTRESSRATE]);
    fread(m->p[i].sigmaDotExt, sizeof(double), m->nExtStressRateSteps, Fps[EXTSTRESSRATE]); 
  }
  
      
  checkDdot(&(m->p[i]), i);
       
  return;
}


void closeFiles(FILE **Fps, int nfps)
{
  int ifps;
  
  for (ifps=0; ifps<nfps; ifps++)
    if (Fps[ifps] != NULL) fclose(Fps[ifps]);
  
  
  return;
}

void xferParamsToM(FaultModel *m, Params *params)
{  
  m->ddotEQ = params->ddotEQ_1;
  m->dMu3 = params->dMu3;
  m->stressOvershootFactor = params->stressOvershootFactor;
  m->stressRateSpecification = params->stressRateSpecification;
  m->nRupture = 0;
  m->areaRupture = 0.0;
  m->npatchWrite = 0;
  m->fA = params->fA;
  m->nTransitions = 0;
  m->minDtWriteCoseismic = params->minDtWrite; /* in case specified the old single minDtWrite */
  m->minDtWriteInterseismic = params->minDtWrite; /* in case specified the old single minDtWrite */
  if (params->minDtWriteCoseismic > 0) m->minDtWriteCoseismic = params->minDtWriteCoseismic;
  if (params->minDtWriteInterseismic > 0) m->minDtWriteInterseismic = params->minDtWriteInterseismic;
  m->minMagWrite = params->minMagWrite;
  m->tLastWrite = params->tStart;
  m->DEBUG = params->DEBUG;
  m->ZBrentUpperBracket = params->ZBrentUpperBracket;
  m->ZBrentTol = params->ZBrentTol;
  m->maxThetaPin = params->maxThetaPin;
  m->lowSigmaAction = params->lowSigmaAction;
  m->slipInState1 = params->slipInState1;
  m->resetStressingRates = params->resetStressingRates;
  m->maxTimeStep = params->maxTimeStep;
  m->unreduceAFix = params->unreduceAFix;
  m->nbuffSnap = params->nbuffSnap;
  m->aReduceNmax = params->aReduceNmax;
  m->aReduceN1 = params->aReduceN1;
  m->aReduceN2 = params->aReduceN2;
  m->fAmin = params->fAmin;
  m->maxDt1Zero = params->maxDt1Zero;

  return;
}

void allocSnapBuffers(FaultModel *m)
{
  int npatches, iwhen;
  
  npatches = m->stopPatch - m->startPatch + 1;
  m->taubuffer = (double *) xmalloc(3 * m->nbuffSnap * npatches * sizeof(double));
  for (iwhen=0; iwhen<3; iwhen++) m->ibuffSnap[iwhen] = 0;
  
  
  return;
}
