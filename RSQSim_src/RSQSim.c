#include "RSQSim.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif

static char svnId[] = {"$Id: RSQSim.c 1065 2017-04-07 21:57:01Z dinger $"};

/* some local utitlities */

void writePED_xx(FaultModel *m, Eqs *eqs, int nslipped, MPI_Comm comm,
        MPI_Offset *disp, int maxbuffer, int *icount1, int *ibuff, int *nrec,
        int *tpmap, int *blocksize, int *nelms, int isOnetime,
        double *buffer_d,
        double *buffer_tau,
        double *buffer_sigma,
        double *buffer_mu,
        int    *buffer_e,
        int    *buffer_p,
        double *buffer_t,
        double *buffer_taup);

void setup(FaultModel *m,int maxbuffer, MPI_Datatype datatype,
                int* nelements, MPI_Datatype *filetype,int* nrec, int* index);
void write1D_static(double t, MPI_File fh, void* array,int nelements,int nrec,
        int maxbuffer,int* icount, void* buffer,
        MPI_Offset *disp, MPI_Datatype datatype,MPI_Datatype filetype,int isOnetime);

void copy_double(FaultModel *m,double t, int nelements_trans, int * index, double* array_d,int which_var);
void copy_char(FaultModel *m,double t, int nelements_trans, int * index, char* array_i,int which_var);
void rollback_x(MPI_Offset *Disp, MPI_Offset Disp_save, int *iount, int icount_save,int nbyte,int maxbuffer,int nrec);


/* restart */
void    readwrite_dumps(Eqs *eqs, FaultModel *m,double *t,int *currentStep,char* buffer,int rw);


Eqs RSQSim(Params *params, FaultModel *m)
{
  /* ---extras ---*/  
  /* 	maxbuffer0: for EQS info
	  maxbuffer1: for PED
	  maxbuffer2: for writeXXX
  */

  int maxbuffer0=params->nbuffEQS;
  int maxbuffer1=params->nbuffPED;
  int maxbuffer2=params->nbuffTRANS;

  eqs_info *buffer;
  buffer = malloc(maxbuffer0*sizeof(eqs_info));
  int icount=0;

  buffer[0].maxx = maxbuffer0;

  /* PED */
  
 MPI_Offset disp=0;
  int icount1=0, ibuff=0, nrec=0, nelms=0;
  int *tpmap = (int*) malloc(maxbuffer1 * sizeof(int));
  int *blocksize = (int*) malloc(maxbuffer1 * sizeof(int));

  double *buffer_d, *buffer_tau,*buffer_sigma,*buffer_mu,*buffer_t, *buffer_taup;
  int *buffer_e,*buffer_p;

  int ntotal_PED = maxbuffer1*(m->stopPatch - m->startPatch+1);

  buffer_d =     (double *)malloc(sizeof(double)*ntotal_PED);
  buffer_tau =   (double *)malloc(sizeof(double)*ntotal_PED);
  buffer_sigma = (double *)malloc(sizeof(double)*ntotal_PED);
  buffer_mu =    (double *)malloc(sizeof(double)*ntotal_PED);
  buffer_t =     (double *)malloc(sizeof(double)*ntotal_PED);
  buffer_taup =  (double *) malloc(sizeof(double)*ntotal_PED);
  buffer_e =     (int *)malloc(sizeof(int)*ntotal_PED);
  buffer_p =     (int *)malloc(sizeof(int)*ntotal_PED);


  Transition mytr, tr;

  int currentStep=0;    // For debugging - Efe
  double startTime=0, endTime=0;
  double startTime2=0,endTime2=0;

  double t;
  Eqs eqs;

#ifdef USE_SERIAL_IO
	if (m->master) printf(" Use original output (serial) \n");
#else
	if (m->master) printf(" Use parallel output  \n");
#endif	
   m->t = params->tStart;
   
  initFaultModel(m, params);  
  if (m->master) printf("initFaultModel() is done.\n");

  initEqs(&eqs, m, params);
  if (m->master) printf("assignPatches() is done.\n");

   m->memorySize+= ntotal_PED *sizeof(double)*7 + maxbuffer1 * sizeof(int)*2 + maxbuffer0*sizeof(eqs_info);

  
  m->t = params->tStart;

	t=m->t;
 /* ---------------------restart -----*/

	FILE * readTimes;
	int num_restart_times=0;
	int i,kkw=0,kkr;
	int *restart_times;
	char restart_buffer[50];
	if ((readTimes=fopen("restart_times","r"))!= (FILE*)NULL){
		if (m->master) fprintf(stdout," Reading restart times file...\n");

		fscanf(readTimes,"%d",&num_restart_times);
		restart_times=(int *) xmalloc(num_restart_times*sizeof(int));
		for (i=0;i<num_restart_times;i++){
			fscanf(readTimes,"%d",&restart_times[i]);
			if (m->master) fprintf(stdout,"%d  %d \n",i,restart_times[i]);
		}
	} 
	
	if ((readTimes=fopen("restart_from","r"))!=(FILE*)NULL){
		/* --read in files-*/
		fscanf(readTimes,"%d",&kkr);
		sprintf(restart_buffer,"CHK%-d",kkr);

		readwrite_dumps(&eqs,m,&t,&currentStep,restart_buffer,0);   /* read in dumps  */
                if (m->master) fprintf(stdout," This is a restart run from eqs #%s, time= %e currentStep= %d  **\n ",restart_buffer,t,currentStep);
 
	}
	else {
	 	if (m->master) fprintf(stdout," This is a fresh run \n");
	}
/*----------------------end restart--*/

  
/* -----transition ------*/
        int nelements_trans,nrec_trans;
        MPI_Datatype filetype_d,filetype_i;
        int index[m->stopPatch-m->startPatch+2];

        MPI_Barrier(MPI_COMM_WORLD);

  if (m->master){
                fprintf(stdout,"before  set up\n");
		fprintf(stdout,"maxbuffer0(EQS)=%d,maxbuffer1 (writePED)=%d,maxbuffer2 (writeXXX)=%d\n",maxbuffer0,maxbuffer1,maxbuffer2);
	}
        setup(m,maxbuffer2,  MPI_DOUBLE,
                &nelements_trans, &filetype_d,&nrec_trans,index);

   	setup(m,maxbuffer2,  MPI_CHAR,
                &nelements_trans, &filetype_i,&nrec_trans,index);

//rollback
	MPI_Offset Disp_trans_save[6][3];
	int icount_trans_save[6][3];
//total number of elements
	int ntotal = 	nelements_trans*maxbuffer2;
	
#ifdef USE_SERIAL_IO
	ntotal=1 ;  /* keep it alive */
#endif	

	
//theta
        int icount_theta_trans0=0,icount_theta_trans1=0,icount_theta_trans2=0;
        double *buffer_theta_trans0,*buffer_theta_trans1,*buffer_theta_trans2;
        MPI_Offset thetaDisp_trans0=0,thetaDisp_trans1=0,thetaDisp_trans2=0;
 
 	if ((params->writeTheta & TRANSITION) > 0)
          buffer_theta_trans0 = (double *) malloc(sizeof(double)*ntotal);
 	if ((params->writeTheta & EVENT) > 0)
          buffer_theta_trans1 = (double *) malloc(sizeof(double)*ntotal);
        if ((params->writeTheta & EVENT_TRANSITION) > 0)
	  buffer_theta_trans2 = (double *) malloc(sizeof(double)*ntotal);

//tau
        int icount_tau_trans0=0, icount_tau_trans1=0, icount_tau_trans2=0;
        double *buffer_tau_trans0, *buffer_tau_trans1, *buffer_tau_trans2;
        MPI_Offset tauDisp_trans0=0, tauDisp_trans1=0, tauDisp_trans2=0;
	if ((params->writeTau & TRANSITION) > 0)
          buffer_tau_trans0 = (double *) malloc(sizeof(double)*ntotal);
        if ((params->writeTau & EVENT) > 0)
	  buffer_tau_trans1 = (double *) malloc(sizeof(double)*ntotal);
	if ((params->writeTau & EVENT_TRANSITION) > 0)
          buffer_tau_trans2 = (double *) malloc(sizeof(double)*ntotal);

//sigma
        int icount_sigma_trans0=0, icount_sigma_trans1=0, icount_sigma_trans2=0;
        double *buffer_sigma_trans0, *buffer_sigma_trans1, *buffer_sigma_trans2;
        MPI_Offset sigmaDisp_trans0=0, sigmaDisp_trans1=0, sigmaDisp_trans2=0;
	if ((params->writeSigma & TRANSITION) > 0)
          buffer_sigma_trans0 = (double *) malloc(sizeof(double)*ntotal);
        if ((params->writeSigma & EVENT) > 0)
	  buffer_sigma_trans1 = (double *) malloc(sizeof(double)*ntotal);
        if ((params->writeSigma & EVENT_TRANSITION) > 0)
	  buffer_sigma_trans2 = (double *) malloc(sizeof(double)*ntotal);
//slip
        int icount_slip_trans0=0, icount_slip_trans1=0, icount_slip_trans2=0;
        double *buffer_slip_trans0, *buffer_slip_trans1, *buffer_slip_trans2;
        MPI_Offset slipDisp_trans0=0, slipDisp_trans1=0, slipDisp_trans2=0;

	if ((params->writeSlip & TRANSITION) > 0)
          buffer_slip_trans0 = (double *) malloc(sizeof(double)*ntotal);
	if ((params->writeSlip & EVENT) > 0)
          buffer_slip_trans1 = (double *) malloc(sizeof(double)*ntotal);
	if ((params->writeSlip & EVENT_TRANSITION) > 0)
          buffer_slip_trans2 = (double *) malloc(sizeof(double)*ntotal);

//ddot
        int icount_ddot_trans0=0,icount_ddot_trans1=0,icount_ddot_trans2=0;
        double *buffer_ddot_trans0,*buffer_ddot_trans1,*buffer_ddot_trans2;
        MPI_Offset ddotDisp_trans0=0,ddotDisp_trans1=0,ddotDisp_trans2=0;

	if ((params->writeSlipSpeed & TRANSITION) > 0)
          buffer_ddot_trans0=(double*)malloc(sizeof(double)*ntotal);
	if ((params->writeSlipSpeed & EVENT) > 0)
          buffer_ddot_trans1=(double*)malloc(sizeof(double)*ntotal);
	if ((params->writeSlipSpeed & EVENT_TRANSITION) > 0)
          buffer_ddot_trans2=(double*)malloc(sizeof(double)*ntotal);
//state
        int icount_state_trans0=0,icount_state_trans1=0,icount_state_trans2=0;
        char *buffer_state_trans0,*buffer_state_trans1,*buffer_state_trans2;
        MPI_Offset stateDisp_trans0=0,stateDisp_trans1=0,stateDisp_trans2=0;

	if ((params->writeState & TRANSITION) > 0)
          buffer_state_trans0=(char*)malloc(sizeof(char)*ntotal);
	if ((params->writeState & EVENT) > 0)
          buffer_state_trans1=(char*)malloc(sizeof(char)*ntotal);
	if ((params->writeState & EVENT_TRANSITION) > 0)
          buffer_state_trans2=(char*)malloc(sizeof(char)*ntotal);

        MPI_Barrier(MPI_COMM_WORLD);
        if (m->master){
                fprintf(stdout,"Nwrite= %d  N=%d done with set up\n",m->npatchWrite,m->N);
                fflush(stdout);
        }


        double *array_d=(double*)malloc(sizeof(double)*nelements_trans);
        char *array_i=(char*)malloc(sizeof(char)*nelements_trans);

#ifdef USE_SERIAL_IO

#else

/* place the initial output in file here  */
//	int ii,jj;
		
	int isOnetime=FALSE;
	
	write_trans(m, params,m->t,   nelements_trans,  nrec_trans,     index, 
	 array_d,    array_i,  filetype_d,  filetype_i,   maxbuffer2,  isOnetime,
	     &icount_tau_trans0,      buffer_tau_trans0, &tauDisp_trans0,
	     &icount_tau_trans1,      buffer_tau_trans1, &tauDisp_trans1,
	     &icount_tau_trans2,      buffer_tau_trans2, &tauDisp_trans2,
	 
	     &icount_sigma_trans0,      buffer_sigma_trans0, &sigmaDisp_trans0,
	     &icount_sigma_trans1,      buffer_sigma_trans1, &sigmaDisp_trans1,
	     &icount_sigma_trans2,      buffer_sigma_trans2, &sigmaDisp_trans2,
	 
	     &icount_slip_trans0,      buffer_slip_trans0, &slipDisp_trans0,
	     &icount_slip_trans1,      buffer_slip_trans1, &slipDisp_trans1,
	     &icount_slip_trans2,      buffer_slip_trans2, &slipDisp_trans2,
	
 	     &icount_ddot_trans0,      buffer_ddot_trans0, &ddotDisp_trans0,
	     &icount_ddot_trans1,      buffer_ddot_trans1, &ddotDisp_trans1,
	     &icount_ddot_trans2,      buffer_ddot_trans2, &ddotDisp_trans2,

	     &icount_theta_trans0,      buffer_theta_trans0, &thetaDisp_trans0,
	     &icount_theta_trans1,      buffer_theta_trans1, &thetaDisp_trans1,
	     &icount_theta_trans2,      buffer_theta_trans2, &thetaDisp_trans2,
	 
	     &icount_state_trans0,    buffer_state_trans0, &stateDisp_trans0,
	     &icount_state_trans1,    buffer_state_trans1, &stateDisp_trans1,
	     &icount_state_trans2,    buffer_state_trans2, &stateDisp_trans2);

#endif

 
  if (m->master) {
        fprintf(stdout,"FaultPatch size: %zd\n",sizeof(FaultPatch));
        fprintf(stdout,"nEq size: %lf\n",params->nEq);
        fprintf(stdout,"double size: %zd\n",sizeof(double));
        fprintf(stdout,"int size: %zd\n",sizeof(int));
        fprintf(stdout,"memory size per partition: %ld MBytes\n", m->memorySize/1024/1024);
        fprintf(stdout,"Number of elements = %d, number of Processes = %d \n",m->N,m->nProc);
        fflush(stdout);
        }
#ifdef USE_OPENMP
  if (m->master) fprintf(stdout,"OMP size: %d\n",omp_get_num_procs());
#endif
#ifdef USE_MPI
  m->initializationTime=MPI_Wtime()-m->initializationTime;
	if (m->master) {
	fprintf(stdout," Initialization time : %lfs \n",m->initializationTime);
	fflush(stdout);
	}
#endif


  while ((eqs.ieq < params->nEq && m->t <= params->maxT && m->nTransitions <= params->maxTrans /* &&
          access("stop", F_OK) == -1 */ ) || m->nRupture > 0)
  {
    if(m->master)
      if(m->nTransitions > currentStep*1000){printf("nTransitions = %d current step (x1000) = %d eqs = %ld\n",
		                                     m->nTransitions,currentStep*1000,eqs.ieq); currentStep++;}   // For debugging - Efe

#ifdef USE_MPI
    startTime=MPI_Wtime();
#endif

    mytr = findMinTransitionTime(m);

#ifdef USE_MPI
    mytr.patch = mytr.nextState*m->N + mytr.patch;  /* pack nextState and patch into patch */

    startTime2=MPI_Wtime();

    MPI_Allreduce(&mytr, &tr, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    endTime2=MPI_Wtime();
    m->allReduceTime+=endTime2-startTime2;

    tr.nextState = tr.patch/m->N;  /* unpack nextState */
    tr.patch = tr.patch % m->N; /* unpack patch */

      endTime=MPI_Wtime();
     //     m->findTransitionTime+=endTime-startTime;
     	 m->findTransitionTime+=startTime2-startTime;

#else
      tr = mytr; 
#endif    

    if (tr.dt < 0)
    {
      if (m->master)
        fprintf(stderr, "Warning: time step negative: dt = %.3f, patch = %d, current state = %d; next state = %d\n",
                tr.dt, tr.patch, m->p[tr.patch].state, tr.nextState);

      tr.dt = 0;
    }
    
    if (tr.patch < 0)
    {
      if (m->master)
        fprintf(stderr, "Warning: tr.patch < 0; all patches locked?  In any case, finishing up\n");

      tr.dt = params->maxT - m->t + 1;
    }
      
    
    if (m->DEBUG > 2 && m->master) 
    {
      fprintf(stderr, "%15.2f %15.4f %4d %1d -> %1d\n", 
              m->t, tr.dt, tr.patch, m->p[tr.patch].state, tr.nextState);
      fflush(stderr);  
    }

#ifdef USE_MPI
    startTime=MPI_Wtime();
#endif

    timeStep(m, tr.dt);

#ifdef USE_MPI
//      Use this barrier to test the load balance in timeStep. Result: completely in-balance.
//      MPI_Barrier(MPI_COMM_WORLD);
    endTime=MPI_Wtime();
    m->timeStepTime+=endTime-startTime;
#endif

    m->t += tr.dt;

#ifdef USE_MPI
    startTime=MPI_Wtime();
#endif

    makeTransition(m, tr);

  //  m->nTransitions++;
#ifdef USE_MPI
    endTime=MPI_Wtime();
    m->transitionTime+=endTime-startTime;
    startTime=MPI_Wtime();
#endif

//    updateEqs(&eqs, m, m->t, tr, params);

    updateEqs(&eqs, m, m->t, tr, params, buffer, &icount,
        &disp, maxbuffer0,maxbuffer1,maxbuffer2, &icount1, &ibuff, &nrec,
        tpmap, blocksize, &nelms,
        buffer_d,
        buffer_tau,
        buffer_sigma,
        buffer_mu,
        buffer_e,
        buffer_p,
        buffer_t,
        buffer_taup,

	nelements_trans,nrec_trans,filetype_d,index,
	&icount_tau_trans0,buffer_tau_trans0, &tauDisp_trans0,
        &icount_tau_trans1,buffer_tau_trans1, &tauDisp_trans1,
        &icount_tau_trans2,buffer_tau_trans2, &tauDisp_trans2,

	&icount_sigma_trans0,buffer_sigma_trans0, &sigmaDisp_trans0,
        &icount_sigma_trans1,buffer_sigma_trans1, &sigmaDisp_trans1,
        &icount_sigma_trans2,buffer_sigma_trans2, &sigmaDisp_trans2,

	&icount_slip_trans0,buffer_slip_trans0, &slipDisp_trans0,
        &icount_slip_trans1,buffer_slip_trans1, &slipDisp_trans1,
        &icount_slip_trans2,buffer_slip_trans2, &slipDisp_trans2,

	&icount_ddot_trans0,buffer_ddot_trans0, &ddotDisp_trans0,
        &icount_ddot_trans1,buffer_ddot_trans1, &ddotDisp_trans1,
        &icount_ddot_trans2,buffer_ddot_trans2, &ddotDisp_trans2,

	&icount_theta_trans0,buffer_theta_trans0, &thetaDisp_trans0,
        &icount_theta_trans1,buffer_theta_trans1, &thetaDisp_trans1,
        &icount_theta_trans2,buffer_theta_trans2, &thetaDisp_trans2,
	filetype_i,
	&icount_state_trans0,buffer_state_trans0, &stateDisp_trans0,
        &icount_state_trans1,buffer_state_trans1, &stateDisp_trans1,
        &icount_state_trans2,buffer_state_trans2, &stateDisp_trans2,

	Disp_trans_save,icount_trans_save,
        array_d,array_i

	
);

#ifdef USE_MPI
    endTime=MPI_Wtime();
    m->updatingEqsTime+=endTime-startTime;
#endif
/*------------write out  ----*/
	if (kkw < num_restart_times && eqs.ieq == restart_times[kkw]){

	/* time to write  */
	if (m->master){
		sprintf(restart_buffer,"mkdir CHK%-d",restart_times[kkw]);
		system(restart_buffer);
		fprintf(stdout,"mkdir for eqs #%ld at time= %lf \n",eqs.ieq,t);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	restart_buffer[0]='\0';
        sprintf(restart_buffer,"CHK%-d",restart_times[kkw]);

	readwrite_dumps(&eqs, m, &t,&currentStep, restart_buffer,1);

	kkw=kkw+1;
	}
/*------------write out  ----*/


/*    if (access("stopnow", F_OK) == 0) break; */

  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (m->master) fprintf(stderr," Finishing up\n");

  //all needed
  if (m->master & icount !=0)
  {
    for (i=0;i<icount;i++) 
    {
      fprintf(m->eqsfp, "%.20lf %lf %lf %lf %lf %lf %lf %.20lf %d\n",buffer[i].t0,buffer[i].M0,buffer[i].Mw,
      buffer[i].x,buffer[i].y,buffer[i].z,buffer[i].area,buffer[i].dt,buffer[i].hypercenter);
    }
    fflush(m->eqsfp);
  }

  #ifdef USE_SERIAL_IO
  finishUp(m, params);
  #else


  if (m->master) fprintf(stderr," Finished EQS\n");

  isOnetime=TRUE;
  int nslipped=1;
  if (params->writePED & icount1 !=0)
  {
    maxbuffer1=icount1;

    writePED_xx(m, &eqs, nslipped, MPI_COMM_WORLD,
                &disp, maxbuffer1, &icount1, &ibuff, &nrec,
                tpmap, blocksize, &nelms, TRUE,
                buffer_d,
                buffer_tau,
                buffer_sigma,
                buffer_mu,
                buffer_e,
                buffer_p,
                buffer_t,
                buffer_taup);

  }	

  /* close PED files */
  if (params->writePED)
  {
    MPI_File_close(&(m->eListfp));
    MPI_File_close(&(m->pListfp));
    MPI_File_close(&(m->dListfp));
    MPI_File_close(&(m->dtauListfp));
    MPI_File_close(&(m->dsigmaListfp));
    MPI_File_close(&(m->dmuListfp));
    MPI_File_close(&(m->tListfp)); 
    MPI_File_close(&(m->taupListfp));
  }
  if (m->master) fprintf(stderr," Finished PED\n");

  write_trans(m, params,m->t,   nelements_trans,  nrec_trans,     index, 
   array_d,    array_i,  filetype_d,  filetype_i,   maxbuffer2,  isOnetime,
       &icount_tau_trans0,      buffer_tau_trans0, &tauDisp_trans0,
       &icount_tau_trans1,      buffer_tau_trans1, &tauDisp_trans1,
       &icount_tau_trans2,      buffer_tau_trans2, &tauDisp_trans2,

       &icount_sigma_trans0,      buffer_sigma_trans0, &sigmaDisp_trans0,
       &icount_sigma_trans1,      buffer_sigma_trans1, &sigmaDisp_trans1,
       &icount_sigma_trans2,      buffer_sigma_trans2, &sigmaDisp_trans2,

       &icount_slip_trans0,      buffer_slip_trans0, &slipDisp_trans0,
       &icount_slip_trans1,      buffer_slip_trans1, &slipDisp_trans1,
       &icount_slip_trans2,      buffer_slip_trans2, &slipDisp_trans2,

       &icount_ddot_trans0,      buffer_ddot_trans0, &ddotDisp_trans0,
       &icount_ddot_trans1,      buffer_ddot_trans1, &ddotDisp_trans1,
       &icount_ddot_trans2,      buffer_ddot_trans2, &ddotDisp_trans2,

       &icount_theta_trans0,      buffer_theta_trans0, &thetaDisp_trans0,
       &icount_theta_trans1,      buffer_theta_trans1, &thetaDisp_trans1,
       &icount_theta_trans2,      buffer_theta_trans2, &thetaDisp_trans2,

       &icount_state_trans0,    buffer_state_trans0, &stateDisp_trans0,
       &icount_state_trans1,    buffer_state_trans1, &stateDisp_trans1,
       &icount_state_trans2,    buffer_state_trans2, &stateDisp_trans2);



   if ((params->writeTau & EVENT) > 0) MPI_File_set_size(m->taufp[1],tauDisp_trans1);
   if ((params->writeSigma & EVENT) > 0) MPI_File_set_size(m->sigmafp[1],sigmaDisp_trans1);
   if ((params->writeSlip & EVENT) > 0) MPI_File_set_size(m->slipfp[1],slipDisp_trans1);
   if ((params->writeSlipSpeed & EVENT) > 0) MPI_File_set_size(m->slipSpeedfp[1],ddotDisp_trans1);
   if ((params->writeTheta & EVENT) > 0) MPI_File_set_size(m->thetafp[1],thetaDisp_trans1);
   if ((params->writeState & EVENT) > 0) MPI_File_set_size(m->statefp[1],stateDisp_trans1);

   if ((params->writeTau & EVENT_TRANSITION) > 0) MPI_File_set_size(m->taufp[2],tauDisp_trans2);
   if ((params->writeSigma & EVENT_TRANSITION) > 0) MPI_File_set_size(m->sigmafp[2],sigmaDisp_trans2);
   if ((params->writeSlip & EVENT_TRANSITION) > 0) MPI_File_set_size(m->slipfp[2],slipDisp_trans2);
   if ((params->writeSlipSpeed & EVENT_TRANSITION) > 0) MPI_File_set_size(m->slipSpeedfp[2],ddotDisp_trans2);
   if ((params->writeTheta & EVENT_TRANSITION) > 0) MPI_File_set_size(m->thetafp[2],thetaDisp_trans2);
   if ((params->writeState & EVENT_TRANSITION) > 0) MPI_File_set_size(m->statefp[2],stateDisp_trans2);


  /* close trans files */
  int iwhen,jwhen;
  for (iwhen=0; iwhen<3; iwhen++)
  {
    jwhen = rint(pow(2, iwhen));
    if ((params->writeTau & jwhen) > 0) MPI_File_close(&(m->taufp[iwhen]));
    if ((params->writeSigma & jwhen) > 0) MPI_File_close(&(m->sigmafp[iwhen]));
    if ((params->writeSlip & jwhen) > 0)  MPI_File_close(&(m->slipfp[iwhen]));
    if ((params->writeSlipSpeed & jwhen) > 0)  MPI_File_close(&(m->slipSpeedfp[iwhen]));
    if ((params->writeTheta & jwhen) > 0)  MPI_File_close(&(m->thetafp[iwhen]));
    if ((params->writeState & jwhen) > 0)  MPI_File_close(&(m->statefp[iwhen]));

  }


  if (m->master) fprintf(stderr," Finished TRANS\n");


   /*  final output always */

   int nbuffer=1;
   int icount_final;
   double * buffer_double_final=malloc(sizeof(double)*nelements_trans);
   char * buffer_char_final=malloc(sizeof(char)*nelements_trans);

   MPI_Offset Disp_final=0;

  setup(m,nbuffer,  MPI_DOUBLE,
          &nelements_trans, &filetype_d,&nrec_trans,index);

  setup(m,nbuffer,  MPI_CHAR,
          &nelements_trans, &filetype_i,&nrec_trans,index);

  MPI_File file_trans[6]=
  {(m->taufp[3]),(m->sigmafp[3]),(m->slipfp[3]),(m->slipSpeedfp[3]),(m->thetafp[3]),(m->statefp[3])};

  //tau, sigma,slip,ddot,theta,state 
  int icomp;
  for (icomp=0;icomp<5;icomp++)
  {
    icount_final=0;Disp_final=0;
    copy_double(m,m->t, nelements_trans,  index,  array_d,icomp+1);
    write1D_static(m->t, file_trans[icomp],  array_d, nelements_trans,nrec_trans,
            nbuffer,&icount_final,  buffer_double_final,
            &Disp_final, MPI_DOUBLE,filetype_d,FALSE);
	    MPI_File_set_size(file_trans[icomp],Disp_final);

	    MPI_File_close(&file_trans[icomp]);
  }
  
  icount_final=0;Disp_final=0;
  copy_char(m,m->t, nelements_trans,  index,  array_i,6);
  write1D_static(m->t, (m->statefp[3]),  array_i, nelements_trans,nrec_trans,
          nbuffer,&icount_final,  buffer_char_final,
          &Disp_final, MPI_CHAR,filetype_i,FALSE);

 	  MPI_File_set_size(file_trans[5],Disp_final);
	  MPI_File_close(&file_trans[5]);

  if (m->master) fprintf(stderr," Finished final\n");

	 
#endif


  free(array_d);
  free(array_i);

 free(tpmap);
  free(blocksize);

  return(eqs);
}   



//void updateEqs(Eqs *eqs, FaultModel *m, double t, Transition tr, Params *params)

void updateEqs(Eqs *eqs, FaultModel *m, double t, Transition tr, Params *params,
        eqs_info *buffer,int *icount,
        MPI_Offset *disp, int maxbuffer0, int maxbuffer1, int maxbuffer2,
        int *icount1, int *ibuff, int *nrec, int *tpmap, int *blocksize, int *nelms,
        double *buffer_d,
        double *buffer_tau,
        double *buffer_sigma,
        double *buffer_mu,
        int    *buffer_e,
        int    *buffer_p,
        double *buffer_t,
        double *buffer_taup,

	int nelements_trans,int nrec_trans,MPI_Datatype filetype_d,int *index,
	int *icount_tau_trans0,double* buffer_tau_trans0, MPI_Offset *tauDisp_trans0,
        int *icount_tau_trans1,double* buffer_tau_trans1, MPI_Offset *tauDisp_trans1,
        int *icount_tau_trans2,double* buffer_tau_trans2, MPI_Offset *tauDisp_trans2,
	
	int *icount_sigma_trans0,double* buffer_sigma_trans0, MPI_Offset *sigmaDisp_trans0,
        int *icount_sigma_trans1,double* buffer_sigma_trans1, MPI_Offset *sigmaDisp_trans1,
        int *icount_sigma_trans2,double* buffer_sigma_trans2, MPI_Offset *sigmaDisp_trans2,
	
	int *icount_slip_trans0,double* buffer_slip_trans0, MPI_Offset *slipDisp_trans0,
        int *icount_slip_trans1,double* buffer_slip_trans1, MPI_Offset *slipDisp_trans1,
        int *icount_slip_trans2,double* buffer_slip_trans2, MPI_Offset *slipDisp_trans2,
	
	int *icount_ddot_trans0,double* buffer_ddot_trans0, MPI_Offset *ddotDisp_trans0,
        int *icount_ddot_trans1,double* buffer_ddot_trans1, MPI_Offset *ddotDisp_trans1,
        int *icount_ddot_trans2,double* buffer_ddot_trans2, MPI_Offset *ddotDisp_trans2,
	
	int *icount_theta_trans0,double* buffer_theta_trans0, MPI_Offset *thetaDisp_trans0,
        int *icount_theta_trans1,double* buffer_theta_trans1, MPI_Offset *thetaDisp_trans1,
        int *icount_theta_trans2,double* buffer_theta_trans2, MPI_Offset *thetaDisp_trans2,
	
	MPI_Datatype filetype_i,
	int *icount_state_trans0,char* buffer_state_trans0, MPI_Offset *stateDisp_trans0,
        int *icount_state_trans1,char* buffer_state_trans1, MPI_Offset *stateDisp_trans1,
        int *icount_state_trans2,char* buffer_state_trans2, MPI_Offset *stateDisp_trans2,
	
	MPI_Offset Disp_trans_save[6][3],int icount_trans_save[6][3],
        double* array_d,char*array_i


)
{
  int i, j, nslipped, root;
  double M0, area;
  double tLastWriteTmp;

#ifdef USE_SERIAL_IO

#else 
/*      double *array_d =( double *) malloc(sizeof(double)*nelements_trans);
      char *array_i = (char *) malloc(sizeof(char)*nelements_trans);
*/
#endif  

#ifndef USE_MPI
  int ip;
#endif
  
  tLastWriteTmp = m->tLastWrite;
  
  if (!(eqs->eq) && m->nRupture > 0) /* an eq just started */
  {
#ifdef USE_OPENMP
#pragma omp parallel for shared(m,eqs) private(i)
#endif

    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      eqs->dsv[i] = m->p[i].d;
      eqs->tausv[i] = m->p[i].tau;
      eqs->sigmasv[i] = m->p[i].sigma;
      eqs->musv[i] = m->p[i].tau/m->p[i].sigma;

    }
    
    eqs->eq = TRUE;

    if (m->master)
    {
      fprintf(stderr, "event number %ld just started (t = %.3f yrs)\n", eqs->ieq, t/86400/365.25);
      fflush(stderr);
      eqs->t0 = t;
      eqs->hypocenter = tr.patch;
    }

#ifdef USE_SERIAL_IO
    if ((params->writeTau & EVENT) > 0) writeTau(m, t, 1, TRUE);
    if ((params->writeSigma & EVENT) > 0) writeSigma(m, t, 1, TRUE);
    if ((params->writeSlip & EVENT) > 0) writeSlip(m, t, 1, TRUE);
    if ((params->writeSlipSpeed & EVENT) > 0) writeSlipSpeed(m, t, 1, TRUE);
    if ((params->writeState & EVENT) > 0) writeState(m, t, 1, TRUE);
    if ((params->writeTheta & EVENT) > 0) writeTheta(m, t, 1, TRUE);
#else
    Disp_trans_save[0][1] = *tauDisp_trans1;
    Disp_trans_save[0][2] = *tauDisp_trans2;
    Disp_trans_save[1][1] = *sigmaDisp_trans1;
    Disp_trans_save[1][2] = *sigmaDisp_trans2;
    Disp_trans_save[2][1] = *slipDisp_trans1;
    Disp_trans_save[2][2] = *slipDisp_trans2;
    Disp_trans_save[3][1] = *ddotDisp_trans1;
    Disp_trans_save[3][2] = *ddotDisp_trans2;
    Disp_trans_save[4][1] = *thetaDisp_trans1;
    Disp_trans_save[4][2] = *thetaDisp_trans2;
    Disp_trans_save[5][1] = *stateDisp_trans1;
    Disp_trans_save[5][2] = *stateDisp_trans2;

    icount_trans_save[0][1] = *icount_tau_trans1;
    icount_trans_save[0][2] = *icount_tau_trans2;
    icount_trans_save[1][1] = *icount_sigma_trans1;
    icount_trans_save[1][2] = *icount_sigma_trans2;
    icount_trans_save[2][1] = *icount_slip_trans1;
    icount_trans_save[2][2] = *icount_slip_trans2;
    icount_trans_save[3][1] = *icount_ddot_trans1;
    icount_trans_save[3][2] = *icount_ddot_trans2;
    icount_trans_save[4][1] = *icount_theta_trans1;
    icount_trans_save[4][2] = *icount_theta_trans2;
    icount_trans_save[5][1] = *icount_state_trans1;
    icount_trans_save[5][2] = *icount_state_trans2;

    if ((params->writeTau & EVENT) > 0) 
    {
      copy_double(m,t, nelements_trans,  index,  array_d,1);

  //	        fprintf(stdout,"Pass step 2 %d %d %d-- \n",*icount_tau_trans0,*icount_tau_trans1,*icount_tau_trans2);

      write1D_static(t, (m->taufp[1]),  array_d, nelements_trans,nrec_trans,
      maxbuffer2,icount_tau_trans1,  buffer_tau_trans1,
      tauDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
  //       fprintf(stdout,"Pass step 3 %d %d %d-- \n",*icount_tau_trans0,*icount_tau_trans1,*icount_tau_trans2);

    }

    if ((params->writeSigma & EVENT) > 0) 
    {
      copy_double(m,t, nelements_trans,  index,  array_d,2);

      write1D_static(t, (m->sigmafp[1]),  array_d, nelements_trans,nrec_trans,
      maxbuffer2,icount_sigma_trans1,  buffer_sigma_trans1,
      sigmaDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
    }
    
    if ((params->writeSlip & EVENT) > 0) 
    {
      copy_double(m,t, nelements_trans,  index,  array_d,3);

      write1D_static(t, (m->slipfp[1]),  array_d, nelements_trans,nrec_trans,
      maxbuffer2,icount_slip_trans1,  buffer_slip_trans1,
      slipDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
    }
    
    if ((params->writeSlipSpeed & EVENT) > 0) 
    {
      copy_double(m,t, nelements_trans,  index,  array_d,4);

      write1D_static(t, (m->slipSpeedfp[1]),  array_d, nelements_trans,nrec_trans,
      maxbuffer2,icount_ddot_trans1,  buffer_ddot_trans1,
      ddotDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
    }
    
    if ((params->writeTheta & EVENT) > 0) 
    {
      copy_double(m,t, nelements_trans,  index,  array_d,5);

      write1D_static(t, (m->thetafp[1]),  array_d, nelements_trans,nrec_trans,
      maxbuffer2,icount_theta_trans1,  buffer_theta_trans1,
      thetaDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
    }
//	
    if ((params->writeState & EVENT) > 0) 
    {
      copy_char(m,t, nelements_trans,  index,  array_i,6);

      write1D_static(t, (m->statefp[1]),  array_i, nelements_trans,nrec_trans,
      maxbuffer2,icount_state_trans1,  buffer_state_trans1,
      stateDisp_trans1, MPI_CHAR,filetype_i,FALSE);
    }
#endif    

#ifdef USE_SERIAL_IO
    if (m->master)
    {
      if ((params->writeTau & EVENT_TRANSITION) > 0) m->tauOff[2] = ftello(m->taufp[2]);
      if ((params->writeSigma & EVENT_TRANSITION) > 0) m->sigmaOff[2] = ftello(m->sigmafp[2]);
      if ((params->writeSlip & EVENT_TRANSITION) > 0) m->slipOff[2] = ftello(m->slipfp[2]);
      if ((params->writeSlipSpeed & EVENT_TRANSITION) > 0) m->slipSpeedOff[2] = ftello(m->slipSpeedfp[2]);
      if ((params->writeState & EVENT_TRANSITION) > 0) m->stateOff[2] = ftello(m->statefp[2]);
      if ((params->writeTheta & EVENT_TRANSITION) > 0) m->thetaOff[2] = ftello(m->thetafp[2]);
    } /* if (m->master) */
#endif    
    tLastWriteTmp = t;
  } /* if (!(eqs->eq) && m->nRupture > 0) an eq just started */
  
  if (eqs->eq && m->nRupture == 0) /* an eq just stopped */
  {
    eqs->M0 = 0;
    eqs->area = 0;
    nslipped = 0;

/* think over */

    for (i = m->startPatch; i <= m->stopPatch; i++)
    {
      if (eqs->slipped[i] == TRUE)
      {
        eqs->dd[nslipped] = m->p[i].d - eqs->dsv[i];
        eqs->dtau[nslipped] = m->p[i].tau - eqs->tausv[i];
        eqs->dsigma[nslipped] = m->p[i].sigma - eqs->sigmasv[i];
        eqs->dmu[nslipped] = m->p[i].tau/m->p[i].sigma - eqs->musv[i];
        eqs->pNumber[nslipped] = i;
        eqs->tauFirstRupture[nslipped] = eqs->tauFirstRupture[i]; /* shouldn't be a problem with overwriting... */
        eqs->sigmaFirstRupture[nslipped] = eqs->sigmaFirstRupture[i]; /* shouldn't be a problem with overwriting... */
        eqs->M0 += m->p[i].op.area*eqs->dd[nslipped]; 
        eqs->area += m->p[i].op.area; /* L and W are in meters */
        eqs->slipped[i] = FALSE;
        nslipped++;
      }
    }
    for (i = 0; i < m->N; i++) eqs->slipped[i] = FALSE;  /* reset slipped for all patches, not just local ones */

#ifdef USE_MPI    
    MPI_Reduce(&(eqs->M0), &M0, 1, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
    MPI_Reduce(&(eqs->area), &area, 1, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
#else
    M0 = eqs->M0;
    area = eqs->area;
#endif

    
    eqs->dt = t - eqs->t0;

    eqs->eq = FALSE;
    m->areaRupture = 0.0;
    
    if (m->master) fprintf(stderr, "event number %ld just finished\n", eqs->ieq); 

#ifdef USE_SERIAL_IO
    if ((params->writeTau & EVENT) > 0) writeTau(m, t, 1, FALSE);
    if ((params->writeSigma & EVENT) > 0) writeSigma(m, t, 1, FALSE);
    if ((params->writeSlip & EVENT) > 0) writeSlip(m, t, 1, FALSE);
    if ((params->writeSlipSpeed & EVENT) > 0) writeSlipSpeed(m, t, 1, FALSE);
    if ((params->writeState & EVENT) > 0) writeState(m, t, 1, FALSE);
    if ((params->writeTheta & EVENT) > 0) writeTheta(m, t, 1, FALSE);
#else
      if ((params->writeTau & EVENT) > 0) {

                copy_double(m,t, nelements_trans,  index,  array_d,1);

                write1D_static(t, (m->taufp[1]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_tau_trans1,  buffer_tau_trans1,
                tauDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
        }
       if ((params->writeSigma & EVENT) > 0) {

                copy_double(m,t, nelements_trans,  index,  array_d,2);

                write1D_static(t, (m->sigmafp[1]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_sigma_trans1,  buffer_sigma_trans1,
                sigmaDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
        }
       if ((params->writeSlip & EVENT) > 0) {

                copy_double(m,t, nelements_trans,  index,  array_d,3);

                write1D_static(t, (m->slipfp[1]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_slip_trans1,  buffer_slip_trans1,
                slipDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
        }
       if ((params->writeSlipSpeed & EVENT) > 0) {

                copy_double(m,t, nelements_trans,  index,  array_d,4);

                write1D_static(t, (m->slipSpeedfp[1]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_ddot_trans1,  buffer_ddot_trans1,
                ddotDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
        }
       if ((params->writeTheta & EVENT) > 0) {

                copy_double(m,t, nelements_trans,  index,  array_d,5);

                write1D_static(t, (m->thetafp[1]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_theta_trans1,  buffer_theta_trans1,
                thetaDisp_trans1, MPI_DOUBLE,filetype_d,FALSE);
        }

//state
       if ((params->writeState & EVENT) > 0) {

                copy_char(m,t, nelements_trans,  index,  array_i,6);

                write1D_static(t, (m->statefp[1]),  array_i, nelements_trans,nrec_trans,
                maxbuffer2,icount_state_trans1,  buffer_state_trans1,
                stateDisp_trans1, MPI_CHAR,filetype_i,FALSE);
        }

#endif    

#ifdef USE_SERIAL_IO
    if (params->writePED) writePED(m, eqs, nslipped);
#else
   if (params->writePED)  writePED_xx(m, eqs, nslipped, MPI_COMM_WORLD,
                disp, maxbuffer1, icount1, ibuff, nrec,
                tpmap, blocksize, nelms, FALSE,
                buffer_d,
                buffer_tau,
                buffer_sigma,
                buffer_mu,
                buffer_e,
                buffer_p,
                buffer_t,
                buffer_taup);

#endif

    if (m->master)
    {
      /* write out some basic info on this event */
      M0 = M0*params->lameMu*1e6; /* 1e6 because lameMu is in MPa */

     if (*icount < buffer[0].maxx){
        buffer[*icount].t0= eqs->t0;
        buffer[*icount].M0= M0;
        buffer[*icount].Mw= (log10(M0) - 9.1)/1.5;
        buffer[*icount].x = m->p[eqs->hypocenter].op.center[0];
        buffer[*icount].y = m->p[eqs->hypocenter].op.center[1];
        buffer[*icount].z = m->p[eqs->hypocenter].op.center[2];
        buffer[*icount].area = area;
        buffer[*icount].dt = eqs->dt;
        buffer[*icount].hypercenter = eqs->hypocenter;
      }

/*	printf(" icount= %d maxx= %d \n",*icount,buffer[0].maxx); */
      if (*icount == (buffer[0].maxx-1)) {

/* 	fprintf(m->eqsfp, "%.20lf %lf %lf %lf %lf %lf %lf %.20lf %d\n", eqs->t0[eqs->ieq], M0, (log10(M0) - 9.1)/1.5,
              m->p[eqs->hypocenter[eqs->ieq]].op.center[0], m->p[eqs->hypocenter[eqs->ieq]].op.center[1],
              m->p[eqs->hypocenter[eqs->ieq]].op.center[2], area, eqs->dt[eqs->ieq], eqs->hypocenter[eqs->ieq]);
      fflush(m->eqsfp);
*/

     for (i=0;i<buffer[0].maxx;i++) {
          fprintf(m->eqsfp, "%.20lf %lf %lf %lf %lf %lf %lf %.20lf %d\n",buffer[i].t0,buffer[i].M0,buffer[i].Mw,
          buffer[i].x,buffer[i].y,buffer[i].z,buffer[i].area,buffer[i].dt,buffer[i].hypercenter);
        }
      fflush(m->eqsfp);
      (*icount)=0;	/* restart saving*/
      }else
      {
      (*icount)++;   /* keep saving */
      }
      
      

      /* if event is too small, roll back the file pointers */
#ifdef USE_SERIAL_IO
      if ((log10(M0) - 9.1)/1.5 < m->minMagWrite) rollBackFilePointers(params, m);
#endif
    } /* if (m->master) */


#ifdef USE_SERIAL_IO

#else
	MPI_Bcast(&M0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if ((log10(M0) - 9.1)/1.5 < m->minMagWrite){

//void rollback_x(MPI_Offset *Disp, MPI_Offset Disp_save, int *iount, int icount_save,int nbyte,int maxbuffer,int nrec);

//tau-double
	int nbyte_d,nbyte_c;
	MPI_Type_size(MPI_DOUBLE,&nbyte_d);
        MPI_Type_size(MPI_CHAR,&nbyte_c);
/*
		if (m->iProc==1){
	 fprintf(stdout,"eqs %d: disp = %ld (disp0= %ld) icount =%d (icout0=%d) nbyte_d=%d maxbuffer2=%d nrec_trans=%d \n", eqs->ieq,
		*tauDisp_trans1,Disp_trans_save[0][1],*icount_tau_trans1,icount_trans_save[0][1],nbyte_d,maxbuffer2,nrec_trans);
	fflush(stdout);
		}
*/

	rollback_x(tauDisp_trans1,Disp_trans_save[0][1],icount_tau_trans1,icount_trans_save[0][1],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(sigmaDisp_trans1,Disp_trans_save[1][1],icount_sigma_trans1,icount_trans_save[1][1],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(slipDisp_trans1,Disp_trans_save[2][1],icount_slip_trans1,icount_trans_save[2][1],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(ddotDisp_trans1,Disp_trans_save[3][1],icount_ddot_trans1,icount_trans_save[3][1],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(thetaDisp_trans1,Disp_trans_save[4][1],icount_theta_trans1,icount_trans_save[4][1],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(stateDisp_trans1,Disp_trans_save[5][1],icount_state_trans1,icount_trans_save[5][1],nbyte_c,maxbuffer2,nrec_trans);

        rollback_x(tauDisp_trans2,Disp_trans_save[0][2],icount_tau_trans2,icount_trans_save[0][2],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(sigmaDisp_trans2,Disp_trans_save[1][2],icount_sigma_trans2,icount_trans_save[1][2],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(slipDisp_trans2,Disp_trans_save[2][2],icount_slip_trans2,icount_trans_save[2][2],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(ddotDisp_trans2,Disp_trans_save[3][2],icount_ddot_trans2,icount_trans_save[3][2],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(thetaDisp_trans2,Disp_trans_save[4][2],icount_theta_trans2,icount_trans_save[4][2],nbyte_d,maxbuffer2,nrec_trans);
        rollback_x(stateDisp_trans2,Disp_trans_save[5][2],icount_state_trans2,icount_trans_save[5][2],nbyte_c,maxbuffer2,nrec_trans);

// finish all roll back
	}

#endif	    
    
    
    (eqs->ieq)++;

    /* reset neighborInRuptureState */
    for (i = m->startPatch; i <= m->stopPatch; i++) m->p[i].neighborInRuptureState = 0;

    if (m->resetStressingRates > 0)
    {
      if (eqs->ieq % m->resetStressingRates == 0)
      {
      /* reset stressing rates, to avoid drift due to roundoff. */

    #ifdef USE_OPENMP
    #pragma omp parallel for shared(m) private(i)
    #endif

        for (i = m->startPatch; i <= m->stopPatch; i++)
        {
          m->p[i].tauDot = m->p[i].tauDotDrive;
          if (!m->p[i].sigmaClipped) m->p[i].sigmaDot = m->p[i].sigmaDotDrive;
        }

        /* update slip speeds on state 3 patches */
    #ifdef USE_OPENMP
    #pragma omp parallel for shared(m) private(i)
    #endif

        for (i = m->startPatch; i <= m->stopPatch; i++) 
        {
          if (m->p[i].state == 3) setDdot3(&(m->p[i]));
        }


        if (m->state3patches) /* for now at least don't do this for models that contain no state 3 patches at all */
        {
          /* update stressing rates */
          /* need to send all slip speeds to all processes 
             FIXME: see note in initFaultModel.c */
      #ifdef USE_MPI
        /* can't get MPI_Allgatherv() to work so doing a series of MPI_Gatherv()s for now, ok with 
           updated version (1.6.4) of OpenMPI,  MPI_Allgatherv() is now working  */
        MPI_Allgatherv(&(m->p[m->startPatch].ddot), 1, m->SendFaultDouble,
                       m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
                       MPI_COMM_WORLD);
    /*    for (root = 0; root < m->nProc; root++)
          MPI_Gatherv(&(m->p[m->startPatch].ddot), 1, m->SendFaultDouble,
                         m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
                         root, MPI_COMM_WORLD); */
      #else
          for (ip = 0;  ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].ddot;
      #endif

          /* adjust tauDot and sigmaDot for state 3 creeping patches */
    #ifdef USE_OPENMP
    #pragma omp parallel for shared(m,eqs) private(i,j)
    #endif
          for (i = m->startPatch; i <= m->stopPatch; i++)
          {
            for (j = 0; j < m->N; j++)
            {
              if (m->p[j].state == 3)
              {
                m->p[i].tauDot += m->p[i].op.Ktau[j]*m->allFaultPropDouble[j];
                if (!m->p[i].sigmaClipped) m->p[i].sigmaDot += m->p[i].op.Ksigma[j]*m->allFaultPropDouble[j];
              }
            }
          }
        } /* if (m->state3patches) */
        
        /* adjust stressing rates for external stressing sources */
        if (m->nExtStressRateSteps > 0)
        {
          if (m->iExtStressRateStep > 0)
          {
            for (i = m->startPatch; i <= m->stopPatch; i++)
            {
              m->p[i].tauDot += m->p[i].tauDotExt[m->iExtStressRateStep - 1];
              m->p[i].sigmaDot += m->p[i].sigmaDotExt[m->iExtStressRateStep - 1];
            }
          }
        }
        
      } /* if (eqs->ieq % m->resetStressingRates == 0) */

    } /* if (m->resetStressingRates > 0) */
    
    
    tLastWriteTmp = t;
  } /* if (eqs->eq && m->nRupture == 0)  an eq just stopped */
  
  if ((m->p[tr.patch].state == 2 || m->p[tr.patch].state == SLOWSLIP_2A) && !eqs->slipped[tr.patch]) 
  {
    eqs->tFirstRupture[tr.patch] = t;
    eqs->tauFirstRupture[tr.patch] = m->p[tr.patch].tau; /* this will be wrong for other procs' patches */
    eqs->sigmaFirstRupture[tr.patch] = m->p[tr.patch].sigma; /* this will be wrong for other procs' patches */
    eqs->slipped[tr.patch] = TRUE;
    m->areaRupture += m->p[tr.patch].op.area;
  }
  
  if ((params->writeTau & TRANSITION) > 0 && ( (eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic) ||
                                               (!eqs->eq && t - m->tLastWrite > m->minDtWriteInterseismic) ))
  {
#ifdef USE_SERIAL_IO
    writeTau(m, t, 0, FALSE);
#else
 	copy_double(m,t, nelements_trans,  index,  array_d,1);

                write1D_static(t, (m->taufp[0]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_tau_trans0,  buffer_tau_trans0,
                tauDisp_trans0, MPI_DOUBLE,filetype_d,FALSE);
#endif
     tLastWriteTmp = t;
  }
  if ((params->writeSigma & TRANSITION) > 0 && ( (eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic) ||
                                               (!eqs->eq && t - m->tLastWrite > m->minDtWriteInterseismic) ))
  {
#ifdef USE_SERIAL_IO
    writeSigma(m, t, 0, FALSE);
#else 
 	copy_double(m,t, nelements_trans,  index,  array_d,2);

                write1D_static(t, (m->sigmafp[0]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_sigma_trans0,  buffer_sigma_trans0,
                sigmaDisp_trans0, MPI_DOUBLE,filetype_d,FALSE);  
#endif
    tLastWriteTmp = t;
  }
  if ((params->writeSlip & TRANSITION) > 0 && ( (eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic) ||
                                               (!eqs->eq && t - m->tLastWrite > m->minDtWriteInterseismic) ))
  {
#ifdef USE_SERIAL_IO
   writeSlip(m, t, 0, FALSE);
#else
	copy_double(m,t, nelements_trans,  index,  array_d,3);

                write1D_static(t, (m->slipfp[0]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_slip_trans0,  buffer_slip_trans0,
                slipDisp_trans0, MPI_DOUBLE,filetype_d,FALSE);
#endif   
    tLastWriteTmp = t;
  }
  if ((params->writeSlipSpeed & TRANSITION) > 0 && ( (eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic) ||
                                               (!eqs->eq && t - m->tLastWrite > m->minDtWriteInterseismic) ))
  {
#ifdef USE_SERIAL_IO
  writeSlipSpeed(m, t, 0, FALSE);
#else
	copy_double(m,t, nelements_trans,  index,  array_d,4);

                write1D_static(t, (m->slipSpeedfp[0]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_ddot_trans0,  buffer_ddot_trans0,
                ddotDisp_trans0, MPI_DOUBLE,filetype_d,FALSE);
#endif   
    tLastWriteTmp = t;
  }
  if ((params->writeState & TRANSITION) > 0 && ( (eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic) ||
                                               (!eqs->eq && t - m->tLastWrite > m->minDtWriteInterseismic) ))
  {
#ifdef USE_SERIAL_IO
   writeState(m, t, 0, FALSE);
#else
	copy_char(m,t, nelements_trans,  index,  array_i,6);

                write1D_static(t, (m->statefp[0]),  array_i, nelements_trans,nrec_trans,
                maxbuffer2,icount_state_trans0,  buffer_state_trans0,
                stateDisp_trans0, MPI_CHAR,filetype_i,FALSE);   
#endif  
    tLastWriteTmp = t;
  }
  if ((params->writeTheta & TRANSITION) > 0 && ( (eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic) ||
                                               (!eqs->eq && t - m->tLastWrite > m->minDtWriteInterseismic) ))
  {
#ifdef USE_SERIAL_IO
    writeTheta(m, t, 0, FALSE);
#else
	copy_double(m,t, nelements_trans,  index,  array_d,5);

                write1D_static(t, (m->thetafp[0]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_theta_trans0,  buffer_theta_trans0,
                thetaDisp_trans0, MPI_DOUBLE,filetype_d,FALSE);
#endif   
    tLastWriteTmp = t;
  }
/*****************************/  

  if (((params->writeTau & EVENT_TRANSITION) > 0) && eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic)
  {
#ifdef USE_SERIAL_IO
   writeTau(m, t, 2, FALSE);
#else
		copy_double(m,t, nelements_trans,  index,  array_d,1);

                write1D_static(t, (m->taufp[2]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_tau_trans2,  buffer_tau_trans2,
                tauDisp_trans2, MPI_DOUBLE,filetype_d,FALSE);
#endif        
    tLastWriteTmp = t;
  }
  if (((params->writeSigma & EVENT_TRANSITION) > 0) && eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic)
  {
#ifdef USE_SERIAL_IO
   writeSigma(m, t, 2, FALSE);
#else
		copy_double(m,t, nelements_trans,  index,  array_d,2);

                write1D_static(t, (m->sigmafp[2]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_sigma_trans2,  buffer_sigma_trans2,
                sigmaDisp_trans2, MPI_DOUBLE,filetype_d,FALSE);
#endif     
    tLastWriteTmp = t;
  }
  if (((params->writeSlip & EVENT_TRANSITION) > 0) && eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic)
  {
#ifdef USE_SERIAL_IO
   writeSlip(m, t, 2, FALSE);
#else
 		copy_double(m,t, nelements_trans,  index,  array_d,3);

                write1D_static(t, (m->slipfp[2]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_slip_trans2,  buffer_slip_trans2,
                slipDisp_trans2, MPI_DOUBLE,filetype_d,FALSE);
#endif    
    tLastWriteTmp = t;
  }
  if (((params->writeSlipSpeed & EVENT_TRANSITION) > 0) && eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic)
  {
#ifdef USE_SERIAL_IO
  writeSlipSpeed(m, t, 2, FALSE);
#else
 		copy_double(m,t, nelements_trans,  index,  array_d,4);

                write1D_static(t, (m->slipSpeedfp[2]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_ddot_trans2,  buffer_ddot_trans2,
                ddotDisp_trans2, MPI_DOUBLE,filetype_d,FALSE);
#endif 
    tLastWriteTmp = t;
  }
  if (((params->writeState & EVENT_TRANSITION) > 0) && eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic)
  {
#ifdef USE_SERIAL_IO
  writeState(m, t, 2, FALSE);
#else
 		copy_char(m,t, nelements_trans,  index,  array_i,6);

                write1D_static(t, (m->statefp[2]),  array_i, nelements_trans,nrec_trans,
                maxbuffer2,icount_state_trans2,  buffer_state_trans2,
                stateDisp_trans2, MPI_CHAR,filetype_i,FALSE);
#endif  
    tLastWriteTmp = t;
  }
  if (((params->writeTheta & EVENT_TRANSITION) > 0) && eqs->eq && t - m->tLastWrite > m->minDtWriteCoseismic)
  {
#ifdef USE_SERIAL_IO
  writeTheta(m, t, 2, FALSE);
#else
		copy_double(m,t, nelements_trans,  index,  array_d,5);

                write1D_static(t, (m->thetafp[2]),  array_d, nelements_trans,nrec_trans,
                maxbuffer2,icount_theta_trans2,  buffer_theta_trans2,
                thetaDisp_trans2, MPI_DOUBLE,filetype_d,FALSE);
#endif   
    tLastWriteTmp = t;
  }
  
  if (m->master)
  {  
    if (params->writeTransitions)
    {
      fwrite(&t, sizeof(t), 1, m->transfp);
      fwrite(&(tr.patch), sizeof(tr.patch), 1, m->transfp);
      fwrite(&(m->p[tr.patch].state), sizeof(m->p[tr.patch].state), 1, m->transfp);
    }
  }
     
  m->tLastWrite = tLastWriteTmp;

#ifdef USE_SERIAL_IO

#else
/*	free(array_d);
	free(array_i);
*/
#endif  
  return;
}
  
void initEqs(Eqs *eqs, FaultModel *m, Params *params)
{
  int i;
  
  eqs->eq = 0;
  eqs->ieq = 0;
  eqs->nEq = params->nEq;
  eqs->slipped = xmalloc(m->N * sizeof(int));
  for (i=0; i<m->N; i++) eqs->slipped[i] = 0;
  eqs->dsv = xmalloc(m->N * sizeof(double));
  eqs->tausv = xmalloc(m->N * sizeof(double));
  eqs->sigmasv = xmalloc(m->N * sizeof(double));
  eqs->musv = xmalloc(m->N * sizeof(double));
  eqs->dd = xmalloc(m->N * sizeof(double));
  eqs->dtau = xmalloc(m->N * sizeof(double));
  eqs->dsigma = xmalloc(m->N * sizeof(double));
  eqs->dmu = xmalloc(m->N * sizeof(double));
  eqs->tFirstRupture = xmalloc(m->N * sizeof(double));
  eqs->tauFirstRupture = xmalloc(m->N * sizeof(double));
  eqs->sigmaFirstRupture = xmalloc(m->N * sizeof(double));
  eqs->pNumber = xmalloc(m->N * sizeof(int));
  eqs->nSlipped = xmalloc(m->nProc * sizeof(int));
  eqs->start = xmalloc(m->nProc * sizeof(int));

 m->memorySize += 2*m->N*sizeof(int) + 5*m->N*sizeof(double) + 2*m->nProc*sizeof(int);
 
  return;
} /* void initEqs(Eqs *eqs, FaultModel *m, Params *params) */

void finishUp(FaultModel *m, Params *params) /* truncate files if necessary */
{  
  if (m->master)
  {
#ifdef USE_SERIAL_IO
    if ((params->writeTau & EVENT) > 0)
    {
      fflush(m->taufp[1]);
      ftruncate(fileno(m->taufp[1]), ftello(m->taufp[1]));
    }
    if ((params->writeSigma & EVENT) > 0)
    {
      fflush(m->sigmafp[1]);
      ftruncate(fileno(m->sigmafp[1]), ftello(m->sigmafp[1]));
    }
    if ((params->writeSlip & EVENT) > 0)
    {
      fflush(m->slipfp[1]);
      ftruncate(fileno(m->slipfp[1]), ftello(m->slipfp[1]));
    }
    if ((params->writeSlipSpeed & EVENT) > 0)
    {
      fflush(m->slipSpeedfp[1]);
      ftruncate(fileno(m->slipSpeedfp[1]), ftello(m->slipSpeedfp[1]));
    }
    if ((params->writeState & EVENT) > 0) 
    {
      fflush(m->statefp[1]);
      ftruncate(fileno(m->statefp[1]), ftello(m->statefp[1]));
    }
    if ((params->writeTheta & EVENT) > 0) 
    {
      fflush(m->thetafp[1]);
      ftruncate(fileno(m->thetafp[1]), ftello(m->thetafp[1]));
    }

    if ((params->writeTau & EVENT_TRANSITION) > 0)
    {
      fflush(m->taufp[2]);
      ftruncate(fileno(m->taufp[2]), ftello(m->taufp[2]));
    }
    if ((params->writeSigma & EVENT_TRANSITION) > 0)
    {
      fflush(m->sigmafp[2]);
      ftruncate(fileno(m->sigmafp[2]), ftello(m->sigmafp[2]));
    }
    if ((params->writeSlip & EVENT_TRANSITION) > 0)
    {
      fflush(m->slipfp[2]);
      ftruncate(fileno(m->slipfp[2]), ftello(m->slipfp[2]));
    }
    if ((params->writeSlipSpeed & EVENT_TRANSITION) > 0)
    {
      fflush(m->slipSpeedfp[2]);
      ftruncate(fileno(m->slipSpeedfp[2]), ftello(m->slipSpeedfp[2]));
    }
    if ((params->writeState & EVENT_TRANSITION) > 0) 
    {
      fflush(m->statefp[2]);
      ftruncate(fileno(m->statefp[2]), ftello(m->statefp[2]));
    }
    if ((params->writeTheta & EVENT_TRANSITION) > 0) 
    {
      fflush(m->thetafp[2]);
      ftruncate(fileno(m->thetafp[2]), ftello(m->thetafp[2]));
    }
#endif  
    
    fprintf(stdout, "Total number of transitions: %d\n", m->nTransitions);
    
  } /* if (m->master) */ 
   
  /* always write out final model state for possible restart */
  writeTau(m, m->t, 3, FALSE);
  writeSigma(m, m->t, 3, FALSE);
  writeSlip(m, m->t, 3, FALSE);
  writeSlipSpeed(m, m->t, 3, FALSE);
  writeState(m, m->t, 3, FALSE);
  writeTheta(m, m->t, 3, FALSE);
  
} /* void finishUp(FaultModel *m, Params *params) truncate files if necessary */


void rollBackFilePointers(Params *params, FaultModel *m)
{
#ifdef USE_SERIAL_IO
  if ((params->writeTau & EVENT) > 0) fseeko(m->taufp[1], m->tauOff[1], SEEK_SET);
  if ((params->writeSigma & EVENT) > 0) fseeko(m->sigmafp[1], m->sigmaOff[1], SEEK_SET);
  if ((params->writeSlip & EVENT) > 0) fseeko(m->slipfp[1], m->slipOff[1], SEEK_SET);
  if ((params->writeSlipSpeed & EVENT) > 0) fseeko(m->slipSpeedfp[1], m->slipSpeedOff[1], SEEK_SET);
  if ((params->writeState & EVENT) > 0) fseeko(m->statefp[1], m->stateOff[1], SEEK_SET);
  if ((params->writeTheta & EVENT) > 0) fseeko(m->thetafp[1], m->thetaOff[1], SEEK_SET);
  if ((params->writeTau & EVENT_TRANSITION) > 0) fseeko(m->taufp[2], m->tauOff[2], SEEK_SET);
  if ((params->writeSigma & EVENT_TRANSITION) > 0) fseeko(m->sigmafp[2], m->sigmaOff[2], SEEK_SET);
  if ((params->writeSlip & EVENT_TRANSITION) > 0) fseeko(m->slipfp[2], m->slipOff[2], SEEK_SET);
  if ((params->writeSlipSpeed & EVENT_TRANSITION) > 0) fseeko(m->slipSpeedfp[2], m->slipSpeedOff[2], SEEK_SET);
  if ((params->writeState & EVENT_TRANSITION) > 0) fseeko(m->statefp[2], m->stateOff[2], SEEK_SET);
  if ((params->writeTheta & EVENT_TRANSITION) > 0) fseeko(m->thetafp[2], m->thetaOff[2], SEEK_SET);
#endif
  return;
} /* void rollBackFilePointers(Params *params, FaultModel *m) */



void writePED(FaultModel *m, Eqs *eqs, int nslipped)
{
  int ip, nTotalSlipped;
  double *dd, *dtau, *dsigma, *dmu, *tauFirstRupture, *sigmaFirstRupture;

#ifdef USE_SERIAL_IO  
  if (m->master)
  {
    dd = xmalloc(m->N * sizeof(double));
    dtau = xmalloc(m->N * sizeof(double));
    dsigma = xmalloc(m->N * sizeof(double));
    dmu = xmalloc(m->N * sizeof(double));
    tauFirstRupture = xmalloc(m->N * sizeof(double));
    sigmaFirstRupture = xmalloc(m->N * sizeof(double));

  }

#ifdef USE_MPI 
  MPI_Gather(&nslipped, 1, MPI_INT, eqs->nSlipped, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (m->master)
  { 
    eqs->start[0] = 0;
    for (ip=1; ip < m->nProc; ip++) eqs->start[ip] = eqs->start[ip-1] + eqs->nSlipped[ip-1];
    nTotalSlipped = eqs->start[m->nProc - 1] + eqs->nSlipped[m->nProc - 1];
  }
  
  MPI_Gatherv(eqs->dd, nslipped, MPI_DOUBLE, 
              dd, eqs->nSlipped, eqs->start, MPI_DOUBLE, 
              0, MPI_COMM_WORLD); 
  MPI_Gatherv(eqs->dtau, nslipped, MPI_DOUBLE, 
              dtau, eqs->nSlipped, eqs->start, MPI_DOUBLE, 
              0, MPI_COMM_WORLD); 
  MPI_Gatherv(eqs->dsigma, nslipped, MPI_DOUBLE,
              dsigma, eqs->nSlipped, eqs->start, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
  MPI_Gatherv(eqs->dmu, nslipped, MPI_DOUBLE,
              dmu, eqs->nSlipped, eqs->start, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
  MPI_Gatherv(eqs->pNumber, nslipped, MPI_INT,
              m->allFaultPropInt, eqs->nSlipped, eqs->start, MPI_INT,
              0, MPI_COMM_WORLD);
  MPI_Gatherv(eqs->tauFirstRupture, nslipped, MPI_DOUBLE,
              tauFirstRupture, eqs->nSlipped, eqs->start, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
  MPI_Gatherv(eqs->sigmaFirstRupture, nslipped, MPI_DOUBLE,
              sigmaFirstRupture, eqs->nSlipped, eqs->start, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

#else
  nTotalSlipped = nslipped;
  eqs->nSlipped[0] = nslipped;
  for (ip = 0; ip < nslipped; ip++) 
  {
    m->allFaultPropInt[ip] = eqs->pNumber[ip];
    dd[ip] = eqs->dd[ip];
    dtau[ip] = eqs->dtau[ip];
    dsigma[ip] = eqs->dsigma[ip];
    dmu[ip] = eqs->dmu[ip];
    tauFirstRupture[ip] = eqs->tauFirstRupture[ip];
    sigmaFirstRupture[ip] = eqs->sigmaFirstRupture[ip];
  }
#endif
              
  if (m->master)
  {
    for (ip=0; ip < nTotalSlipped; ip++)
    {
      fprintf(m->eListfp, "%d\n", eqs->ieq + 1);  /* + 1 is so that event numbers match R's numbers */
      fprintf(m->pListfp, "%d\n", m->allFaultPropInt[ip] + 1); /* + 1 is so that patch numbers match R's numbers */
      fprintf(m->dListfp, "%.5lf\n", dd[ip]);
      fprintf(m->dtauListfp, "%.5lf\n", dtau[ip]);
      fprintf(m->dsigmaListfp, "%.5lf\n", dsigma[ip]);
      fprintf(m->dmuListfp, "%.5lf\n", dmu[ip]);
      fprintf(m->tListfp, "%.5lf\n", eqs->tFirstRupture[m->allFaultPropInt[ip]]);
      fprintf(m->taupListfp, "%.5lf\n", tauFirstRupture[ip]);
      fprintf(m->sigmapListfp, "%.5lf\n", sigmaFirstRupture[ip]);
    }
    fflush(m->eListfp);
    fflush(m->pListfp);
    fflush(m->dListfp);
    fflush(m->dtauListfp);
    fflush(m->dsigmaListfp);
    fflush(m->dmuListfp);
    fflush(m->tListfp);
    fflush(m->taupListfp);
    fflush(m->sigmapListfp);
  }
  
  if (m->master)
  {
    free(dd);
    free(dtau);
    free(dsigma);
    free(dmu);
    free(tauFirstRupture);
    free(sigmaFirstRupture);
  }
#endif
  return;
}
