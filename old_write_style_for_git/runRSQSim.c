#include "RSQSim.h"
#ifdef USE_OPENMP
#include <omp.h>
#include <math.h>
#endif

static char svnId[] = {"$Id: runRSQSim.c 796 2014-07-25 00:44:36Z dinger $"};

int main(int argc, char **argv)
{
  Eqs eqs;
  FaultModel m;
	int k;
/*
 *
 */
  FILE *screenout,*times_saved;
  m.memorySize=0;
  char tmpfname[RSQSIM_MAXLINE];
/*
 *
 */

  Params params;
  
  m.memorySize += sizeof(Eqs) + sizeof(FaultModel) + sizeof(Params);
 
  
  void *parambuffer;
  int bufsize;
  char comment[] = {"#"};  /* comment char in faultFname */
  
  time_t date;
  char hostname[RSQSIM_MAXLINE], workingDir[RSQSIM_MAXLINE];
  int iarg;
  
  void setUpDataTypes(FaultModel *m);
  setUpMPI(&m); /* if compiled with -DUSE_MPI this will start up MPI */

 int tid,omp_num_procs,itemp;

#ifdef USE_OPENMP

/* ---llnl openMP
 *   Fork a team of threads with each thread having a private tid variable */
#pragma omp parallel private(tid)
  {

  /* Obtain and print thread id */
  tid = omp_get_thread_num();

  /* Only master thread does this */
  if (tid == 0)
    {
    omp_num_procs = omp_get_num_threads();
    itemp = omp_get_num_procs();

//         printf("Number of threads = %d\n",omp_num_procs);
     }
   }  /* All threads join master thread and terminate */
#endif

  getParams(argc, argv, &params, &m);
  if (m.master) 
  {

//    getParams(argc, argv, &params, &m);

    sprintf(tmpfname, "%s.Screen.Time.out", params.outFnameInfix);
    screenout = xfopen(tmpfname, "w");
    sprintf(tmpfname, "%s.times_saved.out", params.outFnameInfix);
    times_saved = xfopen(tmpfname, "w");
   
    fprintf(screenout,"\tThere are %d MPI processes\n", m.nProc);
#ifdef USE_OPENMP
    fprintf(screenout,"\tThere are %d processors per node\n", itemp);
    fprintf(screenout,"\tThere are %d OpenMP THREADS   per process\n", omp_num_procs);
#endif
    fflush(screenout);
    
    m.N = nlines(params.faultFname, comment);
    
#ifdef USE_MPI_XX
    packParams(&params, &parambuffer, &bufsize);
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(parambuffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);
    free(parambuffer);
#endif

    /* write out date, command line, working directory, and hostname to main output file */
    date = time(NULL);
    fprintf(m.eqsfp, "StartTime = %s", ctime(&date));
    gethostname(hostname, RSQSIM_MAXLINE);
    fprintf(m.eqsfp, "Hostname = %s\n", hostname);
    getcwd(workingDir, RSQSIM_MAXLINE);
    fprintf(m.eqsfp, "RunDirectory = %s\n", workingDir);
    fprintf(m.eqsfp, "CommandLine = ");
    for (iarg=0; iarg<argc; iarg++) fprintf(m.eqsfp, "%s ", argv[iarg]);
    fprintf(m.eqsfp, "\n");
  }
  else
  {
#ifdef USE_MPI_XX
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    parambuffer = xmalloc(bufsize);
    MPI_Bcast(parambuffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);
    unpackParams(&params, parambuffer, bufsize); /* unpack params */
    free(parambuffer);
#endif 
  }
  
#ifdef USE_MPI
  MPI_Bcast(&(m.N), 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  assignPatches(&m);

  if (m.master) 
  {
    printf("assignPatches() is done.\n");
    fprintf(screenout,"assignPatches() is done.\n");
    for (k=0;k<m.nProc;k++)
      fprintf(screenout," k= %d  from start patch [%d] to stop patch [%d] # of patches %d\n",
	      k,m.startPatches[k],m.startPatches[k]+m.nPatches[k]-1,m.nPatches[k]);
  }
  
  if (params.DEBUG)
  {
    char hostname[RSQSIM_MAXLINE];
    
    gethostname(hostname, RSQSIM_MAXLINE);
    
    fprintf(stdout, "process %d (pid %d on host %s) of %d has patches %d through %d\n", 
            m.iProc, getpid(), hostname, m.nProc, m.startPatch, m.stopPatch);
  }
  
  if (m.master)
  {
    if (params.DEBUG > 10 && params.DEBUG <=15) 
    {
      fprintf(stdout, "DEBUG > 10 [pid = %d]; any key to continue: ", getpid());
      fflush(stdout);
      getc(stdin);
    }
  }
  
  if (params.DEBUG > 15)
  {
    int iwait = 0;
    fprintf(stdout, "process %d (pid %d on host %s) ready for attaching\n", m.iProc, getpid(), hostname);
    while (iwait == 0) 
      sleep(5);
  }
  
  setUpDataTypes(&m);
  if (m.master) {
		printf("setUpDataTypes() is done.\n");
  		fprintf(screenout,"setUpDataTypes() is done.\n");
		fflush(screenout);
	}

#ifdef USE_MPI
	//m.initializationTime=MPI_Wtime()-m.initializationTime;

#endif


  eqs = RSQSim(&params, &m);

#ifdef USE_MPI  
  m.endTime = MPI_Wtime();
  if (m.master) {
		printf("Total wall clock time: %lf second\n", m.endTime - m.startTime);
		fprintf(screenout,"Total wall clock time: %lf second\n", m.endTime - m.startTime);
	}

/*
 *	profiling and uing OPENMP
 */

 double initializationTime=0;
  double findTransitionTime=0;
  double timeStepTime=0;
  double transitionTime=0;
  double updatingEqsTime=0;
  double allReduceTime=0;
  long memorySize=0;


  MPI_Allreduce(&(m.initializationTime),&initializationTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&(m.findTransitionTime),&findTransitionTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&(m.timeStepTime),&timeStepTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&(m.transitionTime),&transitionTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&(m.updatingEqsTime),&updatingEqsTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&(m.allReduceTime),&allReduceTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&(m.memorySize),&memorySize,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  if (m.master) {
	int OMP=1;
#ifdef USE_OPENMP
	    if (m.iProc == 0) printf("\tThere are %d processors per node\n", omp_get_num_procs());
	OMP=omp_get_num_procs();
#endif
		printf("\tInitialization time: %lf second\n", initializationTime);
  		printf("\tTotal find min transition time: %lf second\n", findTransitionTime);
  		printf("\t\tTotal allreduce time: %lf second\n", allReduceTime);
  		printf("\tTotal time step time: %lf second\n", timeStepTime);
  		printf("\tTotal state transition time: %lf second\n", transitionTime);
  		printf("\tTotal updating equation time: %lf second\n", updatingEqsTime);
  		printf("Total memory: %ld MBytes\n", memorySize/(1024*1024));

		fprintf(screenout,"\tInitialization time: %lf second\n", initializationTime);
                fprintf(screenout,"\tTotal find min transition time: %lf second\n", findTransitionTime);
                fprintf(screenout,"\t\tTotal allreduce time: %lf second\n", allReduceTime);
                fprintf(screenout,"\tTotal time step time: %lf second\n", timeStepTime);
                fprintf(screenout,"\tTotal state transition time: %lf second\n", transitionTime);
                fprintf(screenout,"\tTotal updating equation time: %lf second\n", updatingEqsTime);
                fprintf(screenout,"Total memory: %ld MBytes\n", memorySize/(1024*1024));
		fprintf(screenout,"\tThere are %d OpenMP processes per node\n", OMP);

		fprintf(screenout,"\n\tEarthquake information : \n");
		 fprintf(screenout,"\t# of EQS =  %ld \n",eqs.ieq);
		 fprintf(screenout,"\t# nEq = %lf \n",eqs.nEq);

	 fprintf(times_saved,"OMP #Proc #N initializationTime findTransitionTime allReduceTime timeStepTime transitionTime updatingEqsTime memorySize(MB) OMP_pernode TotalTime \n");  
	fprintf(times_saved,"%d %d %d %lf %lf %lf %lf %lf %lf %ld %d %lf\n",OMP,m.nProc,m.N,initializationTime,
			findTransitionTime,allReduceTime,timeStepTime,transitionTime,updatingEqsTime,memorySize/(1024*1024),OMP,m.endTime - m.startTime);
	fflush(times_saved);
	}

/*---------all partition times-----------*/

	double sendarray[6]={	m.initializationTime, m.findTransitionTime,m.allReduceTime,
				m.timeStepTime,m.transitionTime,m.updatingEqsTime};
	double *rbuf;
	FILE * times_details;

	if (m.master) {
		rbuf=(double*)xmalloc(m.nProc*6*sizeof(double));
                sprintf(tmpfname, "%s.TIMES_PARTITION.out", params.outFnameInfix);
		times_details=xfopen(tmpfname, "w");
	}
	MPI_Gather(sendarray,6,MPI_DOUBLE,rbuf,6, MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (m.master){
		for (k=0;k<m.nProc;k++){
	fprintf(times_details,"%d %lf %lf %lf %lf %lf %lf \n",k,rbuf[k*6],rbuf[k*6+1],rbuf[k*6+2],
		rbuf[k*6+3],rbuf[k*6+4],rbuf[k*6+5]);

		}
	}
/*
 *	end of profiling
 */

  MPI_Finalize();
#endif

  return(0);

}   


/* start up MPI if desired (if compiled with -DUSE_MPI, then the macro MPI
   will be TRUE, else it will be FALSE) */
void setUpMPI(FaultModel *m)
{
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &(m->nProc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(m->iProc));
    m->startTime = MPI_Wtime();
/*
 *
 */
    if (m->iProc == 0) printf("\tThere are %d MPI processes\n", m->nProc);
#ifdef USE_OPENMP
    if (m->iProc == 0) printf("\tThere are %d processors per node\n", omp_get_num_procs());
#endif
/*
 *
 */

#else
    m->nProc = 1;
    m->iProc = 0;
#endif

/*
 *
 */
m->initializationTime=m->startTime;
  m->findTransitionTime=0;
  m->timeStepTime=0;
  m->transitionTime=0;
  m->updatingEqsTime=0;
  m->allReduceTime=0;
/*
 *
 */
  
  m->master = (m->iProc == 0);
}


/* figure out which patches each process is responsible for
   note: might want to randomly assign patches instead of
         assigning them in order, especially if there are
         a bunch of pinned patches */
void assignPatches(FaultModel *m)
{
  int nPatchesPerProc, nExtraPatches, nPatchesThisProc;
  
  nPatchesPerProc = m->N / m->nProc;
  nExtraPatches = m->N % m->nProc;
  
  if (m->iProc < nExtraPatches)
  {
    m->startPatch = m->iProc*(nPatchesPerProc + 1);
    m->stopPatch = m->startPatch + nPatchesPerProc;
  }
  else
  {
    m->startPatch = nExtraPatches*(nPatchesPerProc + 1) + (m->iProc - nExtraPatches)*nPatchesPerProc;
    m->stopPatch = m->startPatch + nPatchesPerProc - 1;
  }
  
  nPatchesThisProc = m->stopPatch - m->startPatch + 1;

  m->startPatches = xmalloc(m->nProc*sizeof(int));
  m->nPatches = xmalloc(m->nProc*sizeof(int));

        m->memorySize += m->nProc*sizeof(int)*2;


#ifdef USE_MPI   
  MPI_Allgather(&(m->startPatch), 1, MPI_INT, m->startPatches, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&nPatchesThisProc, 1, MPI_INT, m->nPatches, 1, MPI_INT, MPI_COMM_WORLD);
#else
  m->startPatches[0] = 0;
  m->nPatches[0] = nPatchesThisProc;
#endif
  
  return;
}


/* set up the data types used for passing properties of each patch
   back and forth from the master to the slave processes and
   buffers for the master to use to receive these properties */
void setUpDataTypes(FaultModel *m)
{
  int err;
  
#ifdef USE_MPI
  err = MPI_Type_create_hvector(m->stopPatch - m->startPatch + 1, 1, sizeof(FaultPatch), MPI_DOUBLE,
                         &(m->SendFaultDouble));
  err = MPI_Type_commit(&(m->SendFaultDouble));
  err = MPI_Type_create_hvector(m->stopPatch - m->startPatch + 1, 1, sizeof(FaultPatch), MPI_INT,
                         &(m->SendFaultInt));
  err = MPI_Type_commit(&(m->SendFaultInt));
  err = MPI_Type_create_hvector(m->stopPatch - m->startPatch + 1, 1, sizeof(FaultPatch), MPI_CHAR,
                         &(m->SendFaultChar));
  err = MPI_Type_commit(&(m->SendFaultChar));
#endif
  
  m->allFaultPropDouble = xmalloc((m->N)*sizeof(double));
  m->allFaultPropInt = xmalloc((m->N)*sizeof(int));
  m->allFaultPropChar = xmalloc((m->N)*sizeof(char));

        m->memorySize += m->N*(sizeof(double)+sizeof(int)+sizeof(char));

  return;
}
