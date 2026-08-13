#include "RSQSim.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif


void writeTau(FaultModel *m, double t, int fpNum, int setOffset)
{
  int ip;
#ifdef USE_SERIAL_IO  
#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].tau), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (ip = 0;  ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].tau;
#endif

  if (m->master)
  {           
    if (setOffset) m->tauOff[fpNum] = ftello(m->taufp[fpNum]);
    fwrite(&t, sizeof(t), 1, m->taufp[fpNum]);
    if (m->npatchWrite == m->N)
      fwrite(m->allFaultPropDouble, sizeof(double), m->N, m->taufp[fpNum]);
    else
    {
      for (ip = 0; ip < m->N; ip++)
        if (m->p[ip].writePatch !=0) fwrite(&(m->allFaultPropDouble[ip]), sizeof(double), 1, m->taufp[fpNum]);
    }
  }
#endif  
  return;
}

void writeSigma(FaultModel *m, double t, int fpNum, int setOffset)
{
  int ip;
#ifdef USE_SERIAL_IO 
#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].sigma), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (ip = 0; ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].sigma;
#endif

  if (m->master)
  {           
    if (setOffset) m->sigmaOff[fpNum] = ftello(m->sigmafp[fpNum]);
    fwrite(&t, sizeof(t), 1, m->sigmafp[fpNum]);
    if (m->npatchWrite == m->N)
      fwrite(m->allFaultPropDouble, sizeof(double), m->N, m->sigmafp[fpNum]);
    else
    {
      for (ip = 0; ip < m->N; ip++)
        if (m->p[ip].writePatch !=0) fwrite(&(m->allFaultPropDouble[ip]), sizeof(double), 1, m->sigmafp[fpNum]);
    }
  }
#endif
  return;
}

void writeSlip(FaultModel *m, double t, int fpNum, int setOffset)
{
  int ip;
#ifdef USE_SERIAL_IO
#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].d), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (ip = 0; ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].d;
#endif

  if (m->master)
  {           
    if (setOffset) m->slipOff[fpNum] = ftello(m->slipfp[fpNum]);
    fwrite(&t, sizeof(t), 1, m->slipfp[fpNum]);
    if (m->npatchWrite == m->N)
      fwrite(m->allFaultPropDouble, sizeof(double), m->N, m->slipfp[fpNum]);
    else
    {
      for (ip = 0; ip < m->N; ip++)
        if (m->p[ip].writePatch !=0) fwrite(&(m->allFaultPropDouble[ip]), sizeof(double), 1, m->slipfp[fpNum]);
    }
  }
#endif
  return;
}

void writeSlipSpeed(FaultModel *m, double t, int fpNum, int setOffset)
{
  int ip;
#ifdef USE_SERIAL_IO  
#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].ddot), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (ip = 0; ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].ddot;
#endif

  if (m->master)
  {           
    if (setOffset) m->slipSpeedOff[fpNum] = ftello(m->slipSpeedfp[fpNum]);
    fwrite(&t, sizeof(t), 1, m->slipSpeedfp[fpNum]);
    if (m->npatchWrite == m->N)
      fwrite(m->allFaultPropDouble, sizeof(double), m->N, m->slipSpeedfp[fpNum]);
    else
    {
      for (ip = 0; ip < m->N; ip++)
        if (m->p[ip].writePatch !=0) fwrite(&(m->allFaultPropDouble[ip]), sizeof(double), 1, m->slipSpeedfp[fpNum]);
    }
  }
#endif 
  return;
}

void writeState(FaultModel *m, double t, int fpNum, int setOffset)
{
  int ip;
#ifdef USE_SERIAL_IO 
#ifdef USE_MPI
  int err;
  err = MPI_Gatherv(&(m->p[m->startPatch].state), 1, m->SendFaultChar,
                     m->allFaultPropChar, m->nPatches, m->startPatches, MPI_CHAR,
                     0, MPI_COMM_WORLD);
#else
  for (ip = 0; ip < m->N; ip++) m->allFaultPropChar[ip] = m->p[ip].state;
#endif

  if (m->master)
  {           
    if (setOffset) m->stateOff[fpNum] = ftello(m->statefp[fpNum]);
    fwrite(&t, sizeof(t), 1, m->statefp[fpNum]);
    if (m->npatchWrite == m->N)
      fwrite(m->allFaultPropChar, sizeof(char), m->N, m->statefp[fpNum]);
    else
    {
      for (ip = 0; ip < m->N; ip++)
        if (m->p[ip].writePatch !=0) fwrite(&(m->allFaultPropChar[ip]), sizeof(char), 1, m->statefp[fpNum]);
    }
  }
#endif  
  return;
}

void writeTheta(FaultModel *m, double t, int fpNum, int setOffset)
{
  int ip;
#ifdef USE_SERIAL_IO 
#ifdef USE_MPI
  MPI_Gatherv(&(m->p[m->startPatch].theta), 1, m->SendFaultDouble,
              m->allFaultPropDouble, m->nPatches, m->startPatches, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
#else
  for (ip = 0; ip < m->N; ip++) m->allFaultPropDouble[ip] = m->p[ip].theta;
#endif

  if (m->master)
  {           
    if (setOffset) m->thetaOff[fpNum] = ftello(m->thetafp[fpNum]);
    fwrite(&t, sizeof(t), 1, m->thetafp[fpNum]);
    if (m->npatchWrite == m->N)
      fwrite(m->allFaultPropDouble, sizeof(double), m->N, m->thetafp[fpNum]);
    else
    {
      for (ip = 0; ip < m->N; ip++)
        if (m->p[ip].writePatch !=0) fwrite(&(m->allFaultPropDouble[ip]), sizeof(double), 1, m->thetafp[fpNum]);
    }
  }
#endif 
  return;
}
