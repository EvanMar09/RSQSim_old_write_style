#include "RSQSim.h"

void packParams(Params *params, void **parambuffer, int *bufsize)
{
#ifdef USE_MPI
  int position = 0;
  
  *bufsize = sizeof(Params);  /* upper bound since we're not packing all of params */
  *parambuffer = xmalloc(*bufsize);
  MPI_Pack(&(params->lameLambda),              1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->lameMu),                  1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->nEq),                     1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->tStart),                  1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->maxT),                    1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->maxTrans),                1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->stressRateSpecification), 1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeTau),                1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeSigma),              1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeSlip),               1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeSlipSpeed),          1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeState),              1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeTheta),              1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writePED),                1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->writeStiffness),          1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->DEBUG),                   1, MPI_INT,    *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->KZeroFrac),               1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(params->maxTimeStep),             1, MPI_DOUBLE, *parambuffer, *bufsize, &position, MPI_COMM_WORLD);
#endif

  return;
}

void unpackParams(Params *params, void *parambuffer, int bufsize)
{
#ifdef USE_MPI
  int position = 0;


  MPI_Unpack(parambuffer, bufsize, &position, &(params->lameLambda),              1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->lameMu),                  1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->nEq),                     1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->tStart),                  1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->maxT),                    1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->maxTrans),                1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->stressRateSpecification), 1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeTau),                1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeSigma),              1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeSlip),               1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeSlipSpeed),          1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeState),              1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeTheta),              1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writePED),                1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->writeStiffness),          1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->DEBUG),                   1, MPI_INT,    MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->KZeroFrac),               1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(parambuffer, bufsize, &position, &(params->maxTimeStep),             1, MPI_DOUBLE, MPI_COMM_WORLD);
#endif

  return;
}
