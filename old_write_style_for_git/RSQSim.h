#define SVNID $Id$

#ifndef RSQSIM_H
#define RSQSIM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include <time.h>
#include "KRDOkada.h"
#include "array.h"
#include "xfile.h"
#include "xmalloc.h"
#ifdef USE_MPI
#include "mpi.h" 
#endif

#define RSQSIM_MAXLINE 131072
#ifndef MAX
#define MAX(A,B) ( (A) > (B) ? (A) : (B) )
#endif
#ifndef MIN
#define MIN(A,B) ( (A) < (B) ? (A) : (B) )
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifdef USE_MPI
#define MPI TRUE
#else
#define MPI FALSE
#endif

/* constants for whether/when to write out various quantities (can sum these
   to get multiple files written */
#define NEVER 0             /* self-explanatory */
#define TRANSITION 1        /* at every transition */
#define EVENT 2             /* before and after every event */
#define EVENT_TRANSITION 4  /* at every transition during an event */

/* constants for writing out stiffness matrices */
#define NONE 0
#define KTAU 1
#define KSIGMA 2
#define BOTH 3

/* constants for how stressing rates are determined */
#define BACKSLIP 1  /* last column of fault file gives desired slip rate on each patch */
#define DIRECT 2    /* ignore last column of fault file; tauDotDrive and/or sigmaDotDrive
                       are supplied in files */
#define TENSOR_RATE 3   /* a stressing rate tensor is supplied (in stressRateFname) */

#define ALLK 0 /* whether or not to store all the stiffness matrix coeff. or just
                  Ksigma and Ktau.  For now, we are only using the latter two, so
                  don't need the general 6 K */

/* constants for APFail (a priori failure time) */
#define NOTYET  0
#define NOW     1
#define ALREADY 2  

/* value to indicate that a valid 0 -> 1 time has not been precalculated */
#define T01_NA_VAL -1.0

/* some of the state numbers */
#define LOCKED         0
#define NUCLEATE       1
#define RUPTURE        2
#define CREEP          3
#define LOW_SIGMA      4
#define LOW_TAU        5
#define HIGH_THETA     6
#define SLOWSLIP_2A    7
#define SLOWSLIP_2B    8
#define SLOWSLIP_2C    9
#define HIGH_TAU      10
#define MAX_TIME_STEP 11
#define STRESS_RATE_STEP 12

/* to avoid instantaneous switching back and forth
   transitions between the various slowslip rupture
   states */
# define MU_EPS 1e-8

/* constants for FILE pointers to the various optional input files
   that specify a value for each fault element */
#define DDOTEQ         0
#define INITTAU        1
#define INITSIGMA      2
#define INITTHETA      3
#define INITSLIPSPEED  4
#define FP_A           5
#define FP_B           6
#define DC             7
#define MU0            8
#define DDOTSTAR       9
#define DDOTAB        10
#define ALPHA         11
#define TFAIL         12
#define TAUFAIL       13
#define PINNED        14
#define NEIGHBOR      15
#define SLOWSLIP      16
#define TAUDOT        17
#define SIGMADOT      18
#define WRITEPATCH    19
#define EXTSTRESSRATE 20
#define NFPS 21
                  
/* options for lowSigmaAction */
#define LOCK 0
#define CLIP 1

/* choices for unreduceAFix - what is held fixed when unreducing a 
   (and the other is adjusted to fit the constitutive relation) */
#define THETA 0
#define SLIPSPEED 1
    
#include "RSQSim_params.h"
#include "RSQSim_defaultParams.h"
#include "RSQSim_dataStructures.h"
#include "RSQSim_functionDec.h"
  

#endif /* ndef RSQSIM_H */
