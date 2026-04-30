
typedef struct
{
  KRDOkadaFaultPatch op;

  double tau; /* current shear stress at center of patch, in prescribed slip direction */
  double sigma;  /* current normal stress at center of patch */ 
  double theta;  /* current value of state variable */
  double d;      /* total accumulated slip */
  double ddot;   /* current slip speed */
  double mu0;   /*  steady-state friction at ddot = ddotStar */
  double ddotStar; /* normalizing constant in constitutive law */
  double alpha; /* normal stress state variable coeff, see Linker and Dieterich, 1992 */
  double ddotEQ; /* slip speed during state 2 */
  double ddotAB; /* slip speed at boundary between state 2a and 2b */
  double ddot2A; /* slip speed that defines the lower boundary of state 2a */
  double A;      /* rate coeff */
  double Asv;    /* nominal rate coeff; A = Asv unless a neighboring patch is in state 2 and this patch is in state 1 */
  double B;      /* state coeff */
  double Dc;      /* characteristic slip distance in state evolution law */
  double iDc;  /* Inverse of Dc*/
  double H;       /* B/Dc + Kself/sigma */
  double tauDot; /* current value of the shear stress rate of change */
  double sigmaDot; /* current value of the normal stress rate of change */
  double sigmaPin; /* if the normal stress drops below this, the patch will be pinned or sigma will be clipped */
  double muPin;   /* if this is greater than zero, pin any element whose tau/sigma exceeds this */
  double dDotDrive;
  double tauDotDrive;   /* the tectonic (shear) stressing rate */
  double sigmaDotDrive; /* the tectonic (normal) stressing rate */
  double Kself; /* Ktau of the element itself, repeated from op.Ktau[] for convenience */
  double tFail;  /* an optional single a priori failure time - the shear stress will be 
                    instantaneously increased to failure stress (i.e. that which
                    makes the slip speed increase to ddotEQ) at this time */
  double tauFail; /* an optional single a priori failure stress - the shear stress will be 
                    instantaneously increased to this failure stress and the
                    state adjusted to make the slip speed increase to ddotEQ at 
                    this time */
  double muTrans2; /* the mu at which this patch will transition 2->0, equal to muSS(ddotEQ) -
                     sOF*(tau0/sigma0 - muSS(ddotEQ)) */
  double muTrans2Low;  /* values of tau/sigma at which state SLOWSLIP_2x patches */
  double muTrans2High; /* will have their slip speeds updated or be transitioned */
  double muTrans2aLow; /* need a separate one for 2a because we save it according to the highest stress
                          reached while a slow-slip patch is in state 2x */
  double muTrans3Low;  /* values of tau/sigma at which state 3 patch creep */
  double muTrans3High; /* speeds will be updated */
  char state;    /* 0, 1, or 2 */
  int APFail;   /* 0 = a priori failure time has not yet arrived
                   1 = a priori failure time has just been been triggered, time to set tau
                   2 = tau has already been set at the a priori failure time */
  int pinned;  /* 0 = patch is not pinned, i.e. allowed to slip
                  1 = patch is pinned, i.e. not allowed to slip */
  int *neighbors;  /* a list of which patches are "neighbors" of the given patch, while the 
                      given patch is in state 1 and any of the neighbors are in state 2, A 
                      will be reduced to fA*A on this patch */
  int nNeighbors;  /* total number of neighbors */
  int neighborInRuptureState; /* non-zero if any of its neighbors have been in state 2 during
                           the current event, otherwise zero */
  int aReduced;  /* has a been reduced (during rupture), 0 if not, 1 if so */
  int slowSlip;  /* is this a slow-slip patch?  0 if not, 1 if so */
  int writePatch; /* do this patch's faultVals get written out? 0 if not, 1 if so */
  int sigmaClipped; /* has this patch's sigma been clipped at sigmaPin? */
  double t01; /* if patch is in state 0 and has already calculated a 0 -> 1 transition time, it is stored here */
  double *tauDotExt; /* if there is an external stressing history supplied, its rates will be stored here */
  double *sigmaDotExt; /* as above */
} FaultPatch;

typedef struct
{
  FaultPatch *p;
  int N;           /* number of fault patches */
  int npatchWrite; /* number of fault patches whose faultVals get written out */
  double ddotEQ;    /* slip speed during an eq */
  double stressOvershootFactor; /* should be > 0 and < 1, Crude way of implementing dynamic overshoot
				.  Elements will stay in state 2 until
				   the shear stress drops to tauSS(ddotEQ) minus this
                                   fraction of the difference between tau/sigma upon entering
                                   state 2 and tauSS(ddotEQ) */
  int stressRateSpecification; /* BACKSLIP or DIRECT */
  int nRupture;    /* number of patches currently in state 2 */
  double areaRupture;  /* total area of elements which have ruptured in the current event */
  double fA;       /* factor by which to multiply A during rupture propagation */
  double dMu3;     /* state 3 (creeping) patches will have their slip speeds updated
                      when tau/sigma has changed by this much */
#ifdef USE_SERIAL_IO
  FILE *taufp[4];     /* files to which to write tau     transitions -> fp[0] */
  FILE *sigmafp[4];   /* files to which to write sigma   events -> fp[1]*/
  FILE *slipfp[4];    /* files to which to write slip    transitions during events -> fp[2] */
  FILE *slipSpeedfp[4];    /* files to which to write slip    transitions during events -> fp[2] */
  FILE *statefp[4];   /* files to which to write state */
  FILE *thetafp[4];   /* files to which to write theta */
#else
  MPI_File taufp[4];     /* files to which to write tau     transitions -> fp[0] */
  MPI_File sigmafp[4];   /* files to which to write sigma   events -> fp[1] */
  MPI_File slipfp[4];    /* files to which to write slip    transitions during events -> fp[2] */
  MPI_File slipSpeedfp[4];    /* files to which to write slip    transitions during events -> fp[2] */
  MPI_File statefp[4];   /* files to which to write state */
  MPI_File thetafp[4];   /* files to which to write theta */
#endif

  FILE *transfp;   /* files to which to write transitions */
  off_t tauOff[3];  /* places to store the file offsets before writing info for */
  off_t sigmaOff[3];/* an event, in case need to roll back the stream because event */
  off_t slipOff[3]; /* ended up being too small (< minMagWrite).  Don't actually use */
  off_t slipSpeedOff[3]; /* ended up being too small (< minMagWrite).  Don't actually use */
  off_t stateOff[3];/* the [0] element, since if you ask for every transition, it will do so */
  off_t thetaOff[3];/*  */
  FILE *eqsfp;      /* file to which to write eq info */

#ifdef USE_SERIAL_IO
  FILE *eListfp; /* the eList, pList, and dList and tList files provide a list of which elements slipped in */
  FILE *pListfp; /* each event and how much and when they first ruptured.  they are ascii files */
  FILE *dListfp; /*  */
  FILE *dtauListfp; /*  */
  FILE *dsigmaListfp; /*  */
  FILE *dmuListfp; /*  */
  FILE *tListfp; /*  */
  FILE *taupListfp; /* peak shear stress during rupture (actually, shear stress at time of first rupture) */
  FILE *sigmapListfp; /* normal stress at time of first rupture */
#else
  MPI_File eListfp; /* the eList, pList, and dList and tList files provide a list of which elements slipped in */
  MPI_File pListfp; /* each event and how much and when they first ruptured.  they are ascii files */
  MPI_File dListfp; /*  */
  MPI_File dtauListfp; /*  */
  MPI_File dsigmaListfp; /*  */
  MPI_File dmuListfp; /*  */
  MPI_File tListfp; /*  */
  MPI_File taupListfp; /* peak shear stress during rupture (actually, shear stress at time of first rupture) */
#endif

  int nTransitions; /* total number of transitions since model start */
  double minDtWriteCoseismic;  /* minimum time (s) between model info being written out during and eq */
  double minDtWriteInterseismic;  /* minimum time (s) between model info being written out in between eqs */
  double tLastWrite;  /* last time model info was written out */
  double minMagWrite;  /* min magnitude to write out info (basic catalog info gets written for all events) */
  double maxThetaPin; /* see explanation in definition of parameters */
  int DEBUG; /* whether to print out debugging info */
  double ZBrentUpperBracket; /* see explanation in RSQSim_params.h */
  double ZBrentTol;
  int master; /* is this process the master process */
  int nProc, iProc;  /* number of processes and process number */
  int startPatch, stopPatch;  /* which patches this process is responsible for */
  double startTime, endTime; /* values from MPI_Wtime() */
  int maxDt1Zero; /* if the total number of T1.dt == 0 conditions on any one process exceeds this, will quit */

/*
 profiling and using OPENMP*
 *
 */

  double initializationTime;
  double findTransitionTime;
  double timeStepTime;
  double transitionTime;
  double updatingEqsTime;
  double allReduceTime;
  double zbrentTime;
  long memorySize;
  long nIter;
  long numCalls;
  long nnei;
  long cnei;
/*
 * end of profiling
 */


#ifdef USE_MPI
  MPI_Datatype SendFaultDouble,  /* MPI constructed data types for sending */
               SendFaultInt,     /*             fault patch properties */
               SendFaultChar;    
#endif
  int *nPatches;  /* number of patches for each process */
  int *startPatches; /* the startPatches for all the processes */
  double *allFaultPropDouble; /* places to temporarily store a given property  */
  int *allFaultPropInt;       /* for every patch */
  char *allFaultPropChar; 
  int maxNeighbors; /* the maximum number of neighbors for any patch */    
  int state3patches; /* do any rate-strengthening patches exist in this model */  
  int lowSigmaAction; /* what to do if a patch's sigma drops below sigmaFracPin*sigma0 */
  int slipInState1; /* do patches slip in state 1 */
  int resetStressingRates; /* are stressing rates reset after every event? */
  double t; /* current simulation time */
  double maxTimeStep; /* if next transition time is longer than this, will take an intermediate step */
  double nExtStressRateSteps; /* the number of external stressing rate changes */
  double *tExtStressRateSteps; /* the times of the external stressing rate changes */
  int iExtStressRateStep; /* which external stress rate step is the next up */
  int unreduceAFix;  /* should be either THETA or SLIPSPEED */
  
  /* buffers for storing snapshots */
  int nbuffSnap, ibuffSnap[3]; /* the number of snapshots the snapshot buffers can hold and the number stored so far */
  double *taubuffer;
  
  double aReduceNmax; /* if greater than zero, the a reduction during rupture is dependent on 
                         the number of elements currently rupturing: a = a0 * sqrt(nRupt/aReduceNmax)
                         or a = a0 if nRupt >= aReduceNmax */
  double fAmin, aReduceN1, aReduceN2; /* alternative to aReduceNmax.  When a neighboring patch is
                                         rupturing, a nucleating patch will have its a parameter
                                         adjusted to f*a, where f is:
                                           fAmin if the number of currently rupturing patches is less than aReduceN1
                                           1 if the number of currently rupturing patches is greater than aReduceN2
                                           fAmin + (1-fAmin)*fAfunc((n - aReduceN1)/(aReduceN2 - aReduceN1)) otherwise */
  
} FaultModel;

typedef struct
{
  double dt;    /* time to transition */
  int patch;   /* patch which undergoes a transition after dt */
  int nextState; /* state to which the transitioning patch is going to transition */
} Transition;

typedef struct
{
  int eq;      /* is an eq underway */
  long int ieq;     /* current event number */
  double nEq;     /* total number of events for which the simulation will be run */
  double t0;   /* origin time of current event */
  double M0;   /* moment of current event */
  int hypocenter; /* patch index of hypocenter of current event */
  double area;    /* area of fault patches which slipped in current event */
  double dt;   /* duration of slipping of current event */
  double *dsv;   /* slip of each element at beginning of an eq */
  double *tausv;   /* shear stress on each element at beginning of an eq */
  double *sigmasv;   /* normal stress on each element at beginning of an eq */
  double *musv;   /* shear to normal stress ratio on each element at beginning of an eq */
  double *dd;   /* event slip of each element that entered state 2 during an event */
  double *dtau;   /* change in shear stress of each element that entered state 2 during an event */
  double *dsigma;   /* change in normal stress of each element that entered state 2 during an event */
  double *dmu;   /* change in shear to normal stress ratio of each element that entered state 2 during current event */
  double *tFirstRupture;   /* time at which each element first entered state 2 during an event */
  double *tauFirstRupture; /* shear stress at time at which each element first entered state 2 during current event */
  double *sigmaFirstRupture; /* normal stress at time at which each element first entered state 2 during current event */
  int *pNumber;   /* patch number to which the slips in dd correspond */
  int *slipped; /* slipped[i] is 0 or 1 according to whether the i^th patch
                   has slipped in the earthquake currently underway */
  int *nSlipped;  /* the number of patches that slipped in the current event on each process (only meaningful on master) */
  int *start; /* the cumulative sum of nSlipped (only meaningful on master) */
} Eqs;


typedef struct{
        double t0;
        double M0;
        double Mw;
        double x;
        double y;
        double z;
        double area;
        double dt;
        int hypercenter;
        int maxx;

}eqs_info;


