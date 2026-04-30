typedef struct
{
  double A_1;       /* rate coeff (dimensionless), if the same for all fault elements */
  double fA;               /* Reduction factor for A during rupture propagation (A' = fA*A) */ 
  double B_1;       /* state coeff (dimensionless), if the same for all fault elements */
  double Dc_1;      /* critical slip distance (m), if the same for all fault elements */
  double mu0_1;      /* steady-state friction at ddot = ddotStar, if the same for all elements */
  double ddotStar_1;   /* normalizing constant in constitutive law (m/s), if the same for all elements */
  double ddotAB_1;     /* 2a -> 2b transition slip speed (m/s), if the same for all elements */
  double alpha_1;      /*  */
  double theta0_1;  /* inital value of state variable (s), if the same for all elements */
  double tau0_1;    /* initial value of shear stress (MPa), if the same for all elements */
  double sigma0_1;   /* initial normal stress (MPa), if the same for all elements */
  double sigmaFracPin;   /* if the normal stress on a patch drops to this times the 
                            initial normal stress, the patch will be pinned */
  double muPin;   /* if this is greater than zero, then any patch whose tau/sigma exceeds
                     this value will be pinned */
  double maxThetaPin; /* if the state variable theta on any patch exceeds this value + 1e3*t, the
                         patch will be locked */
  double ddotEQ_1;   /* earthquake sliding speed (m/s), if the same for all elements */
  char *ddotEQFname; /* file from which to read ddotEQ for each patch */
  double stressOvershootFactor;   /* should be > 0 and < 1 (and closer to zero than 1), 
                                   Elements will stay in state 2 until
                                   the shear stress drops to tauSS(ddotEQ) - sOF*(tau0 - tauSS(ddotEQ) */
  double lameLambda, lameMu;        /* Lame parameters (MPa) */
  int slowSlip_1;    /* if all patches are regular, this should be 0. If all patches are slow-slip
                      this should be 1.  Otherwise a slipSlipFname will need to be specified
                      (and if so, this param is ignored */
  double nEq;         /* for how many events to run the simulation */
  double maxT;     /* maximum amount of simulated time to run (s) */
  double maxTrans;     /* maximum number of transitions to simulate */
  double tStart;   /* start time of simulation (mostly for use when restarting to continue previous runs) */
  char *faultFname;
  char *outFnameInfix;
  int writeTau; /* whether/when to write out tau for each fault patch */
  int writeSigma; /* whether/when to write out sigma for each fault patch */
  int writeSlip; /* whether/when to write out slip for each fault patch */
  int writeSlipSpeed; /* whether/when to write out slip speed for each fault patch */
  int writeState; /* whether/when to write out state for each fault patch */
  int writeTheta; /* whether/when to write out theta for each fault patch */
  int writePED; /* whether to write the pList, eList, and dList files */
  int writeTransitions; /* whether to write out transitions */
  double minDtWrite; /* only write out tau, sigma, slip, and/or state this often (s) both during or between events (below overrides) */
  double minDtWriteCoseismic; /* only write out tau, sigma, slip, and/or state this often (s) during events */
  double minDtWriteInterseismic; /* only write out tau, sigma, slip, and/or state this often (s) in between events */
  double minMagWrite; /* only write tau, sigma, slip, and/or state for events >= minMagWrite */
  int writeStiffness;  /* which stiffness matrices to write out */
  int stressRateSpecification;  /* how stressing rates are specified */
  double dMu3; /* state 3 (creeping) patches will have their slip speeds updated
                  when tau/sigma has changed by this much */
  char *initTauFname; /* optional file from which to read initial shear stresses
                         should contain a list of initial shear stresses in MPa 
                         in the same order as the faultFname file */
  char *initSigmaFname; /* optional file from which to read initial normal stresses
                         should contain a list of initial normal stresses in MPa 
                         in the same order are the faultFname file */
  char *initThetaFname; /* optional file from which to read initial state
                         should contain a list of initial state in seconds 
                         in the same order are the faultFname file */
  char *initSlipSpeedFname; /* optional file from which to read initial slip speed
                         will only be used for state 3 patches (other patches have
                         their slip speed determined by their initial tau, sigma, and theta)
                         should contain a list of initial state in seconds 
                         in the same order are the faultFname file */
  char *AFname; /* optional file from which to read the rate coeff in the R+S constitutive relation
                         should contain a list of rate coeffs 
                         in the same order are the faultFname file */
  char *BFname; /* optional file from which to read the state coeff in the R+S constitutive relation
                         should contain a list of state coeffs 
                         in the same order are the faultFname file */
  char *DcFname; /* optional file from which to read the critical slip distance in the state 
                    evolution law
                         should contain a list of critical slip distances in meters 
                         in the same order are the faultFname file */
  char *mu0Fname; /* optional file from which to read mu0 in the R+S constitutive relation
                         should contain a list of mu0
                         in the same order as the faultFname file */
  char *ddotStarFname;  /* optional file from which to read ddotStar in the R+S constitutive relation
                           should contain a list of ddotStar in the same order as the patches in the 
                           faultFname file */
  char *ddotABFname;  /* optional file from which to read ddotAB 
                           should contain a list of ddotAB in the same order as the patches in the 
                           faultFname file */
  char *alphaFname; /* optional file from which to read alpha
                         should contain a list of alphas
                         in the same order are the faultFname file */
  char *KTauFname; /* optional file from which to read shear stiffness matrix
                      the i^th line should contain the elements for the i^th
                      patch as the receiving patch */
  char *KSigmaFname; /* optional file from which to read normal stiffness matrix
                      the i^th line should contain the elements for the i^th
                      patch as the receiving patch */
  char *tFailFname; /* optional file from which to read a priori failure times;
                       file should contain list of failure times for each patch
                       in the same order as the faultFname file. 
                       shear stress will be instantaneously incresed to the 
                       failure stress (i.e. that which causes slip speed to
                       increase to ddotEQ) at these times 
                       If a tFailFname is not supplied, there will be no a priori
                       failure times. To leave some patches without a priori failure
                       times enter any negative number on the appropriate lines
                       of the tFailFname file */
  char *tauFailFname; /* optional file from which to read failure stresses.  At the
                         failure times in tFailFname, each patch will have its shear
                         stress set to the stress read from this file and the state,
                         theta, will be adjusted so that the patch fails.  If tauFailFname
                         is not supplied, the state will not be adjusted, and the 
                         shear stress alone will be jumped to bring on failure.  To
                         leave some patches without an a priori failure stress, enter
                         any negative number on the appropriate lines of the tauFailFname
                         file */
  char *tauDotFname;    /* shear stressing rate for each patch, in MPa/s */
  char *sigmaDotFname;  /* normal stressing rate for each patch, in MPa/s */
                        
  char *pinnedFname;  /* optional file from which to read whether or not a patch is allowed
                         to slip.  If not specified, all patches are allowed to slip.  Should
                         be a file with one line per patch, either 0 (for not pinned, allowed
                         to slip) or 1 (pinned, not allowed to slip) */
  char *neighborFname; /* optional file from which read a list of "neighboring" patches for 
                          each patch, i^th line should contain a list of which patches are
                          neighbors of the i^th patch (with indexing starting at 0) */
  char *stressRateFname; /* if stressRateSpecification is TENSOR_RATE, then this file specifies
                            a (uniform) stress-rate tensor which is then projected onto
                            each patch to get the shear and normal stressing rates.  There
                            should be 6 numbers in this file each giving a tensor component
                            of the stress-rate in MPa/s.  Ordering should be as in projectStress.c */
  char *slowSlipFname;  /* optional file from which to read whether each patch is a slow slip
                           patch or not.  Should be a one-column ascii file with each containing
                           either a 0 or a 1 indicating whether the correpsonding patch in the
                           faultFname file is a regular or a slow-slip patch, respectively */
  char *writePatchFname;  /* optional file from which to read whether each patch's faultVals 
                             get written out.  Should be a one-column ascii file with each containing
                             either a 0 or a 1 indicating whether the faultVals of the correpsonding 
                             patch in the faultFname file do not or do get written out, respectively  */
  int DEBUG;          /* whether to print out debugging info */
  double ZBrentUpperBracket;  /* if 0.0, then use minDt found so far, otherwise use this value.  If not
                                 0.0 should be some large value like 1e15 or 1e30. */
  double ZBrentTol;
  int lowSigmaAction;    /* what to do if a patch's sigma drops below sigmaFracPin*sigma0 */
  double KZeroFrac;    /* if not 0.0 then elements of the stiffness matrices for interactions whose
                          coulomb interaction is greater than this fraction of the self coulomb interaction
                          (in absolute value) are set to zero */
  char *KZeroFname;  /* Two column ascii file, 1-based fault element indices, if a line
                        'j i' exists then Ktau[i-1][j-1] and Ksigma[i-1][j-1] will be
                        set to zero */
  int slipInState1;  /* if non-zero then patches slip (through their stress is not transferred to other
                        patches) in state 1.  If zero, then patches are locked in state 1 */
  int resetStressingRates;  /* if non-zero then the stressing rates are reset after every event
                               to avoid drift due to roundoff */
  int pointSource;  /* if zero, then uses Okada or tgf finite Greens functions in the stiffness matrix
                       calculations.  If non-zero, just uses point sources at the centroid of the elements */
  double tgfDist1;  /* if center-to-center distance is less than this multiple of the element size, will use full tgf */
  double tgfDist2;  /* if between tgfDist1 and tgfDist2, will use 12 point sources, and if greater than tgfDist2, 6 pt. sources */
  int nbuffEQS; /* size of earthquake buffer */
  int nbuffPED; /* size of list item buffers */
  int nbuffTRANS; /* size of transition buffer */
  int nbuffSnap; /* size of buffer for snapshots */
  double maxTimeStep; /* will never take a time step bigger than this */
  char *extStressRateFname; /* external stressing history file */
  int unreduceAFix; /* should be either THETA (in which case slip speed is updated) or SLIPSPEED 
                       (in which case theta is updated) */
  double aReduceNmax; /* if greater than zero, the a reduction during rupture is dependent on 
                         the number of elements currently rupturing: a = a0 * sqrt(nRupt/aReduceNmax)
                         or a = a0 if nRupt >= aReduceNmax */
  double fAmin, aReduceN1, aReduceN2; /* alternative to aReduceNmax.  When a neighboring patch is
                                         rupturing, a nucleating patch will have its a parameter
                                         adjusted to f*a, where f is:
                                           fAmin if the number of currently rupturing patches is less than aReduceN1
                                           1 if the number of currently rupturing patches is greater than aReduceN2
                                           fAmin + (1-fAmin)*fAfunc((n - aReduceN1)/(aReduceN2 - aReduceN1)) otherwise */
  int maxDt1Zero; /* if the total number of T1.dt == 0 conditions on any one process exceeds this, will quit */
} Params;

typedef struct
{
  char *name;
  char *type;
  char *valueString;
  size_t offset;
} OneParam;
