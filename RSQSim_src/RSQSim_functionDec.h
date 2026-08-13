void rollBackFilePointers(Params *params, FaultModel *m);
void assignPatches(FaultModel *m);
void setUpMPI(FaultModel *m);
Eqs RSQSim(Params *params, FaultModel *m);
void writeTau(FaultModel *m, double t, int fpNum, int setOffset);
void writeSigma(FaultModel *m, double t, int fpNum, int setOffset);
void writeSlip(FaultModel *m, double t, int fpNum, int setOffset);
void writeSlipSpeed(FaultModel *m, double t, int fpNum, int setOffset);
void writeState(FaultModel *m, double t, int fpNum, int setOffset);
void writeTheta(FaultModel *m, double t, int fpNum, int setOffset);
void writePED(FaultModel *m, Eqs *eqs, int nslipped);

Transition findMinTransitionTime(FaultModel *m);
Transition S0TransTime(FaultModel *m, int patch, double minDt);
Transition S1TransTime(FaultModel *m, int patch);
Transition S2TransTime(FaultModel *m, int patch);
Transition S3TransTime(FaultPatch *p);
Transition S2xTransTime(FaultPatch *p);
Transition checkForLowStress(FaultPatch *p);
Transition checkForHighTheta(FaultPatch *p, double t, double maxThetaPin);

void makeTransition(FaultModel *m, Transition tr);
void aPrioriFail(FaultModel *m, int patch);
void transition0to1(FaultModel *m, int patch);
void transition1to0(FaultModel *m, int patch);
void transition1to2(FaultModel *m, int patch);
void transition1to2a(FaultModel *m, int patch);
void transition2to0(FaultModel *m, int patch);
void transition2ato0(FaultModel *m, int patch);
void transition2ato2b(FaultModel *m, int patch);
void transition2bto2a(FaultModel *m, int patch);
void transition2bto2b(FaultModel *m, int patch);
void transition2bto2c(FaultModel *m, int patch);
void transition2cto2b(FaultModel *m, int patch);
void transition2cto0(FaultModel *m, int patch);
void transition3to3(FaultModel *m, int patch);
void transitionLowSigma(FaultModel *m, int patch, int nextState);
void transitionMaxTimeStep(FaultModel *m);
void transitionStressRateStep(FaultModel *m);
void checkA(FaultModel *m, int patch);
void reduceA(FaultModel *m, int patch);
void unreduceA(FaultModel *m, int patch);

void timeStep(FaultModel *m, double dt);
void updateState0(FaultModel *m, int patch, double dt);
void updateState1(FaultModel *m, int patch, double dt);
void updateState2(FaultModel *m, int patch, double dt);
void updateState3(FaultModel *m, int patch, double dt);
void updateState2x(FaultModel *m, int patch, double dt);

void initFaultModel(FaultModel *m, Params *params);
void openFiles(FaultModel *m, Params *params, FILE **Fps);
void initializeFaultPatch(FaultModel *m, Params *params, FILE **Fps, int i, double *stressRate);
void closeFiles(FILE **Fps, int nfps);
void readFault(FaultModel *m, Params *params);
void allocFaultModel(FaultModel *m);
void freeFaultModel(FaultModel *m);
void allocFaultPatch(FaultModel *m, FaultPatch *f, int iPatch);
void initEqs(Eqs *eqs, FaultModel *m, Params *params);
void allocSnapBuffers(FaultModel *m);

void packParams(Params *params, void **parambuffer, int *bufsize);
void unpackParams(Params *params, void *parambuffer, int bufsize);

void packFaultModel(FaultModel *m, void **mbuffer, int *bufsize);
void unpackFaultModel(FaultModel *m, void *mbuffer, int bufsize);
void packFaultPatch(FaultPatch *p, void *buffer, int *position, int bufsize);
void unpackFaultPatch(FaultPatch *p, void *buffer, int *position, int bufsize);
void packKRDOkadaFaultPatch(KRDOkadaFaultPatch *op, void *buffer, int *position, int bufsize);
void unpackKRDOkadaFaultPatch(KRDOkadaFaultPatch *op, void *buffer, int *position, int bufsize);


//void updateEqs(Eqs *eqs, FaultModel *m, double t, Transition tr, Params *params);
void finishUp(FaultModel *m, Params *params);

double fixedPoint(double (*f)(double x), double x0, double acc, double prec,
                  double xmin, double xmax);

void getParams(int argc, char **argv, Params *params, FaultModel *m);
void checkDdot(FaultPatch *p, int ip);
void xferParamsToM(FaultModel *m, Params *params);

double theta_0(double t, double sigmaDot, double sigma0, double alpha,
               double B, double theta0);
double theta_1(double tau, double sigma, double ddot, double mu0, double ddotStar,
               double A, double B, double Dc);
double theta_1a(double theta0, double d, double Dc, double alpha,
                double B, double sigma0, double sigma);
double d1(double t, double A, double H, double ddot0, double sigma0,
          double SDotExt);
double ddot1(double t, double ddot0, double SDotExt, double A,
             double sigma0, double H, double ddotEQ);
double ddotRS(double ddotStar, double tau, double sigma, double mu0,
              double a, double b, double theta, double Dc);
void setDdot3(FaultPatch *p);
void setDdot2b(FaultPatch *p);

double zbrent(double (*func)(double,FaultPatch *), double x1, double x2, double tol,FaultPatch *p);

int nlines(char *fname, char *comment);
int strflds(char *str);

double fAfunc(double);


void updateEqs(Eqs *eqs, FaultModel *m, double t, Transition tr, Params *params,
        eqs_info *buffer,int *icount,
        MPI_Offset *disp,
        int maxbuffer0, int maxbuffer1, int maxbuffer2,
        int* icount2, int *ibuff, int* nrec, int* tpmap, int *blocksize, int *nelms,
        double *buffer_d,
        double *buffer_tau,
        double *buffer_sigma,
        double *buffer_mu,
        int *buffer_e,
        int  *buffer_p,
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
        double* array_d,char* array_i

);

void write_trans(FaultModel* m, Params* params,double t, int nelements_trans,int nrec_trans, int* index,
	 double* array_d,char* array_i,MPI_Datatype filetype_d,MPI_Datatype filetype_i, int maxbuffer2,int isOnetime,
	 int *icount_tau_trans0,  double* buffer_tau_trans0, MPI_Offset* tauDisp_trans0,
	 int *icount_tau_trans1,  double* buffer_tau_trans1, MPI_Offset* tauDisp_trans1,
	 int *icount_tau_trans2,  double* buffer_tau_trans2, MPI_Offset* tauDisp_trans2,

	 int *icount_sigma_trans0,  double* buffer_sigma_trans0, MPI_Offset* sigmaDisp_trans0,
	 int *icount_sigma_trans1,  double* buffer_sigma_trans1, MPI_Offset* sigmaDisp_trans1,
	 int *icount_sigma_trans2,  double* buffer_sigma_trans2, MPI_Offset* sigmaDisp_trans2,

	 int *icount_slip_trans0,  double* buffer_slip_trans0, MPI_Offset* slipDisp_trans0,
	 int *icount_slip_trans1,  double* buffer_slip_trans1, MPI_Offset* slipDisp_trans1,
	 int *icount_slip_trans2,  double* buffer_slip_trans2, MPI_Offset* slipDisp_trans2,

 	 int *icount_ddot_trans0,  double* buffer_ddot_trans0, MPI_Offset* ddotDisp_trans0,
	 int *icount_ddot_trans1,  double* buffer_ddot_trans1, MPI_Offset* ddotDisp_trans1,
	 int *icount_ddot_trans2,  double* buffer_ddot_trans2, MPI_Offset* ddotDisp_trans2,

	 int *icount_theta_trans0,  double* buffer_theta_trans0, MPI_Offset* thetaDisp_trans0,
	 int *icount_theta_trans1,  double* buffer_theta_trans1, MPI_Offset* thetaDisp_trans1,
	 int *icount_theta_trans2,  double* buffer_theta_trans2, MPI_Offset* thetaDisp_trans2,

	 int *icount_state_trans0,  char* buffer_state_trans0, MPI_Offset* stateDisp_trans0,
	 int *icount_state_trans1,  char* buffer_state_trans1, MPI_Offset* stateDisp_trans1,
	 int *icount_state_trans2,  char* buffer_state_trans2, MPI_Offset* stateDisp_trans2
);
