#include "RSQSim.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif

#include <assert.h>



	

void write_out(MPI_File      fh,
               MPI_Offset    disp,
               int           num,
               MPI_Datatype  buftype,
               MPI_Datatype  filetype,
               void         *buffer)
{
    MPI_Status status;

    MPI_File_set_view(fh, disp, buftype, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, buffer, num, buftype, &status);
}

void writePED_xx(FaultModel *m, Eqs *eqs, int nslipped, MPI_Comm comm,
        MPI_Offset *disp, int maxbuffer, int *icount, int *ibuff, int *nrec,
        int *tpmap, int *blocksize, int *nelms, int isOnetime,
        double *buffer_d,
        double *buffer_tau,
        double *buffer_sigma,
        double *buffer_mu,
        int    *buffer_e,
        int    *buffer_p,
        double *buffer_t,
        double *buffer_taup)
{
    int i, nleft, ntotal;

    if (isOnetime == FALSE)
    {
        memcpy(buffer_d    +(*ibuff), eqs->dd,     nslipped*sizeof(double));
        memcpy(buffer_tau  +(*ibuff), eqs->dtau,   nslipped*sizeof(double));
        memcpy(buffer_sigma+(*ibuff), eqs->dsigma, nslipped*sizeof(double));
        memcpy(buffer_mu   +(*ibuff), eqs->dmu,    nslipped*sizeof(double));

	for (i=0; i<nslipped; i++)
        {
            buffer_t   [*ibuff] = eqs->tFirstRupture[eqs->pNumber[i]];
            buffer_taup[*ibuff] = eqs->tauFirstRupture[eqs->pNumber[i]];
            buffer_e   [*ibuff] = eqs->ieq + 1;
            buffer_p   [*ibuff] = eqs->pNumber[i] + 1;

            if (buffer_p[*ibuff]<=0) fprintf(stdout," Weird ---\n");

            (*ibuff)++;
	}

        MPI_Scan(&nslipped, &nleft, 1, MPI_INT, MPI_SUM, comm);
        nleft -= nslipped;

        MPI_Allreduce(&nslipped, &ntotal, 1, MPI_INT, MPI_SUM, comm);

        if (nslipped > 0) /* skip zero-length write */
        { 
            tpmap[*nelms] = *nrec + nleft;
            blocksize[*nelms] = nslipped;
            (*nelms)++;
        }
        *nrec += ntotal;
        (*icount)++;
    }
        /* write out */
        if (*icount == maxbuffer) {
            MPI_Offset offset_d, offset_i;
            MPI_Datatype filetype_d, filetype_i;

            offset_d = *disp * sizeof(double);
            offset_i = *disp * sizeof(int);

	if (*nelms>0){
            MPI_Type_indexed(*nelms, blocksize, tpmap, MPI_DOUBLE, &filetype_d);
            MPI_Type_commit(&filetype_d);
            MPI_Type_indexed(*nelms, blocksize, tpmap, MPI_INT,    &filetype_i);
            MPI_Type_commit(&filetype_i);
	}else{
		filetype_d=MPI_DOUBLE;
		filetype_i=MPI_INT;
	}
            write_out(m->dListfp,      offset_d, *ibuff, MPI_DOUBLE, filetype_d, buffer_d);
            write_out(m->dtauListfp,   offset_d, *ibuff, MPI_DOUBLE, filetype_d, buffer_tau);
            write_out(m->dsigmaListfp, offset_d, *ibuff, MPI_DOUBLE, filetype_d, buffer_sigma);
            write_out(m->dmuListfp,    offset_d, *ibuff, MPI_DOUBLE, filetype_d, buffer_mu);
            write_out(m->tListfp,      offset_d, *ibuff, MPI_DOUBLE, filetype_d, buffer_t);
            write_out(m->taupListfp,   offset_d, *ibuff, MPI_DOUBLE, filetype_d, buffer_taup);
            write_out(m->eListfp,      offset_i, *ibuff, MPI_INT,    filetype_i, buffer_e);
            write_out(m->pListfp,      offset_i, *ibuff, MPI_INT,    filetype_i, buffer_p);

	if (*nelms>0){
            MPI_Type_free(&filetype_d);
            MPI_Type_free(&filetype_i);
	}
            *disp += *nrec;

            /* reset counters */
            *nrec=0;
            *ibuff=0;
            *icount=0;
            *nelms=0;
        }
    
}



/*
	set up *nrec, *nelements,filetype

*/
void setup(FaultModel *m,int maxbuffer, MPI_Datatype datatype,
		int* nelements, MPI_Datatype *filetype,int* nrec,int *index){
	int nleft=0;
	int nbyte,i,j,ip;
	*nrec=0;
	*nelements=0;

	MPI_Type_size(datatype,&nbyte);
	
        int* pmap= (int*)malloc(sizeof(int)*(m->stopPatch-m->startPatch+2));
  	


	if (m->master){
                        pmap[(*nelements)]= 0;
                        (*nelements)++;
        }else{

		nleft++;
	}


	if (m->npatchWrite == m->N){

		for (ip = m->startPatch; ip <= m->stopPatch; ip++){
			pmap[(*nelements)] =ip+1;  //add one 
			index[(*nelements)]=ip;

			(*nelements)++;
		}

		*nrec = m->N+1;
 
	}else{
	/* find total number */

		 (*nrec)++;
		
		for (ip = 0; ip < m->N; ip++)
        		if (m->p[ip].writePatch !=0) (*nrec)++;

		for (ip=0; ip < m->startPatch;ip++)
			if (m->p[ip].writePatch !=0) 
					nleft++;

		for (ip = m->startPatch; ip <= m->stopPatch; ip++)
        		if (m->p[ip].writePatch !=0) {
				 pmap[(*nelements)]=nleft + (*nelements);
				index[(*nelements)]=ip;

				(*nelements)++;
			}
	}

/*
        MPI_Aint* tpmap= (MPI_Aint*)malloc(sizeof(MPI_Aint)*maxbuffer*(*nelements));

        int *blocksize;blocksize=(int *)malloc((maxbuffer*(*nelements))*sizeof(int));



	for (j=0;j<maxbuffer;j++)
	for (i=0;i<*nelements;i++){
		blocksize[i + j* (*nelements)] = 1;
		tpmap[i + j* (*nelements)] = nbyte* (pmap[i] + j* (*nrec));
	}
 	

	MPI_Type_create_hindexed(maxbuffer*(*nelements),blocksize,tpmap,datatype,filetype);
*/

       int* tpmap= (int*)malloc(sizeof(int)*maxbuffer);

        int *blocksize;blocksize=(int *)malloc(maxbuffer*sizeof(int));

	int ncount= (*nelements >0)? maxbuffer:0;	
	for (j=0;j<ncount;j++){
		blocksize[j]= *nelements;
		tpmap[j]= j*(*nrec)+ pmap[0];
	}


	if (ncount >0){
        	MPI_Type_indexed(ncount,blocksize,tpmap,datatype,filetype);
        	MPI_Type_commit(filetype);
	}else{
		*filetype=datatype;
	}

	free(tpmap);
	free(pmap);
	free(blocksize);
        
	//fprintf(stdout," free blocksize----commit %zd   %zd \n",sizeof(MPI_Aint)*maxbuffer*(*nelements),(maxbuffer*(*nelements))*sizeof(int));
        //sleep(2);

        return;
}


void write1D_static(double t, MPI_File fh, void* array,int nelements,int nrec, 
	int maxbuffer,int* icount, void* buffer,
	MPI_Offset *disp, MPI_Datatype datatype,MPI_Datatype filetype,int isOnetime){
 	
	int nleft,nbyte,ibuff;
        MPI_Status status;
	MPI_Offset temp = nrec;


	MPI_Type_size(datatype,&nbyte);
		ibuff= (*icount )*nelements;

	//	fprintf(stdout," static : %d %d %d nelem=%d nbyte=%d\n",ibuff, *icount,nrec,nelements,nbyte);

	if (isOnetime==FALSE){
	
        memcpy(buffer + ibuff*nbyte,array,nelements*nbyte);
	ibuff += nelements;
	(*icount )++;
	}else{
	// write out once
	}
	
	if (*icount ==maxbuffer){
		
		MPI_File_set_view(fh, *disp, datatype, filetype, "native", MPI_INFO_NULL);
		MPI_File_write_all(fh, buffer, ibuff, datatype, &status);

		MPI_File_sync(fh);
        	*disp = *disp + temp*maxbuffer * nbyte;  /* could be a large number */

//		fprintf(stdout,"disp = %ld ibuff= %d max=%d nelems=%d temp=%ld nbyte=%d\n",*disp,ibuff,maxbuffer,nelements, temp,nbyte);
		*icount=0;	
	}
}


void copy_double(FaultModel *m,double t, int nelements_trans, int * index, double* array_d,int which_var){
/* which_var =1,2,3,4,5 --> Tau,Sigma,Slip,SLipSpeed,Theta   */

	int ii,jj,ii0;
	
	if (m->master){
                array_d[0]=t;   //include  time
		ii0=1;
	}else{
		ii0=0;
	}

         for (ii=ii0;ii<nelements_trans;ii++){
                jj=index[ii];
		
                if (which_var==1) array_d[ii]=m->p[jj].tau;
                if (which_var==2) array_d[ii]=m->p[jj].sigma;
                if (which_var==3) array_d[ii]=m->p[jj].d;
                if (which_var==4) array_d[ii]=m->p[jj].ddot;
                if (which_var==5) array_d[ii]=m->p[jj].theta;

        }

}

void copy_char(FaultModel *m,double t, int nelements_trans, int * index, char* array_d,int which_var){
/* which_var =1,2,3,4,5 --> Tau,Sigma,Slip,SLipSpeed,Theta   */

	int ii,jj,ii0;
	
	if (m->master){
                array_d[0]=9;   //include  time
		ii0=1;
	}else{
		ii0=0;
	}

         for (ii=ii0;ii<nelements_trans;ii++){
                jj=index[ii];
		
                array_d[ii]=m->p[jj].state;
               
        }

}

void rollback_x(MPI_Offset *Disp, MPI_Offset Disp_save, int *icount, int icount_save,int nbyte,int maxbuffer,int nrec){

//	fprintf(stdout," move pointer  %ld \n ",(*Disp- Disp_save) - (nbyte*maxbuffer*nrec));
//	fflush(stdout);

	MPI_Offset temp=nrec;


    if ((*Disp- Disp_save) >=  (nbyte*maxbuffer*temp)){
                /* already on disk */
                        *Disp = Disp_save + icount_save*temp*nbyte;
                        *icount=0;

     }else{
     	*icount=icount_save;
     }
}



/* restart */




/* wrappers*/

size_t myfread(void *ptr,size_t s,size_t x, FILE *fp){
	return(fread(ptr,s,x,fp));
}
size_t myfwrite(void *ptr,size_t s,size_t x, FILE *fp){
        return(fwrite(ptr,s,x,fp));
}

void    readwrite_dumps(Eqs *eqs, FaultModel *m,double *t,int *currentStep,char* buffer,int rw)
{
/* rw=0: read;rw=1:write */

/*  buffer : directory name     */
        int i;
        FILE *wt;
        char wtfile[50];
	char *perm[]={"rb","wb"};

        size_t (*FUNC)(void *,size_t,size_t, FILE *);
	FUNC=&myfread;
	if (rw==1)FUNC=&myfwrite;

        sprintf(wtfile,"%s/P%-d",buffer,m->iProc);
	//fprintf(stdout,"open file:--%s---[%s] \n",wtfile,perm[rw]);

        if((wt=fopen(wtfile,perm[rw]))!=(FILE*)NULL){

                (*FUNC)(t,sizeof(double),1,wt);
                (*FUNC)(currentStep,sizeof(int),1,wt);
                (*FUNC)(&(eqs->ieq),sizeof(int),1,wt);
                (*FUNC)(&(m->npatchWrite),sizeof(int),1,wt);
                (*FUNC)(&(m->nTransitions),sizeof(int),1,wt);
                (*FUNC)(&(m->state3patches),sizeof(int),1,wt);
	        (*FUNC)(&(m->tLastWrite),sizeof(double),1,wt);

		//if (m->master)  fprintf(stdout,"time = %lf, nTransitions= %d \n",*t,m->nTransitions);

                for (i=m->startPatch;i<=m->stopPatch;i++)
                {
			(*FUNC)(&(m->p[i].tau     ),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].sigma   ),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].theta   ),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].d       ),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].ddot    ),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].tauDot  ),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].sigmaDot),sizeof(double),1,wt);
			
			(*FUNC)(&(m->p[i].muTrans2),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].muTrans2Low),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].muTrans2High),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].muTrans2aLow),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].muTrans3Low),sizeof(double),1,wt);
			(*FUNC)(&(m->p[i].muTrans3High),sizeof(double),1,wt);
			
			(*FUNC)(&(m->p[i].APFail),sizeof(int),1,wt);
			(*FUNC)(&(m->p[i].pinned),sizeof(int),1,wt);
			(*FUNC)(&(m->p[i].neighborInRuptureState),sizeof(int),1,wt);
			(*FUNC)(&(m->p[i].aReduced),sizeof(int),1,wt);
			(*FUNC)(&(m->p[i].slowSlip),sizeof(int),1,wt);
			(*FUNC)(&(m->p[i].writePatch),sizeof(int),1,wt);
			
		}
                for (i=0;i<m->N;i++)
		{
			(*FUNC)(&(m->p[i].state),sizeof(char),1,wt);
			//fprintf(stdout,"%lf %lf\n",p.tau,p.sigma);
                }
        	fclose(wt);
	}
}


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
)
{
	int nbuffer;

  

//tau
	copy_double(m,t, nelements_trans,  index,  array_d,1);

//        fprintf(stdout,"Pass step 0 %d %d %d-- \n",icount_tau_trans0,icount_tau_trans1,icount_tau_trans2);

if ((params->writeTau & TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_tau_trans0 !=0 ){
	nbuffer=*icount_tau_trans0;
  
  	setup(m, *icount_tau_trans0,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	
	if (!(isOnetime == TRUE & *icount_tau_trans0 ==0) )
 	write1D_static(t, (m->taufp[0]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_tau_trans0,  buffer_tau_trans0,
                tauDisp_trans0, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeTau & EVENT) > 0){
  if (isOnetime == TRUE & *icount_tau_trans1 !=0 ){
	nbuffer=*icount_tau_trans1;
  
  	setup(m, *icount_tau_trans1,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_tau_trans1 ==0) )
        write1D_static(t, (m->taufp[1]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_tau_trans1,  buffer_tau_trans1,
                tauDisp_trans1, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeTau & EVENT_TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_tau_trans2 !=0 ){
	nbuffer=*icount_tau_trans2;
  
  	setup(m, *icount_tau_trans2,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
        if (!(isOnetime == TRUE & *icount_tau_trans2 ==0) )
	write1D_static(t, (m->taufp[2]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_tau_trans2,  buffer_tau_trans2,
                tauDisp_trans2, MPI_DOUBLE,filetype_d,isOnetime);
}	
	
//	fprintf(stdout,"Pass step 1 %d %d %d-- \n",icount_tau_trans0,icount_tau_trans1,icount_tau_trans2);
//sigma		
	copy_double(m,t, nelements_trans,  index,  array_d,2);

if ((params->writeSigma & TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_sigma_trans0 !=0 ){
	nbuffer=*icount_sigma_trans0;
  
  	setup(m, *icount_sigma_trans0,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_sigma_trans0 ==0) )
 	write1D_static(t, (m->sigmafp[0]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_sigma_trans0,  buffer_sigma_trans0,
                sigmaDisp_trans0, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeSigma & EVENT) > 0){
  if (isOnetime == TRUE & *icount_sigma_trans1 !=0 ){
	nbuffer=*icount_sigma_trans1;
  
  	setup(m, *icount_sigma_trans1,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_sigma_trans1 ==0) )
       write1D_static(t, (m->sigmafp[1]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_sigma_trans1,  buffer_sigma_trans1,
                sigmaDisp_trans1, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeSigma & EVENT_TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_sigma_trans2 !=0 ){
	nbuffer=*icount_sigma_trans2;
  
  	setup(m, *icount_sigma_trans2,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_sigma_trans2 ==0) )
       write1D_static(t, (m->sigmafp[2]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_sigma_trans2,  buffer_sigma_trans2,
                sigmaDisp_trans2, MPI_DOUBLE,filetype_d,isOnetime);
}
//slip
	copy_double(m,t, nelements_trans,  index,  array_d,3);

if ((params->writeSlip & TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_slip_trans0 !=0 ){
	nbuffer=*icount_slip_trans0;
  
  	setup(m, *icount_slip_trans0,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_slip_trans0 ==0) )
 	write1D_static(t, (m->slipfp[0]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_slip_trans0,  buffer_slip_trans0,
                slipDisp_trans0, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeSlip & EVENT) > 0){
  if (isOnetime == TRUE & *icount_slip_trans1 !=0 ){
	nbuffer=*icount_slip_trans1;
  
  	setup(m, *icount_slip_trans1,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_slip_trans1 ==0) )
       write1D_static(t, (m->slipfp[1]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_slip_trans1,  buffer_slip_trans1,
                slipDisp_trans1, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeSlip & EVENT_TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_slip_trans2 !=0 ){
	nbuffer=*icount_slip_trans2;
  
  	setup(m, *icount_slip_trans2,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_slip_trans2 ==0) )
       write1D_static(t, (m->slipfp[2]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_slip_trans2,  buffer_slip_trans2,
                slipDisp_trans2, MPI_DOUBLE,filetype_d,isOnetime);
}

//ddot
	copy_double(m,t, nelements_trans,  index,  array_d,4);

if ((params->writeSlipSpeed & TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_ddot_trans0 !=0 ){
	nbuffer=*icount_ddot_trans0;
  
  	setup(m, *icount_ddot_trans0,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_ddot_trans0 ==0) )
 	write1D_static(t, (m->slipSpeedfp[0]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_ddot_trans0,  buffer_ddot_trans0,
                ddotDisp_trans0, MPI_DOUBLE,filetype_d,isOnetime);
}
if ((params->writeSlipSpeed & EVENT) > 0){
  if (isOnetime == TRUE & *icount_ddot_trans1 !=0 ){
	nbuffer=*icount_ddot_trans1;
  
  	setup(m, *icount_ddot_trans1,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_ddot_trans1 ==0) )
       write1D_static(t, (m->slipSpeedfp[1]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_ddot_trans1,  buffer_ddot_trans1,
                ddotDisp_trans1, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeSlipSpeed & EVENT_TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_ddot_trans2 !=0 ){
	nbuffer=*icount_ddot_trans2;
  
  	setup(m, *icount_ddot_trans2,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_ddot_trans2 ==0) )
       write1D_static(t, (m->slipSpeedfp[2]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_ddot_trans2,  buffer_ddot_trans2,
                ddotDisp_trans2, MPI_DOUBLE,filetype_d,isOnetime);
}

//theta
	copy_double(m,t, nelements_trans,  index,  array_d,5);

if ((params->writeTheta & TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_theta_trans0 !=0 ){
	nbuffer=*icount_theta_trans0;
  
  	setup(m, *icount_theta_trans0,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_theta_trans0 ==0) )
 	write1D_static(t, (m->thetafp[0]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_theta_trans0,  buffer_theta_trans0,
                thetaDisp_trans0, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeTheta & EVENT) > 0){
  if (isOnetime == TRUE & *icount_theta_trans1 !=0 ){
	nbuffer=*icount_theta_trans1;
  
  	setup(m, *icount_theta_trans1,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_theta_trans1 ==0) )
       write1D_static(t, (m->thetafp[1]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_theta_trans1,  buffer_theta_trans1,
                thetaDisp_trans1, MPI_DOUBLE,filetype_d,isOnetime);
}

if ((params->writeTheta & EVENT_TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_theta_trans2 !=0 ){
	nbuffer=*icount_theta_trans2;
  
  	setup(m, *icount_theta_trans2,  MPI_DOUBLE, &nelements_trans, &filetype_d,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_theta_trans2 ==0) )
       write1D_static(t, (m->thetafp[2]),  array_d, nelements_trans,nrec_trans,
                nbuffer,icount_theta_trans2,  buffer_theta_trans2,
                thetaDisp_trans2, MPI_DOUBLE,filetype_d,isOnetime);
}

//state
	copy_char(m,t, nelements_trans,  index,  array_i,6);

if ((params->writeState & TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_state_trans0 !=0 ){
	nbuffer=*icount_state_trans0;
  
  	setup(m, *icount_state_trans0,  MPI_CHAR, &nelements_trans, &filetype_i,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_state_trans0 ==0) )
 	write1D_static(t, (m->statefp[0]),  array_i, nelements_trans,nrec_trans,
                nbuffer,icount_state_trans0,  buffer_state_trans0,
                stateDisp_trans0, MPI_CHAR,filetype_i,isOnetime);
}
if ((params->writeState & EVENT) > 0){
  if (isOnetime == TRUE & *icount_state_trans1 !=0 ){
	nbuffer=*icount_state_trans1;
  
  	setup(m, *icount_state_trans1,  MPI_CHAR, &nelements_trans, &filetype_i,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_state_trans1 ==0) )
       write1D_static(t, (m->statefp[1]),  array_i, nelements_trans,nrec_trans,
                nbuffer,icount_state_trans1,  buffer_state_trans1,
                stateDisp_trans1, MPI_CHAR,filetype_i,isOnetime);
}

if ((params->writeState & EVENT_TRANSITION) > 0){
  if (isOnetime == TRUE & *icount_state_trans2 !=0 ){
	nbuffer=*icount_state_trans2;
  
  	setup(m, *icount_state_trans2,  MPI_CHAR, &nelements_trans, &filetype_i,&nrec_trans,index);
	}else{
	nbuffer=maxbuffer2;
	}
	if (!(isOnetime == TRUE & *icount_state_trans2 ==0) )
       write1D_static(t, (m->statefp[2]),  array_i, nelements_trans,nrec_trans,
                nbuffer,icount_state_trans2,  buffer_state_trans2,
                stateDisp_trans2, MPI_CHAR,filetype_i,isOnetime);
}
/* end transition*/
} 

