

#include "RSQSim.h"

static void getDefaults(Params *params);
static void setOneParam(Params *params, int iP, char *valueString);
static int findParam(char *name);
static void getOneParamString(Params *params, int iP, char *valueString);

static char svnId[] = {"$Id: getParams.c.in 717 2014-03-04 06:39:54Z dinger $"};
static char svnVersion[] = {"@svnversion: @"};

void getParams(int argc, char **argv, Params *params, FaultModel *m)
{
  FILE *pfp;
  char line[RSQSIM_MAXLINE], name[RSQSIM_MAXLINE], value[RSQSIM_MAXLINE]; 
  int iarg, iP, nParams;
  
  char outfname[RSQSIM_MAXLINE];
  
  getDefaults(params);
 
  for (iarg=1; iarg<argc; iarg++)  /* read input control files */
  {
    if (strcmp(argv[iarg], "-") == 0) 
    {
      nParams = sizeof(defaultParams)/sizeof(defaultParams[0]);
     if (m->master) fprintf(stderr, "Current parameters:\n");
      for (iP=0; iP<nParams; iP++)
      {
	getOneParamString(params, iP, value);
	if (m->master) fprintf(stderr, "  %s = %s\n", defaultParams[iP].name, value); /* NB. defaultControls used here just to get name */
      }
      if (m->master) fprintf(stderr, "Enter params (name value) one per line and cntl-D when done\n");
      fflush(stderr);
      pfp = stdin;
    }
    else pfp = xfopen(argv[iarg], "r");
    
    while(fgets(line, RSQSIM_MAXLINE, pfp) != NULL)
    {
      if (line[0] == '#') continue;   /* comment line */
      
      if (sscanf(line, "%s", name) <= 0) /* get name of control paramter */
      {
        if (m->master) fprintf(stderr, "getParams: warning: unable to parse line: %s\n", line);
        continue;  
      }
      
      iP = findParam(name);   /* find position in list */
      if (iP == -1) /* not found */
      {
        if (m->master) fprintf(stderr, "getParams: warning: unknown parameter name %s (ignoring)\n", name);
	continue;
      }
      else if (iP == -2)  /* ambiguous */
      {
        if (m->master) fprintf(stderr, "getParams: warning: ambiguous parameter name %s (ignoring)\n", name);
	continue;
      }
	
      value[0] = '\0';  /* clear any value left over from previous lines */
      if (sscanf(MAX(line+strlen(name), strchr(line, '=')+1), "%s", value) <= 0)
      {
        if (m->master) fprintf(stderr, "getParams: warning: unable to parse line: %s\n", line);
        continue;
      }
      setOneParam(params, iP, value);
    }
    
    if (pfp != stdin) fclose(pfp);
  }
  
  /* write out control parameters used to output file and stderr */
  sprintf(outfname, "eqs.%s.out", params->outFnameInfix);

  /* but first check if eqs*out file already exists */
  if (m->master && fopen(outfname, "r") != (FILE *)NULL)
  {
    if (m->master) fprintf(stderr, "Error: output eqs file (%s) already exists\n", outfname);
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
  }

  if (m->master) m->eqsfp = xfopen(outfname, "w");
  nParams = sizeof(defaultParams)/sizeof(defaultParams[0]);
  if (m->master) fprintf(m->eqsfp, "ProgramVersion =  %s\n", svnVersion+13);
  if (m->master) fprintf(stderr, "ProgramVersion =  %s\n", svnVersion+13);
  if (m->master) fprintf(m->eqsfp, "Control parameters used:\n");
  if (m->master) fprintf(stderr, "Control parameters used:\n");
  for (iP=0; iP<nParams; iP++)
  {
    getOneParamString(params, iP, value);
    if (m->master) fprintf(m->eqsfp, "  %s = %s\n", defaultParams[iP].name, value); /* NB. defaultControls used here just to get name */
    if (m->master) fprintf(stderr, "  %s = %s\n", defaultParams[iP].name, value); /* NB. defaultControls used here just to get name */
  }
  if (m->master) fflush(m->eqsfp);
  if (m->master) fflush(stderr);
  
  return;
}
      
static void getDefaults(Params *params)
{
  int nParams, iP;
  
  nParams = sizeof(defaultParams)/sizeof(defaultParams[0]);
  
  /* set up defaults */
  for (iP=0; iP<nParams; iP++)  
    setOneParam(params, iP, defaultParams[iP].valueString);
  
  return;
}

static void setOneParam(Params *params, int iP, char *valueString)
{

  int vi;
  long int vli;
  float vf;
  double vlf;
  char *vs;


  if (strcmp(defaultParams[iP].type, "d") == 0) /* int */
  {
    sscanf(valueString, "%d", &vi);
    *(int *)((char *)params+defaultParams[iP].offset) = vi;
  }  
  else if (strcmp(defaultParams[iP].type, "ld") == 0) /* long int */
  {
    sscanf(valueString, "%ld", &vli);
    *(long int *)((char *)params+defaultParams[iP].offset) = vli;
  }
  else if (strcmp(defaultParams[iP].type, "f") == 0) /* float */
  {
    sscanf(valueString, "%f", &vf);
    *(float *)((char *)params+defaultParams[iP].offset) = vf;
  }
  else if (strcmp(defaultParams[iP].type, "lf") == 0) /* double */
  {
    sscanf(valueString, "%lf", &vlf);
    *(double *)((char *)params+defaultParams[iP].offset) = vlf;
  }
  else if (strcmp(defaultParams[iP].type, "c") == 0) /* char */
  {
    *(char *)((char *)params+defaultParams[iP].offset) = valueString[0];
  }
  else if (strcmp(defaultParams[iP].type, "s") == 0) /* string MEMORY LEAK*/
  {
    vs = xmalloc((strlen(valueString)+1)*sizeof(char));
    strcpy(vs, valueString);
    *(char **)((char *)params+defaultParams[iP].offset) = vs;
  }
  else
  {
    fprintf(stderr, "getParams.c: unknown type %s\n", defaultParams[iP].type);
    exit(-1);
  }
  
  return;
}


static void getOneParamString(Params *params, int iP, char *valueString)
{

  int vi;
  long int vli;
  float vf;
  double vlf;
  char *vs;


  if (strcmp(defaultParams[iP].type, "d") == 0) /* int */
  {
    memcpy(&vi, (char *)params+defaultParams[iP].offset, sizeof(int));
    sprintf(valueString, "%d", vi);
  }  
  else if (strcmp(defaultParams[iP].type, "ld") == 0) /* long int */
  {
    memcpy(&vli, (char *)params+defaultParams[iP].offset, sizeof(long int));
    sprintf(valueString, "%ld", vli);
  }
  else if (strcmp(defaultParams[iP].type, "f") == 0) /* float */
  {
    memcpy(&vf, (char *)params+defaultParams[iP].offset, sizeof(float));
    sprintf(valueString, "%.10g", vf);
  }
  else if (strcmp(defaultParams[iP].type, "lf") == 0) /* double */
  {
    memcpy(&vlf, (char *)params+defaultParams[iP].offset, sizeof(double));
    sprintf(valueString, "%.20lg", vlf);
  }
  else if (strcmp(defaultParams[iP].type, "c") == 0) /* char (NOT TESTED) */
  {
    sprintf(valueString, "%1c", *((char *)params+defaultParams[iP].offset));
  }
  else if (strcmp(defaultParams[iP].type, "s") == 0) /* string */
  {
    sprintf(valueString, "%s", *(char **)((char *)params+defaultParams[iP].offset));
  }
  else
  {
    fprintf(stderr, "getControls.c: unknown type %s\n", defaultParams[iP].type);
    exit(-1);
  }
  
  return;
}

static int findParam(char *name) /* returns -1 if no match, -2 if multiple matches, otherwise */
{                                  /* returns index of matching parameter */
  int nParams, iP, fP, nMatch;
  size_t namelen;
  
  nParams = sizeof(defaultParams)/sizeof(defaultParams[0]);
  
  namelen = strlen(name);
  
  nMatch = 0;
  for (iP=0; iP<nParams; iP++)
  {
    if (strcmp(name, defaultParams[iP].name) == 0)
    {
      nMatch++;
      fP = iP;
    }
  }
  
  if (nMatch == 0) return(-1);
  if (nMatch > 1) return(-2);
  return(fP);
}
  
  
