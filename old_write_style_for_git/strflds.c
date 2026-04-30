/* returns the number of whitespace separated fields in the input string */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xmalloc.h"

int strflds(char *str)
{
  int nf, nc;
  char *strp;
  char *copyOfStr;
  char whitespace[] = {" \t\n"};
  
  copyOfStr = xmalloc((strlen(str)+1)*sizeof(char));
  strcpy(copyOfStr, str);
  nf = 0;
  if (strtok(copyOfStr, whitespace) != NULL)
  {
    nf++;
    while(strtok(NULL, whitespace) != NULL) nf++;
  }
  free(copyOfStr);
  return(nf);
  
  
/*

  Below always seemed to work fine, until I tried on Mac OS X
  
  Despite, OS X supposedly using gcc-3.3, things don't always
  seem to work as the glibc documentation says they should,
  sscanf() below never returns EOF and doesn't return nc 
  properly as 0 when it gets to the end of the string....

  nf = 0;
  strp = str;
  do
  {
    int tmp;
    tmp = sscanf(strp, "%*s %n", &nc);
    if (sscanf(strp, "%*s %n", &nc) == (int)EOF) break;
    if (nc==0) break;
    strp += nc;
    nf++;
  } while (1);

  return(nf);
  */
}
   
      


