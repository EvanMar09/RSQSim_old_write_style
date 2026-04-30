#include <stdio.h>
#include <stdlib.h>

FILE *xfopen(char *fname, char *mode)
{
   FILE *fp;

   if ( (fp = fopen(fname, mode)) == (FILE *)NULL)
   {
      fprintf(stderr, "error opening %s\n", fname);
      exit(-1);
   }
   else
      return(fp);
}


