/* counts the number of lines in the given file 
  leaving out those whose first char is comment[0];  pass in NULL for comment to count all lines */

#include <stdio.h>
#include <string.h>

#define NLINES_LL 128

int nlines(char *fname, char *comment)
{
   char tmp[NLINES_LL];
   int nl = 0;
   FILE *fp;

   if ( (fp = fopen(fname, "r")) == (FILE *)NULL)
   {
      fprintf(stderr, "nlines: error opening %s\n", fname);
      return(0);
   }

   while (fgets(tmp, NLINES_LL, fp) != (char *)NULL)
   {
     if (comment != (char *)NULL)
     {
       if (tmp[0] == comment[0]) continue;
     }
     
     if (tmp[strlen(tmp)-1] == '\n') nl++;
   }

   fclose(fp);

   return(nl);
}
