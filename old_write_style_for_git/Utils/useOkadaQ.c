#include <stdio.h>
#include <string.h>
#include "../dc3dq.h"

#define MAXLEN 2048

static char *usage[] = 
  {"useOkadaQ infile > outfile",
   " ",
   "where infile's first line gives the geometry of the element:",
   "alpha depth dip W L u1 u2 u3",
   " ",
   "and subsequent lines give the coords of the field points",
   " ",
   "specifying an infile of '-' means read from stdin", 
   " ",
   "it writes to stdout one line per field pt:",
   "ux uy uz uxx uxy uxz uyx uyy uyz uzx uzy uzz",
   ""};

int main(int argc, char **argv)
{
  int i;
  FILE *infp;
  char s[MAXLEN];
  
  __float128 alpha, depth, dip, W, L, u1, u2, u3;
  __float128 x, y, z;
  double alpha_d, depth_d, dip_d, W_d, L_d, u1_d, u2_d, u3_d;
  double x_d, y_d, z_d;

  __float128 ux, uy, uz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz;
  double ux_d, uy_d, uz_d, uxx_d, uxy_d, uxz_d, uyx_d, uyy_d, uyz_d, 
         uzx_d, uzy_d, uzz_d;
  
  if (argc == 1)
  {
    i = 0;
    while (usage[i][0] != '\0')
      fprintf(stderr, "%s\n", usage[i++]);
    
    return(1);
  }
  
  if (strlen(argv[1]) == 1 && argv[1][0] == '-')
  {
    infp = stdin;
  } else
  {
    infp = fopen(argv[1], "r");
  }
  
  fscanf(infp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", 
         &alpha_d, &depth_d, &dip_d, &W_d, &L_d, &u1_d, &u2_d, &u3_d);
  alpha = alpha_d;
  depth = depth_d;
  dip = dip_d;
  W = W_d;
  L = L_d;
  u1 = u1_d;
  u2 = u2_d; 
  u3 = u3_d;
  
  while (fgets(s, MAXLEN, infp) != NULL)
  {
    sscanf(s, "%lf %lf %lf", &x_d, &y_d, &z_d);
    x = x_d; y = y_d; z = z_d;
    dc3dq(alpha, x, y, z, depth, dip, 0.0, L, 0.0, W,
          u1, u2, u3, 
          &ux, &uy, &uz,
          &uxx, &uyx, &uzx,
          &uxy, &uyy, &uzy,
          &uxz, &uyz, &uzz);
    ux_d = ux; uy_d = uy; uz_d = uz;
    uxx_d = uxx; uxy_d = uxy; uxz_d = uxz;
    uyx_d = uyx; uyy_d = uyy; uyz_d = uyz;
    uzx_d = uzx; uzy_d = uzy; uzz_d = uzz;
    fprintf(stdout, "%.20le %.20le %.20le %.20le %.20le %.20le %.20le %.20le %.20le %.20le %.20le %.20le\n",
            ux_d, uy_d, uz_d, uxx_d, uxy_d, uxz_d, uyx_d, uyy_d, uyz_d,
            uzx_d, uzy_d, uzz_d);
  }
  return(0);
}
