#include <stdio.h>
#include <math.h>

void swap(double *a, double *b);
void rotateXyVec(double *x, double *y, double alpha);
void advs(double y1, double y2, double y3, double a, double b, double nu, 
          double B1, double B2, double B3, double *e);
void calcTriStrains(double sx, double sy, double sz, 
                    double x[3], double y[3], double z[3], 
                    double pr, double ss, double ts, double ds, 
                    double S[6]);

static double cot(double x)
{
  return(1./tan(x));
}
