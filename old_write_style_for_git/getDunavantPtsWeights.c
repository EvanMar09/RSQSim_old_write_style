
/* evaluation points and weights from Dunavant [1985] */

#include "RSQSim.h"

int getDunavantPtsWeights(int n, double *w, double *alpha, double *beta, 
                           double *gamma)
{
  double abg[3];
  
  switch(n)
  {
    case 1:
      w[0] = 1.0;
      alpha[0] = 1.0/3.0;
      beta[0] = 1.0/3.0;
      gamma[0] = 1.0/3.0;
      break;
    case 3:
      w[0] = w[1] = w[2] = 1.0/3.0;
      alpha[0] = 2.0/3.0; alpha[1] = 1.0/6.0; alpha[2] = 1.0/6.0;
      beta[0]  = 1.0/6.0; beta[1]  = 2.0/3.0; beta[2]  = 1.0/6.0;
      gamma[0] = 1.0/6.0; gamma[1] = 1.0/6.0; gamma[2] = 2.0/3.0;
      break;
    case 4:
      w[0] = -0.5625; w[1] = w[2] = w[3] = 0.520833333333;
      alpha[0] = 1.0/3.0; alpha[1] = 0.6; alpha[2] = 0.2; alpha[3] = 0.2;
      beta[0]  = 1.0/3.0; beta[1]  = 0.2; beta[2]  = 0.6; beta[3]  = 0.2;
      gamma[0] = 1.0/3.0; gamma[1] = 0.2; gamma[2] = 0.2; gamma[3] = 0.6;
      break;
    case 6:
      w[0] = w[1] = w[2] = 0.223381589678011;
      w[3] = w[4] = w[5] = 0.109951743655322;
      alpha[0] = 0.108103018168070; alpha[1] = 0.445948490915965; alpha[2] = 0.445948490915965;
      beta[0]  = 0.445948490915965; beta[1]  = 0.108103018168070; beta[2]  = 0.445948490915965;
      gamma[0] = 0.445948490915965; gamma[1] = 0.445948490915965; gamma[2] = 0.108103018168070;
      alpha[3] = 0.816847572980459; alpha[4] = 0.091576213509771; alpha[5] = 0.091576213509771;
      beta[3]  = 0.091576213509771; beta[4]  = 0.816847572980459; beta[5]  = 0.091576213509771;
      gamma[3] = 0.091576213509771; gamma[4] = 0.091576213509771; gamma[5] = 0.816847572980459;
      break;
    case 7:
      w[0] = 0.225;
      w[1] = w[2] = w[3] = 0.132394152788506;
      w[4] = w[5] = w[6] = 0.125939180544827;
      alpha[0] = beta[0] = gamma[0] = 1.0/3.0;
      alpha[1] = 0.059715871789770; alpha[2] = 0.470142064105115; alpha[3] = 0.470142064105115;
      beta[1]  = 0.470142064105115; beta[2]  = 0.059715871789770; beta[3]  = 0.470142064105115;
      gamma[1] = 0.470142064105115; gamma[2] = 0.470142064105115; gamma[3] = 0.059715871789770;
      alpha[4] = 0.797426985353087; alpha[5] = 0.101286507323456; alpha[6] = 0.101286507323456;
      beta[4]  = 0.101286507323456; beta[5]  = 0.797426985353087; beta[6]  = 0.101286507323456;
      gamma[4] = 0.101286507323456; gamma[5] = 0.101286507323456; gamma[6] = 0.797426985353087;
      break;
    case 12:
      w[0] = w[1] = w[2] = 0.116786275726379;
      w[3] = w[4] = w[5] = 0.050844906370207; 
      w[6] = w[7] = w[8] = w[9] = w[10] = w[11] = 0.082851075618374;
      abg[0] = 0.501426509658179; abg[1] = 0.249286745170910;
      alpha[0] = abg[0]; alpha[1] = abg[1]; alpha[2] = abg[1];
      beta[0]  = abg[1]; beta[1]  = abg[0]; beta[2]  = abg[1];
      gamma[0] = abg[1]; gamma[1] = abg[1]; gamma[2] = abg[0];
      abg[0] = 0.873821971016996; abg[1] = 0.063089014491502;
      alpha[3] = abg[0]; alpha[4] = abg[1]; alpha[5] = abg[1];
      beta[3]  = abg[1]; beta[4]  = abg[0]; beta[5]  = abg[1];
      gamma[3] = abg[1]; gamma[4] = abg[1]; gamma[5] = abg[0];
      abg[0] = 0.053145049844817; abg[1] = 0.310352451033784; abg[2] = 0.636502499121399;
      alpha[6] = abg[0]; alpha[7] = abg[0]; alpha[8] = abg[1]; alpha[9] = abg[1]; alpha[10] = abg[2]; alpha[11] = abg[2]; 
      beta[6]  = abg[1]; beta[7]  = abg[2]; beta[8]  = abg[0]; beta[9]  = abg[2]; beta[10]  = abg[0]; beta[11]  = abg[1];
      gamma[6] = abg[2]; gamma[7] = abg[1]; gamma[8] = abg[2]; gamma[9] = abg[0]; gamma[10] = abg[1]; gamma[11] = abg[0];
      break;
    default:
      fprintf(stderr, "getDunavantPtsWeights(): error: unsupported number of points\n");
      return(-1);
    }
  
  
  return(0);
}
