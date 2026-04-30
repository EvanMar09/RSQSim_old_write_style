
/* some example fAfuncs - they should be monotonically increasing
   from f(0) = 0 to f(1) = 1 */
#include <math.h>
/* with aReduceN1 = 0, fAmin = 0, and fAfunc() = sqrt(), this should reproduce the
   old aReduceNmax behaviour with aReduceNmax = aReduceN1 */
double fAfunc(double x)
{
  return(sqrt(x));
}
