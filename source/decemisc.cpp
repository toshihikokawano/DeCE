/******************************************************************************/
/**     DeCE MISC                                                            **/
/******************************************************************************/

#include <cmath>

using namespace std;

#include "decemisc.h"
#include "constant.h"

/**********************************************************/
/*     Legendre Function                                  */
/**********************************************************/
double legendre(int n, double t)
{
  double x  = cos(PI*t/180.0);
  double p0 = 1.0;
  double p1 = x;
  double p2 = (3*x*x-1)/2;
  double p3 = 0.0;

  double p  = 0.0;
  if(n==0)        p=p0;
  else if(n==1)   p=p1;
  else if(n==2)   p=p2;
  else{
    for(int i=2 ; i<n ; i++){
      double pn = (double)i;
      p3 = ((2*pn+1)*x*p2-pn*p1)/(pn+1);
      p1=p2;  p2=p3;
    }
    p=p3;
  }
  return(p);
}


