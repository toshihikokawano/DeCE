/******************************************************************************/
/**     DeCE Westcott Factor                                                 **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"

const int    subdiv    = 100;
const double BOLTZMANN = 8.6173324e-5; // Boltzman Constant [eV/K]
const double Ethermal  = 0.02526;      // Thermal energy [eV]
const double eps       = 1.0e-99;

static double specaverage (const int, double *, double *, const double, const double);
inline double maxwellian (const double, const double);

static bool firstcall = true;

/**********************************************************/
/*      Calculate Westcott Factor                         */
/**********************************************************/
void DeceWestcottFactor(ENDFDict *dict, ENDF *lib[], const double temperature)
{
  int k0 = dict->getID(3,102);
  if(k0 < 0){ message << "Capture channel not found"; TerminateCode("DeceWestcottFactor"); }

  double t = temperature;
  /*** when negative, the temeprature is in K */
  if(temperature < 0.0) t *= -BOLTZMANN;

  int     n = lib[k0]->rdata[0].n2;
  double *x = new double [n];
  double *y = new double [n];
  for(int i=0 ; i<n ; i++){
    x[i] = lib[k0]->xdata[2*i  ];
    y[i] = lib[k0]->xdata[2*i+1];
  }

  /*** thermal cross section */
  double sig0 = ENDFInterpolation(lib[k0], Ethermal, false, 0);

  cout << setprecision(6) << setiosflags(ios::scientific);

  if(temperature == 0.0){
    if(firstcall) cout << "# Temperature [K] Westcott-g" << endl;
    double t0 = 100;
    for(int i=0 ; i<5 ; i++){
      t = t0 + i * 100;
      double gw = specaverage(n,x,y,t*BOLTZMANN,sig0);
      cout << setw(14) << t << setw(14) << gw << endl;
    }
  }
  else{
    double gw = specaverage(n,x,y,t,sig0);

    if(temperature < 0.0){
      if(firstcall) cout << "# Temperature [K] Westcott-g" << endl;
      cout << setw(14) << -temperature << setw(14) << gw <<endl;
    }
    else{
      if(firstcall) cout << "# Temperature [K] Westcott-g" << endl;
      cout << setw(14) <<  temperature << setw(14) << gw <<endl;
    }
  }
  firstcall = false;

  delete [] x;
  delete [] y;
}


/**********************************************************/
/*      Maxwellian Average                                */
/**********************************************************/
double specaverage(const int n, double *x, double *y, const double t, const double s0)
{
  double s1 = 0.0, s2 = 0.0, macs = 0.0;
  double c = 1.0 / (s0 * sqrt(Ethermal));

  for(int i=0 ; i<n-1 ; i++){
    double d = x[i+1] - x[i]; // energy width
    if(d == 0.0) continue;

    /*** sub-divide each of energy bins */
    double d1 = d/(double)subdiv;
    double y1 = (y[i+1] - y[i])/(double)subdiv;

    double p1 = 0.0; // sum of cross section x spectrum
    double p2 = 0.0; // sum of spectrum

    for(int j=0 ; j<subdiv ; j++){

      double e1 =  x[i] + d1* j;
      double e2 =  x[i] + d1*(j+1);
      double u1 =  y[i] + y1* j   ;
      double u2 =  y[i] + y1*(j+1);

      /*** Maxwellian spectrum at both bin sides */
      double w1 = maxwellian(e1,t);
      double w2 = maxwellian(e2,t);

      /*** sqrt(E) x sigma */
      u1 *= sqrt(e1);
      u2 *= sqrt(e2);

      p1 += (u1*w1 + u2*w2) * d1 * 0.5;
      p2 += (   w1 +    w2) * d1 * 0.5;

      if(w1 == 0.0 && w2 == 0.0) break;
    }
    s1 += p1;
    s2 += p2;

    if(p1 < eps && p2 < eps) break;
/*
    cout << setprecision(5) << setiosflags(ios::scientific);
    cout << setw(5) << i;
    cout << setw(14) << x[i];
    cout << setw(14) << y[i];
    cout << setw(14) << s0 * sqrt(Ethermal / x[i]);
    cout << setw(14) << maxwellian(x[i],t);
    cout << setw(14) << p1;
    cout << setw(14) << p2;
    if(s2 > 0.0) cout << setw(12) << s1/s2 * c;
    cout << endl;
*/
  }

  if(s2 > 0.0) macs = c * s1 / s2;

  return(macs);
}


inline double maxwellian(const double e, const double t)
{
  double w = e * exp(-e/t);
  if(w < eps) w = 0.0;
  return(w);
}
