/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Calculate MACS from Pointwise ENDF Data                 **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

#include "../source/endflib.h"

const int    subdiv = 100;
const double PI =  3.14159265358979323846;

int main(int, char *[]);
static double maxwell_average(int, double *, double *, double, double);
inline double maxwell        (double, double);


int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: decemacs temperature(keV) PENDF_file" << endl;  exit(-1);
  }

  string   libname = "";
  ifstream fpin;
  ENDF     lib(L);
  double   *x,*y,t = 0.0, awr = 0.0;
  int      n = 0;

  t = atof(argv[1]) * 1000.0;
  libname = argv[2];

  if(t == 0.0){
    cerr << "zero temperature" << endl;
    exit(-1);
  }

  /*** read in MF3 MT102 */
  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;
    exit(-1);
  }
  ENDFReadMF3(&fpin,&lib,102);
  fpin.close();

  Record head = lib.getENDFhead();
  awr = head.c2;

  n = lib.rdata[0].n2;
  x = new double [n];
  y = new double [n];

  for(int i=0 ; i<n ; i++){
    x[i] = lib.xdata[2*i  ];
    y[i] = lib.xdata[2*i+1];
  }

  double macs = maxwell_average(n,x,y,t,awr);
  cout << setprecision(5) << setiosflags(ios::scientific);
  cout << setw(12) << macs * 1000.0 << endl;

  delete [] x;
  delete [] y;
  return(0);
}


double maxwell_average(int n, double *x, double *y, double t, double awr)
{
  double s1 = 0.0, s2 = 0.0, macs = 0.0;

  double c = awr/(awr+1.0);

  for(int i=0 ; i<n-1 ; i++){
    double d = c*(x[i+1] - x[i]);
    if(d == 0.0) continue;

    double d1 = d/(double)subdiv;
    double y1 = (y[i+1] - y[i])/(double)subdiv;

    double p1 = 0.0, p2 = 0.0;
    for(int j=0 ; j<subdiv ; j++){
      double e1 = (x[i] + d1* j   )*c;
      double e2 = (x[i] + d1*(j+1))*c;
      double u1 =  y[i] + y1* j   ;
      double u2 =  y[i] + y1*(j+1);

      double w1 = maxwell(e1,t);
      double w2 = maxwell(e2,t);

      p1 += (u1*w1 + u2*w2)*d1/2.0;
      p2 += (   w1 +    w2)*d1/2.0;

    }
    s1 += p1;
    s2 += p2;
/*
    cout << setw(5) << i;
    cout << setw(12) << x[i]*c;
    cout << setw(12) << y[i];
    cout << setw(12) << maxwell(x[i]*c,t) << endl;
 */
  }

  macs = s1/s2 * 2.0/sqrt(PI);

  return(macs);
}


inline double maxwell(double e, double t)
{
  double w = e * exp(-e/t);
  if(w < 1e-99) w = 0.0;
  return(w);
}
