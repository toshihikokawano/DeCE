/******************************************************************************/
/**     Kalbach-Mann Systematics for double-differential cross section       **/
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "constant.h"
#include "kalbach.h"

static ZAnumber comp(0,0), proj(0,0);
static double  ma = 0.0, mb[7]; // reduced masses
static double  sa = 0.0, sb[7]; // separation energies
static int     incid = 0;       // incident particle ID
static int     ejcid = 0;       // ejectile particle ID

static double ddxKalbachAfac(const int, const int, const double, const double, const double, const double);

/**********************************************************/
/*      Kalbach Systematics for DDX                       */
/**********************************************************/
void ddxKalbach(const int nddx, const double ein, const double eout, const double f, const double p, double *ang, double *ddx)
{
  double ecmsa = ein  * ma * 1e-6;
  double ecmsb = eout * mb[ejcid] * 1e-6;

  /*** Kalbach systematics a-parameter */
  double a = ddxKalbachAfac(incid,ejcid,ecmsb,sa,sb[ejcid],ecmsa);
  double x = p * a/(2.0*sinh(a));

  /*** calculate angular distribution */
  for(int i=0 ; i<nddx ; i++){
    double q = ang[i]/180.0 * PI;
    ddx[i] = x * (cosh(a*cos(q)) + f*sinh(a*cos(q)));
  }
}


/**********************************************************/
/*      Separation Energies Defined in Kalbach Syst.      */
/**********************************************************/
void ddxKalbachSetParm(const double zap, const double zat, const int idp)
{
  unsigned int znum, anum;
  znum = (unsigned int)(zap / 1000.0);
  anum = (unsigned int)(zap - znum*1000.0);
  proj.setZA(znum,anum);

  znum = (unsigned int)(zat / 1000.0);
  anum = (unsigned int)(zat - znum*1000.0);
  comp.setZA(znum,anum);

  ejcid = idp;

  double   ib = 0.0;
  ZAnumber ejec(0,0),resd(0,0);

  sb[0] = mb[0] = 0.0;

  for(int c=1 ; c<7 ; c++){

    switch(c){
    case  1: ejec.setZA(0,1); ib =  0.0  ;  break;
    case  2: ejec.setZA(1,1); ib =  0.0  ;  break;
    case  3: ejec.setZA(2,4); ib = 28.296;  break;
    case  4: ejec.setZA(1,2); ib =  2.225;  break;
    case  5: ejec.setZA(1,3); ib =  7.718;  break;
    case  6: ejec.setZA(2,3); ib =  8.482;  break;
    default: ejec.setZA(0,0); ib =  0.0  ;  break;
    }

    resd = comp - ejec;
    if( (proj.getZ() == ejec.getZ()) && (proj.getA() == ejec.getA()) ) incid = c;

    double xc = (double)(comp.getN()-comp.getZ())*(double)(comp.getN()-comp.getZ());
    double xb = (double)(resd.getN()-resd.getZ())*(double)(resd.getN()-resd.getZ());
    double yc = pow((double)comp.getA(),1.0/3.0);
    double yb = pow((double)resd.getA(),1.0/3.0);
    double zc = (double)comp.getZ()*(double)comp.getZ();
    double zb = (double)resd.getZ()*(double)resd.getZ();

    sb[c] = 15.68 * (comp.getA()    - resd.getA()   )
          - 28.07 * (xc/comp.getA() - xb/resd.getA())
          - 18.56 * (yc*yc - yb*yb)
          + 33.22 * (xc/pow(yc,4.0) - xb/pow(yb,4.0))
          - 0.717 * (zc/yc - zb/yb)
          + 1.211 * (zc/comp.getA() - zb/resd.getA())
          - ib;

    mb[c] = (double)resd.getA() / ((double)(resd.getA() + ejec.getA()));
  }

  sa = sb[incid];
  ma = mb[incid];
}


/**********************************************************/
/*      a(E) Factor in Kalbach Systematics                */
/**********************************************************/
double ddxKalbachAfac(const int a, const int b, const double eout, const double sa, const double sb, const double ecms)
{
  const double et1 = 130.0;
  const double et3 =  41.0;
  const double c1 = 0.04, c2 = 1.8e-06, c3 = 6.7e-07;

  double ea = ecms + sa;
  double eb = eout + sb;

  double e1 = min(ea,et1);
  double e3 = min(ea,et3);

  double x1 = e1*eb/ea;
  double x3 = e3*eb/ea;

  double ma = (a==3) ? 0.0 : 1.0;
  double mb = (b==3) ? 2.0 : ((b==1) ? 0.5 : 1.0);

  return( c1*x1 + c2*pow(x1,3.0) + c3*pow(x3,4.0)*ma*mb );
}
