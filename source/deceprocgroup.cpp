/******************************************************************************/
/**     DeCE Proccessing: Group Cross Section                                **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>

using namespace std;

#include "dece.h"
#include "gfr.h"
#include "terminate.h"
#include "groupstructure.h"

static const int Ndiv = 20;
static const int ncx = 15;
static int mfr[ncx] = {3, 3, 3,  3,  3,  3,  3,  3,   3,   3,   3,   3,   3,   3,   1};
static int mtr[ncx] = {1, 2, 4, 16, 17, 18, 22, 28, 102, 103, 104, 105, 106, 107, 452};

static void DeceGroupAverage (ENDF *, const int, const int, double *, double *);

/**********************************************************/
/*      Generate Group Cross Section                      */
/*      --------                                          */
/*      The group structures are:                         */
/*          0: SAND-IIa 640 group                         */
/*          1: LANL 70 group                              */
/*          2: VITAMINE-J 175 group                       */
/*          3: SAND-IIa 725 group                         */
/*      The weight of averaging is                        */
/*          0: constant                                   */
/*          1: 1/E                                        */
/**********************************************************/
void DeceGenerateGroup(ENDFDict *dict, ENDF *lib[], const int group, const int weight)
{
  double *xdat = NULL, **ydat;

  /*** determine group structure */
  int ng = 0;
  switch(group){
  case  0: ng = grpEnergyPoint0;  xdat = grpEnergyGrid0; break;
  case  1: ng = grpEnergyPoint1;  xdat = grpEnergyGrid1; break;
  case  2: ng = grpEnergyPoint2;  xdat = grpEnergyGrid2; break;
  case  3: ng = grpEnergyPoint3;  xdat = grpEnergyGrid3; break;
  default: break;
  }

  if(ng == 0){
    message << "group number " << group << " not defined";
    WarningMessage();
    return;
  }

  ydat = new double * [ncx];
  for(int j=0 ; j<ncx ; j++){
    ydat[j] = new double [ng];
    for(int k=0 ; k<ng ; k++) ydat[j][k] = 0.0;
  }

  /*** check if total is given */
  int id = dict->getID(3,1);
  if(id < 0){
    message << "MF3 MT1 should exist for processing";
    WarningMessage();
    return;
  }

  /*** excluded channels' MT number negative */
  for(int j=1 ; j<ncx ; j++){
    if(dict->getID(mfr[j],mtr[j]) < 0) mtr[j] *= -1;
  }

  /*** group average */
  for(int j=0 ; j<ncx ; j++){
    id = dict->getID(mfr[j],mtr[j]);
    if(id >= 0) DeceGroupAverage(lib[id],weight,ng,xdat,ydat[j]);
  }

  /*** print heading */
  cout << setprecision(6);
  cout <<"# Emin          Emax        ";
  for(int j=0 ; j<ncx ; j++){
    id = dict->getID(mfr[j],mtr[j]);
    if(id >= 0) cout << "   [" << setw(3) << mfr[j] << setw(5) << mtr[j] << "] ";
  }
  cout << endl;

  /*** print average cross section */
  for(int k=0 ; k<ng-1 ; k++){
    cout << setw(14) << xdat[k];
    cout << setw(14) << xdat[k+1];

    for(int j=0 ; j<ncx ; j++){
      id = dict->getID(mfr[j],mtr[j]);
      if(id >= 0) cout << setw(14) << ydat[j][k];
    }
    cout << endl;
  }

  for(int j=0 ; j<ncx ; j++){
    delete [] ydat[j];
  }
  delete [] ydat;
}


/**********************************************************/
/*      Calculate Average Cross Section                   */
/**********************************************************/
/* Linear Interpolation */
//static inline double linear(double x1, double y1, double x2, double y2, double x)
//{
//  double s = (y2-y1)/(x2-x1) * (x-x1) + y1;
//  if(s < 0.0) s = 0.0;
// return s;
//}

/* Weighting Functions */
static inline double wfunc(const int k, double x)
{
  double w = 1.0;
  if(k == 1) w = 1.0/x;
  return w;
}

/* Composite Simpspon's Rule */
static inline double integ_interval(ENDF *lib, const int n, const double e0, const double e1)
{
  double de = (e1 - e0) / Ndiv;
  double f0 = wfunc(n,e0);
  double f1 = wfunc(n,e1);

  double z0 = ENDFInterpolation(lib,e0,true,0);
  double z1 = ENDFInterpolation(lib,e1,false,0);

  double w =      f0 +      f1;
  double s = z0 * f0 + z1 * f1;

  for(int k=1 ; k<Ndiv ; k++){
    double c = (k%2 == 0) ? 2.0 : 4.0;
    double e = e0 + de * k;
    double f = wfunc(n,e);
    double z = ENDFInterpolation(lib,e,true,0);

    w += c *     f;
    s += c * z * f;
  }

  if(w > 0.0) s = s / w * (e1 - e0); // Simpson's h/3 factor cancels
  else s = 0.0;
  
  return s;
}


void DeceGroupAverage(ENDF *lib, const int weight, const int ng, double *energy, double *sigma)
{
  int np  = lib->rdata[0].n2;

  for(int i=0 ; i<ng-1 ; i++){

    sigma[i] = 0.0;

    double x0 = 0.0, x1 = 0.0;
    double y0 = 0.0, y1 = 0.0;

    /*** boundary energies for the i-th energy group */
    double e0 = energy[i];
    double e1 = energy[i+1];

    /*** (x,y) data just above E0 */
    int j0 = 0;
    for(int j=1 ; j<np ; j++){
      int j2 = j*2;
      if(lib->xdata[j2] >= e0){
        x0 = lib->xdata[j2  ];
        y0 = lib->xdata[j2+1];
        j0 = j; break;
      }
    }

    /*** (x,y) data just below E1 */
    int j1 = 0;
    for(int j=j0 ; j<np ; j++){
      int j2 = j*2;
      if(lib->xdata[j2] >= e1){
        x1 = lib->xdata[j2-2];
        y1 = lib->xdata[j2-1];
        j1 = j; break;
      }
    }

    if((j0 == 0) && (j1 == 0)) continue;

    double s = 0.0;

    /*** each interval is sub-devided by Ndiv */
    if(j0 == j1){
      s = integ_interval(lib,weight,e0,e1);
    }
    else{
      /*** from E0 to the first point (x0,y0) */
      if(x0 > e0) s += integ_interval(lib,weight,e0,x0);

      /*** from the last point (x1,y1) to E1 */
      if(e1 > x1) s += integ_interval(lib,weight,x1,e1);

      /*** all the intervals */
      for(int j=j0 ; j<j1-1 ; j++){
        int j2 = j*2;
        if(lib->xdata[j2+2] > lib->xdata[j2]){
          s += integ_interval(lib,weight,lib->xdata[j2],lib->xdata[j2+2]);
        }
      }
    }

    /*** becase S is integral in [E0,E1], need to divide by the interval width */
    sigma[i] = s / (e1 - e0);
  }
}


