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
    if(id >= 0) cout << setw(7) << mfr[j] << setw(7) << mtr[j];
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
static inline double linear(double x1, double y1, double x2, double y2, double x)
{
  double s = (y2-y1)/(x2-x1) * (x-x1) + y1;
  if(s < 0.0) s = 0.0;
 return s;
}

/* Weighting Functions */
static inline double wfunc(const int k, double x)
{
  double w = 1.0;
  if(k == 1) w = 1.0/x;
  return w;
}

/* Composite Simpspon's Rule */
static inline double sinteg(const int n, const double e0, const double z0, const double e1, const double z1)
{
  double de = (e1 - e0) / Ndiv;
  double f0 = wfunc(n,e0);
  double f1 = wfunc(n,e1);
 
  double w =      f0 +      f1;
  double s = z0 * f0 + z1 * f1;

  for(int k=1 ; k<Ndiv ; k++){
    double c = (k%2 == 0) ? 2.0 : 4.0;
    double e = e0 + de * k;
    double f = wfunc(n,e);
    double z = linear(e0,z0,e1,z1,e);

    w += c *     f;
    s += c * z * f;
  }

  if(w > 0.0) s = s / w * (e1 - e0); // Simpson's h/3 factor cancels
  else s = 0.0;
  
  return s;
}


void DeceGroupAverage(ENDF *lib, const int weight, const int ng, double *energy, double *sigma)
{
  int idx = 0;
  int np  = lib->rdata[idx].n2;

  for(int i=0 ; i<ng-1 ; i++){

    sigma[i] = 0.0;

    double x0a = 0.0, x0b = 0.0, x1a = 0.0, x1b = 0.0;
    double y0a = 0.0, y0b = 0.0, y1a = 0.0, y1b = 0.0;
    double z0 = 0.0, z1 = 0.0;

    /*** boundary energies for the i-th energy group */
    double e0 = energy[i];
    double e1 = energy[i+1];

    int j0 = 0;
    int j1 = 0;
    for(int j=1 ; j<np ; j++){
      int j2 = j*2;
      if(lib->xptr[idx][j2] >= e0){
        x0a = lib->xptr[idx][j2-2];
        y0a = lib->xptr[idx][j2-1];
        x0b = lib->xptr[idx][j2  ];
        y0b = lib->xptr[idx][j2+1];

        z0 = linear(x0a,y0a,x0b,y0b,e0);
        j0 = j; break;
      }
    }

    for(int j=j0 ; j<np ; j++){
      int j2 = j*2;
      if(lib->xptr[idx][2*j] >= e1){
        x1a = lib->xptr[idx][j2-2];
        y1a = lib->xptr[idx][j2-1];
        x1b = lib->xptr[idx][j2  ];
        y1b = lib->xptr[idx][j2+1];
        z1 = linear(x1a,y1a,x1b,y1b,e1);
        j1 = j; break;
      }
    }

    if((j0 == 0) && (j1 == 0)) continue;

    double s = 0.0;

    /*** constant weight, use trapezoid integration */
    if(weight == 0){
      /*** when no node exists inside [E0,E1] */
      if(j0 == j1) s = (e1 - e0) * (z0 + z1) * 0.5;
      /*** general case */
      else{
        s = ((x0b - e0) * (z0 + y0b) + (e1 - x1a) * (y1a + z1)) * 0.5; // both boundaries
        for(int j=j0 ; j<j1-1 ; j++){
          int j2 = j*2;
          s += (lib->xptr[idx][j2+2] - lib->xptr[idx][j2]) * (lib->xptr[idx][j2+3] + lib->xptr[idx][j2+1]) * 0.5;
        }
      }
    }

    /*** 1/E weigth, each interval is sub-devided by Ndiv */
    else if(weight == 1){
      if(j0 == j1){
        s = sinteg(weight,e0,z0,e1,z1);
      }
      else{
        /*** from E0 to the first point (x0b,y0b) */
        if(x0b > e0) s += sinteg(weight,e0,z0,x0b,y0b);

        /*** from the last point (x1a,y1a) to E1 */
        if(e1 > x1a) s += sinteg(weight,x1a,y1a,e1,z1);

        /*** all the intervals */
        for(int j=j0 ; j<j1-1 ; j++){
          int j2 = j*2;
          if(lib->xptr[idx][j2+2] > lib->xptr[idx][j2]){
            s += sinteg(weight,lib->xptr[idx][j2],lib->xptr[idx][j2+1],lib->xptr[idx][j2+2],lib->xptr[idx][j2+3]);
          }
        }
      }
    }

    /*** becase S is integral in [E0,E1], need to divide by the interval width */
    sigma[i] = s / (e1 - e0);
/*
    cout << setprecision(4);
    cout << setw(8) << i;
    cout << setw(12) << e0;
    cout << setw(12) << z0;
    cout << setw(8) << j0;
    cout << setw(12) << x0b;
    cout << setw(12) << y0b;
    cout << setw(8) << j1;
    cout << setw(12) << x1a;
    cout << setw(12) << y1a;
    cout << setw(8) << i+1;
    cout << setw(12) << e1;
    cout << setw(12) << z1;
    cout << setw(12) << sigma[i] << endl;
*/
  }
}


