/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate Elastic Angular Distributions from MF3/4 MT2   **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

int    main (int, char *[]);
void   dataout (ENDF *);
void   processMF4 (ENDF *);
int    processMF4LEG (ENDF *, int);
int    processMF4TAB (ENDF *, int);
double legendre (int, double);

static const int ANGSTEP    =    1;  // calculation angle increment
static const int MAX_ANGLE  =  360;  // max number of angular points
static const int MAX_ENERGY = 1000;  // max number of incident energies

static double **xdat, **ydat, *edat;
static int    *ndat, dptr = 0;

int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: deceangdist ENDF_file" << endl;  exit(-1);
  }

  ifstream fpin;
  string   libname = "";
  ENDF     lib3, lib4;

  libname = argv[1];

  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;  exit(-1);
  }

  xdat = new double * [MAX_ENERGY];
  ydat = new double * [MAX_ENERGY];
  for(int i=0 ; i<MAX_ENERGY ; i++){
    xdat[i] = new double [MAX_ANGLE];
    ydat[i] = new double [MAX_ANGLE];
  }
  edat = new double [MAX_ENERGY];
  ndat = new int [MAX_ENERGY];

  ENDFReadMF3(&fpin,&lib3,2);
  ENDFReadMF4(&fpin,&lib4,2);
  fpin.close();

  processMF4(&lib4);
  dataout(&lib3);

  for(int i=0 ; i<MAX_ENERGY ; i++){
    delete [] xdat[i];
    delete [] ydat[i];
  }
  delete [] xdat;
  delete [] ydat;
  delete [] edat;
  delete [] ndat;

  return(0);
}


/**********************************************************/
/*      Differential Scattering Cross Section Output      */
/**********************************************************/
void dataout(ENDF *lib3)
{
  cout << setprecision(6) << setiosflags(ios::scientific);

  for(int i=0 ; i<dptr ; i++){
    double sig = ENDFInterpolation(lib3,edat[i],false,0);
    if(sig < 1.0e-10) continue; // this could be in the resonance region

    cout << "# Energy[eV]    Angle[deg]    Probability   dsig[b/st]" << endl;
    for(int j=0 ; j<ndat[i] ; j++){
      cout << setw(14) << edat[i];
      cout << setw(14) << xdat[i][j];
      cout << setw(14) << ydat[i][j];
      cout << setw(14) << ydat[i][j] * sig / (2.0*M_PI) << endl;
    }
    cout << endl;
    cout << endl;
  }
}


/**********************************************************/
/*      Process MF=4                                      */
/**********************************************************/
void processMF4(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    ltt  = head.l2;   // 0: isotropic, 1: Legendre, 2: tabulated
  int    idx  = 0;
  Record cont = lib->rdata[idx];
  int    li   = cont.l1;   // 0: non-isotropic, 1: isotropic

  if(li == 1){
    cout << "# isotropic angular distribution" << endl;
  }
  else{
    idx++;
    if(ltt == 1){
      idx = processMF4LEG(lib,idx);
    }
    else if(ltt == 2){
      idx = processMF4TAB(lib,idx);
    }
    else if(ltt == 3){
      idx = processMF4LEG(lib,idx);
      idx = processMF4TAB(lib,idx);
    }
  }
}


/**********************************************************/
/*      Angular Distribution in Legendre Coefficient      */
/**********************************************************/
int processMF4LEG(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];
  int    ne   = cont.n2;
  idx++;

  for(int i=0 ; i<ne ; i++){
    double e  = lib->rdata[idx].c2;
    int    nl = lib->rdata[idx].n1;

    int np = 180/ANGSTEP;
    if( (180 % ANGSTEP) == 0 ) np++;

    if(dptr < MAX_ENERGY){
      ndat[dptr] = np;
      edat[dptr] = e;
    }

    int k = 0;
    for(int t=0 ; t<=180 ; t+=ANGSTEP){
      double f = 0.5;
      for(int j=0 ; j<nl ; j++){
        f += (j+1.5)*lib->xptr[idx][j]*legendre(j+1,(double)t);
      }

      if((dptr < MAX_ENERGY) && (k < MAX_ANGLE)){
        xdat[dptr][k] = t;
        ydat[dptr][k] = f;
      }
      k++;
    }
    idx++;
    dptr++;
  }

  return(idx);
}


/**********************************************************/
/*      Angular Distribution in Table                     */
/**********************************************************/
int processMF4TAB(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];
  int    ne   = cont.n2;
  idx++;

  for(int i=0 ; i<ne ; i++){
    double e  = lib->rdata[idx].c2;
    int    np = lib->rdata[idx].n2;

    if(dptr < MAX_ENERGY){
      ndat[dptr] = np;
      edat[dptr] = e;
    }

    int k = 0;
    for(int j=0 ; j<np ; j++){
      double t = acos(lib->xptr[idx][2*j])/M_PI * 180.0;
      double f = lib->xptr[idx][2*j+1];

      if((dptr < MAX_ENERGY) && (k < MAX_ANGLE)){
        xdat[dptr][k] = t;
        ydat[dptr][k] = f;
      }
      k++;
    }
    idx++;
    dptr++;
  }

  return(idx);
}


/**********************************************************/
/*     Legendre Function                                  */
/**********************************************************/
double legendre(int n, double t)
{
  double x  = cos(M_PI*t/180.0);
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


