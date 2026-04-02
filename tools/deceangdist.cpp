/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate Elastic and Inelastic  Angular Distributions   **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <unistd.h>

using namespace std;

#include "../source/endflib.h"

int    main (int, char *[]);
void   tabulateAngDist (const double, const double, const double, ENDF *, ENDF *);
void   findEnergyMF4 (const double, ENDF *);
void   findEnergyMF6 (const double, ENDF *);
int    calcAngDist (const double, const double, ENDF *);
void   calcAngDistTAB (const double, const int, double *, ENDF *);
void   calcAngDistLEG (const double, const int, double *, ENDF *);
void   printAngDist (const int, const double, const double);
double legendre (int, double);
void   usage ();

inline double interpolation(double x1, double x2, double y1, double y2, double x)
{ return( (y2-y1)/(x2-x1) * (x-x1) + y1 ); }

static const double DEFAULT_ANGSTEP  =  1.0;  // default calculation angle increment
static const int    MAX_ANGLE        =  360;  // max number of angular points
static const int    MAX_EGRID        =   61;
static double default_egrid[MAX_EGRID] = {
  0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,
 1.25,  1.5, 1.75,  2.0, 2.25,  2.5, 2.75,  3.0,  3.5,  4.0,
  4.5,  5.0,  5.5,  6.0,  6.5,  7.0,  7.5,  8.0,  8.5,  9.0,
  9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0,
 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0,
 19.5, 20.0};

static double *xdat, *ydat;
static int dataptr = 0, mf = 4, mt = 2;

/**********************************************************/
/*      Differential Scattering Angular Distribution      */
/**********************************************************/
int main(int argc, char *argv[])
{
  ifstream fpin;
  string   libname = "";
  ENDF     libcs, libad;
  double   angstep = DEFAULT_ANGSTEP;
  double   ein = 0.0;
  bool     autoenergy = false;

  //----------------------------------- Command Line Options
  int p = 0;
  while((p = getopt(argc,argv,"e:s:t:ah")) != -1){
    switch(p){
    case 'e':
      ein = atof(optarg);      break;
    case 's':
      angstep = atof(optarg);  break;
    case 't':
      mt = atoi(optarg);       break;
    case 'a':
      autoenergy = true;
      ein = default_egrid[0];  break;
    case 'h':
      usage();                 break;
    default:                   break;
    }
  }
  if(optind < argc) libname = argv[optind];

  /*** check input data */
  if(libname == "" || ein == 0.0){ usage(); }
  if(angstep < 0.0 || angstep > 180.0){
    cerr << "invalid angle step " << angstep << endl;  exit(-1);
  }
  if(mt != 2){
    if( (mt < 51) || (mt > 90) ){
      cerr << "invalid MT number " << mt << endl; exit(-1);
    }
  }

  /*** scan ENDF file and store data in dict */
  ENDFDict dict;
  if(ENDFScanLibrary(libname,&dict) < 0){
    cerr << "ENDF file cannot open " << libname << endl;  exit(-1);
  }


  //--------------------------------------- Main Calculation
  /*** open ENDF file */
  fpin.open(libname.c_str());

  /*** look for the resonance boundary */
  ENDFReadMF2(&fpin,&libcs,151);
  ENDFMF2boundary(&dict,&libcs);

  /*** allocate data arrays */
  xdat = new double [MAX_ANGLE]; // angular points
  ydat = new double [MAX_ANGLE]; // differential scattering probabilities

  /*** read ENDF tape */
  libcs.resetPOS(); // reuse libcs object
  ENDFReadMF3(&fpin,&libcs,mt);

  /*** determine MF number by MT */
  if(mt == 2){
    mf = 4;
    ENDFRead(&fpin,&libad,mf,mt);
  }
  else{
    int cnd = ENDFRead(&fpin,&libad,4,mt); // try if MF4 is given
    if(cnd >= 0) mf = 4;
    else{
      cnd = ENDFRead(&fpin,&libad,6,mt);   // if not, see MF6
      if(cnd >= 0) mf = 6;
      else{
        cerr << "MF:MT " << mf << ":" << mt << " not given in " << libname << endl; exit(-1);
      }
    }
  }

  /*** close file */
  fpin.close();

  /*** process at each incident energy */
  if(autoenergy){
    for(int i=0 ; i<MAX_EGRID ; i++){
      ein = default_egrid[i] * 1e+6;
      tabulateAngDist(ein,dict.emaxRe,angstep,&libcs,&libad);
    }
  }
  else tabulateAngDist(ein,dict.emaxRe,angstep,&libcs,&libad);


  //------------------------------------------- Post Process
  delete [] xdat;
  delete [] ydat;

  return(0);
}


/**********************************************************/
/*      Tabulate Scattering Angular Distribution          */
/**********************************************************/
void tabulateAngDist(const double ein, const double eres, const double step, ENDF *libcs, ENDF *libad)
{
  if(mf == 4) findEnergyMF4(ein,libad);
  else        findEnergyMF6(ein,libad);

  int nt = calcAngDist(ein,step,libad);

  /*** when inside resonance, take the cross section at 1 keV above the boundary */
  double ex = ein;
  if(ein < eres) ex = eres + 1000.0;

  double sigma = ENDFInterpolation(libcs,ex,false,0);

  /*** print results */
  printAngDist(nt,ein,sigma);
}


/**********************************************************/
/*      Locate Angular Distribution Data at Given Energy  */
/**********************************************************/
void findEnergyMF4(const double ein, ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    ltt  = head.l2;   // 0: isotropic, 1: Legendre, 2: tabulated
  int    idx  = 0;
  Record cont = lib->rdata[idx ++];
  int    li   = cont.l1;   // 0: non-isotropic, 1: isotropic

  bool found  = false;
  int  dtype  = 0;

  if(li == 1){
    dataptr = 0; // isotropic case
    return;
  }
  else{
    /*** find two energy points Ein lies in between  */
    Record cont = lib->rdata[idx ++];
    int    ne = cont.n2;

    for(int i=0 ; i<ne-1 ; i++){
      double e1  = lib->rdata[idx  ].c2;
      double e2  = lib->rdata[idx+1].c2;
      if(e1 < ein && ein <= e2 ){
        dtype = (ltt == 3) ? 1 : ltt;;
        found = true;
        break;
      }
      idx ++;
    }

    /*** when LTT = 3, and Ein is higher than the first energy region */
    if(!found && (ltt == 3)){
      idx ++;
      cont = lib->rdata[idx ++];
      ne = cont.n2;
      for(int i=0 ; i<ne-1 ; i++){
        double e1  = lib->rdata[idx  ].c2;
        double e2  = lib->rdata[idx+1].c2;
        if(e1 < ein && ein <= e2 ){
          dtype = 2;
          found = true;
          break;
        }
        idx++;
      }
    }
  }

  /*** set pointer where angular distribution data are given */
  if(!found) dataptr = 0;
  else if(dtype == 1) dataptr =  idx;
  else                dataptr = -idx; // for tabulated data case, set negative

  return;
}


void findEnergyMF6(const double ein, ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;

  if(nk == 0) return;

  Record cont = lib->rdata[idx ++];
  int    law  = cont.l2;   // LAW, should be 2 or 3

  bool found  = false;
  int  dtype  = 0;

  if(law == 3){
    dataptr = 0; // isotropic case
    return;
  }
  else if(law == 2){
    Record cont = lib->rdata[idx ++];
    double e0  = lib->rdata[idx].c2;
    int    ne  = cont.n2;

    /*** check the lowest energy first */
    if(ein < e0){
      dataptr = 0; // set to isotropic
      return;
    }

    for(int i=0 ; i<ne-1 ; i++){
      double e1  = lib->rdata[idx  ].c2;
      double e2  = lib->rdata[idx+1].c2;
      if(e1 < ein && ein <= e2 ){
        dtype = (lib->rdata[idx].l1 == 0) ? 1 : 2;
        found = true;
        break;
      }
      idx ++;
    }
  }
  else{
    cerr << "cannot handle LAW = " << law << " in MF6" << endl; exit(-1);
  }

  if(!found) dataptr = 0;
  else if(dtype == 1) dataptr =  idx;
  else{
    cerr << "tabulated angular distribution in MF6 not yet implemented" << endl; exit(-1);
  }

  return;
}


/**********************************************************/
/*      Calculate Scattering Angular Distribution         */
/**********************************************************/
int calcAngDist(const double ein, const double step, ENDF *lib)
{
  /*** angle points */
  int nt = 0;
  double theta = 0.0;
  while(theta <= 180.0){
    if(nt >= MAX_ANGLE) break;
    xdat[nt] = theta;
    theta += step;
    nt ++;
  }

  /*** isotropic angular distribution */
  if(dataptr == 0){
    for(int n=0 ; n<nt ; n++) ydat[n] = 0.5;
  }

  /*** data given by table */
  else if(dataptr < 0){
    dataptr *= -1;
    double *x = new double [nt];
    for(int n=0 ; n<nt ; n++) x[n] = cos(M_PI*xdat[n]/180.0);

    calcAngDistTAB(ein,nt,x,lib);
    delete [] x;
  }

  /*** data given by Legendre coefficients */
  else{
    calcAngDistLEG(ein,nt,xdat,lib);
  }

  return nt;
}


void calcAngDistTAB(const double ein, const int nt, double *x, ENDF *lib)
{
  double e1  = lib->rdata[dataptr  ].c2; // energy grid lower than Ein
  double e2  = lib->rdata[dataptr+1].c2; // higher Ein
  int    np1 = lib->rdata[dataptr  ].n2;
  int    np2 = lib->rdata[dataptr+1].n2;

  double py1 = 0.0, py2 = 0.0;

  for(int n=0 ; n<nt ; n++){

    /*** find angle point at E1, and interpolate Y(E1) */
    for(int ip=0 ; ip<np1-1; ip++){
      double pxa = lib->xptr[dataptr][2*ip  ];
      double pxb = lib->xptr[dataptr][2*ip+2];
      if(pxa < x[n] && x[n] <= pxb){
        double pya = lib->xptr[dataptr][2*ip+1];
        double pyb = lib->xptr[dataptr][2*ip+3];
        py1 = interpolation(pxa,pxb,pya,pyb,x[n]);
        break;
      }
    }

    /*** find angle point at E2, and interpolate Y(E2) */
    for(int ip=0 ; ip<np2-1; ip++){
      double pxa = lib->xptr[dataptr+1][2*ip  ];
      double pxb = lib->xptr[dataptr+1][2*ip+2];
      if(pxa < x[n] && x[n] <= pxb){
        double pya = lib->xptr[dataptr+1][2*ip+1];
        double pyb = lib->xptr[dataptr+1][2*ip+3];
        py2 = interpolation(pxa,pxb,pya,pyb,x[n]);
        break;
      }
    }

    /*** interpolate Y(E1) and Y(E2) at Ein */
    ydat[n] = interpolation(e1,e2,py1,py2,ein);
  }
}


void calcAngDistLEG(const double ein, const int nt, double *x, ENDF *lib)
{
  double e1  = lib->rdata[dataptr  ].c2;
  double e2  = lib->rdata[dataptr+1].c2;
  int    nl1 = lib->rdata[dataptr  ].n1;
  int    nl2 = lib->rdata[dataptr+1].n1;

  /*** interpolate Legendre coefficieints at Ein */
  double *pl  = new double [nl2]; // hope nl2 >= nl1
  for(int j=0 ; j<nl2 ; j++){
    double pl1 = 0.0;
    if(j < nl1) pl1 = lib->xptr[dataptr][j];
    double pl2 = lib->xptr[dataptr+1][j];
    pl[j]  = interpolation(e1,e2,pl1,pl2,ein);
  }

  /*** construct angular distribution probabilities */
  for(int n=0 ; n<nt ; n++){
    ydat[n] = 0.5;
    for(int j=0 ; j<nl2 ; j++){
      ydat[n] += (j+1.5) * pl[j] * legendre(j+1,x[n]);
    }
  }

  delete [] pl;
}


/**********************************************************/
/*      Print Differential Scattering Cross Section       */
/**********************************************************/
void printAngDist(const int nt, const double ein, const double sigma)
{
  cout << setprecision(6) << setiosflags(ios::scientific);

  cout << "# Energy[eV]    Angle[deg]    Probability   dsig[b/st]" << endl;
  for(int j=0 ; j<nt ; j++){
    cout << setw(14) << ein;
    cout << setw(14) << xdat[j];
    cout << setw(14) << ydat[j];
    cout << setw(14) << ydat[j] * sigma / (2.0*M_PI) << endl;
  }
  cout << endl;
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


/**********************************************************/
/*     Help                                               */
/**********************************************************/
void usage()
{
  cout <<
    "usage\n"
    " % deceangdist [-e energy] ... ENDF_file\n"
    "      ENDF_file : input ENDF-6 formatted file\n"
    "      -e energy : calculate angular distribution at energy [eV]\n"
    "      -a        : angular distributions at built-in energy grid\n"
    "                  -e option will be ignored\n"
    "      -t MT     : MTnumber for inelastic scattering, 51 - 90\n"
    "                  when not given, MT=4 (elastic scattering) assumed\n"
    "      -h        : this help\n"
    "\n";
  cout << endl;
  exit(0);
}
