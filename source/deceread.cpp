/******************************************************************************/
/**     DeCE READ                                                            **/
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"
#include "masstable.h"

static int    readCSdata    (char *, int, const int, double *, double *);
static int    readISdata    (char *, int, const int, double *, double *, double *);
static int    readNUdata    (char *, int,            double *, double *);
static int    geneCSdata    (int, double *, double *, double, double, double *);
static int    mergeCSdata   (int, double *, double *, double, double *, double *);
static double findBoundary  (ENDF *);
static double loginterpol   (int, double, double *, double *, int *);

static double eps = 1.0e-20;


/**********************************************************/
/*      Read in External Data from a File                 */
/**********************************************************/
void DeceRead(ENDFDict *dict, ENDF *lib, const int mf, const int mt, char *datafile, int ofset, bool mflag)
{
  int      nc = 0, np = 0;
  double   *cx, *cy, *xdat, elev, qm = 0.0, qi = 0.0, et = 0.0;

  if((mf!=1) && (mf!=3)) return;
  if(mt<=0 || mt>=1000) return;

  /*** allocate data array and open data file */
  cx   = new double [MAX_DBLDATA];
  cy   = new double [MAX_DBLDATA];
  xdat = new double [MAX_DBLDATA*2];

  /*** ZA, AWT, and MAT number from Dictionary */
  double za   = dict->head.c1;    //  1000*Z + A number
  double awt  = dict->head.c2;    //  mass ratio to neutron 1.008665
  double elis = dict->cont[0].c1; //  target excitation energy
  int    mat  = dict->mat;

  /*** in the case of cross sections in MF3 */
  if(mf == 3){
    if( (51<=mt && mt <=91) || (600<=mt) )
      nc = readISdata(datafile,ofset,mt,cx,cy,&elev);
    else
      nc = readCSdata(datafile,ofset,mt,cx,cy);
    if(nc==0) TerminateCode("no data to be added from a file for MT = ",mt);

    /*** find Q-values */
    /*** for MT=4, there is no way to get QI */
    qm = qvalue(dict->getProj(), (int)za,mt) + elis;
    qi = qm;
    et = 0.0;

    if(51<=mt && mt<=91){
      qi = elis - elev;
      et = threshold((int)za,qi);
    }
    else if(600<=mt && mt<=849){
      qi = qm - elev;
      if(qi<0) et = threshold((int)za,qi);
    }
    else{
      if(qm<0) et = threshold((int)za,qm);
    }

    /*** check resonance boundary */
    if(mflag){
      double ebtest = findBoundary(lib);
      if(ebtest < dict->emaxRR && ebtest != 1.0e-05){
        cerr << "maybe background cross sections given for MT = "<< mt << " at E1 = " << dict->emaxRR;
        cerr << "  E2 = " << ebtest << endl;
      }
    }

    /*** generate floating point data */
    np = nc;
    if(mflag) np = mergeCSdata(nc,cx,cy,dict->emaxRR,xdat,lib->xptr[0]);
    else      np = geneCSdata(nc,cx,cy,et,dict->emaxRR,xdat);
  }
  /*** for number of prompt/delayed neutrons, MT=455, 456 */
  else{
    nc = readNUdata(datafile,ofset,cx,cy);
    if(nc==0) TerminateCode("no data to be added from a file for MT = ",mt);

    np = nc;
    for(int i=0 ; i<np ; i++){
      xdat[2*i  ] = cx[i];
      xdat[2*i+1] = cy[i];
    }
  }

  /*** make TAB1 */
  Record cont;
  int    idat[4];

  /*** Make HEAD and CONT */
  int lnu = (mt == 3) ? 0 : 2; // tabulated nu case
  lib->setENDFhead(za,awt,0,lnu,0,0);
  lib->setENDFmat(mat);
  lib->setENDFmf(mf);
  lib->setENDFmt(mt);

  if(mflag){
    /*** keep INT in the first range (assume there is only one INT range for the resonance)*/
    if( lib->idata[1] != 2 ){
      cont.setRecord(qm,qi,0,0,2,np);
      idat[0] = lib->idata[0];
      idat[1] = lib->idata[1];
      idat[2] = np;
      idat[3] = 2;
    }
    else{
      cont.setRecord(qm,qi,0,0,1,np);
      idat[0] = np;
      idat[1] = 2;
    }
  }
  else{
    cont.setRecord(qm,qi,0,0,1,np);
    idat[0] = np;
    idat[1] = 2;
  }

  ENDFPackTAB1(cont,idat,xdat,lib);

  //  ENDFWriteHEAD(lib);
  //  ENDFWriteTAB1(lib);
  //  ENDFWriteSEND(lib);

  /*** Clean all */
  delete [] cx;
  delete [] cy;
  delete [] xdat;

  return;
}


/**********************************************************/
/*      Read in Cross Section Data                        */
/**********************************************************/
int readCSdata(char *file, int ofset, const int mt, double *x, double *y)
{
  ifstream fp;
  string   line;

  fp.open(file);
  if(!fp) TerminateCode("cannot open data file",file);

  /*** default CoH3 output file structure in CrossSection.dat */
  if(ofset == 0){
    switch(mt){
    case   1: ofset =  1; break;
    case   2: ofset =  2; break;
    case 102: ofset =  4; break;
    case   4: ofset =  5; break;
    case  18: ofset =  6; break;

    case 103: ofset =  2; break;
    case 107: ofset =  3; break;
    case 104: ofset =  4; break;
    case 105: ofset =  5; break;
    case 106: ofset =  6; break;

    case  16: ofset =  7; break;
    case  28: ofset =  8; break;
    case  22: ofset =  9; break;
    case  32: ofset = 10; break;
    case  33: ofset = 11; break;
    case  34: ofset = 12; break;

    case 111: ofset = 13; break;
    case 112: ofset = 14; break;
    case 115: ofset = 15; break;
    case 108: ofset = 16; break;
    case 117: ofset = 17; break;

    case  17: ofset = 18; break;
    case  41: ofset = 19; break;
    case  24: ofset = 20; break;
    case  45: ofset = 21; break;
    case  44: ofset = 22; break;

    case  37: ofset = 23; break;

    default : break;
    }
  }

  int nc=0;
  while(getline(fp,line)){
    if(line[0] == '#') continue;

    istringstream ss(line);
    ss >>x[nc];
    for(int i=0 ; i<ofset ; i++) ss >> y[nc];

    /*** in case blank line is given, skip it */
    if(x[nc] == 0.0) continue;


#ifdef ENERGY_UNIT_MEV
    x[nc] *= 1e+06;
#endif

#ifdef CROSS_SECTION_UNIT_MB
    y[nc] *= 0.001;
#endif

    nc++;
    if(nc>=MAX_DBLDATA) TerminateCode("too many energy points");
  }

  fp.close();

  return nc;
}


/**********************************************************/
/*      Read in Inelastic Scattering Cross Section Data   */
/**********************************************************/
int readISdata(char *file, int ofset, const int mt, double *x, double *y, double *elev)
{
  ifstream fp;
  string   line;
  double   eth = 0.0;

  fp.open(file);
  if(!fp) TerminateCode("cannot open data file",file);

  if(ofset==0){
    if( (51<=mt) && (mt<=91) )        ofset = mt -50+1;
    else if( (600<=mt) && (mt<=648) ) ofset = mt-600+1;
    else if( (650<=mt) && (mt<=698) ) ofset = mt-650+1;
    else if( (700<=mt) && (mt<=748) ) ofset = mt-700+1;
    else if( (750<=mt) && (mt<=798) ) ofset = mt-750+1;
    else if( (800<=mt) && (mt<=848) ) ofset = mt-800+1;
    else if( (mt==649) || (mt==699) || (mt==749) || (mt==799) || (mt==849) ) ofset = 42;
  }

  getline(fp,line);
  istringstream s1(&line[1]);  // skip comment #
  for(int i=0 ; i<ofset ; i++) s1 >> eth;

#ifdef ENERGY_UNIT_MEV
  *elev = eth * 1e+6;
#else
  *elev = eth;
#endif

  int nc = 0;
  while(getline(fp,line)){
    istringstream s2(line);
    s2 >>x[nc];
    for(int i=0 ; i<ofset ; i++) s2 >> y[nc];

#ifdef ENERGY_UNIT_MEV
    x[nc] *= 1e+06;
#endif
#ifdef CROSS_SECTION_UNIT_MB
    y[nc] *= 0.001;
#endif

    if( (mt >= 600) || (y[nc] > 0.0) ) nc++;
    if(nc>=MAX_DBLDATA) TerminateCode("too many energy points");
  }

  fp.close();

  /*** check non-zero data */
  bool zero = true;
  for(int i=0 ; i<nc ; i++){
    if(y[i]>0.0){
      zero = false;
      break;
    }
  }
  if(zero) nc = 0;

  return nc;
}


/**********************************************************/
/*      Read in Nu-p, Nu-d Data                           */
/**********************************************************/
int readNUdata(char *file, int ofset, double *x, double *y)
{
  ifstream fp;
  string   line;

  fp.open(file);
  if(!fp) TerminateCode("cannot open data file",file);

  if(ofset == 0) ofset = 1;

  int nc=0;
  while(getline(fp,line)){
    if(line[0] == '#') continue;

    istringstream ss(line);
    ss >>x[nc];
    for(int i=0 ; i<ofset ; i++) ss >> y[nc];

    /*** in case blank line is given, skip it */
    if(x[nc] == 0.0) continue;

#ifdef ENERGY_UNIT_MEV
    x[nc] *= 1e+06;
#endif

    nc++;
    if(nc>=MAX_DBLDATA) TerminateCode("too many energy points");
  }

  fp.close();

  return nc;
}


/**********************************************************/
/*      Replace XY Data by Inputs                         */
/**********************************************************/
int geneCSdata(int n, double *x, double *y, double eth, double eres, double *xdat)
{
  int i = 0;

  /*** for the threshold reaction */
  if(eth > 0.0){
    xdat[i++] = eth;
    xdat[i++] = 0.0;

    /*** if E(threshold) is inside the resonance range
         place zeros from Eth to Eres */
    if(eth < eres){
      xdat[i++] = eres;
      xdat[i++] = 0.0;

      int skip = 0;
      double yint = loginterpol(n,eres,x,y,&skip);
      xdat[i++] = eres;
      xdat[i++] = yint;

      for(int j=skip ; j<n ; j++){
        if(y[j] >= eps){
          xdat[i++] = x[j];
          xdat[i++] = y[j];
        }
      }
    }
    /*** if Eth is larger than Eres, omit resonance range */
    else{
      for(int j=0 ; j<n ; j++){
        if(y[j] >= eps){
          xdat[i++] = x[j];
          xdat[i++] = y[j];
        }
      }
    }
  }

  /*** for the non-threshold reaction */
  else{
    /*** Insert thermal point, start at the resonance boundary */
    xdat[i++] = 1.0000e-05;
    xdat[i++] = 0.0;
    xdat[i++] = 2.5300e-02;
    xdat[i++] = 0.0;
    xdat[i++] = eres;
    xdat[i++] = 0.0;

    int skip = 0;
    double yint = loginterpol(n,eres,x,y,&skip);

    if(yint>0.0){
      xdat[i++] = eres;
      xdat[i++] = yint;
    }

    bool onetrip = false;
    for(int j=skip ; j<n ; j++){
      if(y[j] >= eps){
        /*** insert 1-point before data start */
        if(skip>1 && !onetrip && j>=1 && y[j-1]==0.0){
          xdat[i++] = x[j-1];
          xdat[i++] = 0.0;
          onetrip = true;
        }
        xdat[i++] = x[j];
        xdat[i++] = y[j];
      }
    }
  }

  return(i/2);
}


/**********************************************************/
/*      Resonance Background + New XY Data                */
/**********************************************************/
int mergeCSdata(int n, double *x, double *y, double eres, double *xdat, double *xbak)
{
  int i = 0;

  /*** copy old data up to Eres */
  for(int j=0 ; ; j++){
    xdat[i++] = xbak[2*j  ];
    xdat[i++] = xbak[2*j+1];
    if(xbak[2*j] >= eres) break;
  }

  /*** start at the resonance boundary */
  int skip = 0;
  double yint = loginterpol(n,eres,x,y,&skip);

  if(yint>0.0){
    xdat[i++] = eres;
    xdat[i++] = yint;
  }

  bool onetrip = false;
  for(int j=skip ; j<n ; j++){
    if(y[j] >= eps){
      /*** insert 1-point before data start */
      if(skip>1 && !onetrip && j>=1 && y[j-1]==0.0){
        xdat[i++] = x[j-1];
        xdat[i++] = 0.0;
        onetrip = true;
      }
      xdat[i++] = x[j];
      xdat[i++] = y[j];
    }
  }

  return(i/2);
}


/**********************************************************/
/*      Interpolation                                     */
/**********************************************************/
double loginterpol(int n, double e, double *x, double *y, int *idx)
{
  double z = 0.0;

  /*** if same energy point is given */
  if(e == x[0]){
    z = -1;
    *idx = 0;
  }
  /*** if E < X[], use 1/v dependence */
  else if(e < x[0]){
    z = y[0] * sqrt(x[0]/e);
    *idx = 0;
  }
  /*** if outsize, use the same highest value */
  else if(e > x[n-1]){
    z = y[n-1];
    *idx = n;
  }
  /*** otherwise, linear interpolate in log-log */
  else{
    int m = 0;
    for(int i=0 ; i<n-1 ; i++){
      if(x[i] <= e && e < x[i+1]){ m = i; break; }
    }
    z = (y[m+1]-y[m])/(log(x[m+1])-log(x[m]))*(log(e)-log(x[m])) + y[m];
    *idx = m+1;
  }
  return (z);
}


/**********************************************************/
/*      Estimate High-Side of Resonance Range             */
/**********************************************************/
double findBoundary(ENDF *lib)
{
  const double e0 = 1.0e-05;
  double e    = lib->xdata[0];
  Record cont = lib->rdata[0];
  int    np   = cont.n2;

  /*** if x[0] = 0, maybe new data */
  if(e == 0.0) return(e0);

  /*** for threshold reaction, x[0] > 1.0e-5 */
  if(e > e0) return(e);

  /*** in many cases, X-val is duplicated at the boundary of 
       resonance range */
  for(int i=2 ; i<np*2 ; i+=2){
    if( (lib->xdata[i-2] == lib->xdata[i]) &&
        (lib->xdata[i-1] == 0.0) && (lib->xdata[i+1] > 0.0) ){
      e = lib->xdata[i];
    }
  }

  return(e);
}
