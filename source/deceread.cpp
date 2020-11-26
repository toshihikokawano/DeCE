/******************************************************************************/
/**     DeCE READ                                                            **/
/******************************************************************************/

#include <iostream>
#include <ostream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "global.h"
#include "terminate.h"
#include "masstable.h"

static int    readCSdata    (char *, int, const int, double *, double *);
static int    readISdata    (char *, int, const int, double *, double *, double *);
static int    readNUdata    (char *, int,            double *, double *);
static int    geneCSdata    (int, double *, double *, double, double, double *);
static int    geneCSdata1   (int, double *, double *, double, double, double *);
static int    geneCSdata2   (int, double *, double *, double, double *);
static int    mergeCSdata   (int, double *, double *, double, double *, double *);
static double findBoundary  (ENDF *);
static double loginterpol   (int, double, double *, double *, int *);

static bool   charged_particle_file = false;


/**********************************************************/
/*      Read in External Data from a File                 */
/**********************************************************/
void DeceRead(ENDFDict *dict, ENDF *lib, const int mf, const int mt, char *datafile, int ofset, bool mflag)
{
  int      nc = 0, np = 0;
  double   *cx, *cy, *xdat, elev = 0.0, qm = 0.0, qi = 0.0, et = 0.0;
  ostringstream os;

  if((mf != 1) && (mf != 3)){
    message << "MF" << mf << " different from MF1 or MF3";
    WarningMessage();
    return;
  }
  if((mt <= 0) || (mt >= 1000)){
    message << "MT" << mt << " out of range";
    WarningMessage();
    return;
  }

  /*** allocate data array and open data file */
  cx   = new double [MAX_DBLDATA];
  cy   = new double [MAX_DBLDATA];
  xdat = new double [MAX_DBLDATA*2];

  /*** ZA, AWR, and MAT number from Dictionary */
  double za   = dict->getZA();    //  1000*Z + A number
  double awr  = dict->getAWR();   //  mass ratio to neutron 1.008665
  double elis = dict->getELIS();  //  target excitation energy
  int    mat  = dict->getMAT();

  /*** determine if the file is for charged particle incident reactions */
  if(dict->getProj() > 1) charged_particle_file = true;

  /*** in the case of cross sections in MF3 */
  if(mf == 3){
    /*** cross section to discrete levels */
    if( (51 <= mt && mt <= 91) || (600 <= mt && mt <= 849) ){
      nc = readISdata(datafile,ofset,mt,cx,cy,&elev);
    }
    /*** general case */
    else{
      nc = readCSdata(datafile,ofset,mt,cx,cy);
    }
    if(nc == 0){
      message << "no cross section data to be added from " << datafile << " for MT = " << mt;
      WarningMessage();
    }

    /*** find Q-values */
    /*** for MT=4, there is no way to get QI */
    qm = qvalue(dict->getProj(), (int)za,mt) + elis;
    qi = qm;
    et = 0.0;

    /*** determine the threshold energy */
    if( (51 <= mt) && (mt <= 91) ){
      /*** inelastic scattering */
      qi = elis - elev;
      et = threshold((int)za,qi);
    }
    else if( (600 <= mt) && (mt <= 849) ){
      /*** discrete transition by charged particle */
      qi = qm - elev;
      if(qi < 0.0) et = threshold((int)za,qi);
    }
    else{
      /*** all other reactions */
      if(qm < 0.0) et = threshold((int)za,qm);
    }

    /*** check resonance boundary, when data will be merged */
    if(mflag){
      double ebtest = findBoundary(lib);
      if(ebtest < dict->emaxRe && ebtest != 1.0e-05){
        message << "maybe background cross sections given for MT = " << mt << " at E1 = " << dict->emaxRe << "  E2 = " << ebtest;
        WarningMessage();
      }
    }

    message << "Q(mass) " << qm << " Q(level) " << qi << " Threshold Energy " << et << " Resonance Boundary " << dict->emaxRe;
    Notice("DeceRead");

    /*** generate floating point data */
    np = nc;
    if(mflag) np = mergeCSdata(nc,cx,cy,dict->emaxRe,xdat,lib->xptr[0]);
    else      np = geneCSdata(nc,cx,cy,et,dict->emaxRe,xdat);

    message << "number of points added " << np;
    Notice("DeceRead");

  }
  /*** for number of prompt/delayed neutrons, MT=455, 456 */
  else{
    nc = readNUdata(datafile,ofset,cx,cy);
    if(nc == 0){
      message << "no nu data to be added from " << datafile << " for MT = " << mt;
      WarningMessage();
    }

    np = nc;
    for(int i=0 ; i<np ; i++){
      xdat[2*i  ] = cx[i];
      xdat[2*i+1] = cy[i];
    }
  }

  if(np > 1){
    /*** make TAB1 */
    Record cont;
    int    idat[4];

    /*** Make HEAD and CONT */
    int lnu = (mf == 3) ? 0 : 2; // tabulated nu case
    lib->setENDFhead(za,awr,0,lnu,0,0);
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
  }
  else{
    if(!mflag) DeceDelete(dict,mf,mt);
  }

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
    if(line.length() == 0) continue;

    istringstream ss(line);
    ss >> x[nc];
    for(int i=0 ; i<ofset ; i++) ss >> y[nc];

    /*** in case blank line is given, skip it */
    if(x[nc] == 0.0) continue;

    /*** convert energy unit into eV */
    x[nc] *= opt.ReadXdataConversion;

    /*** convert cross section unit into barns */
    y[nc] *= opt.ReadYdataConversion;

    /*** skip data if range is set by options */
    if(DeceCheckReadRange(x[nc])) continue;

    nc++;
    if(nc >= MAX_DBLDATA) TerminateCode("too many energy points");
  }

  fp.close();

  if(nc >= 1){
    message << "MF3:MT" << mt << " " << nc << " points from (" << x[0] << "," << y[0] << ") to (" << x[nc-1] << "," << y[nc-1] << ") imported from " << file;
    Notice("DeceRead:readCSdata");
  }

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

  if(ofset == 0){
    if(      ( 51 <= mt) && (mt <=  91) ) ofset = mt  - 50 + 1;
    else if( (600 <= mt) && (mt <= 640) ) ofset = mt - 600 + 1;
    else if( (650 <= mt) && (mt <= 690) ) ofset = mt - 650 + 1;
    else if( (700 <= mt) && (mt <= 740) ) ofset = mt - 700 + 1;
    else if( (750 <= mt) && (mt <= 790) ) ofset = mt - 750 + 1;
    else if( (800 <= mt) && (mt <= 840) ) ofset = mt - 800 + 1;
    else if( (mt == 649) || (mt == 699) || (mt == 749) || (mt == 799) || (mt == 849) ) ofset = 42;
  }

  getline(fp,line);
  istringstream s1(&line[1]);  // skip comment #
  for(int i=0 ; i<ofset ; i++) s1 >> eth;

  *elev = eth * opt.ReadXdataConversion;

  int nc = 0;
  while(getline(fp,line)){
    if(line.length() == 0) continue;

    istringstream s2(line);
    s2 >> x[nc];
    for(int i=0 ; i<ofset ; i++) s2 >> y[nc];

    /*** skip blank line */
    if(x[nc] == 0.0) continue;

    x[nc] *= opt.ReadXdataConversion;
    y[nc] *= opt.ReadYdataConversion;

    /*** skip data if range is set by options */
    if(DeceCheckReadRange(x[nc])) continue;

    if( (mt >= 600) || (y[nc] > 0.0) ) nc++;
    if(nc >= MAX_DBLDATA) TerminateCode("too many energy points");
  }

  fp.close();

  /*** check non-zero data */
  bool zero = true;
  for(int i=0 ; i<nc ; i++){
    if(y[i] > 0.0){
      zero = false;
      break;
    }
  }
  if(zero) nc = 0;

  if(nc >= 1){
    message << "MF3:MT" << mt << " " << nc << " points from (" << x[0] << "," << y[0] << ") to (" << x[nc-1] << "," << y[nc-1] << ") imported from " << file;
    Notice("DeceRead:readISdata");
  }

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
    if(line.length() == 0) continue;

    istringstream ss(line);
    ss >> x[nc];
    for(int i=0 ; i<ofset ; i++) ss >> y[nc];

    /*** in case blank line is given, skip it */
    if(x[nc] == 0.0) continue;

    x[nc] *= opt.ReadXdataConversion;

    /*** skip data if range is set by options */
    if(DeceCheckReadRange(x[nc])) continue;

    nc++;
    if(nc >= MAX_DBLDATA) TerminateCode("too many energy points");
  }

  fp.close();

  if(nc >= 1){
    message << "MF3:MT455(6) " << nc << " points from (" << x[0] << "," << y[0] << ") to (" << x[nc-1] << "," << y[nc-1] << ") imported from " << file;
    Notice("DeceRead:readNUdata");
  }

  return nc;
}


/**********************************************************/
/*      Replace XY Data by Inputs                         */
/**********************************************************/
int geneCSdata(int n, double *x, double *y, double eth, double eres, double *xdat)
{
  if(n == 0) return 0;

  int i = 0;

  /*** for the threshold reaction */
  if(eth > 0.0){
    i = geneCSdata1(n,x,y,eth,eres,xdat);
  }
  /*** for the non-threshold reaction */
  else{
    i = geneCSdata2(n,x,y,eres,xdat);
  }

  return(i/2);
}


/**********************************************************/
/*      Threshold Reaction Case                           */
/**********************************************************/
int geneCSdata1(int n, double *x, double *y, double eth, double eres, double *xdat)
{
  int i = 0;

  xdat[i++] = eth;
  xdat[i++] = 0.0;

  /*** if E(threshold) is inside the resonance range place zeros from Eth to Eres */
  if(eth < eres){
    xdat[i++] = eres;
    xdat[i++] = 0.0;

    /*** duplicated point at Eres */
    int skip = 0;
    double yint = loginterpol(n,eres,x,y,&skip);
    xdat[i++] = eres;
    xdat[i++] = yint;

    /*** copy all the rest */
    for(int j=skip ; j<n ; j++){
      xdat[i++] = x[j];
      xdat[i++] = y[j];
    }
  }
  /*** if Eth is larger than Eres, omit entire resonance range */
  else{
    /*** remove points below Eth if given */
    int skip = 0;
    for(int j=0 ; j<n ; j++){
      if(x[j] > eth){ skip = j; break; }
    }

    for(int j=skip ; j<n ; j++){
      xdat[i++] = x[j];
      xdat[i++] = y[j];
    }
  }

  return i;
}


/**********************************************************/
/*      Non-Threshold Reaction Case                       */
/**********************************************************/
int geneCSdata2(int n, double *x, double *y, double eres, double *xdat)
{
  const double e0 = 1.0000e-05; // lowest energy
  const double e1 = 2.5300e-02; // thermal point

  int i = 0;
  int skip = 0;

  /*** when resonance region exists */
  if(eres > 0.0){
    /*** Insert thermal point, start at the resonance boundary */
    xdat[i++] = e0;
    xdat[i++] = 0.0;
    if(!charged_particle_file){
      xdat[i++] = e1;
      xdat[i++] = 0.0;
    }

    /*** duplicated point at Eres, we hope Eres is bigger than 0.0253 */
    xdat[i++] = eres;
    xdat[i++] = 0.0;

    /*** find Y-value at Eres */
    double yint = loginterpol(n,eres,x,y,&skip);
    xdat[i++] = eres;
    xdat[i++] = yint;
  }

  /*** no resonance case */
  else{
    xdat[i++] = e0;
    xdat[i++] = loginterpol(n,e0,x,y,&skip);

    if(!charged_particle_file){
      /*** insert data points when energies below thermal are given */
      bool thermal = false;
      for(skip=0 ; skip<n ; skip++){

        double eps = fabs(x[skip] / e1 -1.0);

        if( (e0 < x[skip]) && (x[skip] < e1) ){
          xdat[i++] = x[skip];
          xdat[i++] = y[skip];
        }
        /*** if thermal is already given */
        else if(eps < 1e-10){
          xdat[i++] = x[skip];
          xdat[i++] = y[skip];
          thermal = true;
          break;
        }
      }

      if(!thermal){
        xdat[i++] = e1;
        xdat[i++] = loginterpol(n,e1,x,y,&skip);
      }
      else skip ++;

      /*** avoid long linear interpolation by duplicating the first data point */
      if( (xdat[i-1] == 0.0) && (y[skip] > 0.0) ){
        xdat[i++] = x[skip];
        xdat[i++] = 0.0;
      }
    }
  }

  /*** copy all the rest */
  for(int j=skip ; j<n ; j++){
    xdat[i++] = x[j];
    xdat[i++] = y[j];
  }

  return i;
}


/**********************************************************/
/*      Resonance Background + New XY Data                */
/**********************************************************/
int mergeCSdata(int n, double *x, double *y, double eres, double *xdat, double *xbak)
{
  if(n == 0) return 0;

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
  xdat[i++] = eres;
  xdat[i++] = yint;

  bool onetrip = false;
  for(int j=skip ; j<n ; j++){
    /*** insert 1-point before data start */
    if(skip>1 && !onetrip && j>=1 && y[j-1]==0.0){
      xdat[i++] = x[j-1];
      xdat[i++] = 0.0;
      onetrip = true;
    }
    xdat[i++] = x[j];
    xdat[i++] = y[j];
  }

  return(i/2);
}


/**********************************************************/
/*      Interpolation                                     */
/**********************************************************/
double loginterpol(int n, double e, double *x, double *y, int *idx)
{
  double z = 0.0;

  /*** if E is lower than the given data range, set zero  */
  if(e < x[0]){
    z = 0.0;
    *idx = 0;
  }
  /*** if higher than the range, use the same highest value */
  else if(e > x[n-1]){
    z = y[n-1];
    *idx = n;
  }
  /*** otherwise interpolate */
  else{
    int m = 0;
    for(int i=0 ; i<n-1 ; i++){
      if(x[i] <= e && e < x[i+1]){ m = i; break; }
    }
    /*** when both points are positive, linear interpolation in log-log space */
    if( (y[m] > 0.0) && (y[m+1] > 0.0) ){
      z = (log(y[m+1]) - log(y[m])) / (log(x[m+1]) - log(x[m])) * (log(e) - log(x[m])) + log(y[m]);
      z = exp(z);
    }
    /*** when one point is zero, linear interpolation */
    else{
      z = (y[m+1] - y[m]) / (x[m+1] - x[m]) * (e - x[m]) + y[m];
    }
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
