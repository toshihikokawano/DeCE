/******************************************************************************/
/**     DeCE READ FILE                                                       **/
/******************************************************************************/

#include <iostream>
#include <ostream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "deceread.h"
#include "global.h"
#include "terminate.h"

static int    geneCSdata1 (int, double *, double *, double, double, double *);
static int    geneCSdata2 (int, double *, double *, double, double *);
static double loginterpol   (int, double, double *, double *, int *);

static bool   charged_particle_file = false;


/**********************************************************/
/*      Set Charged Particle File                         */
/**********************************************************/
void readSetCharged()
{
  charged_particle_file = true;
}


/**********************************************************/
/*      Read in Cross Section Data                        */
/**********************************************************/
int readCSdata(char *file, int ofset, const int mt, double *x, double *y)
{
  ifstream fp;
  string   line;

  /*** default CoH3 output file structure in CoHParticleProduction.dat */
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

  if(ofset == 0){
    message << "MF3:MT" << mt << " cannot read from " << file << " because ofset is not given";
    Notice("DeceRead:readCSdata");
    return 0;
  }

  fp.open(file);
  if(!fp){ message << "cannot open data file " << file; TerminateCode("readCSdata"); }

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
    if(nc >= MAX_DBLDATA){ message << "too many energy points, " << nc; TerminateCode("readCSdata"); }
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

  if(ofset == 0){
    if(      ( 51 <= mt) && (mt <=  91) ) ofset = mt  - 50 + 1;
    else if( (600 <= mt) && (mt <= 640) ) ofset = mt - 600 + 1;
    else if( (650 <= mt) && (mt <= 690) ) ofset = mt - 650 + 1;
    else if( (700 <= mt) && (mt <= 740) ) ofset = mt - 700 + 1;
    else if( (750 <= mt) && (mt <= 790) ) ofset = mt - 750 + 1;
    else if( (800 <= mt) && (mt <= 840) ) ofset = mt - 800 + 1;
    else if( (mt == 649) || (mt == 699) || (mt == 749) || (mt == 799) || (mt == 849) ) ofset = 42;
  }

  if(ofset == 0){
    message << "MF3:MT" << mt << " cannot read from " << file << " because ofset is not given";
    Notice("DeceRead:readISdata");
    return 0;
  }

  fp.open(file);
  if(!fp){ message << "cannot open data file " << file; TerminateCode("readISdata"); }

  getline(fp,line);
  istringstream s1(&line[1]);  // skip comment #
  double eth = 0.0;
  for(int i=0 ; i<ofset ; i++) s1 >> eth;

  if(eth == 0.0 && (mt != 600) && (mt != 650) && (mt != 700) && (mt !=750) && (mt != 800)){
    message << "MF3:MT" << mt << " from " << file << " has zero excitation energy, cannot read";
    Notice("DeceRead:readISdata");
    return 0;
  }

  *elev = eth * opt.ReadXdataConversion;

  int nc = 0;
  while(getline(fp,line)){
    if(line.length() == 0) continue;

    istringstream s2(line);
    s2 >> x[nc];
    for(int i=0 ; i<ofset ; i++) s2 >> y[nc];

    if(x[nc] == 0.0) continue;
    x[nc] *= opt.ReadXdataConversion;
    y[nc] *= opt.ReadYdataConversion;
    if(DeceCheckReadRange(x[nc])) continue;

    if( (mt >= 600) || (y[nc] > 0.0) ) nc++;
    if(nc >= MAX_DBLDATA){ message << "too many energy points, " << nc; TerminateCode("readISdata"); }
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
  if(!fp){ message << "cannot open data file " << file; TerminateCode("readNUdata"); }

  if(ofset == 0) ofset = 1;

  int nc = 0;
  while(getline(fp,line)){
    if(line[0] == '#') continue;
    if(line.length() == 0) continue;

    istringstream ss(line);
    ss >> x[nc];
    for(int i=0 ; i<ofset ; i++) ss >> y[nc];

    if(x[nc] == 0.0) continue;
    x[nc] *= opt.ReadXdataConversion;
    if(DeceCheckReadRange(x[nc])) continue;

    nc++;
    if(nc >= MAX_DBLDATA){ message << "too many energy points, " << nc; TerminateCode("readNUdata"); }
  }
  fp.close();

  if(nc >= 1){
    message << "MF3:MT455(6) " << nc << " points from (" << x[0] << "," << y[0] << ") to (" << x[nc-1] << "," << y[nc-1] << ") imported from " << file;
    Notice("DeceRead:readNUdata");
  }

  return nc;
}


/**********************************************************/
/*      Read in Radioactive Nuclide Data                  */
/**********************************************************/
struct Prod readMShead(char *file, const int ofset)
{
  ifstream fp;
  string   line;

  fp.open(file);
  if(!fp){ message << "cannot open data file " << file; TerminateCode("readMShead"); }

  getline(fp,line); // Nuclide
  istringstream s0(&line[13]);  // skip heading "# Nuclide Ex "
  int za = 0;
  double ex[3];
  for(int i=0 ; i<ofset ; i++){
    s0 >> za;
    s0 >> ex[1];
    s0 >> ex[2];
  }
  ex[0] = 0.0;
  ex[1] *= opt.ReadXdataConversion;
  ex[2] *= opt.ReadXdataConversion;

  getline(fp,line); // Particle
  istringstream s1(&line[13]);  // skip heading "# Particle L "
  int pid = 0, nx[3];
  for(int i=0 ; i<ofset ; i++){
    s1 >> pid;
    s1 >> nx[1];
    s1 >> nx[2];
  }
  nx[0] = 0;

  getline(fp,line); // HalfLife
  istringstream s2(&line[13]);  // skip heading "# HalfLife   "
  double th[3];
  for(int i=0 ; i<ofset ; i++){
    s2 >> th[0];
    s2 >> th[1];
    s2 >> th[2];
  }

  struct Prod prd;
  prd.za = za;
  for(int i=0 ; i<3 ; i++){
    prd.nx[i] = nx[i];
    prd.ex[i] = ex[i];
    prd.th[i] = th[i];
  }

  return prd;
}


/**********************************************************/
/*      Read in Radioactive Production Data               */
/**********************************************************/
int readMSdata(char *file, const int mf, const int ofset, double *x, double *y)
{
  ifstream fp;
  string   line;

  fp.open(file);
  if(!fp){ message << "cannot open data file " << file; TerminateCode("readMShead"); }

  int nc = 0;
  while(getline(fp,line)){
    if(line[0] == '#') continue;
    if(line.length() == 0) continue;

    istringstream ss(line);
    ss >> x[nc];
    for(int i=0 ; i<ofset ; i++) ss >> y[nc];

    if(x[nc] == 0.0) continue;
    x[nc] *= opt.ReadXdataConversion;
    y[nc] *= opt.ReadYdataConversion;
    if(DeceCheckReadRange(x[nc])) continue;

    nc++;
    if(nc >= MAX_DBLDATA){ message << "too many energy points, " << nc; TerminateCode("readMSdata"); }
  }
  fp.close();

  if(nc >= 1){
    message << "MF:" << mf << " " << nc << " points from (" << x[0] << "," << y[0] << ") to (" << x[nc-1] << "," << y[nc-1] << ") imported from " << file;
    Notice("DeceRead:readMSdata");
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
  /*** ignore resonance region and replace all the data points by input */
  else if(eth < 0.0){
    i = geneCSdata2(n,x,y,0.0,xdat);
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
int mergeCSdata(int n, double *x, double *y, double eres, double *xdat, const int m, double *xbak)
{
  if(n == 0) return 0;

  int i = 0;

  /*** copy old data up to Eres */
  for(int j=0 ; j<m ; j++){
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


