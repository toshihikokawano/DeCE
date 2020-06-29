/******************************************************************************/
/**     DeCE Proccessing: Pointwise Cross Section at Union Energy Grid       **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>

using namespace std;

#include "dece.h"
#include "gfr.h"
#include "terminate.h"


static int  DeceReconstructResonance (ENDFDict *, ENDF **, double **);
static int  DeceCopyHighEnergyCrossSection (ENDFDict *, ENDF **, int, double **);
static void DeceCheckNegativeCrossSection (const int, double **);

static const int ncx = 4;
static int mtr[ncx] = {1, 2, 102, 18}; // MT numbers for reconstructed cross sections

/**********************************************************/
/*      Generate Pointwise Cross Section                  */
/*      --------                                          */
/*      Cross sections of total, elastic, capture, and    */
/*      fission are calculated on a common energy grid    */
/**********************************************************/
void DeceGeneratePointwise(ENDFDict *dict, ENDF *lib[])
{
  double **xdat;
  xdat = new double * [ncx];
  for(int j=0 ; j<ncx ; j++){
    xdat[j] = new double [MAX_DBLDATA_LARGE];
  }

  /*** total should always exist */
  int id = dict->getID(3,1);
  if(id < 0){
    message << "MF3 MT1 should exit for processing";
    WarningMessage();
    return;
  }

  /*** excluded channels' MT number negative */
  for(int j=1 ; j<ncx ; j++){
    if(dict->getID(3,mtr[j]) < 0) mtr[j] *= -1;
  }

  /*** when resonance range exists */
  int np = 0;
  if(dict->emaxRe > 0.0){
    np = DeceReconstructResonance(dict,lib,xdat);
    np = DeceCopyHighEnergyCrossSection(dict,lib,np,xdat);
    DeceCheckNegativeCrossSection(np,xdat);

    /*** do not use FILE 2 anymore */
    dict->setLRP(2);
    DeceDelete(dict,2,151);
  }
  /*** otherwise use the total cross section as the common grid */
  else{
    DeceCopyHighEnergyCrossSection(dict,lib,np,xdat);
  }

  /*** create TAB1 record, and replace MT1 etc */
  int nr = 1;
  Record head(dict->getZA(),dict->getAWR(),0,0,0,0);
  Record cont(0.0, 0.0, 0, 0, nr, np);
  int idat[2] = {np,2};

  id = dict->getID(3,mtr[0]);
  lib[id]->memresize(L);
  lib[id]->setENDFhead(head);
  lib[id]->resetPOS();

  ENDFPackTAB1(cont,idat,xdat[0],lib[id]);
  message << "MF3 MT" << mtr[0] << " replaced by resonance contribution + background";
  Notice("DeceGeneratePointwise");
//ENDFWrite(lib[id]);

  for(int j=1 ; j<ncx ; j++){
    if(mtr[j] > 0){
      id = dict->getID(3,mtr[j]);
      lib[id]->memresize(L);
      lib[id]->setENDFhead(head);
      lib[id]->resetPOS();
      ENDFPackTAB1(cont,idat,xdat[j],lib[id]);
      message << "MF3 MT" << mtr[j] << " replaced by resonance contribution + background";
      Notice("DeceGeneratePointwise");
//    ENDFWrite(lib[id]);
    }
  }

  for(int j=0 ; j<ncx ; j++){
    delete [] xdat[j];
  }
  delete [] xdat;
}


/**********************************************************/
/*      Reconstructing Pointwise Cross Section in RRR     */
/**********************************************************/
int DeceReconstructResonance(ENDFDict *dict, ENDF *lib[], double **xdat)
{
  System sys;
  Pcross crs;
  double *edat;

  edat = new double [MAX_DBLDATA_LARGE/2];

  int kres = dict->getID(2,151);

  gfrReadHEADData(&sys,lib[kres]);

  /*** resolved resonance cross sections and background cross sections */
  int np1 = gfrAutoEnergyRRR(&sys,lib[kres],edat,dict->emaxRR);

  for(int i=0 ; i<np1 ; i++){
    int i2 = i*2;
    crs = gfrCrossSection(1,edat[i],&sys,lib[kres]);

    for(int j=0 ; j<ncx ; j++) xdat[j][i2] = edat[i];

    xdat[0][i2+1] = crs.total;
    xdat[1][i2+1] = crs.elastic;
    xdat[2][i2+1] = crs.capture;
    if(dict->isFission()) xdat[3][i2+1] = crs.fission;
  }

  /*** unresolved resonance cross sections */
  bool lssf = false;
  if(dict->emaxRe == dict->emaxUR) lssf = true;

  int np2 = 0;
  if(lssf){
    np2 = gfrAutoEnergyURR(edat,dict->emaxRR,dict->emaxUR);

    for(int i=np1 ; i<np1+np2 ; i++){
      int i2 = i*2;

      crs = gfrCrossSection(2,edat[i-np1],&sys,lib[kres]);

      for(int j=0 ; j<ncx ; j++) xdat[j][i2] = edat[i-np1];

      xdat[0][i2+1] = crs.total;
      xdat[1][i2+1] = crs.elastic;
      xdat[2][i2+1] = crs.capture;
      if(dict->isFission()) xdat[3][i2+1] = crs.fission;
    }
  }

  int np = np1 + np2;

  /*** add background */
  for(int i=0 ; i<np ; i++){
    int i2 = i*2;
    double e = xdat[0][i2];
    for(int j=0 ; j<ncx ; j++){
      if(mtr[j] > 0){
        xdat[j][i2+1] += ENDFInterpolation(lib[dict->getID(3,mtr[j])],e,true,0);
      }
    }
  }

  delete [] edat;

  return(np);
}


/**********************************************************/
/*      Add Smooth Part in MT3                            */
/**********************************************************/
int DeceCopyHighEnergyCrossSection(ENDFDict *dict, ENDF *lib[], int np0, double **xdat)
{
  int k = 0;
  double x0 = 0.0;

  if(np0 > 0){
    /*** duplicate point at resonance boundary */
    k  = 2 * (np0-1);
    x0 = xdat[0][k];

    for(int j=0 ; j<ncx ; j++){
      if(mtr[j] > 0){
        xdat[j][k+2  ] = x0;
        xdat[j][k+2+1] = ENDFInterpolation(lib[dict->getID(3,mtr[j])],x0,false,0);
      }
    }
    k += 4;
  }

  /*** use energy grid of total as the common grid */
  int tid = dict->getID(3,mtr[0]);
  int i=0;
  int nr = lib[tid]->rdata[0].n1;
  for(int ir=0 ; ir<nr ; ir++){
    int np = lib[tid]->idata[2*ir];

    for(int ip=i ; ip<np ; ip++){
      double x = lib[tid]->xdata[2*ip  ];
      double y = lib[tid]->xdata[2*ip+1];

      if(x > x0){
        xdat[0][k++] = x;
        xdat[0][k++] = y;
      }
    }
    i = lib[tid]->idata[2*ir]-1;
  }

  /*** get other reactions at the common energy grid */
  int np1 = k/2;
  for(int i=np0 ; i<np1 ; i++){
    int i2 = i*2;
    double x1 = xdat[0][i2];

    for(int j=1 ; j<ncx ; j++){
      if(mtr[j] > 0){
        xdat[j][i2  ] = x1;
        xdat[j][i2+1] = ENDFInterpolation(lib[dict->getID(3,mtr[j])],x1,false,0);
      }
    }
  }

  return(np1);
}


/**********************************************************/
/*      Check If Negative Cross Section Happened          */
/**********************************************************/
void DeceCheckNegativeCrossSection(const int np, double **xdat)
{
  double x, y;

  for(int i=0 ; i<np ; i++){
    int i2 = i*2;

    for(int j=0 ; j<ncx ; j++){
      if(mtr[j] > 0){
        x = xdat[j][i2  ];
        y = xdat[j][i2+1];

        if(y < 0.0){
          message << "negative section " << y << " detected at " << x << " in MT " << mtr[j];
          WarningMessage();
        }
      }
    }
  }
}


