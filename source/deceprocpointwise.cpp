/******************************************************************************/
/**     DeCE Proccessing Point                                               **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>

using namespace std;

#include "dece.h"
#include "gfr.h"


static void DeceReconstructResonance      (ENDFDict *, ENDF **, int *);
static void DeceAddHighEnergyCrossSection (ENDFDict *, ENDF **, int *);
static void DeceCheckNegativeCrossSection (ENDFDict *, ENDF **, int *);

static const int n_cross_sec = 4;

/**********************************************************/
/*      Generate Pointwise Cross Section                  */
/**********************************************************/
void DeceGeneratePointwise(ENDFDict *dict, ENDF *lib[])
{
  int kid[n_cross_sec];

  /*** 901 - 904 are temporary MTs */
  for(int i=0 ; i<n_cross_sec ; i++){
    kid[i] = dict->getID(3,901 + i);
  }

  DeceReconstructResonance(dict,lib,kid);
  DeceAddHighEnergyCrossSection(dict,lib,kid);
  DeceCheckNegativeCrossSection(dict,lib,kid);

  DeceDelete(dict,3,1);
  if(dict->getID(3,  2) >= 0) DeceDelete(dict,3,  2);
  if(dict->getID(3,102) >= 0) DeceDelete(dict,3,102);
  if(dict->getID(3, 18) >= 0) DeceDelete(dict,3, 18);
  else                        DeceDelete(dict,3,904);

  dict->head.l1 = 2;

//ENDFWrite(lib[kid[2]]);
}


/**********************************************************/
/*      Reconstructing Pointwise Cross Section in RRR     */
/**********************************************************/
void DeceReconstructResonance(ENDFDict *dict, ENDF *lib[], int *kid)
{
  System sys;
  Pcross crs;
  double *elab;

  elab = new double [MAX_DBLDATA_LARGE/2];

  int kres = dict->getID(2,151);

  gfrReadHEADData(&sys,lib[kres]);

  /*** resolved resonance cross sections and background cross sections */
  int np1 = gfrAutoEnergyRRR(&sys,lib[kres],elab,dict->emaxRR);
  for(int i=0 ; i<np1 ; i++){
    int i2 = i*2;
    crs = gfrCrossSection(1,elab[i],&sys,lib[kres]);

    for(int k=0 ; k<n_cross_sec; k++) lib[kid[k]]->xdata[i2  ] = elab[i];

    lib[kid[0]]->xdata[i2+1] = crs.total;
    lib[kid[1]]->xdata[i2+1] = crs.elastic;
    lib[kid[2]]->xdata[i2+1] = crs.capture;
    if(dict->isFission()) lib[kid[3]]->xdata[i2+1] = crs.fission;

    if(dict->getID(3,  1) >= 0)
      lib[kid[0]]->xdata[i2+1] += ENDFInterpolation(lib[dict->getID(3,  1)],elab[i],true,0);

    if(dict->getID(3,  2) >= 0)
      lib[kid[1]]->xdata[i2+1] += ENDFInterpolation(lib[dict->getID(3,  2)],elab[i],true,0);

    if(dict->getID(3,102) >= 0)
      lib[kid[2]]->xdata[i2+1] += ENDFInterpolation(lib[dict->getID(3,102)],elab[i],true,0);

    if( dict->isFission() && (dict->getID(3, 18) >= 0) )
      lib[kid[3]]->xdata[i2+1] += ENDFInterpolation(lib[dict->getID(3, 18)],elab[i],true,0);
  }

  /*** unresolved resonance cross sections */
  int np2 = gfrAutoEnergyURR(elab,dict->emaxRR,dict->emaxUR);

  for(int i=np1 ; i<np1+np2 ; i++){
    int i2 = i*2;

    crs = gfrCrossSection(2,elab[i-np1],&sys,lib[kres]);

    for(int k=0 ; k<n_cross_sec; k++) lib[kid[k]]->xdata[i2  ] = elab[i-np1];

    lib[kid[0]]->xdata[i2+1] = crs.total;
    lib[kid[1]]->xdata[i2+1] = crs.elastic;
    lib[kid[2]]->xdata[i2+1] = crs.capture;
    if(dict->isFission()) lib[kid[3]]->xdata[i2+1] = crs.fission;
  }
  int np = np1 + np2;
  
  Record head(dict->head.c1,dict->head.c2,0,0,0,0);

  lib[kid[0]]->setENDFhead(head);
  lib[kid[0]]->rdata[0].setRecord(0.0, 0.0, 0, 0, 1, np);
  lib[kid[0]]->idata[0] = np;
  lib[kid[0]]->idata[1] = 2;

  lib[kid[1]]->setENDFhead(head);
  lib[kid[1]]->rdata[0].setRecord(0.0, 0.0, 0, 0, 1, np);
  lib[kid[1]]->idata[0] = np;
  lib[kid[1]]->idata[1] = 2;

  lib[kid[2]]->setENDFhead(head);
  lib[kid[2]]->rdata[0].setRecord(0.0, 0.0, 0, 0, 1, np);
  lib[kid[2]]->idata[0] = np;
  lib[kid[2]]->idata[1] = 2;

  if(dict->isFission()){
    lib[kid[3]]->setENDFhead(head);
    lib[kid[3]]->rdata[0].setRecord(0.0, 0.0, 0, 0, 1, np);
    lib[kid[3]]->idata[0] = np;
    lib[kid[3]]->idata[1] = 2;
  }

  delete [] elab;
}


/**********************************************************/
/*      Add Smooth Part in MT3                            */
/**********************************************************/
void DeceAddHighEnergyCrossSection(ENDFDict *dict, ENDF *lib[], int *kid)
{
  int tid = dict->getID(3,1);
  int np0 = lib[kid[0]]->rdata[0].n2;

  /*** duplicate point at resonance boundary */
  int k = 2 * (np0-1);
  double x = lib[kid[0]]->xdata[k];
  double x0 = x;

  k += 2;
  double y = ENDFInterpolation(lib[tid],x,false,0);

  lib[kid[0]]->xdata[k++] = x;
  lib[kid[0]]->xdata[k++] = y;


  int i=0;
  int nr = lib[tid]->rdata[0].n1;
  for(int ir=0 ; ir<nr ; ir++){
    int np = lib[tid]->idata[2*ir  ];

    for(int ip=i ; ip<np ; ip++){
      x = lib[tid]->xdata[2*ip  ];
      y = lib[tid]->xdata[2*ip+1];

      if(x > x0){
        lib[kid[0]]->xdata[k++] = x;
        lib[kid[0]]->xdata[k++] = y;
      }
    }
    i = lib[tid]->idata[2*ir]-1;
  }

  int np1 = k/2;
  lib[kid[0]]->rdata[0].n2 = np1;
  lib[kid[0]]->idata[0]    = np1;

  if(dict->getID(3,2) >= 0){
    for(int i=np0 ; i<np1 ; i++){
      int i2 = i*2;
      double e = lib[kid[0]]->xdata[i2];
      lib[kid[1]]->xdata[i2  ] = e;
      lib[kid[1]]->xdata[i2+1] = ENDFInterpolation(lib[dict->getID(3,2)],e,false,0);
    }
    lib[kid[1]]->rdata[0].n2 = np1;
    lib[kid[1]]->idata[0]    = np1;
  }

  if(dict->getID(3,102) >= 0){
    for(int i=np0 ; i<np1 ; i++){
      int i2 = i*2;
      double e = lib[kid[0]]->xdata[i2];
      lib[kid[2]]->xdata[i2  ] = e;
      lib[kid[2]]->xdata[i2+1] = ENDFInterpolation(lib[dict->getID(3,102)],e,false,0);
    }
    lib[kid[2]]->rdata[0].n2 = np1;
    lib[kid[2]]->idata[0]    = np1;
  }

  if(dict->isFission() && (dict->getID(3,18) >= 0)){
    for(int i=np0 ; i<np1 ; i++){
      int i2 = i*2;
      double e = lib[kid[0]]->xdata[i2];
      lib[kid[3]]->xdata[i2  ] = e;
      lib[kid[3]]->xdata[i2+1] = ENDFInterpolation(lib[dict->getID(3,18)],e,false,0);
    }
    lib[kid[3]]->rdata[0].n2 = np1;
    lib[kid[3]]->idata[0]    = np1;
  }
}


/**********************************************************/
/*      Check If Negative Cross Section Happened          */
/**********************************************************/
void DeceCheckNegativeCrossSection(ENDFDict *dict, ENDF *lib[], int *kid)
{
  int np = lib[kid[0]]->rdata[0].n2;
  double x, y;

  for(int i=0 ; i<np ; i++){
    int i2 = i*2;
    x = lib[kid[0]]->xdata[i2];

    y = lib[kid[0]]->xdata[i2+1];
    if(y < 0.0) cerr << "negative total cross section " << y << " detected at " << x << endl;

    if(dict->getID(3,2) >= 0){
      y = lib[kid[1]]->xdata[i2+1];
      if(y < 0.0) cerr << "negative elastic cross section " << y << " detected at " << x << endl;
    }

    if(dict->getID(3,102) >= 0){
      y = lib[kid[2]]->xdata[i2+1];
      if(y < 0.0) cerr << "negative capture cross section " << y << " detected at " << x << endl;
    }

    if(dict->getID(3,18) >= 0){
      y = lib[kid[3]]->xdata[i2+1];
      if(y < 0.0) cerr << "negative fission cross section " << y << " detected at " << x << endl;
    }
  }
}


