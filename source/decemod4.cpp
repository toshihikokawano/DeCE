/******************************************************************************/
/**     DeCE MANIPULATE MF4                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decemisc.h"
#include "gfr.h"
#include "terminate.h"

static const int defaultenergypoints = 100;

static void smoothedlegendre (ENDFDict *, ENDF *[], const int, const int, double *, double **);
static void updateMF4 (const int, const int, const int, int **, double **, Record *, ENDF *);


/**********************************************************/
/*      Scattering Angular Distribution from Resonances   */
/**********************************************************/
void DeceResonanceAngularDistribution(ENDFDict *dict, ENDF *lib[], int np)
{
  const int mf = 4;
  const int mt = 2;
  int k4 = dict->getID(mf,mt);
  if(k4 < 0){
    message << "MF/MT = " << mf << "/" << mt << " not found";
    TerminateCode("DeceResonanceAngulardistribution");
  }

  int kres = dict->getID(2,151);
  if(kres < 0){
    message << "resonance parameters not found";
    TerminateCode("DeceResonanceAngulardistribution");
  }

  /*** energy grid for smoothed angular distributions */
  if(np == 0) np = defaultenergypoints + 1;
  else np ++; // add boundary point
  double de = dict->emaxRR / (np - 1);

  /*** L-max of Legendre coefficients */
  int psize = 2 * gfrLMax(dict,lib) + 1;

  /*** allocate data */
  double *xres = new double [np];
  double **pres = new double * [np];
  for(int i=0 ; i<np ; i++){
    xres[i] = de * i;
    pres[i] = new double [psize];
    for(int l=0 ; l<psize ; l++) pres[i][l] = 0.0;
  }

  /*** smoothed Legendre coefficients */
  smoothedlegendre(dict,lib,np,psize,xres,pres);


  Record head = lib[k4]->getENDFhead();
  int    ltt  = head.l2;   // 0: isotropic, 1: Legendre, 2: tabulated
  int    idx  = 0;
  int    li   = lib[k4]->rdata[idx++].l1;   // 0: non-isotropic, 1: isotropic

  if(li != 0){
    message << "MF/MT = " << mf << "/" << mt << " is isotropic";
    TerminateCode("DeceResonanceAngulardistribution");
  }

  int ne = lib[k4]->rdata[idx++].n2;

  /*** total number of energy points */
  int npmax = 0;  // total number of energy points in the new file
  int ne0   = 0;  // number of energy points inside RRR given in the original file
  for(int i=0 ; i<ne ; i++){
    if(lib[k4]->rdata[idx + i].c2 >= xres[np-1]){
      ne0 = i;
      npmax = np + ne - ne0;
      break;
    }
  }

  /*** allocate double-data block for formatting */
  double **xdat = new double * [npmax];
  int **idat = new int * [npmax];
  Record *cont  = new Record [npmax];

  /*** when all the data are Legendre coefficients */
  if(ltt == 1){
    /*** highest L for memory allocation */
    int nlmax = 0;
    for(int i=0 ; i<ne ; i++){
      int nl = lib[k4]->rdata[idx + i].n1;
      if(nl > nlmax) nlmax = nl;
    }

    /*** Npmax energy points x Nlmax L-values */
    for(int i=0 ; i<npmax ; i++){
      xdat[i] = new double [nlmax];
      idat[i] = new int [2]; // not used
      for(int j=0 ; j<nlmax ; j++) xdat[i][j] = 0.0;
    }

    /*** copy energies and Legendre coefficients */
    int n = 0;
    for(int i=0 ; i<np ; i++){
      cont[i].setRecord(0.0,xres[i],0,0,psize-1,0);
      for(int l=0 ; l<psize-1 ; l++) xdat[i][l] = pres[i][l+1];
    }
    n += np;

    /*** copy the rest */
    for(int i=ne0 ; i<ne ; i++){
      int nl = lib[k4]->rdata[idx + i].n1;
      cont[n].setRecord(0.0,lib[k4]->rdata[idx + i].c2,0,0,nl,0);

      for(int j=0 ; j<nl ; j++){
        xdat[n][j] = lib[k4]->xptr[idx + i][j];
      }
      n ++;
    }
  }

  /*** when tabulated data */
  else if(ltt == 2){
    /*** look for the max number of angles */
    int namax = 0;
    for(int i=0 ; i<ne ; i++){
      int na = lib[k4]->rdata[idx + i].n2;
      if(na > namax) namax = na;
    }

    /*** Npmax energy points x Namax angles x 2 (angle,prob) sets */
    for(int i=0 ; i<npmax ; i++){
      xdat[i] = new double [2 * namax];
      idat[i] = new int [2];
      for(int j=0 ; j<2*namax ; j++) xdat[i][j] = 0.0;
    }

    /*** prepare angle, cos(t), from backward to front */
    double step = 2.0 / (namax - 1.0);
    double *tdat = new double [namax];

    for(int j=0 ; j<namax ; j++){
      tdat[j] = -1.0 + step * j;
      if(fabs(tdat[j]) < 1e-10) tdat[j] = 0.0;
    }

    /*** copy energies and tabulated angular distributions */
    int n = 0;
    for(int i=0 ; i<np ; i++){
      for(int j=0 ; j<namax ; j++){
        double t = acos(tdat[j]) / M_PI * 180.0;
        double f = 0.5;
        for(int l=1 ; l<psize ; l++){
          double p = (pres[i][0] == 0.0) ? 0.0 : pres[i][l] / pres[i][0];
          f += (l+0.5) * p * legendre(l,t);
        }
        xdat[n][2*j  ] = tdat[j];
        xdat[n][2*j+1] = f;
      }
      cont[n].setRecord(0.0,xres[i],0,0,1,namax);
      idat[n][0] = namax;
      idat[n][1] = 2; // linear interpolation between angles
      n ++;
    }

    /*** copy the rest */
    for(int i=ne0 ; i<ne ; i++){
      cont[n] = lib[k4]->rdata[idx + i];
      for(int j=0 ; j<2*cont[n].n2 ; j++){
        xdat[n][j] = lib[k4]->xptr[idx + i][j];
      }
      idat[n][0] = lib[k4]->iptr[idx + i][0];
      idat[n][1] = lib[k4]->iptr[idx + i][1];

      n ++;
    }

    delete [] tdat;
  }

  /*** replace MF4/MT2 */
  updateMF4(mf,mt,npmax,idat,xdat,cont,lib[k4]);

  message << "angular distributions in MF4/MT2 below ";
  message << setw(13) << setprecision(6) << xres[np-1];
  message << " replaced by resonance-reconstructed data, ";
  message << setw(4) << np-1 << " energy points";
  Notice("DeceResonanceAngularDistribution");

  for(int i=0 ; i<npmax ; i++){
    delete [] xdat[i];
    delete [] idat[i];
  }
  delete [] xdat;
  delete [] idat;
  delete [] cont;

  for(int i=0 ; i<np ; i++) delete [] pres[i];
  delete [] pres;
  delete [] xres;
}


/**********************************************************/
/*      Calculate Smoothed Legendre Coefficients          */
/**********************************************************/
void smoothedlegendre(ENDFDict *dict, ENDF *lib[], const int np, const int psize, double *eave, double **pave)
{
  /*** allocate coefficients memory */
  double *edat = new double [MAX_DBLDATA/2];
  double **pdat = new double * [MAX_DBLDATA/2];
  for(int i=0 ; i<MAX_DBLDATA/2 ; i++){
    pdat[i] = new double [psize];
    edat[i] = 0.0;
    for(int l=0 ; l<psize ; l++) pdat[i][l] = 0.0;
  }
  
  /*** calculate Legendre coefficients */
  int ne = gfrAngDistSimple(dict,lib,psize,edat,pdat);
  // for(int n=0 ; n<ne ; n++){
  //   cout << edat[n];
  //   for(int l=0 ; l<psize ; l++) cout << " " << pdat[n][l];
  //   cout << endl;
  // }

  double width = dict->emaxRR / (np - 1); // Gaussian averaging width = Delta E
  double y[psize];

  /*** at each energy point, calculate Gaussian average */
  for(int i=1 ; i<np ; i++){
    double e0 = eave[i] - 2 * width;
    double e1 = eave[i] + 2 * width;

    double x = 0.0;
    for(int l=0 ; l<psize ; l++) y[l] = 0.0;

    for(int j=0 ; j<ne-1 ; j++){
      if( (edat[j] < e0) || (e1 < edat[j]) ) continue;

      double d = eave[i] - edat[j];
      double w = exp(-d*d / (width*width)) * (edat[j+1] - edat[j]);

      x += w;
      for(int l=0 ; l<psize ; l++) y[l] += w * pdat[j][l];
    }
    for(int l=0 ; l<psize ; l++) pave[i][l] = (x > 0.0) ? y[l] / x : 0.0;
  }

  /*** isotropic at zero energy */
  eave[0] = edat[0];
  pave[0][0] = 1.0;
  for(int l=1 ; l<psize ; l++) pave[0][l] = 0.0;

  /*** normalize Legendre coefficients */
  for(int n=0 ; n<np ; n++){
    for(int l=1 ; l<psize ; l++) if(pave[n][0] != 0.0) pave[n][l] = pave[n][l]/pave[n][0];
    pave[n][0] = 1.0;
  }

  delete [] edat;
  for(int i=0 ; i<MAX_DBLDATA/2 ; i++) delete [] pdat[i];
  delete [] pdat;
/*
  for(int n=0 ; n<np ; n++){
    cout << eave[n];
    for(int l=1 ; l<psize ; l++){
      if(pave[n][0] == 0.0) cout << " " << 0.0;
      else cout << " " << pave[n][l]/pave[n][0];
    }
    cout << endl;
  }
*/
}


/**********************************************************/
/*      Store Legendre Coefficients in ENDF lib (MF4)     */
/**********************************************************/
void updateMF4(const int mf, const int mt, const int ne, int **itab, double **xtab, Record *xcont, ENDF *lib4)
{
  ENDF   lib;
  Record cont;
  int    idat[2];

  /*** Make HEAD and CONT */
  int    lvt  = 0; // transformation matrix not given
  int    ltt  = lib4->getENDFhead().l2; // Legendre parameters given
  int    li   = 0; // not isotropic
  int    lct  = 2; // center-of-mass system

  Record head = lib4->getENDFhead();
  lib.setENDFmat( lib4->getENDFmat() );
  lib.setENDFmf(mf);
  lib.setENDFmt(mt);

  lib.setENDFhead(head.c1, head.c2, lvt, ltt, 0, 0);

  /*** extra card for MF4 */
  cont.setRecord(0.0, head.c2, li, lct, 0, 0);
  ENDFPackCONT(cont,&lib);

  /*** Make TAB2 (LIST) */
  cont.setRecord(0.0, 0.0, 0, 0, 1, ne);
  idat[0] = ne;     // there are NE incident energies
  idat[1] = 2;      // lin-lin interpolation

  if( (ltt == 1) || (ltt == 3) ) ENDFPackTAB2( cont, xcont, idat, xtab, &lib);
  else ENDFPackTAB21(cont, xcont, idat, itab, xtab, &lib);

  if(ltt == 3){
    /*** copy TAB21 from source */
    int idx = lib4->rdata[1].n2 + 2; // skip first 2 subsections
    ENDFPackCopyTAB21(lib4,&lib,idx);
  }

//ENDFWrite(&lib);
  ENDFLibCopy(&lib,lib4);
}

