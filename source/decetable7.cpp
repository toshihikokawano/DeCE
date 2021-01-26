/******************************************************************************/
/**     DeCE TABLE for MF7                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "global.h"

static void DeceTableMF7MT2c(ENDF *);
static void DeceTableMF7MT2i(ENDF *);
static void DeceTableMF7MT4i(ENDF *);


/**********************************************************/
/*      Process MF=7                                      */
/**********************************************************/
void DeceTableMF7(ENDF *lib)
{
  int mt = lib->getENDFmt();
  Record head = lib->getENDFhead();
  int lthr = head.l1;

  cout << "# Thermal neutron scattering law data" << endl;
  if(mt == 2){
    if(lthr == 1){
      cout << "# coherent elastic scattering" << endl;
      DeceTableMF7MT2c(lib);
    }
    else{
      cout << "# incoherent elastic scattering" << endl;
      DeceTableMF7MT2i(lib);
    }
  }
  else{
    cout << "# incoherent inelastic scattering" << endl;
    DeceTableMF7MT4i(lib);
  }
}


/**********************************************************/
/*      Coherent Elastic Scattering                       */
/**********************************************************/
void DeceTableMF7MT2c(ENDF *lib)
{
  int idx = 0;

  int lt = lib->rdata[idx].l1;
  int nr = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  cout << "#           LT" << setw(14) << lt  << "  temperature dependence flag / number of temp" << endl;
  cout << "#           NP" << setw(14) << np  << "  number of Bragg edges given" << endl;

  double t = lib->rdata[idx].c1;
  cout << "#            T"; outVal(t); cout << "  temperature [K]" << endl;

  ENDF tmp;
  tmp.setENDFmf(lib->getENDFmf());
  tmp.setENDFmt(lib->getENDFmf());
  tmp.setENDFhead(lib->getENDFhead());

  tmp.rdata[0] = lib->rdata[idx];
  tmp.checkDataSize(2*nr,2*np);
  for(int i=0 ; i<2*nr ; i++) tmp.idata[i] = lib->iptr[idx][i];
  for(int i=0 ; i<2*np ; i++) tmp.xdata[i] = lib->xptr[idx][i];

  ENDFPrint1Dim(&tmp,0,"E","S(E,T)");
  idx ++;

  for(int j=1 ; j<=lt ; j++){
    t = lib->rdata[idx].c1;
    cout << "#            T"; outVal(t); cout << "  temperature [K]" << endl;
    for(int i=0 ; i<np ; i++) tmp.xdata[2*i+1] = lib->xptr[idx][i];
    ENDFPrint1Dim(&tmp,0,"E","S(E,T)");
    idx ++;
  }
}


/**********************************************************/
/*      Incoherent Elastic Scattering                     */
/**********************************************************/
void DeceTableMF7MT2i(ENDF *lib)
{
  int idx = 0;

  double sb = lib->rdata[idx].c1;
  int    np = lib->rdata[idx].n2;
  cout << "#           SB"; outVal(sb); cout << "  characteristic bound cross section [b]" << endl;
  cout << "#           NP" << setw(14) << np  << "  number of temperatures" << endl;
  ENDFPrint1Dim(lib,idx,"T","W");
}


/**********************************************************/
/*      Incoherent Inelastic Scattering                   */
/**********************************************************/
void DeceTableMF7MT4i(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int lat   = head.l2;
  int lasym = head.n1;

  cout << "#          LAT" << setw(14) << lat   << "  temperature to be used, 0: actual, 1: const 0.0253 eV" << endl;
  cout << "#        LASYM" << setw(14) << lasym << "  0: symmetric S, 1: asymmetric S" << endl;

  int idx = 0;
  int lln = lib->rdata[idx].l1;
  int ni  = lib->rdata[idx].n1;
  int ns  = lib->rdata[idx].n2;

  int a[3];
  a[0] = a[1] = a[2] = -1;
  if(ns == 1){
    a[0] = lib->xptr[idx][ 6];
  }
  else if(ns == 2){
    a[0] = lib->xptr[idx][ 6];
    a[1] = lib->xptr[idx][12];
  }
  else if(ns == 3){
    a[0] = lib->xptr[idx][ 6];
    a[1] = lib->xptr[idx][12];
    a[2] = lib->xptr[idx][18];
  }

  cout << "#          LLN" << setw(14) << lln   << "  0: S stored, 1: log(S) stored" << endl;
  cout << "#           NS" << setw(14) << ns    << "  number of non-principal scattering atom types" << endl;

  cout << "# N             B(N)" << endl;
  for(int i=0 ; i<ni ; i++){
    cout << setw(14) << i+1;
    outVal(lib->xptr[idx][i]);
    cout << endl;
  }
  cout << endl;
  cout << endl;
  idx ++;

  int nb = lib->rdata[idx].n2;
  cout << "#           NB" << setw(14) << nb    << "  total number of beta values" << endl;
  idx ++;

  ENDF tmp;
  tmp.setENDFmf(lib->getENDFmf());
  tmp.setENDFmt(lib->getENDFmf());
  tmp.setENDFhead(lib->getENDFhead());

  for(int ib=0 ; ib<nb ; ib++){
    double temp = lib->rdata[idx].c1;
    double beta = lib->rdata[idx].c2;
    int lt = lib->rdata[idx].l1;
    int nr = lib->rdata[idx].n1;
    int np = lib->rdata[idx].n2;
    cout << "#           NT" << setw(14) << lt+1   << "  total number of temperatures" << endl;
    cout << "#           NP" << setw(14) << np     << "  number of alpha values" << endl;

    cout << "#         BETA"; outVal(beta);  cout << "  beta" << endl;
    cout << "#         TEMP"; outVal(temp);  cout << "  temperature [K]" << endl;

    tmp.rdata[0] = lib->rdata[idx];
    tmp.checkDataSize(2*nr,2*np);
    for(int i=0 ; i<2*nr ; i++) tmp.idata[i] = lib->iptr[idx][i];
    for(int i=0 ; i<2*np ; i++) tmp.xdata[i] = lib->xptr[idx][i];
    ENDFPrint1Dim(&tmp,0,"alpha","S(a,b,T)");
    idx ++;

    for(int j=1 ; j<=lt ; j++){
      temp = lib->rdata[idx].c1;
      beta = lib->rdata[idx].c2;
      cout << "#         BETA"; outVal(beta);  cout << "  beta" << endl;
      cout << "#         TEMP"; outVal(temp);  cout << "  temperature [K]" << endl;

      for(int i=0 ; i<np ; i++) tmp.xdata[2*i+1] = lib->xptr[idx][i];
      ENDFPrint1Dim(&tmp,0,"alpha","S(a,b,T)");
      idx ++;
    }
  }

  cout << "# Teff" << 0 << endl;
  ENDFPrint1Dim(lib,idx,"T","Teff");

  for(int is=0 ; is<ns ; is++){
    if(a[is] == 0){
      cout << "# Teff" << is+1 << endl;
      ENDFPrint1Dim(lib,idx,"T","Teff");
      idx ++;
    }
  }
}


