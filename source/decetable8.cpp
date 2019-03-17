/******************************************************************************/
/**     DeCE TABLE for MF8                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"

static void DeceTableMF8Other(ENDF *);
static void DeceTableMF8MT454(ENDF *);
static void DeceTableMF8MT457(ENDF *);


/**********************************************************/
/*      Process MF=8                                      */
/**********************************************************/
void DeceTableMF8(ENDF *lib)
{
  int mt = lib->getENDFmt();

  cout << "# Radioactive decay and fission production yield data" << endl;
  if(mt == 454)      cout << "# independent yield" << endl;
  else if(mt == 459) cout << "# cumulative yield" << endl;
  else if(mt == 457) cout << "# radioactive decay data" << endl;
  else               cout << "# radioactive nuclide production" << endl;

  switch(mt){
  case 454:
  case 459: DeceTableMF8MT454(lib); break;
  case 457: DeceTableMF8MT457(lib); break;
  default:  DeceTableMF8Other(lib); break;
  }
}


/**********************************************************/
/*      Fission Product Yields                            */
/**********************************************************/
void DeceTableMF8MT454(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    le1  = head.l1; // energy dependence flag

  cout << "#         LE+1" << setw(14) << le1 << "  number of neutron energies" << endl;

  for(int ie=0 ; ie<le1 ; ie++){
    double e    = lib->rdata[ie].c1;
    int    i    = lib->rdata[ie].l1;
    int    nfp  = lib->rdata[ie].n2;
    cout << "#            E" << setw(14) << e   << "  neutron energy" << endl;
    cout << "#          NFP" << setw(14) << nfp << "  number of FPs" << endl;
    if(ie == 0){
      cout << "#           LE" << setw(14) << i   << "  energy-dependent flag" << endl;
    }else{
      cout << "#            I" << setw(14) << i   << "  interpolation scheme" << endl;
    }

    cout << "#    Z      A   Level         FPY           delta-FPY   " << endl;

    for(int n=0 ; n<nfp ; n++){
        int j = 4*n;
        double za = lib->xptr[ie][j];
        int    z = (int)(za/1000.0);
        int    a = (int)(za - 1000*z);
        
        cout << setw(7) << z << setw(7) << a;
        outVal(lib->xptr[ie][j+1]);
        outVal(lib->xptr[ie][j+2]);
        outVal(lib->xptr[ie][j+3]); 
        cout << endl;
    }
    cout << endl;
    cout << endl;
  }
}



/**********************************************************/
/*      Radioactive Decay Data                            */
/**********************************************************/
void DeceTableMF8MT457(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    lis  = head.l1;   // level number indicator
  int    liso = head.l2;   // level number indicator
  int    nsp  = head.n2;   // number of radiation types

  cout << "#          LIS" << setw(14) << lis  << "  target level" << endl;
  cout << "#         LISO" << setw(14) << liso << "  isomeric state number" << endl;
  cout << "#          NSP" << setw(14) << nsp  << "  total number of STYP for which spectra are given" << endl;

  int idx = 0;
  double t12  = lib->rdata[idx].c1;
  double dt12 = lib->rdata[idx].c2;

  outVal(t12); outVal(dt12); cout << "  half-life [sec]" << endl;
  outVal(lib->xptr[idx][0]);  outVal(lib->xptr[idx][1]); cout << "   beta energy [eV]" << endl;
  outVal(lib->xptr[idx][2]);  outVal(lib->xptr[idx][3]); cout << "  gamma energy [eV]" << endl;
  outVal(lib->xptr[idx][4]);  outVal(lib->xptr[idx][5]); cout << "  alpha energy [eV]" << endl;
  idx++;

  double spi = lib->rdata[idx].c1;
  double par = lib->rdata[idx].c2;
  int    nkd = lib->rdata[idx].n2;

  cout << "#          SPI"; outVal(spi); cout << "  spin of the nuclide in its LIS state" << endl;
  cout << "#          PAR"; outVal(par); cout << "  parity of the nuclide" << endl;
  cout << "#          NKD" << setw(14) << nkd << "  total number of decay modes" << endl;

  cout << "# RTYP          RFS           Q             BR"<< endl;
  for(int n=0 ; n<nkd ; n++){
    int j = 6*n;
    outVal(lib->xptr[idx][j  ]); outVal(lib->xptr[idx][j+1]);
    outVal(lib->xptr[idx][j+2]); outVal(lib->xptr[idx][j+4]); cout << endl;
    cout <<"                            ";
    outVal(lib->xptr[idx][j+3]); outVal(lib->xptr[idx][j+5]); cout << endl;
  }
  idx++;

  for(int n=0 ; n<nsp ; n++){
    double styp = lib->rdata[idx].c2;
    int    lcon = lib->rdata[idx].l1;
    int    ner  = lib->rdata[idx].n2;

    cout << "#         STYP" << setw(14) << (int)styp << "  decay radiation type" << endl;
    cout << "#         LCON" << setw(14) << lcon << "  0: no continuum, 1: only continuum, 2: both discrete and continuum" << endl;
    cout << "#          NER" << setw(14) << ner << "  total number of discrete energies" << endl;


    cout << "# FD            ER            FC"<< endl;
    outVal(lib->xptr[idx][0]); outVal(lib->xptr[idx][2]); outVal(lib->xptr[idx][4]); cout << endl;
    outVal(lib->xptr[idx][1]); outVal(lib->xptr[idx][3]); outVal(lib->xptr[idx][5]); cout << endl;
    idx++;

    if( (lcon == 0) || (lcon == 2) ){
      for(int k=0 ; k<ner ; k++){
        outVal(lib->rdata[idx].c1);  outVal(lib->rdata[idx].c2); cout << "   discrete energy [eV]" << endl;

        cout << "# RTYP          TYPE          RI            RIS";
        if(styp == 0.0) cout << "          RICC          RICK          RICL";
        cout << endl;

        outVal(lib->xptr[idx][ 0]); outVal(lib->xptr[idx][ 1]);
        outVal(lib->xptr[idx][ 2]); outVal(lib->xptr[idx][ 4]);
        if(styp == 0.0){
          outVal(lib->xptr[idx][ 6]); outVal(lib->xptr[idx][ 8]); outVal(lib->xptr[idx][10]);
        }
        cout << endl;

        cout <<"                            ";
        outVal(lib->xptr[idx][ 3]); outVal(lib->xptr[idx][ 5]);
        if(styp == 0.0){
          outVal(lib->xptr[idx][ 7]); outVal(lib->xptr[idx][ 9]); outVal(lib->xptr[idx][11]);
        }
        cout << endl;
        cout << endl;
        idx++;
      }
    }

    if( (lcon == 1) || (lcon == 2) ){
      int lcov = lib->rdata[idx].l2;
      cout << "#         RTYP"; outVal(lib->rdata[idx].c1); cout << "  decay mode" << endl;
      cout << "#         LCOV" << setw(14) << lcov << "  covariance data flag" << endl;
      ENDFPrint1Dim(lib,idx);
      idx++;
    }
  }
}



/**********************************************************/
/*      All Other MT SubSections                          */
/**********************************************************/
void DeceTableMF8Other(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    lis  = head.l1;   // level number indicator
  int    ns   = head.n1;   // number of subsections
  int    no   = head.n2;   // 0: complete decay chain, 1: given in MT457 in MATP

  cout << "#           NS" << setw(14) << ns  << "  number of subsections" << endl;
  cout << "#           NO" << setw(14) << no  << "  0: complete decay chain given, 1: in MT457 of MATP" << endl;
  cout << "#          LIS" << setw(14) << lis << "  target level" << endl;

  for(int n=0 ; n<ns ; n++){
    Record cont = lib->rdata[n];
    double zap  = cont.c1;
    double elfs = cont.c2;
    int    lmf  = cont.l1;   // file number
    int    lfs  = cont.l2;   // level number indicator
    int    nd   = cont.n1 / 6;
    int    matp = cont.n2;

    cout << "#          ZAP"; outVal(zap); cout << "  ZA produced" << endl;
    cout << "#         ELFS"; outVal(elfs); cout << "  excitation energy" << endl;
    cout << "#          LMF" << setw(14) << lmf << "  file number for multiplicity or cross section" << endl;
    cout << "#          LFS" << setw(14) << lfs << "  product level" << endl;
    cout << "#         MATP" << setw(14) << matp << "  material number of product" << endl;

    if(no == 0){
      cout << "#           ND" << setw(14) << nd << "  number of branches" << endl;
      cout << "# Half Life     Decay(RTYP)   Next ZA       Branching     End-Point E   CT" << endl;

      for(int i=0 ; i<nd ; i++){
        int j = 6*i;
        outVal(lib->xdata[j  ]); outVal(lib->xdata[j+1]); outVal(lib->xdata[j+2]);
        outVal(lib->xdata[j+3]); outVal(lib->xdata[j+4]); outVal(lib->xdata[j+5]);
        cout << endl;
      }
    }
  }
}
