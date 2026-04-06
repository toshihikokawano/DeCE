/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Print Energy Spectra in MF8 Decay Data                  **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "endflib.h"

int main(int, char *[]);
void tabulateSpectra(ENDF *);

/**********************************************************/
/*      Decay Particle Energy Spectra                     */
/**********************************************************/
int main(int argc, char *argv[])
{
  ifstream fpin;   // file pointer to input library
  ENDF     lib;    // allocate object

  if(argc <= 1){
    cerr << "ENDF file not supplied" << endl;
    exit(-1);
  }

  /*** read in MF8 MT457 */
  fpin.open(argv[1]);
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;
    exit(-1);
  }
  int c = ENDFReadMF8(&fpin,&lib,457);
  fpin.close();
  if(c < 0){
    cerr << "MF8/MT457 not given" << endl;
    exit(-1);
  }

  tabulateSpectra(&lib);
  
  return(0);
}


void tabulateSpectra(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nsp  = head.n2;   // number of radiation types

  int idx = 0;
  double t12   = lib->rdata[idx].c1; // half-life
  idx++;
  idx++;

  for(int n=0 ; n<nsp ; n++){
    double styp = lib->rdata[idx].c2; // radiation type (0:g, 1:b, 2:EC,b+, 4:alpha, 5:n, 6:SF, 7:p, 8:e
    int    lcon = lib->rdata[idx].l1; // 0: no continuum, 1: only continuum, 2: both discrete and continuum
    int    ner  = lib->rdata[idx].n2; // total number of discrete energies

    double fd   = lib->xptr[idx][0];
    double erav = lib->xptr[idx][2];
    double fc   = lib->xptr[idx][4];
    idx++;

    if( (lcon == 0) || (lcon == 2) ){
      for(int k=0 ; k<ner ; k++){
        double er   = lib->rdata[idx].c1;
        double rtyp = lib->xptr[idx][ 0];
        double type = lib->xptr[idx][ 1];
        double ri   = lib->xptr[idx][ 2];
        double ris  = lib->xptr[idx][ 4];
        double ricc = 0.0, rick = 0.0, ricl = 0.0;
        if(styp == 0.0){
          ricc = lib->xptr[idx][ 6];
          rick = lib->xptr[idx][ 8];
          ricl = lib->xptr[idx][10];
        }
        idx++;
      }
    }

    if( (lcon == 1) || (lcon == 2) ){
      int lcov = lib->rdata[idx].l2;
      double rtyp = lib->rdata[idx].c1;

      int nr = lib->rdata[idx].n1;
      int i = 0;
      for(int ir=0 ; ir<nr ; ir++){
        cout.setf(ios::scientific, ios::floatfield);
        for(int ip=i ; ip<lib->iptr[idx][2*ir] ; ip++){
          cout << 

          outVal(lib->xptr[idx][2*ip  ]);
          outVal(lib->xptr[idx][2*ip+1]);
          cout << endl;
        }
        i = lib->idata[2*ir]-1;
        if(ir < nr-1){
          cout << endl;
          cout << endl;
        }
      }
      cout << endl;
      cout << endl;

      idx++;
    }
  }
}
