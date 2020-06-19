/******************************************************************************/
/**     DeCE TABLE for MF13                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"


/**********************************************************/
/*      Process MF=13                                     */
/**********************************************************/
void DeceTableMF13(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nk   = head.n1;   // number of subsections

  cout << "# Photon production cross section" << endl;
  cout << "#           NK" << setw(14) << nk  << "  number of discrete photons (incl. continuum)" << endl;
  cout << endl;

  int n0 = 0;
  if(nk != 1){
    cout << "#             " << "  total photon production cross section" << endl;
    ENDFPrint1Dim(lib,n0++,"Energy","CrossSection");
  }

  for(int n=n0 ; n<nk+n0 ; n++){
    Record cont = lib->rdata[n];
    double eg   = cont.c1;
    double es   = cont.c2;
    int    lp   = cont.l1;
    int    lf   = cont.l1;

    cout << "#   DiscPhoton" << setw(14) << n << endl;
    cout << "#           LP" << setw(14) << lp << "  0: origin unknown, 1: nonprimary photon, 2: primary photons" << endl;
    cout << "#           LF" << setw(14) << lf << "  1: tabulated in File 15, 2: discrete photons" << endl;
    cout << "#           EG"; outVal(eg); cout << "  photon energy (LP=0,1) or binding energy (LP=2)" << endl;
    cout << "#           ES"; outVal(es); cout << "  level energy from which the photon originates" << endl;

    ENDFPrint1Dim(lib,n,"Energy","CrossSection");
  }
}
