/******************************************************************************/
/**     DeCE TABLE for MF12                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"


/**********************************************************/
/*      Process MF=12                                     */
/**********************************************************/
void DeceTableMF12(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    l0   = head.l1;   // 1: multiplicity, 2: transition probability

  cout << "# Photon production multiplicities and transition probability arrays" << endl;
  cout << "#           L0" << setw(14) << l0 << "  1: multiplicity, 2: transition probability" << endl;

  if(l0 == 1){
    int    nk   = head.n1;   // number of subsections
    cout << "#           NK" << setw(14) << nk << "  number of subsections" << endl;
    for(int n=0 ; n<nk ; n++){
      Record cont = lib->rdata[n];
      double eg   = cont.c1;
      double es   = cont.c2;
      int    lp   = cont.l1;
      int    lf   = cont.l1;

      cout << "#           LP" << setw(14) << lp << "  0: origin unknown, 1: nonprimary photon, 2: primary photons" << endl;
      cout << "#           LF" << setw(14) << lf << "  1: tabulated in File 15, 2: discrete photons" << endl;
      cout << "#           EG"; outVal(eg); cout << "  photon energy" << endl;
      cout << "#           ES"; outVal(es); cout << "  level energy from which the photon originates" << endl;

      ENDFPrint1Dim(lib,n);
    }
  }
  else if(l0 == 2){
    int    lg   = head.l2;   // 1: all gamma, 2: complex case
    int    ns   = head.n1;   // number of levels

    cout << "#           NS" << setw(14) << ns << "  number of levels below the current level" << endl;

    Record cont = lib->rdata[0];
    double es   = cont.c1;
    int    lp   = cont.l1;
    int    nt   = cont.n2;

    cout << "#           LG" << setw(14) << lg << "  1: simple case, 2: complex case" << endl;
    cout << "#           LP" << setw(14) << lp << "  0: origin unknown, 1: nonprimary photon, 2: primary photons" << endl;
    cout << "#           NT" << setw(14) << nt << "  number of transitions" << endl;
    cout << "#           ES"; outVal(es); cout << "  level energy" << endl;


    cout.setf(ios::scientific, ios::floatfield);
    if(lg == 1){
        for(int i=0 ; i<nt ; i++){
        cout << setprecision(6) << setw(14) << lib->xptr[0][2*i  ];
        cout << setprecision(6) << setw(14) << lib->xptr[0][2*i+1];
        cout << endl;
      }
    }
    else if(lg == 2){
        for(int i=0 ; i<nt ; i++){
        cout << setprecision(6) << setw(14) << lib->xptr[0][3*i  ];
        cout << setprecision(6) << setw(14) << lib->xptr[0][3*i+1];
        cout << setprecision(6) << setw(14) << lib->xptr[0][3*i+2];
        cout << endl;
      }
    }
  }
}
