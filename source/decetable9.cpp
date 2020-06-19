/******************************************************************************/
/**     DeCE TABLE for MF9                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"


/**********************************************************/
/*      Process MF=9                                      */
/**********************************************************/
void DeceTableMF9(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    lis  = head.l1;   // level number indicator
  int    ns   = head.n1;   // number of subsections

  cout << "# Multiplicities for production of radioactive nuclides" << endl;
  cout << "#           NS" << setw(14) << ns  << "  number of final states" << endl;
  cout << "#          LIS" << setw(14) << lis << "  target level" << endl;
  cout << endl;

  for(int n=0 ; n<ns ; n++){
    Record cont = lib->rdata[n];
    double qm   = cont.c1;
    double qi   = cont.c2;
    int    lfs  = head.l2;   // level number indicator

    cout << "#   FinalState" << setw(14) << n << endl;
    cout << "#           QM"; outVal(qm); cout << "  mass difference Q-value" << endl;
    cout << "#           QI"; outVal(qi); cout << "  reaction Q-value" << endl;
    cout << "#          LFS" << setw(14) << lfs << "  product level" << endl;

    ENDFPrint1Dim(lib,n,"Energy","Multiplicity");
  }
}
