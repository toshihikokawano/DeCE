/******************************************************************************/
/**     DeCE TABLE for MF15                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"


/**********************************************************/
/*      Process MF=15                                     */
/**********************************************************/
void DeceTableMF15(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nc   = head.n1;   // number of subsections

  cout << "# Continuous photon energy spectra" << endl;
  cout << "#           NC" << setw(14) << nc  << "  number of partial distributions" << endl;

  int idx = 0;
  for(int n=0 ; n<nc ; n++){
    Record cont = lib->rdata[idx];
    int lf = cont.l2;
    ENDFPrint1Dim(lib,idx);
    idx ++;

    if(lf == 1){
      int ne = lib->rdata[idx++].n2;
      for(int i=0 ; i<ne ; i++){
        double e  = lib->rdata[idx].c2;
        int    np = lib->rdata[idx].n2;

        cout << "#            E"; outVal(e); cout << endl;
        cout << "#           NP" << setw(13) << np << endl;
        for(int j=0 ; j<np ; j++){
          outVal(lib->xptr[idx][2*j  ]);
          outVal(lib->xptr[idx][2*j+1]);
          cout << endl;
        }
        idx++;

        cout << endl;
        cout << endl;
      }
    }
    else{
      message << "table command cannot process MF = 15, LF = " << lf;
      WarningMessage();
    }
  }
}
