/******************************************************************************/
/**     DeCE TABLE for MF5                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"


/**********************************************************/
/*      Process MF=5                                      */
/**********************************************************/
void DeceTableMF5(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;

  cout << "# Energy spectrum" << endl;
  cout << "#           NK" << setw(14) << nk << "  number of subsections" << endl;

  for(int n=0 ; n<nk ; n++){
    Record cont = lib->rdata[idx++];
    int lf = cont.l2;

    if(lf == 1){
      int ne = lib->rdata[idx++].n2;
      for(int i=0 ; i<ne ; i++){
        double e  = lib->rdata[idx].c2;
        int    nf = lib->rdata[idx].n2;

        cout << "#            E"; outVal(e); cout << endl;
        cout << "#           NF" << setw(14) << nf << endl;
        for(int j=0 ; j<nf ; j++){
          outVal(lib->xptr[idx][2*j  ]);
          outVal(lib->xptr[idx][2*j+1]);
          cout << endl;
        }
        idx++;

        cout << endl;
        cout << endl;
      }
    }
    else if(lf == 12){
      int ne = lib->rdata[idx++].n2;
      for(int i=0 ; i<ne ; i++){
        double e  = lib->rdata[idx].c2;
        int    nf = lib->rdata[idx].n2;

        cout << "#            E"; outVal(e); cout << endl;
        cout << "#           NF" << setw(14) << nf << endl;
        for(int j=0 ; j<nf ; j++){
          outVal(lib->xptr[idx][2*j  ]);
          outVal(lib->xptr[idx][2*j+1]);
          cout << endl;
        }
        idx++;

        cout << endl;
        cout << endl;
      }
    }
    else WarningMessage("table command cannot process MF = 5, LF = ",lf);
  }
}

