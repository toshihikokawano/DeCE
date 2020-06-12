/******************************************************************************/
/**     DeCE TABLE for MF5                                                   **/
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
    Record cont = lib->rdata[idx];
    int lf = cont.l2;

    cout << "#            K" << setw(14) << n << "  subsection" << endl;
    int nek = lib->rdata[idx].n2;
    cout << "#           NE" << setw(14) << nek << endl;
    for(int i=0 ; i<nek ; i++){
      outVal(lib->xptr[idx][2*i  ]);
      outVal(lib->xptr[idx][2*i+1]);
      cout << endl;
    }
    idx ++;

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
    else if((lf == 5) || (lf == 11)){
      int ne1 = lib->rdata[idx].n2;
      cout << "#           NE" << setw(14) << ne1 << endl;
      for(int i=0 ; i<ne1 ; i++){
        outVal(lib->xptr[idx][2*i  ]);
        outVal(lib->xptr[idx][2*i+1]);
        cout << endl;
      }
      idx++;

      cout << endl;
      cout << endl;

      int ne2 = lib->rdata[idx].n2;
      cout << "#           NE" << setw(14) << ne2 << endl;
      for(int i=0 ; i<ne2 ; i++){
        outVal(lib->xptr[idx][2*i  ]);
        outVal(lib->xptr[idx][2*i+1]);
        cout << endl;
      }
      idx++;

      cout << endl;
      cout << endl;
    }
    else if((lf == 7) || (lf == 9) || (lf == 12)){
      int ne = lib->rdata[idx].n2;
      cout << "#           NE" << setw(14) << ne << endl;
      if(lf == 12){
        cout << "#          EFL"; outVal(lib->rdata[idx].c1); cout << endl;
        cout << "#          EFH"; outVal(lib->rdata[idx].c2); cout << endl;
      }
      for(int i=0 ; i<ne ; i++){
        outVal(lib->xptr[idx][2*i  ]);
        outVal(lib->xptr[idx][2*i+1]);
        cout << endl;
      }
      idx++;

      cout << endl;
      cout << endl;
    }
    else{
      message << "table command cannot process MF = 5, LF = " << lf;
      WarningMessage();
    }
  }
}

