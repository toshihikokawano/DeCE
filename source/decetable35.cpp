/******************************************************************************/
/**     DeCE TABLE for MF35                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"


/**********************************************************/
/*      Process MF=35                                     */
/**********************************************************/
void DeceTableMF35(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int nk  = head.n1;
  int idx = 0;

  cout << nk << endl;

  for(int n=0 ; n<nk ; n++){
    Record cont = lib->rdata[idx];

    double e1 = cont.c1;
    double e2 = cont.c2;
//  int    ls = cont.l1; // should be 1
    int    lb = cont.l2; // should be 7
    int    ne = cont.n2;

    cout << "#           LB" << setw(14) << lb << "  absolute covariance" << endl;
    cout << "#           E1" << setw(14) << e1 << "  lowest incident energy for this subsection" << endl;
    cout << "#           E2" << setw(14) << e2 << "  higheset incident energy" << endl;
    cout << "#           NE" << setw(14) << ne << "  number of outgoing energies" << endl;

    double *eptr, *cptr;
    eptr = lib->xptr[idx];
    cptr = &lib->xptr[idx][ne];

    for(int i=0; i<ne-1; i++){
      int ki = i+i*(ne-1)-i*(i+1)/2;

      outVal(eptr[i]);
      outVal(sqrt(cptr[ki]));

      for(int j=0; j<=i; j++){
        int kj = j+j*(ne-1)-j*(j+1)/2;
        int k  = (j < i) ?  i+j*(ne-1)-j*(j+1)/2 : j+i*(ne-1)-i*(i+1)/2;
        int c  = 0;
        if(cptr[ki]*cptr[kj] != 0.0){
          c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
          if((i == j) && (c == 999)) c = 1000;
        }
        cout << setw(5) << c;
      }
      cout << endl;
    }
    outVal(eptr[ne-1]);
    outVal(0.0);
    cout << endl;

    idx++;

    cout << endl;
    cout << endl;
  }
}
