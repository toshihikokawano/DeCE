/******************************************************************************/
/**     DeCE TABLE for MF33                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"


/**********************************************************/
/*      Process MF=33                                     */
/**********************************************************/
void DeceTableMF33(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nl   = head.n2;
  int    idx  = 0;

  cout << "# NL:" << setw(3) << nl << endl;

  for(int n=0 ; n<nl ; n++){
    Record cont = lib->rdata[idx++];
    int nc = cont.n1;
    int ni = cont.n2;

    cout << "# NC:" << setw(3) << nc;
    cout << "  NI:" << setw(3) << ni << endl;

    /*** NC type */
    for(int i=0 ; i<nc ; i++){
      int    lty = lib->rdata[idx++].l2;

      double e1  = lib->rdata[idx].c1;
      double e2  = lib->rdata[idx].c2;
      int    nci = lib->rdata[idx].n2;
      if(lty == 0){
        cout << "# Emin: " << e1 << "  Emax: " << e2 << "  NCI: " << nci << endl;
        cout << "# MT            Const"<< endl;
        for(int j=0; j<nci ; j++){
          outVal((int)lib->xptr[idx][2*j+1]);
          outVal((int)lib->xptr[idx][2*j  ]);
          cout << endl;
        }
      }
      idx++;
    }

    /*** NI type */
    for(int i=0 ; i<ni ; i++){
      int lt = lib->rdata[idx].l1;
      int lb = lib->rdata[idx].l2;
      int nt = lib->rdata[idx].n1;
      int np = lib->rdata[idx].n2;

      cout << "# LT: " << lt << " LB: " << lb << " NT: " << nt << " NP: " << np << endl;
      if(lb == 1 || lb == 2){
        for(int i=0; i<np-1; i++){
          outVal(         lib->xptr[idx][2*i  ]  );
          outVal(sqrt(abs(lib->xptr[idx][2*i+1])));
          cout << endl;
        }
        outVal(lib->xptr[idx][2*(np-1)]);
        outVal(0.0);
        cout << endl;
      }

      else if(lb == 5){
        double *eptr, *cptr;
        eptr = lib->xptr[idx];
        cptr = &lib->xptr[idx][np];

        for(int i=0; i<np-1; i++){
          int ki = i+i*(np-1)-i*(i+1)/2;

          outVal(eptr[i]);
          outVal(sqrt(cptr[ki]));

          for(int j=0; j<=i; j++){
            int kj = j+j*(np-1)-j*(j+1)/2;
            int k  = (j<i) ?  i+j*(np-1)-j*(j+1)/2 : j+i*(np-1)-i*(i+1)/2;
            int c  = 0;
            if(cptr[ki]*cptr[kj] != 0.0){
              c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
            }
            cout << setw(5) << c;
          }
          cout << endl;
        }
        outVal(eptr[np-1]);
        outVal(0.0);
        cout << endl;
      }
      idx++;
    }
  }
}
