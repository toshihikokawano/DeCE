/******************************************************************************/
/**     DeCE TABLE for MF34                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"


/**********************************************************/
/*      Process MF=34                                     */
/**********************************************************/
void DeceTableMF34(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    ltt  = head.l2;
  int    nmt1 = head.n2;
  int    idx  = 0;

  cout << "#  LTT:" << setw(3) << ltt << " NMT1:" << setw(3) << nmt1 << endl;

  for(int m=0 ; m<nmt1 ; m++){
    Record cont = lib->rdata[idx];
    int mat1 = cont.l1;
    int mt1  = cont.l2;
    int nl   = cont.n1;
    int nl1  = cont.n2;
    idx ++;

    cout << "# MAT1:" << setw(3) << mat1 << "  MT1:"  << setw(3) << mt1;
    cout << "   NL:" << setw(3) << nl   << "  NL1:"  << setw(3) << nl1 << endl;

    int l   = lib->rdata[idx].l1;
    int l1  = lib->rdata[idx].l2;
    int lct = lib->rdata[idx].n1;
    int ni  = lib->rdata[idx].n2;
    idx ++;

    cout << "#    L:" << setw(3) << l << "   L1:" << setw(3) << l1;
    cout << "  LCT:" << setw(3) << lct << "   NI:" << setw(3) << ni << endl;

    int ls  = lib->rdata[idx].l1;
    int lb  = lib->rdata[idx].l2;
    int nt  = lib->rdata[idx].n1;
    int ne  = lib->rdata[idx].n2;

    cout << "#   LS:" << setw(3) << ls << "   LB:" << setw(3) << lb;
    cout << "   LT:" << setw(3) << nt << "   NE:" << setw(3) << ne << endl;

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
  }
  idx++;
}
