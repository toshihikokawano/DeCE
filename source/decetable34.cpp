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

  cout << "# Angular Distribution Covanriance Matrix" << endl;
  cout << "#          LTT" << setw(14) << ltt << "  1: data are given as Legendre, 2: Legendre starting with a0, 3: L or L1=0" << endl;
  cout << "#         NMT1" << setw(14) << nmt1 << "  number of MT subsections" <<  endl;
  cout << endl;

  for(int m=0 ; m<nmt1 ; m++){
    Record cont = lib->rdata[idx];
    int mat1 = cont.l1;
    int mt1  = cont.l2;
    int nl   = cont.n1;
    int nl1  = cont.n2;
    idx ++;

    cout << "#   Subsection" << setw(14) << m << endl;
    cout << "#         MAT1" << setw(14) << mat1 << "  other MAT number" << endl;
    cout << "#          MT1" << setw(14) << mt1 << "   other reaction type" << endl;
    cout << "#           NL" << setw(14) << nl << "  number of Legendre coefficient for reaction MT" << endl;
    cout << "#          NL1" << setw(14) << nl1 << "  number of Legendre coefficient for reaction MT1" << endl;

    int l   = lib->rdata[idx].l1;
    int l1  = lib->rdata[idx].l2;
    int lct = lib->rdata[idx].n1;
    idx ++;

    cout << "#            L" << setw(14) << l  << "  Legendre coefficient index for MT" << endl;
    cout << "#           L1" << setw(14) << l1 << "  Legendre coefficient index for MT1" << endl;
    cout << "#          LCT" << setw(14) << lct << "  0:cordinate same as MF4, 1:LAB, 2:CMS" << endl;

    int ls  = lib->rdata[idx].l1;
    int lb  = lib->rdata[idx].l2;
    int ne  = lib->rdata[idx].n2;

    cout << "#           LS" << setw(14) << ls << "  0:asymmetric, 1:symmetric when LB = 5" << endl;
    cout << "#           LB" << setw(14) << lb << "  covariance type" << endl;
    cout << "#           NE" << setw(14) << ne << "  number of energies" << endl;

    double *eptr, *cptr;
    eptr = lib->xptr[idx];
    cptr = &lib->xptr[idx][ne];

    cout << "# Coefficient  Uncertainty" << endl;

    for(int i=0; i<ne-1; i++){
      int ki = i+i*(ne-1)-i*(i+1)/2;

      outVal(eptr[i]);
      outVal(sqrt(cptr[ki]));
/*
      for(int j=0; j<=i; j++){
        int kj = j+j*(ne-1)-j*(j+1)/2;
        int k  = (j < i) ?  i+j*(ne-1)-j*(j+1)/2 : j+i*(ne-1)-i*(i+1)/2;
        int c  = 0;
        if(cptr[ki]*cptr[kj] != 0.0) c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
        if((i == j) && (c == 999)) c = 1000;
        cout << setw(5) << c;
      }
*/
      cout << endl;
    }
    outVal(eptr[ne-1]);
    outVal(0.0);
    cout << endl;
  }
  idx++;
}
