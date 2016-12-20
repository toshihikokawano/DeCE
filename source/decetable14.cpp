/******************************************************************************/
/**     DeCE TABLE for MF14                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"
#include "constant.h"


/**********************************************************/
/*      Process MF=14                                     */
/**********************************************************/
void DeceTableMF14(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    li   = head.l1;   // 0: non-isotropic, 1: isotropic
  int    ltt  = head.l2;   // 1: Legendre, 2: tabulated
  int    nk   = head.n1;   // number of discrete photons including the continuum
  int    ni   = head.n2;   // number of isotropic photon ang. dist.
  int    idx  = 0;

  cout << "# Photon angular distribution" << endl;
  cout << "#           LI" << setw(14) << li  << "  0:non-isotropic, 1:isotropic" << endl;
  cout << "#          LTT" << setw(14) << ltt << "  1:Legendre, 2:tabulated" << endl;

  if(li == 0){
    cout << "#           NK" << setw(14) << nk << "  number of discrete photons including continuum" << endl;
    cout << "#           NI" << setw(14) << ni << "  number of isotropic photon angular distributions" << endl;
    if(ltt == 1){

      cout << "# energy        level from" << endl;
      for(int i=0 ; i<ni ; i++){
        outVal(lib->rdata[idx].c1);
        outVal(lib->rdata[idx].c2);
        cout << "  isotropic" << endl;
        idx ++;
      }
      cout << endl;

      for(int i=0 ; i<nk-ni ; i++){
        double eg = lib->rdata[idx].c1;
        double sg = lib->rdata[idx].c2;
        int ne = lib->rdata[idx].n2;
        idx ++;

        cout << "# energy        level from" << endl;
        outVal(eg);
        outVal(sg);
        cout << endl;

        for(int j=0 ; j<ne ; j++){
          double e  = lib->rdata[idx].c2;
          int    nl = lib->rdata[idx].n1;

          cout << "#           NL" << setw(14) << nl << endl;
          cout << "# Energy       Leg. Order    Coefficient" << endl;
          outVal(e); outVal(0); outVal(1.0); cout << endl;
          for(int l=0 ; l<nl ; l++){
            outVal(e);
            outVal(l+1);
            outVal(lib->xptr[idx][l]);
            cout << endl;
          }
          idx ++;
        }
        cout << endl;
      }
    }

    else if(ltt == 2){
      cout << "# energy        level from" << endl;
      for(int i=0 ; i<ni ; i++){
        outVal(lib->rdata[idx].c1);
        outVal(lib->rdata[idx].c2);
        cout << "  isotropic" << endl;
        idx ++;
      }
      cout << endl;

      for(int i=0 ; i<nk-ni ; i++){
        double eg = lib->rdata[idx].c1;
        double sg = lib->rdata[idx].c2;
        int ne = lib->rdata[idx].n2;
        idx ++;

        cout << "# energy        level from" << endl;
        outVal(eg);
        outVal(sg);
        cout << endl;

        for(int j=0 ; j<ne ; j++){
          double e  = lib->rdata[idx].c2;
          int    np = lib->rdata[idx].n2;

          cout << "#            E"; outVal(e); cout << endl;
          cout << "#           NP" << setw(14) << np << endl;
          cout << "# energy        angle         probability" << endl;
          for(int l=0 ; l<np ; l++){
            outVal(e);
            outVal(acos(lib->xptr[idx][2*l])/PI * 180.0);
            outVal(lib->xptr[idx][2*l+1]);
            cout << endl;
          }
          idx ++;
        }
        cout << endl;
      }
    }
  }
  else{
    cout << "# isotropic angular distribution" << endl;
  }
}
