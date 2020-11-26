/******************************************************************************/
/**     DeCE QVALUE                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Change Q-values in MF3                            */
/**********************************************************/
void DeceChangeQvalue(ENDFDict *dict, ENDF *lib[], const int mt, double qm, double qi)
{
  int k = dict->getID(3,mt);

  if(k < 0) TerminateCode("MT number not found",mt);

  lib[k]->rdata[0].c1 = qm;
  lib[k]->rdata[0].c2 = qi;
}


/**********************************************************/
/*      Check Q-values and Threshold Energies in MF3      */
/**********************************************************/
void DeceCheckThreshold(ENDFDict *dict, ENDF *lib[], bool fix)
{
  double qm, qi, e0, e1, de;
  double tolerance = 1.0; // eV

  cout.setf(ios::fixed, ios::floatfield);
  cerr.setf(ios::fixed, ios::floatfield);

  for(int i=0 ; i<dict->getSEC() ; i++){

    if(dict->mf[i] == 3){
      int k = dict->getID(3,dict->mt[i]);
      if(k < 0) continue;

      Record head = lib[k]->getENDFhead();
      double mass = head.c2;

      qm = lib[k]->rdata[0].c1;
      qi = lib[k]->rdata[0].c2;
      e0 = lib[k]->xptr[0][0];

      e1 = 0.0;
      if((qm == 0.0) && (qi == 0.0)) e1 = 1.0e-5;
      else{
        e1 = -qi * (mass + 1.0) / mass;
        if(e1 <= 0.0) e1 = 1.0e-5;
      }

      if(e0 > 1.0e+7) tolerance = 10.0;
      else tolerance = 1.0;

      de = e0 - e1;
      if(de >= tolerance){
        cout << setw( 5) << dict->mt[i];
        cout << setw(14) << setprecision(2) << e0;
        cout << setw(14) << setprecision(2) << e1;
        cout << setw(14) << setprecision(2) << de << endl;

        if(fix) lib[k]->xptr[0][0] = e1;
      }

      message << "MF" << dict->mf[i] << ":MT" << dict->mt[i];
      message << " given threshold " << setprecision(6) << e0;
      message << "   calc. from mass " << setprecision(6) << e1;
      message << "   diff. " << setprecision(6) << e0 - e1;

      Notice("DeceCheckThreshold");
    }
  }
}

