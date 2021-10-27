/******************************************************************************/
/**     DeCE TABLE for MF3                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "global.h"


/**********************************************************/
/*      Process MF=3                                      */
/**********************************************************/
void DeceTableMF3(ENDF *lib)
{
  Record cont = lib->rdata[0];
  double qm   = cont.c1;
  double qi   = cont.c2;
  int    lr   = cont.l2;

  cout << "# Cross section" << endl;
  cout << "#           QM"; outVal(qm); cout << "  mass difference Q-value" << endl;
  cout << "#           QI"; outVal(qi); cout << "  reaction Q-value" << endl;

  if(lr > 0) cout << "#           LR" << setw(14) << lr << "  inelastic, then decay by this reaction type" << endl;

  if((opt.WriteXdataConversion != 1.0) || (opt.WriteYdataConversion != 1.0)){
    ENDF tmp;
    ENDFLibCopy(lib,&tmp);
    for(int i=0 ; i<tmp.getNX()/2; i++){
      tmp.xdata[2*i  ] *= opt.WriteXdataConversion;
      tmp.xdata[2*i+1] *= opt.WriteYdataConversion;
    }
    ENDFPrint1Dim(&tmp,0,"Energy","CrossSection");
  }
  else ENDFPrint1Dim(lib,0,"Energy","CrossSection");
}


