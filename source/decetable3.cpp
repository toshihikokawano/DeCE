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

  cout << "# Cross section" << endl;
  cout << "#           QM"; outVal(qm); cout << "  mass difference Q-value" << endl;
  cout << "#           QI"; outVal(qi); cout << "  reaction Q-value" << endl;

  if((opt.WriteXdataConversion != 1.0) || (opt.WriteYdataConversion != 1.0)){
    DataSize size = lib->getSIZE();
    ENDF tmp(size);
    ENDFLibCopy(lib,&tmp);
    for(int i=0 ; i<tmp.getNX()/2; i++){
      tmp.xdata[2*i  ] *= opt.WriteXdataConversion;
      tmp.xdata[2*i+1] *= opt.WriteYdataConversion;
    }
    ENDFPrint1Dim(&tmp,0);
  }
  else ENDFPrint1Dim(lib,0);
}


