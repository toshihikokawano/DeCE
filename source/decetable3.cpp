/******************************************************************************/
/**     DeCE TABLE for MF3                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"


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

  ENDFPrint1Dim(lib,0);
}


