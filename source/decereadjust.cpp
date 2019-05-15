/******************************************************************************/
/**     DeCE READJUST                                                        **/
/******************************************************************************/

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Readjust Individual Subsection by Summed Section  */
/**********************************************************/
void DeceReadjust(ENDFDict *dict, ENDF *lib[], const int mt, const int mtmp)
{
  int k0 = dict->getID(3,mt);
  int k1 = dict->getID(3,mtmp);

  if(k0 < 0) TerminateCode("MT number not found",mt);

  /*** MTs for subsections */
  int mt0 = 0, mt1 = 0;
  switch(mt){
  case   4:  mt0 =  51; mt1 =  91; break;
  case 103:  mt0 = 600; mt1 = 649; break;
  case 104:  mt0 = 650; mt1 = 699; break;
  case 105:  mt0 = 700; mt1 = 749; break;
  case 106:  mt0 = 750; mt1 = 799; break;
  case 107:  mt0 = 800; mt1 = 849; break;
  default: break;
  }
  if((mt0 == 0) && (mt1 == 0)) TerminateCode("MT number not processed by READJUST",mt);

  /*** sum subsections from MT0 to MT1 and store in Mtmp = 99 */
  DeceCalc(dict,lib,mtmp,mt0,mt1,':');

  /*** apply ratio MT1 / MT0 to each subsection (MT2) */
  for(int mt=mt0 ; mt<=mt1 ; mt++){
    int k2 = dict->getID(3,mt);
    if(k2 < 0) continue;

    int np = lib[k2]->rdata[0].n2;
    for(int ip=0 ; ip<np ; ip++){
      double x = lib[k2]->xdata[2*ip  ];
      double y = lib[k2]->xdata[2*ip+1];

      double y0 = ENDFInterpolation(lib[k0],x,false,0);
      double y1 = ENDFInterpolation(lib[k1],x,false,0);

      double r  = (y1 == 0.0) ? 0.0 : y0 / y1;

      lib[k2]->xdata[2*ip+1] = y * r;
    }
  }

  message << "MF3:MT" << mt << " used for renormalizing "; 
  message << "from MF3:MT" << mt0;
  message << " to MF3:MT"  << mt1;
  Notice("DeceReadjust");
}

