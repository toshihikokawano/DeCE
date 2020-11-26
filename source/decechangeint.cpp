/******************************************************************************/
/**     DeCE CHANGEINT                                                       **/
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Change Interpolation Law in MF3                   */
/**********************************************************/
void DeceChangeInt(ENDFDict *dict, ENDF *lib[], const int mt, int range, int point, int intlaw)
{
  int k = dict->getID(3,mt);

  if(k < 0) TerminateCode("MT number not found",mt);

  Record r  = lib[k]->rdata[0];
  int    nr = r.n1;
  int    np = r.n2;

  if(point < 0 || point >= np) point = np;

  if(range <= 0) TerminateCode("data range incorrect",range);
  if(range > nr+1) TerminateCode("data range incorrect",range);

  if(range == nr+1) lib[k]->rdata[0].n1 = nr+1;

  lib[k]->idata[(range-1)*2  ] = point;
  lib[k]->idata[(range-1)*2+1] = intlaw;

  message << "interporation range " << range << " is up to " << point << " and INT = " << intlaw;
  Notice("DeceChangeInt");
}


