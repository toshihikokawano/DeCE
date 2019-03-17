/******************************************************************************/
/**     DeCE HEADER FIXER                                                    **/
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"
#include "masstable.h"
#include "constant.h"

/**********************************************************/
/*      Fix AWR in DICT                                   */
/**********************************************************/
void DeceFixAWR(ENDFDict *dict)
{

  /*** ZA and AWR from Dictionary */
  double za   = dict->getZA();    //  1000*Z + A number

  int z = za/1000;
  int a = za - z*1000;

  double mass = mass_excess(z,a);
  double awr  = (mass / AMUNIT + a) / MNEUTRON;

  dict->setAWR(awr);
}


