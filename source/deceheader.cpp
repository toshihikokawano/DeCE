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

  /*** ZA and AWT from Dictionary */
  double za   = dict->head.c1;    //  1000*Z + A number

  int z = za/1000;
  int a = za - z*1000;

  double mass = mass_excess(z,a);
  double awt  = (mass / AMUNIT + a) / MNEUTRON;

  dict->head.c2 = awt;
}


