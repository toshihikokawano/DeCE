/******************************************************************************/
/**     DeCE APPLYFUNC                                                       **/
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      TAB1 Data Multiplied by a Function                */
/**********************************************************/
void DeceApplyFunc(ENDFDict *dict, ENDF *lib[], const int mf, const int mt,
                   int kind, double p1, double p2, double p3)
{
  if(p3 == 0){ message << "P3 Parameter zero"; TerminateCode("DeceApplyfunc"); }
  if(mf != 3) return;

  int k0 = dict->getID(3,mt);
  if(k0 < 0){ message << "MT number " << mt << " not found"; TerminateCode("DeceApplyfunc"); }

  Record r  = lib[k0]->rdata[0];
  int    np = r.n2;

  for(int ip=0 ; ip<np ; ip++){
    double x = lib[k0]->xdata[2*ip];
    double f = 1.0;

    /*** if range range is set */
    if(DeceCheckEditRange(x)) continue;

    /*** energies in MeV */
    x *= 1e-6;

    /*** Fermi function */
    if(kind == 1)      f = p1/( 1.0+ exp( (x-p2)/p3 ) ) + 1.0;
    /*** Gaussian */
    else if(kind == 2) f = p1*exp( -(x-p2)*(x-p2)/p3 )  + 1.0;
    /*** reversed Fermi function */
    else if(kind == 3) f = p1/( 1.0+ exp(-(x-p2)/p3 ) ) + 1.0;

    lib[k0]->xdata[2*ip+1] *= f;
  }

  message << "MF" << mf << ":MT" << mt << " applied ";
  if(     kind == 1) message << "Fermi function p1/{1+exp(x-p2)/p3} + 1";
  else if(kind == 2) message << "Gauss function p1 exp{-(x-p2)(x-p2)/p3} + 1";
  else if(kind == 3) message << "reversed Fermi function p1/{1+exp(-(x-p2)/p3)} + 1";
  message << " where parameters p1 = " << p1 << " p2 = " << p2 << " p3 =" << p3;
  Notice("DeceApplyFunc");
}

