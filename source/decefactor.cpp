/******************************************************************************/
/**     DeCE FACTOR                                                          **/
/******************************************************************************/

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      TAB1 Data Multiplied by Factor                    */
/*      Rescale Y-data if x > 0.0                         */
/**********************************************************/
void DeceFactor(ENDFDict *dict, ENDF *lib[], const int mf, const int mt,
                double x, double y, double xmin, double xmax)
{
  if(mf != 3) return;

  int k0 = dict->getID(3,mt);
  if(k0 < 0) TerminateCode("MT number not found",mt);

  double f  = 0.0;
  /*** all Y data are multiplied by a factor */
  if(x == 0.0){
    f = y;
  }
  /*** renormalize all Y data at given (x,y) */
  else{
    double z = ENDFInterpolation(lib[k0],x,false,0);
    if(z != 0.0) f = y/z;
  }

  Record r  = lib[k0]->rdata[0];
  int    np = r.n2;
  /*** if energy range is given */
  if( (xmin < xmax) && (xmin >= 0.0) && (xmax > 0.0) ){
    for(int ip=0 ; ip<np ; ip++){
      if( (xmin <= lib[k0]->xdata[2*ip]) && (lib[k0]->xdata[2*ip] <= xmax) ){
        lib[k0]->xdata[2*ip+1] *= f;
      }
    }
  }
  /*** manipulate entire energy range */
  else{
    for(int ip=0 ; ip<np ; ip++) lib[k0]->xdata[2*ip+1] *= f;
  }

  ostringstream os;
  os << "MF" << mf << ":MT" << mt << " rescaled by factor " << f;
  if( (xmin < xmax) && (xmin >= 0.0) && (xmax > 0.0) ){
    os << "in the energy range [" << xmin << "," << xmax <<"]";
  }
  Notice("DeceFactor",os.str());
}


