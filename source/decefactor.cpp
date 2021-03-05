/******************************************************************************/
/**     DeCE FACTOR                                                          **/
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"

static void DeceFactorMF3(ENDF *, double, double, double, double);
static void DeceFactorMF4P1(ENDF *, double, double, double);


/**********************************************************/
/*      TAB1 Data Multiplied by Factor                    */
/*      Rescale Y-data if x > 0.0                         */
/**********************************************************/
void DeceFactor(ENDFDict *dict, ENDF *lib[], const int mf, const int mt,
                double x, double y, double xmin, double xmax)
{
  if( (mf != 3) && (mf != 4) ){
    message << "DeceFactor does not work for MF = " << mf <<", MT = " << mt;
    WarningMessage();
    return;
  }

  int k0 = dict->getID(mf,mt);
  if(k0 < 0){ message << "MT number " << mt << " not found"; TerminateCode("DeceFactor"); }

  if(mf == 3){
    DeceFactorMF3(lib[k0],x,y,xmin,xmax);
  }
  else if(mf == 4){
    if(mt != 2){
      message << "DeceFactor does not work for MF = " << mf <<", MT = " << mt;
      WarningMessage();
      return;
    }
    DeceFactorMF4P1(lib[k0],y,xmin,xmax);
  }

  message << "MF" << mf << ":MT" << mt << " rescaled by factor " << y;
  if( (xmin < xmax) && (xmin >= 0.0) && (xmax > 0.0) ){
    message << "in the energy range [" << xmin << "," << xmax <<"]";
  }
  Notice("DeceFactor");
}


/**********************************************************/
/*      Rescaling MF3 Cross Section Data                  */
/**********************************************************/
void DeceFactorMF3(ENDF *lib, double x, double y, double xmin, double xmax)
{
  double f  = 0.0;
  /*** all Y data are multiplied by a factor */
  if(x == 0.0){
    f = y;
  }
  /*** renormalize all Y data at given (x,y) */
  else{
    double z = ENDFInterpolation(lib,x,false,0);
    if(z != 0.0) f = y/z;
  }

  Record r  = lib->rdata[0];
  int    np = r.n2;
  /*** if energy range is given */
  if( (xmin < xmax) && (xmin >= 0.0) && (xmax > 0.0) ){
    for(int ip=0 ; ip<np ; ip++){
      if( (xmin <= lib->xdata[2*ip]) && (lib->xdata[2*ip] <= xmax) ){
        lib->xdata[2*ip+1] *= f;
      }
    }
  }
  /*** manipulate entire energy range */
  else{
    for(int ip=0 ; ip<np ; ip++) lib->xdata[2*ip+1] *= f;
  }
}


/**********************************************************/
/*      Rescaling MF4, P1 Legendre Coefficients           */
/**********************************************************/
void DeceFactorMF4P1(ENDF *lib, double f, double xmin, double xmax)
{
  int idx = 0;
  Record head = lib->getENDFhead();
  Record cont = lib->rdata[idx ++];

  if(head.l2 != 1 && head.l2 != 3) return; // not Legendre
  if(cont.l1 == 1) return; // isotropic

  cont = lib->rdata[idx ++];
  int ne = cont.n2;

  bool range = false;
  if( (xmin < xmax) && (xmin >= 0.0) && (xmax > 0.0) ) range = true;

  int j = 0; // P1
  for(int i=0 ; i<ne ; i++){
    double e  = lib->rdata[idx].c2;

    if(range){
      if( (xmin <= e) && (e <= xmax) ) lib->xptr[idx][j] *= f;
    }
    else lib->xptr[idx][j] *= f;

    idx++;
  }
}

