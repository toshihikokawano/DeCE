/******************************************************************************/
/**     DeCE POINT                                                           **/
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"

static void addpoint(ENDF *, const double, const double);
static void delpoint(ENDF *, const double);


/**********************************************************/
/*      Manipulate One Point in TAB1 Data                 */
/**********************************************************/
void DecePoint(ENDFDict *dict, ENDF *lib[], const int mf, const int mt, double x, double y, string op)
{
  if(mf!=3) return;

  int k0 = dict->getID(3,mt);

  if(k0<0) TerminateCode("MT number not found",mt);

  Record r  = lib[k0]->rdata[0];
  int    np = r.n2;

  /*** insert one point */
  if(op=="addpoint"){
    addpoint(lib[k0],x,y);
  }

  /*** remove one point */
  else if(op=="delpoint"){

    if(x < y){
      /*** when range of delete given */
      ENDF tmp(L);
      ENDFLibCopy(lib[k0],&tmp);
      for(int ip=0 ; ip<np ; ip++){
        double z = tmp.xdata[2*ip];
        if((x <= z) && (z <= y)){
          delpoint(lib[k0],z);
        }
      }
    }else{
      /*** when one point given */
      delpoint(lib[k0],x);
    }
  }

  r  = lib[k0]->rdata[0];
  np = r.n2;

  int nc = np/3 + 4;
  if(np%3==0 && np!=0) nc--;
  dict->nc[k0] = nc;

  //  ENDFWriteHEAD(lib[k0]);
  //  ENDFWriteTAB1(lib[k0]);
  //  ENDFWriteSEND(lib[k0]);
}


/**********************************************************/
/*      Insert One Point                                  */
/**********************************************************/
void addpoint(ENDF *lib, const double x, const double y)
{
  Record r  = lib->rdata[0];
  int    nr = r.n1;
  int    np = r.n2;

  /*** if the end of array */
  if(x >= lib->xdata[2*(np-1)]){
    lib->xdata[2*np  ] = x;
    lib->xdata[2*np+1] = y;
    lib->idata[2*nr-2] ++;
  }
  /*** general case, find the interval */
  else{
    int irp=0, ipp=0;
    int i = 0;
    for(int ir=0 ; ir<nr ; ir++){
      for(int ip=i ; ip<lib->idata[2*ir] ; ip++){
        if(ip==np-1) break;
        if( lib->xdata[2*ip] <= x && x < lib->xdata[2*ip+2]){
          irp = ir;
          ipp = ip+1;
          break;
        }
        i = lib->idata[2*ir]-1;
      }
    }

    /*** shift the array elements to make a space */
    for(int ip=np-1 ; ip>=ipp ; ip--){
      lib->xdata[2*ip+2] = lib->xdata[2*ip  ];
      lib->xdata[2*ip+3] = lib->xdata[2*ip+1];
    }
    /*** add a new point */
    lib->xdata[2*ipp  ] = x;
    lib->xdata[2*ipp+1] = y;
    lib->idata[2*irp] ++;
  }
  np++;
  r.n2 = np;
  lib->rdata[0] = r;
}


/**********************************************************/
/*      Delete One Point                                  */
/**********************************************************/
void delpoint(ENDF *lib, const double x)
{
  Record r  = lib->rdata[0];
  int    nr = r.n1;
  int    np = r.n2;

  /*** find the point to be deleted */
  int irp=0, ipp=0;
  int i = 0;
  bool found = false;
  for(int ir=0 ; ir<nr ; ir++){
    for(int ip=i ; ip<lib->idata[2*ir] ; ip++){
      if( lib->xdata[2*ip] == x ){
        irp = ir;
        ipp = ip;
        found = true;
        break;
      }
      i = lib->idata[2*ir]-1;
    }
  }
  /*** shift the array from high to low */
  if(found){
    for(int ip=ipp ; ip<np-1 ; ip++){
      lib->xdata[2*ip  ] = lib->xdata[2*ip+2];
      lib->xdata[2*ip+1] = lib->xdata[2*ip+3];
    }
    lib->idata[2*irp] --;
  }
  np--;
  r.n2 = np;
  lib->rdata[0] = r;
}


/**********************************************************/
/*      TAB1 Data Multiplied by Factor                    */
/*      Rescale Y-data if x>0.0                           */
/**********************************************************/
void DeceFactor(ENDFDict *dict, ENDF *lib[], const int mf, const int mt,
                double x, double y, double xmin, double xmax)
{
  if(mf!=3) return;

  int k0 = dict->getID(3,mt);
  if(k0<0) TerminateCode("MT number not found",mt);

  double f  = 0.0;
  /*** all Y data are multiplied by a factor */
  if(x==0.0){
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

  //  ENDFWriteHEAD(lib[k0]);
  //  ENDFWriteTAB1(lib[k0]);
  //  ENDFWriteSEND(lib[k0]);
}

/**********************************************************/
/*      TAB1 Data Multiplied by a Function                */
/**********************************************************/
void DeceApplyFunc(ENDFDict *dict, ENDF *lib[], const int mf, const int mt,
                   int kind, double p1, double p2, double p3)
{
  if(p3==0) TerminateCode("P3 Parameter zero");
  if(mf!=3) return;

  int k0 = dict->getID(3,mt);
  if(k0<0) TerminateCode("MT number not found",mt);

  Record r  = lib[k0]->rdata[0];
  int    np = r.n2;

  for(int ip=0 ; ip<np ; ip++){
    double x = lib[k0]->xdata[2*ip] * 1e-6;
    double f = 1.0;

    /*** Fermi function */
    if(kind == 1)      f = p1/( 1.0+ exp( (x-p2)/p3 ) ) + 1.0;
    /*** Gaussian */
    else if(kind == 2) f = p1*exp( -(x-p2)*(x-p2)/p3 )  + 1.0;
    /*** reversed Fermi function */
    else if(kind == 3) f = p1/( 1.0+ exp(-(x-p2)/p3 ) ) + 1.0;

    lib[k0]->xdata[2*ip+1] *= f;
  }
}


/**********************************************************/
/*      Change Interpolation Law in MF3                   */
/**********************************************************/
void DeceChangeInt(ENDFDict *dict, ENDF *lib[], const int mt, int range, int point, int intlaw)
{
  int k = dict->getID(3,mt);

  if(k<0) TerminateCode("MT number not found",mt);

  Record r  = lib[k]->rdata[0];
  int    nr = r.n1;
  int    np = r.n2;

  if(point > np) TerminateCode("data point address exceeds the highest",point);
  if(range <= 0) TerminateCode("data range incorrect",range);
  if(range > nr+1) TerminateCode("data range incorrect",range);

  if(range == nr+1) lib[k]->rdata[0].n1 = nr+1;
  if(point < 0) point = np;

  lib[k]->idata[(range-1)*2  ] = point;
  lib[k]->idata[(range-1)*2+1] = intlaw;
}

