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
  if(mf != 3) return;

  int k0 = dict->getID(3,mt);

  if(k0 < 0) TerminateCode("MT number not found",mt);

  Record r  = lib[k0]->rdata[0];
  int    np = r.n2;

  /*** insert one point */
  if(op == "addpoint"){
    addpoint(lib[k0],x,y);
  }

  /*** remove one point */
  else if(op == "delpoint"){

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

  dict->setLineCount(k0,lib[k0]->rdata[0].n2);

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
        if((lib->xdata[2*ip] <= x) && (x < lib->xdata[2*ip+2])){
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

