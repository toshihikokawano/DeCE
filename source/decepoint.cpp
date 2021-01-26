/******************************************************************************/
/**     DeCE POINT                                                           **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"

static void addpoint(ENDF *, const double, const double);
static void delpoint(ENDF *, const double);
static void modpoint(ENDF *, const double, const double);

static const double eps = 1.0e-7;


/**********************************************************/
/*      Manipulate One Point in TAB1 Data                 */
/**********************************************************/
void DecePoint(ENDFDict *dict, ENDF *lib[], const int mf, const int mt, double x, double y, string op)
{
  if(mf != 3){
    message << "currently "<< op << " command works only for MF3";
    WarningMessage();
    return;
  }

  int k0 = dict->getID(3,mt);

  if(k0 < 0){ message << "MT number " << mt << " not found"; TerminateCode("DecePoint"); }

  Record r  = lib[k0]->rdata[0];
  int    np = r.n2;

  /*** insert one point */
  if(op == "addpoint"){
    /*** check memory size if an extra point can be added */
    if(lib[k0]->checkDataSize(0,1)){
      message << "cannot add more points since the data size reached the maximum of " << MAX_DBLDATA;
      WarningMessage();
      return;
    }

    addpoint(lib[k0],x,y);
  }

  /*** remove one point */
  else if(op == "delpoint"){
    if(x < y){
      /*** when range of delete given */
      ENDF tmp;
      ENDFLibCopy(lib[k0],&tmp);

      for(int ip=0 ; ip<np ; ip++){
        double z = tmp.xdata[2*ip];
        if((x <= z) && (z <= y)){
          delpoint(lib[k0],z);
        }
      }
    }
    else{
      /*** when one point given */
      delpoint(lib[k0],x);
    }
  }

  /*** change one data */
  if(op == "modpoint"){
    modpoint(lib[k0],x,y);
  }

//ENDFWriteHEAD(lib[k0]);
//ENDFWriteTAB1(lib[k0]);
//ENDFWriteSEND(lib[k0]);
}


/**********************************************************/
/*      Insert One Point                                  */
/**********************************************************/
void addpoint(ENDF *lib, const double x, const double y)
{
  Record r  = lib->rdata[0];
  int    nr = r.n1;
  int    np = r.n2;

  int irp = 0, ipp = 0;

  /*** if the end of array */
  if(x >= lib->xdata[2*(np-1)]){
    lib->xdata[2*np  ] = x;
    lib->xdata[2*np+1] = y;
    lib->idata[2*nr-2] ++;
  }
  /*** general case, find the interval */
  else{
    /*** interpolation range (IRP) and position (IPP) the new data will be inserted */
    int i = 0;
    for(int ir=0 ; ir<nr ; ir++){
      for(int ip=i ; ip<lib->idata[2*ir] ; ip++){
        if(ip == np-1) break;
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

    /*** increment interpolation ranges after IRP */
    for(int ir=irp ; ir<nr ; ir++) lib->idata[2*ir] ++;
  }

  /*** increment total number of data points */
  np++;
  r.n2 = np;
  lib->rdata[0] = r;

  /*** update position pointer */ 
  lib->xptr[lib->getPOS()] = &lib->xdata[2*np];

  message << "MF" << lib->getENDFmf() << "MT" <<lib->getENDFmt() << " data (";
  message << setw(13) << setprecision(6) << x << ",";
  message << setw(13) << setprecision(6) << y << ") inserted";
  Notice("DecePoint:addpoint");
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
  int irp = 0, ipp = 0;
  int i = 0;
  bool found = false;
  for(int ir=0 ; ir<nr ; ir++){
    for(int ip=i ; ip<lib->idata[2*ir] ; ip++){

      double dx = 0.0;
      if(x == 0.0) dx = fabs(lib->xdata[2*ip] - x);
      else dx = fabs(lib->xdata[2*ip] / x - 1.0);

      if(dx < eps){
        irp = ir;
        ipp = ip;
        found = true;
        break;
      }
      i = lib->idata[2*ir]-1;
    }
    if(found) break;
  }

  if(!found){
    message << "MF" << lib->getENDFmf() << "MT" <<lib->getENDFmt();
    message << " does not have data a point at " << setw(13) << setprecision(6) << x;
    WarningMessage();
    return;
  }

  /*** shift the array from high to low */
  double xd = lib->xdata[2*ipp  ];
  double yd = lib->xdata[2*ipp+1];

  for(int ip=ipp ; ip<np-1 ; ip++){
    lib->xdata[2*ip  ] = lib->xdata[2*ip+2];
    lib->xdata[2*ip+1] = lib->xdata[2*ip+3];
  }

  /*** decrement interpolation range after IRP */
  for(int ir=irp ; ir<nr ; ir++) lib->idata[2*ir] --;

  np--;
  r.n2 = np;
  lib->rdata[0] = r;

  /*** update position pointer */ 
  lib->xptr[lib->getPOS()] = &lib->xdata[2*np];

  message << "MF" << lib->getENDFmf() << "MT" <<lib->getENDFmt() << " data (";
  message << setw(13) << setprecision(6) << xd << ",";
  message << setw(13) << setprecision(6) << yd << ") removed";
  Notice("DecePoint:delpoint");
}


/**********************************************************/
/*      Change One Data Point                             */
/**********************************************************/
void modpoint(ENDF *lib, const double x, const double y)
{
  Record r  = lib->rdata[0];
  int    nr = r.n1;
  int    np = r.n2;

  int ipp = 0;

  int i = 0;
  bool found = false;
  for(int ir=0 ; ir<nr ; ir++){
    for(int ip=i ; ip<lib->idata[2*ir] ; ip++){

      double dx = 0.0;
      if(x == 0.0) dx = fabs(lib->xdata[2*ip] - x);
      else dx = fabs(lib->xdata[2*ip] / x - 1.0);

      if(dx < eps){
        ipp = ip;
        found = true;
        break;
      }
      if(ip == np-1) break;
      i = lib->idata[2*ir]-1;
    }
    if(found) break;
  }

  if(!found){
    message << "MF" << lib->getENDFmf() << "MT" <<lib->getENDFmt();
    message << " does not have data a point at " << setw(13) << setprecision(6) << x;
    WarningMessage();
    return;
  }

  /*** update data at the given X */
  lib->xdata[2*ipp+1] = y;

  message << "MF" << lib->getENDFmf() << "MT" <<lib->getENDFmt() << " data (";
  message << setw(13) << setprecision(6) << x << ",";
  message << setw(13) << setprecision(6) << y << ") updated";
  Notice("DecePoint:modpoint");
}


