/******************************************************************************/
/**     DeCE TABLE for MF4                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "decemisc.h"
#include "global.h"
#include "terminate.h"
#include "constant.h"

static int decetable4LEG (ENDF *, int);
static int decetable4TAB (ENDF *, int);


/**********************************************************/
/*      Process MF=4                                      */
/**********************************************************/
void DeceTableMF4(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    lvt  = head.l1;   // transformation matrix flag (deprecated)
  int    ltt  = head.l2;   // 0: isotropic, 1: Legendre, 2: tabulated
  int    idx  = 0;
  Record cont = lib->rdata[idx];
  int    li   = cont.l1;   // 0: non-isotropic, 1: isotropic
  int    lct  = cont.l2;   // 1: lab-system, 2: cms-system

  cout << "# Angular distribution" << endl;
  cout << "#          LVT" << setw(14) << lvt << "  0:transformation matrix not give, 1:given" << endl;
  cout << "#          LTT" << setw(14) << ltt << "  0:isotropic, 1:Legendre, 2:tabulated" << endl;
  cout << "#           LI" << setw(14) << li  << "  0:non-isotropic, 1:isotropic" << endl;
  cout << "#          LCT" << setw(14) << lct << "  1:LAB, 2:CMS" << endl;

  if(li == 1){
    cout << "# isotropic angular distribution" << endl;
  }
  else{
    idx++;
    if(ltt == 1){
      idx = decetable4LEG(lib,idx);
    }
    else if(ltt == 2){
      idx = decetable4TAB(lib,idx);
    }
    else if(ltt == 3){
      idx = decetable4LEG(lib,idx);
      idx = decetable4TAB(lib,idx);
    }
  }
}


/**********************************************************/
/*      Angular Distribution in Legendre Coefficient      */
/**********************************************************/
int decetable4LEG(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];
  int    ne   = cont.n2;
  idx++;

  int da = opt.AngleStep;

  cout << "#           NE" << setw(14) << ne << "  number of incident energy points" << endl;

  for(int i=0 ; i<ne ; i++){
    double e  = lib->rdata[idx].c2;
    int    nl = lib->rdata[idx].n1;

    if(da > 0){
      int np = 180/da;
      if( (180%da) == 0 ) np++;
      cout << "#           NP" << setw(14) << np << endl;
      cout << "# energy        angle         probability" << endl;
      for(int t=0 ; t<=180 ; t+=da){
        double f=0.5;
        for(int j=0 ; j<nl ; j++){
          f += (j+1.5)*lib->xptr[idx][j]*legendre(j+1,(double)t);
        }
        outVal(e);
        outVal(t);
        outVal(f);
        cout << endl;
      }
    }else{
      cout << "#           NL" << setw(14) << nl << endl;
      cout << "# Energy       Leg. Order    Coefficient" << endl;
      outVal(e); outVal(0); outVal(1.0); cout << endl;
      for(int j=0 ; j<nl ; j++){
        outVal(e);
        outVal(j+1);
        outVal(lib->xptr[idx][j]);
        cout << endl;
      }
    }
    idx++;

    cout << endl;
    cout << endl;
  }

  return(idx);
}


/**********************************************************/
/*      Angular Distribution in Table                     */
/**********************************************************/
int decetable4TAB(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];
  int    ne   = cont.n2;
  idx++;

  cout << "#           NE" << setw(14) << ne << "  number of incident energy points" << endl;

  for(int i=0 ; i<ne ; i++){
    double e  = lib->rdata[idx].c2;
    int    np = lib->rdata[idx].n2;

    cout << "#            E"; outVal(e); cout << endl;
    cout << "#           NP" << setw(14) << np << endl;
    cout << "# energy        angle         probability" << endl;
    for(int j=0 ; j<np ; j++){
      outVal(e);
      outVal(acos(lib->xptr[idx][2*j])/PI * 180.0);
      outVal(lib->xptr[idx][2*j+1]);
      cout << endl;
    }
    idx++;

    cout << endl;
    cout << endl;
  }

  return(idx);
}
