/******************************************************************************/
/**     DeCE CHECK TOTAL                                                     **/
/******************************************************************************/

#include <iostream>
#include <iomanip>


using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Scan ENDF Tape and Print Index                    */
/**********************************************************/
void DeceCheckTotal(ENDFDict *dict, ENDF *lib[])
{
  const int ndiv = 10;

  int *mt = new int [1000];

  int mf = 3;
  int k0 = dict->getID(mf,1);

  if(k0 < 0 ){
    message << "MF/MT = 3/1 not given";
    TerminateCode("DeceCheckTotal");
  }

  /*** scan all MTs for summation */
  int k = 0;
  for(int i=0 ; i<dict->getSEC() ; i++){

    int mx = dict->mt[i];
  
    if(dict->mf[i] != mf) continue;
    if(dict->getID(mf,mx) < 0) continue;
    if( (mx == 1) || (mx == 3) || (mx == 4) || (mx == 5) || (mx > 118) ) continue;

    mt[k++] = mx;
  }

  /*** summation check */
  cout << setprecision(6);
  int np = lib[k0]->rdata[0].n2;

  for(int n=0 ; n<np ; n++){
    /*** first, at each total cross section  energy */
    double x0 = lib[k0]->xdata[2*n];
    double y0 = lib[k0]->xdata[2*n+1];

    double s0 = 0.0;
    for(int i=0 ; i<k ; i++){
      s0 += ENDFInterpolation(lib[dict->getID(mf,mt[i])],x0,false,0);
    }

    cout << setw( 6) <<  n << "      ";
    cout << setw(14) << x0;
    cout << setw(14) << y0;
    cout << setw(14) << s0;
    if(y0 == 0.0){
      cout << setw(14) << y0 - s0;
    }
    else{
      cout << setw(14) << 1.0 - s0/y0;
    }
    cout << endl;

    if(n == np-1) break;

    /*** second, in between */
    double dx = (lib[k0]->xdata[2*n+2] - x0) / ndiv;
    for(int j=1 ; j<ndiv ; j++){
      double x1 = x0 + j * dx;
      double y1 = ENDFInterpolation(lib[k0],x1,false,0);

      double s1 = 0.0;
      for(int i=0 ; i<k ; i++){
        s1 += ENDFInterpolation(lib[dict->getID(mf,mt[i])],x1,false,0);
      }

      cout << "      " << setw( 6) <<  j;
      cout << setw(14) << x1;
      cout << setw(14) << y1;
      cout << setw(14) << s1;
      if(y1 == 0.0){
        cout << setw(14) << y1 - s1;
      }
      else{
        cout << setw(14) << 1.0 - s1/y1;
      }
      cout << endl;
    }
 }

  delete [] mt;
}
