/******************************************************************************/
/**     DeCE SCAN INDEX                                                      **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <algorithm>


using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Scan ENDF Tape and Print Index                    */
/**********************************************************/
void DeceScanIndex(ENDFDict *dict)
{
  int *mt = new int [1000];

  /*** scan all MF numbers, up to 40 */
  for(int mf=1 ; mf <= 40 ; mf++){

    int k = 0;
    for(int i=0 ; i<dict->getSEC() ; i++){
      if( (dict->mf[i] == mf) && (dict->getID(mf,dict->mt[i]) != -2) ) mt[k++] = dict->mt[i];
    }

    if(k > 0){
      sort(mt,mt+k);

      cout << "MF " << setw(2) << "  " <<  mf << endl;

      int c = 0;
      for(int i=0 ; i<k ; i++){
        cout << setw(4) << mt[i];
        if(c != 0 && (c+1)%20 == 0) cout << endl;
        c++;
      }
      if(c%20 != 0) cout << endl;
    }
  }

  delete [] mt;
}
