/******************************************************************************/
/**     DeCE SCAN INDEX                                                      **/
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Scan ENDF Tape and Print Index                    */
/**********************************************************/
void DeceScanIndex(ENDFDict *dict)
{
  /*** scan all MF numbers, up to 40 */
  for(int mf=1 ; mf <= 40 ; mf++){

    bool given = false;
    for(int i=0 ; i<dict->sec ; i++) if(dict->mf[i] == mf) given = true;

    if(given){
      cout << "MF " << setw(2) << "  " <<  mf << endl;

      int c=0;
      for(int i=0 ; i<dict->sec ; i++){
        if(dict->mf[i] == mf){
          if(dict->getID(mf,dict->mt[i]) > 0){
            cout << setw(4) << dict->mt[i];
            if(c != 0 && (c+1)%20 == 0) cout << endl;
            c++;
          }
        }
      }
      if(c%20 != 0) cout << endl;
    }
  }
}
