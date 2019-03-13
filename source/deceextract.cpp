/******************************************************************************/
/**     DeCE EXTRACT                                                         **/
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Select ENDF Tape or In Memory                     */
/**********************************************************/
void DeceExtract(ENDFDict *dict, ENDF *lib[], ifstream *fp, const int mf, const int mt)
{
  if(mf <= 0) return;

  /*** extract all MT sections in MF */
  if(mt == 0){
    for(int i=0 ; i<dict->sec ; i++){
      if(mf == dict->mf[i]){
        /*** check if in Lib */
        int k = dict->getID(mf,dict->mt[i]);
        if(k == -2) continue;
        else if(k < 0) ENDFExtract(fp,mf,dict->mt[i]);
        else           ENDFWrite(lib[k]);
      }
    }
  }
  /*** when MF and MT specified */
  else{
    for(int i=0 ; i<dict->sec ; i++){
      if( (mf == dict->mf[i]) && (mt == dict->mt[i]) ){
        /*** check if in Lib */
        int k = dict->getID(mf,mt);
        if(k == -2) continue;
        else if(k < 0) ENDFExtract(fp,mf,dict->mt[i]);
        else           ENDFWrite(lib[k]);
      }
    }
  }
}


