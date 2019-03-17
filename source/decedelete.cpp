/******************************************************************************/
/**     DeCE DELETE                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Delete Section from Dictionary                    */
/**********************************************************/
void DeceDelete(ENDFDict *dict, const int mf, const int mt)
{
  if(mf <= 0) return;

  /*** mark id = -2 to be deleted */
  if(mt == 0){
    for(int i=0 ; i<dict->getSEC() ; i++){
      if(mf == dict->mf[i]){
        dict->setID(i,-2);
        ostringstream os;
        os << "MF" << mf << ":MT" << dict->mt[i] << " deleted";
        Notice("DeceDelete",os.str());
      }
    }
  }
  else{
    for(int i=0 ; i<dict->getSEC() ; i++){
      if( (mf == dict->mf[i]) && (mt == dict->mt[i]) ){
        dict->setID(i,-2);
        ostringstream os;
        os << "MF" << mf << ":MT" << mt << " deleted";
        Notice("DeceDelete",os.str());
      }
    }
  }
}


