/******************************************************************************/
/**     DeCE COPY / EXTRACT / SCAN                                           **/
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


/**********************************************************/
/*      Copy One Section from a File                      */
/**********************************************************/
void DeceLibRead(ENDFDict *dict, ENDF *lib, char *file)
{
  ifstream fpin;

  fpin.open(file);
  if(!fpin) TerminateCode("cannot open data file",file);

  ENDFRead(&fpin,lib,lib->getENDFmf(),lib->getENDFmt());
  fpin.close();

  /*** replace ZA and AWR from Dictionary */
  Record head = dict->head;
  double za   = head.c1;
  double awr  = head.c2;

  head = lib->getENDFhead();
  head.c1 = za;
  head.c2 = awr;

  lib->setENDFhead(head);
}


/**********************************************************/
/*      Delete Section from Dictionary                    */
/**********************************************************/
void DeceDelete(ENDFDict *dict, const int mf, const int mt)
{
  if(mf <= 0) return;

  /*** mark id = -2 to be deleted */
  if(mt == 0){
    for(int i=0 ; i<dict->sec ; i++){
      if(mf == dict->mf[i]){
        dict->setID(i,-2);
        ostringstream os;
        os << "MF" << mf << ":MT" << dict->mt[i] << " deleted";
        Notice("DeceDelete",os.str());
      }
    }
  }
  else{
    for(int i=0 ; i<dict->sec ; i++){
      if( (mf == dict->mf[i]) && (mt == dict->mt[i]) ){
        dict->setID(i,-2);
        ostringstream os;
        os << "MF" << mf << ":MT" << mt << " deleted";
        Notice("DeceDelete",os.str());
      }
    }
  }
}


/**********************************************************/
/*      Scan ENDF Tape                                    */
/**********************************************************/
void DeceScan(ENDFDict *dict)
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
