/******************************************************************************/
/**     DeCE TABLE                                                           **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"


/**********************************************************/
/*      Select ENDF Tape or In Memory                     */
/**********************************************************/
void DeceTable(ENDFDict *dict, ENDF *lib[], ifstream *fp, const int mf, const int mt, const int d)
{
  if( (mf <= 5) || (mf == 8) || (mf == 10) || (mf == 12) || (mf == 13) || (mf == 14) || (mf == 15) || (mf == 31) || (mf == 33) || (mf == 34)){
    int k = dict->getID(mf,mt);

    /*** if already in memory */
    if(k>=0)    DeceLibToTable(lib[k],NULL,d);
    /*** otherwise read a data file */
    else if(k==-1) DeceFileToTable(fp,mf,mt,d);
  }
  else if(mf == 6){
    int k3 = dict->getID(3,mt);
    int k6 = dict->getID(6,mt);

    /*** in memory */
    if(k3>=0 && k6>=0) DeceTableMF6(lib[k3],lib[k6]);
    /*** from a file */
    else DeceFileToTable(fp,mf,mt,d);
  }
  else WarningMessage("table command cannot process MF = ",mf);
}


/**********************************************************/
/*      Read ENDF Data into Lib                           */
/**********************************************************/
void DeceFileToTable(ifstream *fp, const int mf, const int mt, const int d)
{
  /*** more sections to be added */
  if( (1<= mf && mf<=10) || mf==12 || mf==13 || mf==14 || mf == 15 || mf==31 || mf==33 || mf==34){
    ENDF lib(L);
    ENDF sup(M);
    int  c = 0;

    switch(mf){
    case  1: c = ENDFReadMF1( fp,&lib,mt);  break;
    case  2: c = ENDFReadMF2( fp,&lib   );  break;
    case  3: c = ENDFReadMF3( fp,&lib,mt);  break;
    case  4: c = ENDFReadMF4( fp,&lib,mt);  break;
    case  5: c = ENDFReadMF5( fp,&lib,mt);  break;
    case  6: c = ENDFReadMF6( fp,&lib,mt);
                 ENDFReadMF3( fp,&sup,mt);  break;
    case  8: c = ENDFReadMF8( fp,&lib,mt);  break;
    case  9: c = ENDFReadMF9( fp,&lib,mt);  break;
    case 10: c = ENDFReadMF10(fp,&lib,mt);  break;
    case 12: c = ENDFReadMF12(fp,&lib,mt);  break;
    case 13: c = ENDFReadMF13(fp,&lib,mt);  break;
    case 14: c = ENDFReadMF14(fp,&lib,mt);  break;
    case 15: c = ENDFReadMF15(fp,&lib,mt);  break;
    case 31: c = ENDFReadMF31(fp,&lib,mt);  break;
    case 33: c = ENDFReadMF33(fp,&lib,mt);  break;
    case 34: c = ENDFReadMF34(fp,&lib,mt);  break;
    default:                                break;
    }

    if(c < 0) WarningMessage("no section given");
    else      DeceLibToTable(&lib,&sup,d);
  }
  else WarningMessage("table command cannot process MF",mf);
}


/**********************************************************/
/*      Read ENDF Data and Interpolate                    */
/**********************************************************/
void DeceDataPoint(ifstream *fp, const int mf, const int mt, const double e)
{
  /*** work only for MF=1, 2 and 3 */
  if((mf != 1) && (mf != 3)){
    WarningMessage("point command works for MF =1 and 3 only");
    return;
  }

  double y = 0.0;
  ENDF lib(L);
  if(mf == 1){
    ENDFReadMF1(fp,&lib,mt);
    Record head = lib.getENDFhead();

    int lnu  = head.l2;
    if(lnu == 2){
      if( (mt == 452) || (mt == 456) ) y = ENDFInterpolation(&lib,e,false,0);
      else if(mt == 455)               y = ENDFInterpolation(&lib,e,false,1);
    }
  }
  else{
    ENDFReadMF3(fp,&lib,mt);
    y = ENDFInterpolation(&lib,e,false,0);
  }

  cout.setf(ios::fixed, ios::floatfield);
  outVal(e);
  outVal(y);
  cout << endl;
}


/**********************************************************/
/*      In Memory Case                                    */
/**********************************************************/
void DeceLibToTable(ENDF *lib, ENDF *sup, int d)
{
  Record head = lib->getENDFhead();
  double za   = head.c1;
  int z       = (int)za/1000;
  int a       = (int)za - z*1000;
  int mat     = lib->getENDFmat();
  int mf      = lib->getENDFmf();
  int mt      = lib->getENDFmt();


  cout << "#    [" << setw(5) <<  mat
       << " : " << setw(2) << mf
       << " : " << setw(3) << mt
       << " ]   "
       << setw(3) << z << " - " << setw(3) << a << endl;

  cout.setf(ios::fixed, ios::floatfield);

  switch(mf){
  case  1: DeceTableMF1(lib);     break;
  case  2: DeceTableMF2(lib);     break;
  case  3: DeceTableMF3(lib);     break;
  case  4: DeceTableMF4(lib,d);   break;
  case  5: DeceTableMF5(lib);     break;
  case  6: DeceTableMF6(sup,lib); break;
  case  8: DeceTableMF8(lib);     break;
  case  9: DeceTableMF9(lib);     break;
  case 10: DeceTableMF10(lib);    break;
  case 12: DeceTableMF12(lib);    break;
  case 13: DeceTableMF13(lib);    break;
  case 14: DeceTableMF14(lib);    break;
  case 15: DeceTableMF15(lib);    break;
  case 31:
  case 33: DeceTableMF33(lib);    break;
  case 34: DeceTableMF34(lib);    break;
  default:                        break;
  }
}

