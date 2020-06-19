/******************************************************************************/
/**     DeCE CALC                                                            **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"

static void    addsection    (char, ENDF *, ENDF *);
static int     qselect       (int, int, double, double);


/**********************************************************/
/*      Arithmetric Operations for TAB1-type Data         */
/**********************************************************/
void DeceCalc(ENDFDict *dict, ENDF *lib[], const int mtdest, const int mtsrc1, const int mtsrc2, const char op)
{
  int k0 = dict->getID(3,mtdest);
  int k1 = dict->getID(3,mtsrc1);
  int k2 = dict->getID(3,mtsrc2);
  int mf = 3;

  if((k1 < 0) && (k2 < 0)){
    ostringstream os;
    os << "both source MT numbers " << mtsrc1 << " and " << mtsrc2 << " not found";
    TerminateCode(os.str());
  }

  /*** copy src1 to dest */
  if(k2 < 0){
    ENDFLibCopy(lib[k1],lib[k0]);
    lib[k0]->setENDFmt(mtdest);
  }
  /*** copy src2 to dest */
  else if(k1 < 0){
    ENDFLibCopy(lib[k2],lib[k0]);
    lib[k0]->setENDFmt(mtdest);
  }
  else{
    /*** deterimine Q-value */
    double qm = 0.0, qi = 0.0;
    if(mtdest >= 4){
      int kq = qselect(k1,k2,lib[k1]->rdata[0].c1+lib[k1]->rdata[0].c2,
                             lib[k2]->rdata[0].c1+lib[k2]->rdata[0].c2);
      qm = lib[kq]->rdata[0].c1;
      qi = lib[kq]->rdata[0].c2;
    }

    /*** copy SRC1 to DEST */
    ENDFLibCopy(lib[k1],lib[k0]);
    lib[k0]->setENDFmt(mtdest);

    /*** add all MTs in the range to DEST */
    if(op == ':'){
      for(int i=0 ; i<dict->getSEC() ; i++){
        if( (dict->mf[i] == mf) && (dict->mt[i] >= mtsrc1+1 && dict->mt[i] <= mtsrc2) ){
          k2 = dict->getID(3,dict->mt[i]);
          if(k2 < 0) continue;
          addsection('+',lib[k2],lib[k0]);
        }
      }
    }
    /*** add SRC2 to DEST */
    else addsection(op,lib[k2],lib[k0]);

    /*** restore Q values */
    lib[k0]->rdata[0].c1 = qm;
    lib[k0]->rdata[0].c2 = qi;
  }

//  ENDFWriteHEAD(lib[k0]);
//  ENDFWriteTAB1(lib[k0]);
//  ENDFWriteSEND(lib[k0]);
}


/**********************************************************/
/*      Sum Prompt and Delayed Neutrons in MF1            */
/**********************************************************/
void DeceCalc452(ENDFDict *dict, ENDF *lib[])
{
  int k0 = dict->getID(1,452);
  int k1 = dict->getID(1,455);
  int k2 = dict->getID(1,456);
  int mf = 1;

  if(k1 < 0) TerminateCode("MT number not found",455);
  if(k2 < 0) TerminateCode("MT number not found",456);

  Record hd = lib[k1]->getENDFhead();
  Record hp = lib[k2]->getENDFhead();

  if( (hd.l2 != 2) || (hp.l2 !=2 ) ){
    TerminateCode("LNU should be 2 to sum MT 455 and 456");
  }

  /*** first, copy 456(prompt) to 452(total) */
  ENDFLibCopy(lib[k2],lib[k0]);
  lib[k0]->setENDFmt(452);

  /*** create a temporaly MF3-like library for delayed-nu */
  ENDF lib455(M);
  lib455.setENDFmf(mf);
  lib455.setENDFmt(455);
  lib455.setENDFhead(hd);
  lib455.rdata[0] = lib[k1]->rdata[1]; // in order to skip LIST

  for(int i=0 ; i<lib455.rdata[0].n1 ; i++){
    lib455.idata[i*2  ] = lib[k1]->iptr[1][i*2  ];
    lib455.idata[i*2+1] = lib[k1]->iptr[1][i*2+1];
  }

  for(int i=0 ; i<lib455.rdata[0].n2 ; i++){
    lib455.xdata[i*2  ] = lib[k1]->xptr[1][i*2  ];
    lib455.xdata[i*2+1] = lib[k1]->xptr[1][i*2+1];
  }

  /*** add 455(delayed) to 452(total) */
  addsection('+',&lib455,lib[k0]);
}


/**********************************************************/
/*      Add SRC Section Data to DEST                      */
/**********************************************************/
void addsection(char pm, ENDF *src, ENDF *dest)
{
  double   *z;


  z = new double [MAX_DBLDATA];

  int n = ENDFMergeXdata(src,dest,z);

  for(int i=0 ; i<n ; i++){
    double y1 = 0.0, y2 = 0.0;

    /*** duplicated point */
    if( (i < n-1) && (z[i*2] == z[i*2+2]) ){
      y1 = ENDFInterpolation(dest,z[i*2],true,0);
      y2 = ENDFInterpolation(src ,z[i*2],true,0);
    }
    else{
      y1 = ENDFInterpolation(dest,z[i*2],false,0);
      y2 = ENDFInterpolation(src ,z[i*2],false,0);
    }

    switch(pm){
    case '+':  z[i*2+1] = y1+y2; break;
    case '-':  z[i*2+1] = y1-y2; break;
    case '*':  z[i*2+1] = y1*y2; break;
    case '/':  z[i*2+1] = (y2!=0.0) ? y1/y2 : 0.0; break;
    default :  break;
    }

    if(z[i*2+1] < 0.0){
      message << "MF3:MT" << dest->getENDFmt();
      message << " negative value " << setw(13) << z[i*2+1];
      message << " detected at "    << setw(13) << z[i*2];
      message << " where Y1 = "     << setw(13) << y1;
      message << " and Y2 = "       << setw(13) << y2;
      WarningMessage();
    }
  }

  Record rs = src->rdata[0];
  double qms = rs.c1;
  double qis = rs.c2;

  Record rd = dest->rdata[0];
  rd.n1 = 1;
  rd.n2 = n;

  /*** ad hoc, take source Q-value */
  rd.c1 = qms;
  rd.c2 = qis;
  dest->rdata[0] = rd;

  /*** force lin-lin interpolation */
  int idat[2];
  idat[0] = n;
  idat[1] = 2;

  dest->resetPOS();
  ENDFPackTAB1(rd,idat,z,dest);

  delete [] z;
}


/**********************************************************/
/*      Select an Appropriate Q-value                     */
/**********************************************************/
int qselect(int k1, int k2, double q1, double q2)
{
  int k = 0;

  /*** If both are the threshold reaction
       we take the higher one.
       For example, Q1 = -10, Q2 = -20, then Q = -10 */
  if(q1 < 0.0 && q2 < 0.0){
    k = (abs(q1) < abs(q2)) ? k1 : k2;
  }
  /*** If one of Q-values is positive, take this one */
  else if( (q1 > 0.0 && q2 <= 0.0) || (q1 <= 0.0 && q2 > 0.0) ){
    k = (q1 > 0.0) ? k1 : k2;
  }
  /*** If one of Q-values is zero, Q is set to zero */
  else if( (q1 == 0.0) || (q2 == 0.0) ){
    k = (q1 == 0.0) ? k1 : k2;
  }
  /*** If both are positive, take the larger one.
       In this case, the Q-value is not so well defined. */
  else if( (q1 > 0.0) && (q2 > 0.0) ){
    k = (q1 > q2) ? k1 : k2;
  }
  else k = k1;

  return(k);
}
