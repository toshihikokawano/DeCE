/******************************************************************************/
/**     DeCE READ                                                            **/
/******************************************************************************/

#include <iostream>
#include <ostream>

using namespace std;

#include "dece.h"
#include "deceread.h"
#include "global.h"
#include "terminate.h"
#include "masstable.h"

static void   DeceReadMF1   (ENDFDict *, ENDF *, const int, char *, const int);
static void   DeceReadMF3   (ENDFDict *, ENDF *, const int, char *, const int, const int);
static void   DeceReadMF8   (ENDFDict *, ENDF *, const int, char *, const int, const int);
static void   DeceReadMF9   (ENDFDict *, ENDF *, const int, char *, const int, const int);
static struct Qval qvalues (const int, const int, const int, const int, const double, const double);
static double findBoundary  (ENDF *);

static double *cx, *cy, *xdat;


/**********************************************************/
/*      Read in External Data from a File                 */
/*      readflag: 0 general case                          */
/*                1 add background in RR (merge)          */
/*                2 force replace all the cross section   */
/**********************************************************/
void DeceRead(ENDFDict *dict, ENDF *lib, const int mf, const int mt, char *datafile, const int ofset, const int readflag)
{
  /*** allocate data array */
  cx   = new double [MAX_DBLDATA];
  cy   = new double [MAX_DBLDATA];
  xdat = new double [MAX_DBLDATA*2];

  /*** for each MF case */
  if     (mf ==  1) DeceReadMF1(dict,lib,mt,datafile,ofset);
  else if(mf ==  3) DeceReadMF3(dict,lib,mt,datafile,ofset,readflag);
  else if(mf ==  8) DeceReadMF8(dict,lib,mt,datafile,ofset,readflag);
  else if(mf ==  9) DeceReadMF9(dict,lib,mt,datafile,ofset,readflag);
  else{
    message << "MF" << mf << " should be 1, 3, or 9";
    WarningMessage();
  }

  ENDFWrite(lib);

  /*** Clean all */
  delete [] cx;
  delete [] cy;
  delete [] xdat;
}


/**********************************************************/
/*      Read External File in MF 1                        */
/**********************************************************/
void DeceReadMF1(ENDFDict *dict, ENDF *lib, const int mt, char *datafile, const int ofset)
{
  const int mf = 1;

  /*** number of prompt/delayed neutrons, MT=455, 456 */
  int np = readNUdata(datafile,ofset,cx,cy);

  if(np == 0){
    message << "no nu data to be added from " << datafile << " for MT = " << mt;
    WarningMessage();

    DeceDelete(dict,mf,mt);
    return;
  }

  for(int i=0 ; i<np ; i++){
    xdat[2*i  ] = cx[i];
    xdat[2*i+1] = cy[i];
  }

  /*** make TAB1 */
  Record cont;
  int    idat[2];

  /*** Make HEAD and CONT */
  int lnu = 2; // tabulated nu case
  lib->setENDFhead(dict->getZA(),dict->getAWR(),0,lnu,0,0);
  lib->setENDFmat(dict->getMAT());
  lib->setENDFmf(1);
  lib->setENDFmt(mt);

  cont.setRecord(0.0,0.0,0,0,1,np);
  idat[0] = np;
  idat[1] = 2;

  ENDFPackTAB1(cont,idat,xdat,lib);

  message << "number of points added " << np << " in MF:" << mf << " MT:" << mt;
  Notice("DeceRead:DeceReadMF1");

  return;
}


/**********************************************************/
/*      Read External File in MF 3                       */
/**********************************************************/
void DeceReadMF3(ENDFDict *dict, ENDF *lib, const int mt, char *datafile, const int ofset, const int readflag)
{
  const int mf = 3;
  double elev = 0.0;

  /*** determine if the file is for charged particle incident reactions */
  if(dict->getProj() > 1) readSetCharged();

  /*** cross section to discrete levels */
  int nc = 0;
  if( (51 <= mt && mt <= 91) || (600 <= mt && mt <= 849) ){
    nc = readISdata(datafile,ofset,mt,cx,cy,&elev);
  }
  /*** general case */
  else{
    nc = readCSdata(datafile,ofset,mt,cx,cy);
  }
  if(nc == 0){
    message << "no cross section data to be added from " << datafile << " for MT = " << mt;
    WarningMessage();
  }

  /*** calculate Qm and Qi */
  struct Qval q = qvalues((int)dict->getZA(),dict->getProj(),(int)dict->getZA(),mt,dict->getELIS(),elev);

  /*** check resonance boundary, when data will be merged */
  if(readflag == 1){
    double ebtest = findBoundary(lib);
    if(ebtest < dict->emaxRe && ebtest != 1.0e-05){
      message << "maybe background cross sections given for MT = " << mt << " at E1 = " << dict->emaxRe << "  E2 = " << ebtest;
      WarningMessage();
    }
  }

  /*** generate floating point data */
  int np = nc;
  if(readflag == 1){
    /*** check if b.g. is given */
    if(lib->rdata[0].n2 == 0){
      message << "no background cross section is given for MT = " << mt;
      WarningMessage();
      np = 0;
    }
    else np = mergeCSdata(nc,cx,cy,dict->emaxRe,xdat,lib->rdata[0].n2,lib->xptr[0]);
  }
  else if(readflag == 2) np = geneCSdata(nc,cx,cy,-1.0,0.0,xdat);
  else                   np = geneCSdata(nc,cx,cy,q.et,dict->emaxRe,xdat);

  if(np > 1){
    /*** make TAB1 */
    Record cont;
    int    idat[4];

    /*** Make HEAD and CONT */
    lib->setENDFhead(dict->getZA(),dict->getAWR(),0,0,0,0);
    lib->setENDFmat(dict->getMAT());
    lib->setENDFmf(mf);
    lib->setENDFmt(mt);

    if(readflag == 1){
      /*** keep INT in the first range (assume there is only one INT range for the resonance)*/
      if( lib->idata[1] != 2 ){
        cont.setRecord(q.qm,q.qi,0,0,2,np);
        idat[0] = lib->idata[0];
        idat[1] = lib->idata[1];
        idat[2] = np;
        idat[3] = 2;
      }
      else{
        cont.setRecord(q.qm,q.qi,0,0,1,np);
        idat[0] = np;
        idat[1] = 2;
      }
    }
    else{
      cont.setRecord(q.qm,q.qi,0,0,1,np);
      idat[0] = np;
      idat[1] = 2;
    }

    ENDFPackTAB1(cont,idat,xdat,lib);
  }
  else{
    if(readflag != 1) DeceDelete(dict,mf,mt);
  }

  message << "number of points added " << np << " in MF:" << mf << " MT:" << mt;
  Notice("DeceRead:DeceReadMF3");

  return;
}


/**********************************************************/
/*      Read External File in MF 8                        */
/**********************************************************/
void DeceReadMF8(ENDFDict *dict, ENDF *lib, const int mt, char *datafile, const int ofset1, const int ofset2)
{
  const int mf = 8;
  double elfs = 0.0;
  int lfs = 0, izap = 0, lmf = 9;

  /*** number of subsections, limited up to meta2 */
  int ns = 1;
  if(ofset2 != 0) ns = 2;

  int no = 1; // decay chain not given here

  /*** Make HEAD and CONT */
  lib->setENDFhead(dict->getZA(),dict->getAWR(),dict->getLIS(),dict->getLISO(),ns,no);
  lib->setENDFmat(dict->getMAT());
  lib->setENDFmf(mf);
  lib->setENDFmt(mt);

  for(int is=0 ; is<ns ; is++){

    /*** isomeric ratio data */
    int ofset = (is == 0) ? ofset1 : ofset2;
    readMSdata(datafile,ofset,cx,cy,&lfs,&izap,&elfs);

    /*** make TAB1 */
    Record cont;
    cont.setRecord((double)izap,elfs,lmf,lfs,0,0);

    ENDFPackCONT(cont,lib);

    message << "isomer state in " << izap << " excitation energy:" << elfs << " LFS:" << lfs;
    Notice("DeceRead:DeceReadMF8");
  }

  return;
}


/**********************************************************/
/*      Read External File in MF 9                        */
/**********************************************************/
void DeceReadMF9(ENDFDict *dict, ENDF *lib, const int mt, char *datafile, const int ofset1, const int ofset2)
{
  const int mf = 9;
  double elev = 0.0;
  int lfs = 0, izap = 0;

  /*** number of subsections, limited up to meta2 */
  int ns = 1;
  if(ofset2 != 0) ns = 2;

  /*** Make HEAD and CONT */
  lib->setENDFhead(dict->getZA(),dict->getAWR(),dict->getLIS(),0,ns,0);
  lib->setENDFmat(dict->getMAT());
  lib->setENDFmf(mf);
  lib->setENDFmt(mt);

  int ntotal = 0;
  for(int is=0 ; is<ns ; is++){

    /*** isomeric ratio data */
    int ofset = (is == 0) ? ofset1 : ofset2;
    int nc = readMSdata(datafile,ofset,cx,cy,&lfs,&izap,&elev);

    /*** calculate Qm and Qi */
    struct Qval q = qvalues((int)dict->getZA(),dict->getProj(),(int)dict->getZA(),mt,dict->getELIS(),elev);

    /*** ignore resonance range, and include Ethresh only */
    int np = nc;
    np = geneCSdata(nc,cx,cy,q.et,0.0,xdat);

    /*** make TAB1 */
    Record cont;
    int    idat[2];

    cont.setRecord(q.qm,q.qi,izap,lfs,1,np);
    idat[0] = np;
    idat[1] = 2;

    ENDFPackTAB1(cont,idat,xdat,lib);

    message << "number of points added " << np << " in MF:" << mf << " MT:" << mt << " LFS:" << lfs;
    Notice("DeceRead:DeceReadMF9");

    ntotal += np;
  }

  if(ntotal == 0){
    message << "no isomeric ratio data to be added from " << datafile << " for MT = " << mt;
    WarningMessage();
    DeceDelete(dict,mf,mt);
  }

  return;
}


/**********************************************************/
/*      Calculate Qm, Qi                                  */
/**********************************************************/
struct Qval qvalues(const int za, const int proj, const int targ, const int mt, const double elis, const double el)
{
  /*** find Q-values */
  double qm = mass_qvalue(proj,targ,mt) + elis;
  double qi = qm;
  double et = 0.0;

  /*** determine the threshold energy */
  /*** for MT=4, there is no way to get QI unless Elev is given */
  if( (mt == 4) || ((51 <= mt) && (mt <= 91)) ){
    /*** inelastic scattering */
    qi = elis - el;
    et = mass_threshold(za,qi);
  }
  else if( (600 <= mt) && (mt <= 849) ){
    /*** discrete transition by charged particle */
    qi = qm - el;
    if(qi < 0.0) et = mass_threshold((int)za,qi);
  }
  else{
    /*** all other reactions */
    if(qm < 0.0) et = mass_threshold((int)za,qm);
  }

  struct Qval q = {qm,qi,et};

  message << "Q(mass) " << q.qm << " Q(level) " << q.qi << " Threshold Energy " << q.et;
  Notice("DeceRead:qvalues");

  return q;
}


/**********************************************************/
/*      Estimate High-Side of Resonance Range             */
/**********************************************************/
double findBoundary(ENDF *lib)
{
  const double e0 = 1.0e-05;
  double e    = lib->xdata[0];
  Record cont = lib->rdata[0];
  int    np   = cont.n2;

  /*** if x[0] = 0, maybe new data */
  if(e == 0.0) return(e0);

  /*** for threshold reaction, x[0] > 1.0e-5 */
  if(e > e0) return(e);

  /*** in many cases, X-val is duplicated at the boundary of 
       resonance range */
  for(int i=2 ; i<np*2 ; i+=2){
    if( (lib->xdata[i-2] == lib->xdata[i]) &&
        (lib->xdata[i-1] == 0.0) && (lib->xdata[i+1] > 0.0) ){
      e = lib->xdata[i];
    }
  }

  return(e);
}
