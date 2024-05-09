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
static void   DeceReadMF8   (ENDFDict *, ENDF *, const int, const int, char *, const int);
static void   DeceReadMF9   (ENDFDict *, ENDF *, const int, const int, char *, const int);
static struct Qval qvalues  (const int, const int, const int, const int, const double, const double);
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
  if( (mf != 1) && (mf != 3)){
    message << "MF" << mf << " should be 1 or 3";
    WarningMessage();
    return;
  }

  /*** allocate data array */
  cx   = new double [MAX_DBLDATA];
  cy   = new double [MAX_DBLDATA];
  xdat = new double [MAX_DBLDATA*2];

  /*** for each MF case */
  if(mf == 1) DeceReadMF1(dict,lib,mt,datafile,ofset);
  else        DeceReadMF3(dict,lib,mt,datafile,ofset,readflag);
//ENDFWrite(lib);

  /*** Clean all */
  delete [] cx;
  delete [] cy;
  delete [] xdat;
}


/**********************************************************/
/*      Read in Radioactive Data from a File              */
/*      up to 3 levels for each nuclide                   */
/*      specfied by ofset[0], [1], and [2]                */
/**********************************************************/
void DeceReadRadioactive(ENDFDict *dict, ENDF *lib8, ENDF *lib9, const int mf, const int mt, char *datafile, const int ofset)
{
  if( (mf != 9) && (mf != 10)){
    message << "MF" << mf << " should be 9 or 10";
    WarningMessage();
    return;
  }

  /*** allocate data array */
  cx   = new double [MAX_DBLDATA];
  cy   = new double [MAX_DBLDATA];
  xdat = new double [MAX_DBLDATA*2];

  /*** create MF8 first */
  DeceReadMF8(dict,lib8,mf,mt,datafile,ofset);
//ENDFWrite(lib8);

  /*** create MF9 or 10 */
  DeceReadMF9(dict,lib9,mf,mt,datafile,ofset);
//ENDFWrite(lib9);

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
  else{
    np = geneCSdata(nc,cx,cy,q.et,dict->emaxRe,xdat);
  }

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
void DeceReadMF8(ENDFDict *dict, ENDF *lib, const int lmf, const int mt, char *datafile, int ofset)
{
  const int mf = 8;
  const int no = 1; // decay chain not given here
  const int maxlevel = 3; // default number of max levels, gs, m1, and m2

  struct Prod prd = readMShead(datafile,mt,ofset);

  /*** number of states to be included (for T1/2 > 0) */
  int ns = 0;
  for(int i=0 ; i<maxlevel ; i++) if(prd.th[i] > 0.0) ns ++;

  /*** Make HEAD and CONT */
  lib->setENDFhead(dict->getZA(),dict->getAWR(),dict->getLIS(),dict->getLISO(),ns,no);
  lib->setENDFmat(dict->getMAT());
  lib->setENDFmf(mf);
  lib->setENDFmt(mt);

  for(int is=0 ; is<maxlevel ; is++){
    if(prd.th[is] == 0.0) continue;

    /*** radioactive data */
    /*** make TAB1 */
    Record cont;
    cont.setRecord((double)prd.za,prd.ex[is],lmf,prd.nx[is],0,0);
    ENDFPackCONT(cont,lib);
    message << "radioactive state in " << prd.za << " excitation energy:" << prd.ex[is] << " LFS:" << prd.nx[is];
    Notice("DeceRead:DeceReadMF8");
  }
}


/**********************************************************/
/*      Read External File in MF 9                        */
/**********************************************************/
void DeceReadMF9(ENDFDict *dict, ENDF *lib, const int mf, const int mt, char *datafile, const int ofset)
{
  const int maxlevel = 3;

  struct Prod prd = readMShead(datafile,mt,ofset);
  double *cy0  = new double [MAX_DBLDATA];
  double *cy1  = new double [MAX_DBLDATA];
  double *cy2  = new double [MAX_DBLDATA];

  int ns = 0;
  for(int i=0 ; i<maxlevel ; i++) if(prd.th[i] > 0.0) ns ++;

  /*** Make HEAD and CONT */
  lib->setENDFhead(dict->getZA(),dict->getAWR(),dict->getLIS(),0,ns,0);
  lib->setENDFmat(dict->getMAT());
  lib->setENDFmf(mf);
  lib->setENDFmt(mt);

  int nc = readMSdata(datafile,mf,3*prd.col-2,cx,cy0); // ground state
           readMSdata(datafile,mf,3*prd.col-1,cx,cy1); // first meta
           readMSdata(datafile,mf,3*prd.col  ,cx,cy2); // second meta

  int ntotal = 0;
  for(int is=0 ; is<maxlevel ; is++){
    if(prd.th[is] == 0.0) continue;

    /*** isomeric ratio data */
    if(mf == 9){
      for(int j=0 ; j<nc ; j++){
        double ctot = cy0[j] + cy1[j] + cy2[j];
        if(ctot > 0.0){
          if     (is == 0) cy[j] = cy0[j] / ctot;
          else if(is == 1) cy[j] = cy1[j] / ctot;
          else             cy[j] = cy2[j] / ctot;
        }
        else cy[j] = 0.0;
      }
    }
    /*** production cross section data */
    else{
      for(int j=0 ; j<nc ; j++){
        if     (is == 0) cy[j] = cy0[j];
        else if(is == 1) cy[j] = cy1[j];
        else             cy[j] = cy2[j];
      }
    }

    /*** calculate Qm and Qi */
    struct Qval q = qvalues((int)dict->getZA(),dict->getProj(),(int)dict->getZA(),mt,dict->getELIS(),prd.ex[is]);

    int np = geneCSdata(nc,cx,cy,q.et,0.0,xdat); // ignore resonance range

    /*** non-threshold reaction */
    if(q.et <= 0.0){
      /*** look for the first non-zero data */
      double cz = 0.0;
      int    nz = 0;
      for(int j=0 ; j<np ; j++){
        if(xdat[2*j+1] != 0.0){ cz = xdat[2*j+1]; nz = j; break; }
      }
      for(int j=0 ; j<nz ; j++) xdat[2*j+1] = cz;

    }

    /*** make TAB1 */
    Record cont;
    int    idat[2];

    cont.setRecord(q.qm,q.qi,prd.za,prd.nx[is],1,np);
    idat[0] = np;
    idat[1] = 2;

    ENDFPackTAB1(cont,idat,xdat,lib);

    message << "number of points added " << np << " in MF:" << mf << " MT:" << mt << " LFS:" << prd.nx[is] << " Ex:" << prd.ex[is] << " T(1/2):" << prd.th[is] << " sec";
    Notice("DeceRead:DeceReadMF9");

    ntotal += np;
  }

  if(ntotal == 0){
    message << "no radiactive production data to be added from " << datafile << " for MT = " << mt;
    WarningMessage();
    DeceDelete(dict,mf,mt);
    DeceDelete(dict,8,mt);
  }

  delete [] cy0;
  delete [] cy1;
  delete [] cy2;

  return;
}


/**********************************************************/
/*      Calculate Qm, Qi                                  */
/**********************************************************/
struct Qval qvalues(const int za, const int proj, const int targ, const int mt, const double elis, const double el)
{
  /*** find Q-values */
  double qm = mass_qvalue(proj,targ,mt);
  double qi = qm + elis - el;
  double et = mass_threshold((int)za,qi);

  /*** for exothermic reaction, there is no threshold */
  if(et < 0.0) et = 0.0;

  struct Qval q = {qm,qi,et};

  message << "Q-values for MT = " << mt << "  Ex(target): " << elis << "  Ex(residual): " << el << "  Q(mass): " << q.qm << "  Q(level): " << q.qi << "  Threshold: " << q.et;
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
