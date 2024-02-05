/******************************************************************************/
/**     DeCE TABLE for MF1                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"

static void   DeceTableMF1MT452 (ENDF *);
static void   DeceTableMF1MT455 (ENDF *);
static void   DeceTableMF1MT458 (ENDF *);
static void   DeceTableMF1MT460 (ENDF *);


/**********************************************************/
/*      Process MF=1                                      */
/**********************************************************/
void  DeceTableMF1(ENDF *lib)
{
  int mt = lib->getENDFmt();
  if(     mt == 452) DeceTableMF1MT452(lib);
  else if(mt == 456) DeceTableMF1MT452(lib);
  else if(mt == 455) DeceTableMF1MT455(lib);
  else if(mt == 458) DeceTableMF1MT458(lib);
  else if(mt == 460) DeceTableMF1MT460(lib);
  else{
    message << "invalid MT number " << mt;
    WarningMessage();
  }
}


/**********************************************************/
/*      MT452: Prompt Neutron or Total Neutron Yields     */
/**********************************************************/
void  DeceTableMF1MT452(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    mt   = lib->getENDFmt();

  if(mt == 452) cout << "# Total neutron yields" << endl;
  else          cout << "# Prompt neutron yields" << endl;

  /*** LNU : nu given by polynomial (1) or table (2) */
  int lnu  = head.l2;

  /*** given by polynomials */
  if(lnu == 1){
    int nc = lib->rdata[1].n1;
    cout << "#           NC" << setw(14) << nc << "  number of polynomial terms" << endl;
    ENDFPrintLIST(lib,0,"Order","Coeff");
  }
  /*** tabulated */
  else if(lnu == 2){
    if(mt == 452) ENDFPrint1Dim(lib,0,"Energy","Nu_total");
    else          ENDFPrint1Dim(lib,0,"Energy","Nu_prompt");
  }
}


/**********************************************************/
/*      MT455: Delaye Neutron Yields                      */
/**********************************************************/
void  DeceTableMF1MT455(ENDF *lib)
{
  Record head = lib->getENDFhead();

  /*** LDG : energy-dependent delayed-group constants */
  int ldg  = head.l1;

  /*** LNU : nu given by polynomial (1) or table (2) */
  int lnu  = head.l2;

  /*** precursor decay constants */
  int nnf = lib->rdata[0].n1;
  cout << "# Delayed neutron decay constants" << endl;
  cout << "#          LDG" << setw(14) << ldg << "  energy dependent delayed-groups" << endl;
  cout << "#          NNF" << setw(14) << nnf << "  number of precursor families"  << endl;
  cout << "#   lambda grp   decay const" << endl;
  ENDFPrintLIST(lib,0,"Lambda Group","Decay Const");
  cout << endl;
  cout << endl;

  /*** given by polynomials */
  cout << "# Delayed neutron yields" << endl;
  if(lnu == 1){
    int nc = lib->rdata[1].n1;
    cout << "#           NC" << setw(14) << nc << "  number of polynomial terms" << endl;
    ENDFPrintLIST(lib,1);
  }
  /*** tabulated */
  else if(lnu == 2){
    ENDFPrint1Dim(lib,1,"Energy","Nu_delayed");
  }
}


/**********************************************************/
/*      MT455: Total Energy Release by Fission            */
/**********************************************************/
void  DeceTableMF1MT458(ENDF *lib)
{
  Record head = lib->getENDFhead();

  /*** LFC : energy-dependence given by polynomial (1) or table (2) */
  int lfc  = head.l2;
  int nfc  = head.n2;

  cout << "# Energy release due to fission" << endl;

  if(lfc == 0){
    int nply = lib->rdata[0].l2;

    cout << "#         LPLY" << setw(14) << nply << "  order of polynomial expansion" << endl;

    int i = 0;
    for(int j=0 ; j<=nply ; j++){

      cout << "#   L   EFR           dEFR        fragment kinetic energy" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   ENP           dENP        kinetic energy of prompt fission neutrons" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   END           dEND        kinetic energy of delayed fission neutrons" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   EGP           dEGP        total energy release by prompt fission gammas" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   EGD           dEGD        total energy release by delayed fission gammas" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   EB            dEB         total energy release by delayed betas" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   ENU           dENU        energy carried away by neutrinos" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   ER            dER         ET - ENU = pseudo-Q value" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;

      cout << "#   L   ET            dET         total energy" << endl;
      cout << setw(5) << j; outVal(lib->xdata[i++]); outVal(lib->xdata[i++]); cout << endl;
      cout << endl;
    }
  }
  else{
    int id =  0;
    cout << "#  EFR           dEFR        fragment kinetic energy" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;


    cout << "#  ENP           dENP        kinetic energy of prompt fission neutrons" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  END           dEND        kinetic energy of delayed fission neutrons" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  EGP           dEGP        total energy release by prompt fission gammas" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  EGD           dEGD        total energy release by delayed fission gammas" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  EB            dEB         total energy release by delayed betas" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  ENU           dENU        energy carried away by neutrinos" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  ER            dER         ET - ENU = pseudo-Q value" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;

    cout << "#  ET            dET         total energy" << endl;
    outVal(lib->xdata[2*id]); outVal(lib->xdata[2*id + 1]); cout << endl;
    for(int k=1 ; k<=nfc ; k++){
      int ifc = lib->rdata[k].l2 - 1;
      if(id == ifc){ ENDFPrint1Dim(lib,k); break; }
    }
    id ++;
    cout << endl;
  }
}


/**********************************************************/
/*      MT460: Delayed Photon Data                        */
/**********************************************************/
void  DeceTableMF1MT460(ENDF *lib)
{
  Record head = lib->getENDFhead();

  cout << "# Delayed photon data" << endl;

  /*** LO : discrete representation (0) or continuous (1) */
  int lo = head.l1;

  /*** discrete */
  if(lo == 0){
    int ng = head.n1;
    cout << "#           NG" << setw(14) << ng << "  number of discrete photons" << endl;
    for(int ig=0 ; ig<ng ; ig++) ENDFPrint1Dim(lib,ig,"Time","Photon Multiplicity");
  }
  /*** continuous */
  else{
    ENDFPrintLIST(lib,0,"Lambda Group","Decay Const");
  }
}
