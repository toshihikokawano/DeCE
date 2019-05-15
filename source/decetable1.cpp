/******************************************************************************/
/**     DeCE TABLE for MF1                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"

static void   DeceTableMF1MT452 (ENDF *);
static void   DeceTableMF1MT455 (ENDF *);
static void   DeceTableMF1MT458 (ENDF *);


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
    ENDFPrintLIST(lib,0);
  }
  /*** tabulated */
  else if(lnu == 2){
    ENDFPrint1Dim(lib,0);
  }
}


/**********************************************************/
/*      MT455: Delaye Neutron Yields                      */
/**********************************************************/
void  DeceTableMF1MT455(ENDF *lib)
{
  Record head = lib->getENDFhead();

  /*** LNU : nu given by polynomial (1) or table (2) */
  int lnu  = head.l2;

  /*** precursor decay constants */
  int nnf = lib->rdata[0].n1;
  cout << "# Delayed neutron decay constants" << endl;
  cout << "#          NNF" << setw(14) << nnf << "  number of precursor families"  << endl;
  cout << "#   lambda grp   decay const" << endl;
  for(int i=0 ; i<nnf ; i++){
    outVal(i); outVal(lib->xptr[0][i]);
    cout << endl;
  }
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
    ENDFPrint1Dim(lib,1);
  }
}


/**********************************************************/
/*      MT455: Total Energy Release by Fission            */
/**********************************************************/
void  DeceTableMF1MT458(ENDF *lib)
{
  cout << "# Energy release due to fission" << endl;
  cout << "# EFR / dEFR  ";  outVal(lib->xdata[0]); outVal(lib->xdata[1]);
  cout << " fragment kinetic energy" << endl;

  cout << "# ENP / dENP  ";  outVal(lib->xdata[2]); outVal(lib->xdata[3]);
  cout << " kinetic energy of prompt fission neutrons" << endl;

  cout << "# END / dEND  ";  outVal(lib->xdata[4]); outVal(lib->xdata[5]);
  cout << " kinetic energy of delayed fission neutrons" << endl;

  cout << "# EGP / dEGP  ";  outVal(lib->xdata[6]); outVal(lib->xdata[7]);
  cout << " total energy release by prompt fission gammas" << endl;

  cout << "# EGD / dEGD  ";  outVal(lib->xdata[8]); outVal(lib->xdata[9]);
  cout << " total energy release by delayed fission gammas" << endl;

  cout << "# EB  / dEB   ";  outVal(lib->xdata[10]); outVal(lib->xdata[11]);
  cout << " total energy release by delayed betas" << endl;

  cout << "# ENU / dENU  ";  outVal(lib->xdata[12]); outVal(lib->xdata[13]);
  cout << " energy carried away by neutrinos" << endl;

  cout << "# ER  / dER   ";  outVal(lib->xdata[14]); outVal(lib->xdata[15]);
  cout << " ET - ENU = pseudo-Q value" << endl;

  cout << "# ET  / dET   ";  outVal(lib->xdata[16]); outVal(lib->xdata[17]);
  cout << " total energy" << endl;
}

