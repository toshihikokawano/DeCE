/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Merge Angular Distributions in Two Energy Regions       **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>

using namespace std;

#include "../source/endflib.h"

int  main (int, char *[]);
void mergeMF4 (ENDF *, ENDF *, ENDF *);


/**********************************************************/
/*      Merge 2 Angular Distribution Data from Files      */
/**********************************************************/
int main(int argc, char *argv[])
{
  if(argc < 2){
    cout << "usage: decemergemf4 ENDF_low_energy ENDF_high_energy\n" << endl;
    exit(-1);
  }

  ifstream fpin;
  ENDF     lib0, lib1;

  string libname0 = argv[1];
  string libname1 = argv[2];

  /*** read MF4 MT2 from low energy data */
  fpin.open(libname0.c_str());
  if(!fpin){
    cerr << "ENDF file " << libname0 << " cannot open" << endl;
    exit(-1);
  }
  ENDFReadMF4(&fpin,&lib0,2);
  fpin.close();

  /*** read MF4 MT2 from high energy data */
  fpin.open(libname1.c_str());
  if(!fpin){
    cerr << "ENDF file " << libname1 << " cannot open" << endl;
    exit(-1);
  }
  ENDFReadMF4(&fpin,&lib1,2);
  fpin.close();

  /*** process MF4 */
  ENDF lib;
  mergeMF4(&lib0,&lib1,&lib);

  /*** print results */
  ENDFWriteHEAD(&lib);
  ENDFWriteCONT(&lib);
  ENDFWriteTAB2(&lib);
  ENDFWriteTAB21(&lib);
  ENDFWriteSEND(&lib);

  return(0);
}


/**********************************************************/
/*      Tabulate Scattering Angular Distribution          */
/**********************************************************/
void mergeMF4(ENDF *lib0, ENDF *lib1, ENDF *lib)
{
  Record head = lib0->getENDFhead();

  if(head.l2 != 1){
    cerr << "low energy region should be Legendre coefficient data" << endl;
    exit(-1);
  }

  int idx0 = 0;
  Record cont = lib0->rdata[idx0++];

  /*** Make HEAD and CONT */
  head.l2 = 3; // Legendre coefficients and tabulated data

  lib->setENDFmat(lib0->getENDFmat());
  lib->setENDFmf(4);
  lib->setENDFmt(2);
  lib->setENDFhead(head);
  ENDFPackCONT(cont,lib);

  /*** low energy region, copy first CONT record, determine highest energy */
  cont = lib0->rdata[idx0++];
  int    ne0 = cont.n2;
  double emax = lib0->rdata[ne0+1].c2;

  /*** copy TAB2 */
  ENDFPackTAB2(cont,&lib0->rdata[idx0],lib0->iptr[idx0],&lib0->xptr[idx0],lib);

  /*** high energy region */
  head = lib1->getENDFhead();
  int ltt = head.l2;   // 0: isotropic, 1: Legendre, 2: tabulated, 3: Leg + table

  int idx1 = 0;
  /*** when high energy part is all tabulated */
  if(ltt == 2){
    idx1 ++;
    cont = lib1->rdata[idx1++];
    int ne1 = cont.n2;

    /*** find energy point that is higher than Emax */
    int i0 = 0;
    for(int i=0 ; i<ne1 ; i++){
      if(lib1->rdata[idx1].c2 > emax){
        i0 = idx1;
        cont.n2 = ne1 - i;
        break;
      }
      idx1++;
    }

    int idat[2];
    idat[0] = cont.n2;
    idat[1] = 2;

    ENDFPackTAB21(cont,&lib1->rdata[i0],idat,&lib1->iptr[i0],&lib1->xptr[i0],lib);
  }
  else if(ltt == 3){
    idx1++;
    idx1 += lib1->rdata[idx1].n2 + 1; // skip low energy part
    cont = lib1->rdata[idx1];
    int i0 = idx1;
    int i1 = idx1 + 1;

    double etest = lib1->rdata[i1].c2;
    if(etest < emax){
      cerr << "low energy boundary " << emax << " higher than the low-side of the high energy region, " << etest << endl;
    }

    ENDFPackTAB21(cont,&lib1->rdata[i1],lib1->iptr[i0],&lib1->iptr[i1],&lib1->xptr[i1],lib);
  }
  else{
    cerr << "high energy region should be tabulated data" << endl;
    exit(-1);
  }
}
