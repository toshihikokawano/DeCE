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
void mergeMF4 (ENDF *, ENDF *);


/**********************************************************/
/*      Merge 2 Angular Distribution Data from Files      */
/**********************************************************/
int main(int argc, char *argv[])
{
  if(argc < 2){
    cout <<
    "usage\n"
    " % decemergemf4 ENDF_low_energy ENDF_high_energy\n" << endl;
    exit(-1);
  }

  ifstream fpin;
  string   libname0 = "", libname1 = "";
  ENDF     lib0, lib1;

  libname0 = argv[1];
  libname1 = argv[2];

  /*** read MF4 MT2 from low energy data */
  fpin.open(libname0.c_str());
  ENDFReadMF4(&fpin,&lib0,2);
  fpin.close();

  /*** read MF4 MT2 from high energy data */
  fpin.open(libname1.c_str());
  ENDFReadMF4(&fpin,&lib1,2);
  fpin.close();

  mergeMF4(&lib0,&lib1);

  return(0);
}


/**********************************************************/
/*      Tabulate Scattering Angular Distribution          */
/**********************************************************/
void mergeMF4(ENDF *lib0, ENDF *lib1)
{
  ENDF lib;
  Record head = lib0->getENDFhead();

  if(head.l2 != 1){
    cerr << "low energy region should be Legendre coefficients" << endl;
    exit(-1);
  }

  int idx0 = 0;
  Record cont = lib0->rdata[idx0++];

  /*** Make HEAD and CONT */
  head.l2 = 3; // Legendre parameters and tabulated data

  lib.setENDFmat(lib0->getENDFmat());
  lib.setENDFmf(4);
  lib.setENDFmt(2);
  lib.setENDFhead(head);
  ENDFPackCONT(cont,&lib);

  /*** low energy region, copy fist CONT record */
  cont = lib0->rdata[idx0++];
  int    ne0 = cont.n2;
  double emax = lib0->rdata[ne0+1].c2;

  /*** copy TAB2 */
  ENDFPackTAB2(cont,&lib0->rdata[idx0],lib0->iptr[idx0],&lib0->xptr[idx0],&lib);

  /*** high energy region */
  head = lib1->getENDFhead();
  int ltt = head.l2;   // 0: isotropic, 1: Legendre, 2: tabulated, 3: Leg + table

  int idx1 = 0;
  /*** when high energy part is all tabulated */
  if(ltt == 2){
    idx1 ++;
    cont = lib1->rdata[idx1++];
    int ne1 = cont.n2;

    int i0 = 0;
    int i1 = ne1;
    for(int i=0 ; i<ne1 ; i++){
      if(lib1->rdata[idx1].c2 > emax){
        i0 = idx1;
        i1 = ne1 - i;
        break;
      }
      idx1++;
    }
    cont.n2 = i1;

    int idat[2];
    idat[0] = i1;
    idat[1] =  2;

    ENDFPackTAB21(cont,&lib1->rdata[i0],idat,&lib1->iptr[i0],&lib1->xptr[i0],&lib);
  }
  else if(ltt == 3){
    idx1++;
    idx1 += lib1->rdata[idx1].n2 + 1;
    cont = lib1->rdata[idx1];
    int i0 = idx1;
    int i1 = idx1 + 1;

    double etest = lib1->rdata[i1].c2;
    if(etest < emax){
      cerr << "low energy boundary " << emax << " higher than the low-side of the high energy region, " << etest << endl;
    }

    ENDFPackTAB21(cont,&lib1->rdata[i1],lib1->iptr[i0],&lib1->iptr[i1],&lib1->xptr[i1],&lib);
  }
  else{
    cerr << "high energy region should be tabulated data" << endl;
    exit(-1);
  }

  ENDFWriteHEAD(&lib);
  ENDFWriteCONT(&lib);
  ENDFWriteTAB2(&lib);
  ENDFWriteTAB21(&lib);
  ENDFWriteSEND(&lib);
}
