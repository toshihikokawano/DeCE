/******************************************************************************/
/**                                                                          **/
/**     ENDF I/O : ENDF MF Data Read / Write Library                         **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "endflib.h"

static string intscheme(int);


/**********************************************************/
/* Read / Write Wrapper Functions                         */
/**********************************************************/
int ENDFRead(ifstream *fp, ENDF *lib, const int mf, const int mt)
{
  int idx = 0;

  switch(mf){
  case  1: idx = ENDFReadMF1( fp,lib,mt); break;
  case  2: idx = ENDFReadMF2( fp,lib   ); break;
  case  3: idx = ENDFReadMF3( fp,lib,mt); break;
  case  4: idx = ENDFReadMF4( fp,lib,mt); break;
  case  5: idx = ENDFReadMF5( fp,lib,mt); break;
  case  6: idx = ENDFReadMF6( fp,lib,mt); break;
  case  8: idx = ENDFReadMF8( fp,lib,mt); break;
  case  9: idx = ENDFReadMF9( fp,lib,mt); break;
  case 10: idx = ENDFReadMF10(fp,lib,mt); break;
  case 12: idx = ENDFReadMF12(fp,lib,mt); break;
  case 13: idx = ENDFReadMF13(fp,lib,mt); break;
  case 14: idx = ENDFReadMF14(fp,lib,mt); break;
  case 15: idx = ENDFReadMF15(fp,lib,mt); break;
  case 31: idx = ENDFReadMF31(fp,lib,mt); break;
  case 32: idx = ENDFReadMF32(fp,lib   ); break;
  case 33: idx = ENDFReadMF33(fp,lib,mt); break;
  default:                                break;
  }

  return(idx);
}


void ENDFWrite(ENDF *lib)
{
  int mf = lib->getENDFmf();
  switch(mf){
  case  1: ENDFWriteMF1( lib); break;
  case  2: ENDFWriteMF2( lib); break;
  case  3: ENDFWriteMF3( lib); break;
  case  4: ENDFWriteMF4( lib); break;
  case  5: ENDFWriteMF5( lib); break;
  case  6: ENDFWriteMF6( lib); break;
  case  8: ENDFWriteMF8( lib); break;
  case  9: ENDFWriteMF9( lib); break;
  case 10: ENDFWriteMF10(lib); break;
  case 12: ENDFWriteMF12(lib); break;
  case 13: ENDFWriteMF13(lib); break;
  case 14: ENDFWriteMF14(lib); break;
  case 15: ENDFWriteMF15(lib); break;
  case 31: ENDFWriteMF31(lib); break;
  case 32: ENDFWriteMF32(lib); break;
  case 33: ENDFWriteMF33(lib); break;
  default:                     break;
  }
}


/**********************************************************/
/* MF1 :                                                  */
/*       General Information                              */
/**********************************************************/
int ENDFReadMF1(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,1,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    lnu  = head.l2;

  /*** for delayed neutron, read precursor decay constants */
  if(mt == 455){
    ENDFReadLIST(fp,lib);
    /*** given by polynomials */
    if(lnu==1) ENDFReadLIST(fp,lib);
    /*** tabulated */
    else       ENDFReadTAB1(fp,lib);
  }
  /*** prompt neutron or total neutron */
  else if( (mt == 452) || (mt == 456) ){
    if(lnu==1) ENDFReadLIST(fp,lib);
    else       ENDFReadTAB1(fp,lib);
  }
  /*** fission energy release */
  else if(mt == 458){
    ENDFReadLIST(fp,lib);
  }
  /*** delayed photon data */
  else if(mt == 460){
    int lo = head.l1;
    int ng = head.n1;
    if(lo == 1){
      for(int ing=0 ; ing<ng ; ing++) ENDFReadTAB1(fp,lib);
    }
    else ENDFReadLIST(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF1(ENDF *lib)
{
  int    mt   = lib->getENDFmt();
  Record head = lib->getENDFhead();
  int    lnu  = head.l2;

  ENDFWriteHEAD(lib);

  if(mt==455){
    ENDFWriteLIST(lib);
    if(lnu==1)  ENDFWriteLIST(lib);
    else        ENDFWriteTAB1(lib);
  }
  else if( (mt == 452) || (mt == 456) ){
    if(lnu==1)  ENDFWriteLIST(lib);
    else        ENDFWriteTAB1(lib);
  }
  else if(mt == 458){
    ENDFWriteLIST(lib);
  }
  else if(mt == 460){
    int lo = head.l1;
    int ng = head.n1;
    if(lo == 1){
      for(int ing=0 ; ing<ng ; ing++) ENDFWriteTAB1(lib);
    }
    else ENDFWriteLIST(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF2 :                                                  */
/*       Resonance Parameters                             */
/**********************************************************/
int ENDFReadMF2(ifstream *fp, ENDF *lib)
{
  if( (ENDFSeekHead(fp,lib,2,151))< 0 ) return(-1);

  /*** second card */
  lib->resetPOS();
  Record cont = ENDFReadCONT(fp,lib);
  int lfw = cont.l2;
  int ner = cont.n1;

  /*** for each energy range */
  for(int i=0 ; i<ner ; i++){

    /*** Erange, LRU, LRF */
    cont = ENDFReadCONT(fp,lib);
    int lru  = cont.l1;
    int lrf  = cont.l2;
    int nro  = cont.n1;

    /*** energy dependent scattering radius */
    if(nro == 1) ENDFReadTAB1(fp,lib);

    /*** Scattering radius only */
    if(lru == 0) ENDFReadCONT(fp,lib);

    /*** Resolved Resonance Region */
    else if (lru == 1){
      /*** SLBW, MLBW, Reich-Moore */
      if( (lrf == 1) || (lrf == 2)  || (lrf == 3) ){
        cont = ENDFReadCONT(fp,lib);
        int nls  = cont.n1;
        for(int inls=0 ; inls<nls ; inls++) ENDFReadLIST(fp,lib);
      }
      /*** R-Matrix Limited */
      else if(lrf == 7){
        cont = ENDFReadCONT(fp,lib);
        int krm = cont.l2;
        int njs = cont.n1;
        if( (krm == 1) || (krm == 2) || (krm == 3) ) ENDFReadLIST(fp,lib);

        for(int injs=0 ; injs<njs ; injs++){
          ENDFReadLIST(fp,lib);
          ENDFReadLIST(fp,lib);
        }
      }
      else{
        cout << "LRU = " << lru << "  LRF = " << lrf
             << "  has not been implemented" << endl;
      }
    }

    /*** Unresolved Resonance Region */
    else if(lru == 2){
      /*** all URR parameters are energy independent */
      if(lrf == 1){
        if(lfw == 0){ // Case A
          cont = ENDFReadCONT(fp,lib);
          int nls  = cont.n1;
          for(int inls=0 ; inls<nls ; inls++) ENDFReadLIST(fp,lib);
        }
        else{ // case B
          cont = ENDFReadLIST(fp,lib);
          int nls  = cont.n2;
          for(int inls=0 ; inls<nls ; inls++){
            cont = ENDFReadCONT(fp,lib);
            int njs  = cont.n1;
            for(int injs=0 ; injs<njs ; injs++) ENDFReadLIST(fp,lib);
          }
        }
      }
      /*** all URR parameters given */
      if(lrf == 2){ // case C
        cont = ENDFReadCONT(fp,lib);
        int nls  = cont.n1;
        for(int inls=0 ; inls<nls ; inls++){
          cont = ENDFReadCONT(fp,lib);
          int njs  = cont.n1;
          for(int injs=0 ; injs<njs ; injs++) ENDFReadLIST(fp,lib);
        }
      }
    }
  }

  return(lib->getPOS());
}


void ENDFWriteMF2(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record cont = ENDFWriteCONT(lib);
  int lfw = cont.l2;
  int ner = cont.n1;

  for(int i=0 ; i<ner ; i++){
    cont = ENDFWriteCONT(lib);
    int lru  = cont.l1;
    int lrf  = cont.l2;
    int nro  = cont.n1;

    if(nro == 1) ENDFWriteTAB1(lib);

    if(lru == 0) ENDFWriteCONT(lib);
    else if (lru == 1){
      if( (lrf == 1) || (lrf == 2)  || (lrf == 3) ){
        cont = ENDFWriteCONT(lib);
        int nls  = cont.n1;
        for(int inls=0 ; inls<nls ; inls++) ENDFWriteLIST(lib);
      }
      else if(lrf == 7){
        cont = ENDFWriteCONT(lib);
        int krm = cont.l2;
        int njs = cont.n1;
        if( (krm == 1) || (krm == 2) || (krm == 3) ) ENDFWriteLIST(lib);

        for(int injs=0 ; injs<njs ; injs++){
          ENDFWriteLIST(lib);
          ENDFWriteLIST(lib);
        }
      }
    }
    else if(lru == 2){
      if(lrf == 1){
        if(lfw == 0){
          cont = ENDFWriteCONT(lib);
          int nls  = cont.n1;
          for(int inls=0 ; inls<nls ; inls++) ENDFWriteLIST(lib);
        }
        else{
          cont = ENDFWriteLIST(lib);
          int nls  = cont.n2;
          for(int inls=0 ; inls<nls ; inls++){
            cont = ENDFWriteCONT(lib);
            int njs  = cont.n1;
            for(int injs=0 ; injs<njs ; injs++) ENDFWriteLIST(lib);
          }
        }
      }
      if(lrf == 2){
        cont = ENDFWriteCONT(lib);
        int nls  = cont.n1;
        for(int inls=0 ; inls<nls ; inls++){
          cont = ENDFWriteCONT(lib);
          int njs  = cont.n1;
          for(int injs=0 ; injs<njs ; injs++) ENDFWriteLIST(lib);
        }
      }
    }
  }
}


/**********************************************************/
/* MF3 :                                                  */
/*       Neutron Cross Sections                           */
/**********************************************************/
int ENDFReadMF3(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,3,mt))< 0 ) return(-1);

  lib->resetPOS();
  ENDFReadTAB1(fp,lib);

  return(lib->getPOS());
}


void ENDFWriteMF3(ENDF *lib)
{
  ENDFWriteHEAD(lib);
  ENDFWriteTAB1(lib);
  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF4 :                                                  */
/*       Angular Distributions of Secondary Particles     */
/**********************************************************/
int ENDFReadMF4(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,4,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    lvt  = head.l1; // transformation matrix, for backward compatibility
  int    ltt  = head.l2; // representation, Legendre or tabulated
  int    li   = 0;       // isotropic flag

  if(lvt == 1) ENDFReadLIST(fp,lib);
  else{
    Record cont = ENDFReadCONT(fp,lib);
    li = cont.l1;
  }

  if(li == 0){
    /*** Legendre coefficients, LTT=1 */
    if(ltt == 1){
      ENDFReadTAB2(fp,lib);
    }
    /*** Tabulated probabilities, LTT=2 */
    else if(ltt == 2){
      ENDFReadTAB21(fp,lib);
    }
    /*** Angular dstribution over two energy range, LTT=3 */
    else if(ltt == 3){
      ENDFReadTAB2(fp,lib);
      ENDFReadTAB21(fp,lib);
    }
  }
  else{
    /*** isotropic angular distribution */
  }

  return(lib->getPOS());
}


void ENDFWriteMF4(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int lvt = head.l1;
  int ltt = head.l2;
  int li  = lib->rdata[0].l1;

  ENDFWriteHEAD(lib);
  if(lvt == 1) ENDFWriteLIST(lib);
  else         ENDFWriteCONT(lib);

  if(li == 0){
    if(ltt == 1){
      ENDFWriteTAB2(lib);
    }
    else if(ltt == 2){
      ENDFWriteTAB21(lib);
    }
    else if(ltt == 3){
      ENDFWriteTAB2(lib);
      ENDFWriteTAB21(lib);
    }
  }
  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF5 :                                                  */
/*       Energy Distributions of Secondary Particles      */
/**********************************************************/
int ENDFReadMF5(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,5,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nk   = head.n1;  // number of subsections

  for(int n=0 ; n<nk ; n++){
    /*** (E,pk(E)) for each subsection */
    Record cont = ENDFReadTAB1(fp,lib);
    int lf = cont.l2;

    /*** FL=1, arbitrary tabulated */
    if(lf == 1)
      ENDFReadTAB21(fp,lib);
    /*** LF=5 (general evaporation), 11 (Watt) */
    else if( (lf == 5) || (lf == 11) ){
      ENDFReadTAB1(fp,lib);
      ENDFReadTAB1(fp,lib);
    }
    /*** LF=7 (Maxwellian), 9 (Evaporation), 12 (Madland-Nix) */
    else
      ENDFReadTAB1(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF5(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nk   = head.n1;

  for(int n=0 ; n<nk ; n++){
    Record cont = ENDFWriteTAB1(lib);
    int lf = cont.l2;

    if(lf == 1)
      ENDFWriteTAB21(lib);
    else if( (lf == 5) || (lf == 11) ){
      ENDFWriteTAB1(lib);
      ENDFWriteTAB1(lib);
    }
    else
      ENDFWriteTAB1(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF6 :                                                  */
/*       Product Energy - Angle Distributions             */
/**********************************************************/
int ENDFReadMF6(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,6,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nk   = head.n1;  // number of subsections

  for(int n=0 ; n<nk ; n++){
    Record cont = ENDFReadTAB1(fp,lib);
    int law = cont.l2;

    /*** LAW=0, unknown distribution */
    if(law == 0) continue;
    /*** LAW=1, continuum energy-angle distribution */
    else if(law == 1) ENDFReadTAB2(fp,lib);
    /*** LAW=2, discrete two-body scattering */
    else if(law == 2) ENDFReadTAB2(fp,lib);
    /*** LAW=3, isotropic discrete emission */
    else if(law == 3) continue;
    /*** LAW=4, discrete two-body recoils */
    else if(law == 4) continue;
    /*** LAW=5, charged-particle elastic scattering */
    else if(law == 5) ENDFReadTAB2(fp,lib);
    /*** LAW=6, N-body phase-space discribution */
    else if(law == 6) ENDFReadCONT(fp,lib);
    /*** LAW=7, laboratory angle-energy law */
    else if(law == 7) ENDFReadTAB22(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF6(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nk   = head.n1;

  for(int n=0 ; n<nk ; n++){
    Record cont = ENDFWriteTAB1(lib);
    int law = cont.l2;

    if( (law == 0) || (law == 3) || (law == 4) ) continue;
    else if( (law == 1) || (law == 2) || (law == 5) ) ENDFWriteTAB2(lib);
    else if(law == 6) ENDFWriteCONT(lib);
    else if(law == 7) ENDFWriteTAB22(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF7 : not implemented                                  */
/*       Thermal Neutron Scattering Law Data              */
/**********************************************************/


/**********************************************************/
/* MF8 :                                                  */
/*       Radioactive Decay and Fission Product Data       */
/**********************************************************/
int ENDFReadMF8(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,8,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();

  if( (mt == 454) || (mt == 459) ){
    int  le1  = head.l1;
    for(int l=0 ; l<le1 ; l++) ENDFReadLIST(fp,lib);
  }
  else if( mt == 457 ){
    int  nsp  = head.n2;  // number of radiation types
    ENDFReadLIST(fp,lib);
    ENDFReadLIST(fp,lib);
    for(int n=0 ; n<nsp ; n++){
      Record cont0 = ENDFReadLIST(fp,lib);
      int    lcon  = cont0.l1; // continuum flag

      if( (lcon == 0) || (lcon == 2) ){
        int  ner  = cont0.n2; // total number of tabulated discrete energies
        for(int k=0 ; k<ner ; k++) ENDFReadLIST(fp,lib);
      }
      if(lcon != 0){
        Record cont1 = ENDFReadTAB1(fp,lib);
        int    lcov  = cont1.l2; // covariance flag
        if(lcov == 1) ENDFReadLIST(fp,lib);
      }
    }
  }
  else{
    int  ns = head.n1;  // number of states
    int  no = head.n2;  // =0 complete decay chain
                        // =1 decay chain given in MT457
    if(no == 0){
      for(int n=0 ; n<ns ; n++) ENDFReadLIST(fp,lib);
    }
    else if(no == 1){
      for(int n=0 ; n<ns ; n++) ENDFReadCONT(fp,lib);
    }
  }

  return(lib->getPOS());
}


void ENDFWriteMF8(ENDF *lib)
{
  ENDFWriteHEAD(lib);
  Record head = lib->getENDFhead();

  int mt = lib->getENDFmt();
  if( (mt == 454) || (mt == 459) ){
    int  le1  = head.l1;
    for(int l=0 ; l<le1 ; l++) ENDFWriteLIST(lib);
  }
  else if( mt == 457 ){
    int  nsp  = head.n2;  // number of radiation types
    ENDFWriteLIST(lib);
    ENDFWriteLIST(lib);
    for(int n=0 ; n<nsp ; n++){
      Record cont0 = ENDFWriteLIST(lib);
      int    lcon  = cont0.l1; // continuum flag

      if( (lcon == 0) || (lcon == 2) ){
        int  ner  = cont0.n2; // total number of tabulated discrete energies
        for(int k=0 ; k<ner ; k++) ENDFWriteLIST(lib);
      }
      if(lcon != 0){
        Record cont1 = ENDFWriteTAB1(lib);
        int    lcov  = cont1.l2; // covariance flag
        if(lcov == 1) ENDFWriteLIST(lib);
      }
    }
  }
  else{
    int  ns = head.n1;  // number of states
    int  no = head.n2;  // =0 complete decay chain
                        // =1 decay chain given in MT457
    if(no == 0){
      for(int n=0 ; n<ns ; n++) ENDFWriteLIST(lib);
    }
    else if(no == 1){
      for(int n=0 ; n<ns ; n++) ENDFWriteCONT(lib);
    }
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF9 :                                                  */
/*       Multiplicities for Production of                 */
/*       Radioactive Nuclides                             */
/**********************************************************/
int ENDFReadMF9(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,9,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int ns = head.n1;  // number of final states

  for(int n=0 ; n<ns ; n++) ENDFReadTAB1(fp,lib);

  return(lib->getPOS());
}


void ENDFWriteMF9(ENDF *lib)
{
  ENDFWriteHEAD(lib);
  Record head = lib->getENDFhead();
  int ns = head.n1;

  for(int n=0 ; n<ns ; n++) ENDFWriteTAB1(lib);
  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF10 :                                                 */
/*        Cross Sections for Production of                */
/*        Radioactive Nuclides                            */
/**********************************************************/
int ENDFReadMF10(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,10,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    ns   = head.n1;   // number of subsections

  for(int n=0 ; n<ns ; n++) ENDFReadTAB1(fp,lib);

  return(lib->getPOS());
}


void ENDFWriteMF10(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    ns   = head.n1;   // number of subsections

  for(int n=0 ; n<ns ; n++) ENDFWriteTAB1(lib);

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF12 :                                                 */
/*        Photon Production Multiplicities                */
/*        and Transition probability arrays               */
/**********************************************************/
int ENDFReadMF12(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,12,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    l0   = head.l1;  // option

  /*** L0=1: multiplicities */
  if(l0 == 1){
    int nk  = head.n1;  // number of subsections
    if(nk > 1) ENDFReadTAB1(fp,lib);
    for(int n=0 ; n<nk ; n++) ENDFReadTAB1(fp,lib);
  }
  /*** L0=2: transition probability array */
  else if(l0 == 2) ENDFReadLIST(fp,lib);

  return(lib->getPOS());
}


void ENDFWriteMF12(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    l0   = head.l1;

  if(l0 == 1){
    int nk  = head.n1;
    if(nk > 1) ENDFWriteTAB1(lib);
    for(int n=0 ; n<nk ; n++) ENDFWriteTAB1(lib);
  }
  else if(l0 == 2) ENDFWriteLIST(lib);

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF13 :                                                 */
/*        Photon Production Cross Sections                */
/**********************************************************/
int ENDFReadMF13(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,13,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nk   = head.n1;  // number of subsections

  if(nk > 1) ENDFReadTAB1(fp,lib);
  for(int n=0 ; n<nk ; n++) ENDFReadTAB1(fp,lib);

  return(lib->getPOS());
}


void ENDFWriteMF13(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nk   = head.n1;

  if(nk > 1) ENDFWriteTAB1(lib);
  for(int n=0 ; n<nk ; n++) ENDFWriteTAB1(lib);

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF14 :                                                 */
/*        Photon Angular Distributions                    */
/**********************************************************/
int ENDFReadMF14(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,14,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    li   = head.l1;
  int    ltt  = head.l2;

  /*** anisotropic distribution */
  if(li == 0){
    int nk  = head.n1;
    int ni  = head.n2;
    for(int n=0 ; n<ni ; n++) ENDFReadCONT(fp,lib);
    for(int n=ni ; n<nk ; n++){
      if(ltt == 1)      ENDFReadTAB2(fp,lib);
      else if(ltt == 2) ENDFReadTAB21(fp,lib);
    }
  }
  /*** isotropic distribution */
  // nothing to do

  return(lib->getPOS());
}


void ENDFWriteMF14(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    li   = head.l1;
  int    ltt  = head.l2;

  if(li == 0){
    int nk  = head.n1;
    int ni  = head.n2;
    for(int n=0 ; n<ni ; n++) ENDFWriteCONT(lib);
    for(int n=ni ; n<nk ; n++){
      if(ltt == 1)      ENDFWriteTAB2(lib);
      else if(ltt == 2) ENDFWriteTAB21(lib);
    }
  }
  /*** isotropic distribution */
  // nothing to do

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF15 :                                                 */
/*        Continuous Photon Energy Spectra                */
/**********************************************************/
int ENDFReadMF15(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,15,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nc   = head.n1;  // number of subsections

  for(int n=0 ; n<nc ; n++){
    ENDFReadTAB1(fp,lib);
    ENDFReadTAB21(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF15(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nc   = head.n1;

  for(int n=0 ; n<nc ; n++){
    ENDFWriteTAB1(lib);
    ENDFWriteTAB21(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF31 :                                                 */
/*        Covariances of Number of Neutrons               */
/**********************************************************/
int ENDFReadMF31(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,31,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nl   = head.n2;  // number of subsections

  for(int n=0 ; n<nl ; n++){
    Record cont = ENDFReadCONT(fp,lib);
    int nc = cont.n1;
    int ni = cont.n2;

    /*** NC-type sub-subsections */
    for(int i=0 ; i<nc ; i++){
      ENDFReadCONT(fp,lib);
      ENDFReadLIST(fp,lib);
    }

    /*** NI-type sub-subsections */
    for(int i=0 ; i<ni ; i++) ENDFReadLIST(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF31(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nl   = head.n2;

  for(int n=0 ; n<nl ; n++){
    Record cont = ENDFWriteCONT(lib);
    int nc = cont.n1;
    int ni = cont.n2;
    for(int i=0 ; i<nc ; i++){
      ENDFWriteCONT(lib);
      ENDFWriteLIST(lib);
    }
    for(int i=0 ; i<ni ; i++) ENDFWriteLIST(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF32 :                                                 */
/*        Covariances of Resonance Parameters             */
/**********************************************************/
int ENDFReadMF32(ifstream *fp, ENDF *lib)
{
  if( (ENDFSeekHead(fp,lib,32,151))< 0 ) return(-1);

  /*** second card */
  lib->resetPOS();
  Record cont = ENDFReadCONT(fp,lib);
  int ner  = cont.n1;

  /*** for each energy range */
  for(int i=0 ; i<ner ; i++){

    cont = ENDFReadCONT(fp,lib);
    int lru  = cont.l1;
    int nro  = cont.n1;

    if(nro !=0 ){
      /*** Energy Dependent Scattering Radius */
      cont = ENDFReadCONT(fp,lib);
      int ni = cont.n2;

      /*** NI-type sub-subsections */
      for(int i=0 ; i<ni ; i++) ENDFReadLIST(fp,lib);
    }

    /*** SPI, AP, NLS */
    cont = ENDFReadCONT(fp,lib);
    int lcomp = cont.l2;
    int nls   = cont.n1;

    /*** Compatible Resolved Resonance Subsection Format */
    if(lcomp == 0){
      for(int inls=0 ; inls<nls ; inls++) ENDFReadLIST(fp,lib);
    }

    /*** General Resolved Resonance Subsection Format */
    else{
      if(lru == 1){
        cont = ENDFReadCONT(fp,lib);
        int nsrs = cont.n1;  // within a resonance
        int nlrs = cont.n2;  // long-range
        for(int insrs=0 ; insrs<nsrs ; insrs++) ENDFReadLIST(fp,lib);
        for(int inlrs=0 ; inlrs<nlrs ; inlrs++) ENDFReadLIST(fp,lib);
      }
      else{
        for(int inls=0 ; inls<nls ; inls++) ENDFReadLIST(fp,lib);
        ENDFReadLIST(fp,lib);
      }
    }
  }

  return(lib->getPOS());
}


void ENDFWriteMF32(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record cont = ENDFWriteCONT(lib);
  int ner = cont.n1;

  for(int i=0 ; i<ner ; i++){
    cont = ENDFWriteCONT(lib);
    int lru  = cont.l1;
    int nro  = cont.n1;

    if(nro !=0 ){
      cont = ENDFWriteCONT(lib);
      int ni = cont.n2;
      for(int i=0 ; i<ni ; i++) ENDFWriteLIST(lib);
    }

    cont = ENDFWriteCONT(lib);
    int lcomp = cont.l2;
    int nls   = cont.n1;

    if(lcomp == 0){
      for(int inls=0 ; inls<nls ; inls++) ENDFWriteLIST(lib);
    }
    else{
      if(lru == 1){
        cont = ENDFWriteCONT(lib);
        int nsrs = cont.n1;
        int nlrs = cont.n2;

        for(int insrs=0 ; insrs<nsrs ; insrs++) ENDFWriteLIST(lib);
        for(int inlrs=0 ; inlrs<nlrs ; inlrs++) ENDFWriteLIST(lib);
      }
      else{
        for(int inls=0 ; inls<nls ; inls++) ENDFWriteLIST(lib);
        ENDFWriteLIST(lib);
      }
    }
  }
}


/**********************************************************/
/* MF33 :                                                 */
/*        Covariances of Neutron Cross Sections           */
/**********************************************************/
int ENDFReadMF33(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,33,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nl   = head.n2;  // number of subsections

  for(int n=0 ; n<nl ; n++){
    Record cont = ENDFReadCONT(fp,lib);
    int nc = cont.n1;
    int ni = cont.n2;

    /*** NC-type sub-subsections */
    for(int i=0 ; i<nc ; i++){
      ENDFReadCONT(fp,lib);
      ENDFReadLIST(fp,lib);
    }

    /*** NI-type sub-subsections */
    for(int i=0 ; i<ni ; i++) ENDFReadLIST(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF33(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nl   = head.n2;

  for(int n=0 ; n<nl ; n++){
    Record cont = ENDFWriteCONT(lib);
    int nc = cont.n1;
    int ni = cont.n2;
    for(int i=0 ; i<nc ; i++){
      ENDFWriteCONT(lib);
      ENDFWriteLIST(lib);
    }
    for(int i=0 ; i<ni ; i++) ENDFWriteLIST(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/* MF34 :                                                 */
/*        Covariances of Scattering Angular Distributions */
/**********************************************************/
int ENDFReadMF34(ifstream *fp, ENDF *lib, const int mt)
{
  if( (ENDFSeekHead(fp,lib,34,mt))< 0 ) return(-1);

  lib->resetPOS();
  Record head = lib->getENDFhead();
  int    nmt1 = head.n2;  // number of subsections

  for(int n=0 ; n<nmt1 ; n++){
    ENDFReadCONT(fp,lib);
    Record cont = ENDFReadCONT(fp,lib);
    int ni = cont.n2;
    for(int i=0 ; i<ni ; i++) ENDFReadLIST(fp,lib);
  }

  return(lib->getPOS());
}


void ENDFWriteMF34(ENDF *lib)
{
  ENDFWriteHEAD(lib);

  Record head = lib->getENDFhead();
  int    nmt  = head.n2;

  for(int n=0 ; n<nmt ; n++){
    ENDFWriteCONT(lib);
    Record cont = ENDFWriteCONT(lib);
    int ni = cont.n2;
    for(int i=0 ; i<ni ; i++) ENDFWriteLIST(lib);
  }

  ENDFWriteSEND(lib);
}


/**********************************************************/
/*      Print LIST Data                                   */
/**********************************************************/
void ENDFPrintLIST(ENDF *lib, const int idx)
{
  cout.setf(ios::scientific, ios::floatfield);

  int n = lib->rdata[idx].n1;
  for(int i=0 ; i<n ; i++){
    cout << setw(11) << i;
    cout << setprecision(6) << setw(13) << lib->xptr[idx][i];
    cout << endl;
  }
}


/**********************************************************/
/*      Print 1-Dimensional Data                          */
/**********************************************************/
void ENDFPrint1Dim(ENDF *lib, const int idx)
{
  int nr = lib->rdata[idx].n1;
  cout << "#           NR" << setw(14) << nr << "  number of interpolation range" << endl;

  int i=0;
  for(int ir=0 ; ir<nr ; ir++){
    cout << "#           NP" << setw(14) << lib->iptr[idx][2*ir];
    cout << setw(9) << intscheme(lib->iptr[idx][2*ir+1]).c_str();
    cout << "  interpolation" << endl;

    cout.setf(ios::scientific, ios::floatfield);

    for(int ip=i ; ip<lib->iptr[idx][2*ir] ; ip++){
      cout << setprecision(6) << setw(14) << lib->xptr[idx][2*ip  ];
      cout << setprecision(6) << setw(14) << lib->xptr[idx][2*ip+1];
      cout << endl;
    }
    i = lib->idata[2*ir]-1;
    if(ir<nr-1){
      cout << endl;
      cout << endl;
    }
  }
  cout << endl;
  cout << endl;
}


inline string intscheme(int i)
{
  string interpol;
  switch(i){
  case  1: interpol = "constant"; break;
  case  2: interpol = " lin-lin"; break;
  case  3: interpol = " lin-log"; break;
  case  4: interpol = " log-lin"; break;
  case  5: interpol = " log-log"; break;
  default: interpol = "        "; break;
  }
  return (interpol);
}


