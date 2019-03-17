/******************************************************************************/
/**     DeCE MANIPULATE MF6                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "dece.h"
#include "terminate.h"
#include "masstable.h"


/**********************************************************/
/*      Correct Boundary Energy in MF6                    */
/**********************************************************/
void DeceBoundCorrect(ENDFDict *dict, ENDF *lib[], const int mt)
{
  int k3 = dict->getID(3,mt);
  int k6 = dict->getID(6,mt);

  if(k3 < 0) TerminateCode("MT number in MF3 not found",mt);
  if(k6 < 0) TerminateCode("MT number in MF6 not found",mt);

  double x0   = lib[k3]->xdata[0];
  Record head = lib[k6]->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;

  for(int n=0 ; n<nk ; n++){
    int law = lib[k6]->rdata[idx].l2;
    if(law != 1) return;

    lib[k6]->xptr[idx][0] = x0;
    idx++;

    int np = lib[k6]->rdata[idx].n2;
    idx++;

    for(int n=0 ; n<np ; n++){
      if(n == 0){
        Record r = lib[k6]->rdata[idx];
        lib[k6]->rdata[idx].setRecord(r.c1, x0, r.l1, r.l2, r.n1, r.n2);
      }
      idx++;
    }
  }
}


/**********************************************************/
/*      Duplicate Highest Points                          */
/**********************************************************/
void DeceDuplicatePoint(ENDFDict *dict, ENDF *lib0[], const int mt, double x)
{
  int k = dict->getID(6,mt);
  if(k < 0) TerminateCode("MT number in MF6 not found",mt);

  ENDF lib1(L);
  const int MAXDAT = 100;
  double **xdat;
  Record *cont;
  int    idat[MAXDAT];
  double zdat[MAXDAT];

  xdat = new double * [MAXDAT];
  cont = new Record [MAXDAT];

  Record head = lib0[k]->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;

  lib1.setENDFmat( lib0[k]->getENDFmat() );
  lib1.setENDFmf(6);
  lib1.setENDFmt(mt);
  lib1.setENDFhead(head);

  for(int ik=0 ; ik<nk ; ik++){
    int law  = lib0[k]->rdata[idx].l2;
    int nr1  = lib0[k]->rdata[idx].n1;
    int np1  = lib0[k]->rdata[idx].n2;

    Record ctab1 = lib0[k]->rdata[idx];
    for(int i=0 ; i<nr1*2 ; i++){
      idat[i] = lib0[k]->iptr[idx][i];
    }
    for(int i=0 ; i<np1 ; i++){
       zdat[i*2  ] = lib0[k]->xptr[idx][i*2  ];
       zdat[i*2+1] = lib0[k]->xptr[idx][i*2+1];
    }
    /*** add one point in TAB1 */
    ctab1.n2 = np1 + 1;
    idat[0]  = ctab1.n2;
    zdat[np1*2  ] = x;
    zdat[np1*2+1] = zdat[np1*2-1];
    ENDFPackTAB1(ctab1, idat, zdat, &lib1);

    /*** increment index */
    idx++;
    if( (law == 1) || (law == 2) || (law == 5) ){
      int    nr2  = lib0[k]->rdata[idx].n1;
      int    np2  = lib0[k]->rdata[idx].n2;

      Record ctab2 = lib0[k]->rdata[idx];
      for(int i=0 ; i<nr2*2 ; i++){
        idat[i] = lib0[k]->iptr[idx][i];
      }
      idx++;
 
      for(int i0=0 ; i0<np2 ; i0++){
        cont[i0] = lib0[k]->rdata[idx];
        xdat[i0] = lib0[k]->xptr[idx];
        idx++;
      }

      /*** repleat the last LIST */
      ctab2.n2 = np2 + 1;
      idat[0]  = ctab2.n2;
      cont[np2] = lib0[k]->rdata[idx-1];
      cont[np2].c2 = x;
      xdat[np2] = lib0[k]->xptr[idx-1];

      ENDFPackTAB2(ctab2, cont, idat, xdat, &lib1);
    }
  }

  ENDFLibCopy(&lib1,lib0[k]);
//ENDFWrite(lib0[k]);

  delete [] xdat;
  delete [] cont;
}


/**********************************************************/
/*      Generate Yield x Cross Section in MF=3            */
/**********************************************************/
void DeceGenProdCS(ENDFDict *dict, ENDF *lib[], const int mt1, const int zap1)
{
  int mt0 = 5;

  int k3 = dict->getID(3,mt0);
  int k6 = dict->getID(6,mt0);

  if(k3 < 0) TerminateCode("MT number in MF3 not found",mt0);
  if(k6 < 0) TerminateCode("MT number in MF6 not found",mt0);

  double *xdat;
  xdat = new double [MAX_DBLDATA];

  Record head = lib[k6]->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;
  int    ndat = 0;

  /*** for each sub block, make lib for (E,yield) */
  for(int ik=0 ; ik<nk ; ik++){
    Record cont = lib[k6]->rdata[idx];
    int    zap  = (int)cont.c1;
    int    law  = cont.l2;
    int    np   = cont.n2;

    if(zap == zap1){
      for(int i=0 ; i<np ; i++){
        xdat[i*2  ]  = lib[k6]->xptr[idx][i*2  ];
        xdat[i*2+1] = ENDFInterpolation(lib[k3],xdat[i*2  ],false,0) * lib[k6]->xptr[idx][i*2+1];
      }
      ndat = np;
      break;
    }

    /*** increment index */
    idx++;
    if( (law == 1) || (law == 2) || (law == 5) ){
      int ne = lib[k6]->rdata[idx].n2; idx++;
      for(int i0=0 ; i0<ne ; i0++) idx++;
    }
    else if(law == 6) idx++;
    else if(law == 7){
      int ne = lib[k6]->rdata[idx].n2; idx++;
      for(int ine=0 ; ine<ne ; ine++)  idx += lib[k6]->rdata[idx].n2 + 1;
    }
  }


  int mat = lib[k3]->getENDFmat();

  k3 = dict->getID(3,mt1);

  lib[k3]->setENDFmat(mat);
  lib[k3]->setENDFmf(3);
  lib[k3]->setENDFmt(mt1);
  lib[k3]->setENDFhead(head);

  Record cont;
  int    idat[2];
  int    za = (int)dict->getZA();

  double qm = qvalue(dict->getProj(),za,mt1);
  cont.setRecord(qm,qm,0,0,1,ndat);

  idat[0] = ndat;
  idat[1] = 2;
  ENDFPackTAB1(cont,idat,xdat,lib[k3]);

  ENDFWriteHEAD(lib[k3]);
  ENDFWriteTAB1(lib[k3]);
  ENDFWriteSEND(lib[k3]);

  delete [] xdat;
}
