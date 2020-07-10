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
/*      --------                                          */
/*      adjust the threshold energy given in MF6          */
/*      to be the same as it in MF3                       */
/**********************************************************/
void DeceBoundCorrect(ENDFDict *dict, ENDF *lib0[], const int mt)
{
  int k3 = dict->getID(3,mt);
  int k6 = dict->getID(6,mt);

  if(k3 < 0) TerminateCode("MT number in MF3 not found",mt);
  if(k6 < 0) TerminateCode("MT number in MF6 not found",mt);

  /*** allocate temporal arrays */
  int npmax = lib0[k6]->getXSIZE();
  int nrmax = lib0[k6]->getISIZE();
  int nbmax = lib0[k6]->getPOS();

  int *idat1, *idat2;
  double **xdat, *zdat;
  Record *cont;

  idat1 = new int [nrmax];
  idat2 = new int [nrmax];
  zdat = new double [npmax];
  xdat = new double * [nbmax];
  cont = new Record [nbmax];


  /*** get the energy range in MF3, the first and last elements in the array */
  int    np0 = lib0[k3]->rdata[0].n2;
  double x0  = lib0[k3]->xdata[0];
  double x1  = lib0[k3]->xdata[2*np0-2];

  message << "MF" << lib0[k3]->getENDFmf() << "MT" <<lib0[k3]->getENDFmt() << " has the energy range of [";
  message << setw(13) << setprecision(6) << x0 << ",";
  message << setw(13) << setprecision(6) << x1 << "]";
  Notice("DeceMod6:DeceBoundCorrect");


  ENDF lib1(L);

  Record head = lib0[k6]->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;

  lib1.setENDFmat( lib0[k6]->getENDFmat() );
  lib1.setENDFmf(6);
  lib1.setENDFmt(mt);
  lib1.setENDFhead(head);

  for(int n=0 ; n<nk ; n++){
    Record ctab1, ctab2;
    int nr1 = 0, np1 = 0, nr2 = 0 , np2 = 0;

    /*** first TAB1 in MF6 */
    ctab1 = lib0[k6]->rdata[idx];
    int zap = (int)ctab1.c1;
    int law = ctab1.l2;
    nr1 = ctab1.n1;
    np1 = ctab1.n2;

    /*** copy data in TAB1 */
    for(int i=0 ; i<nr1*2 ; i++){ idat1[i] = lib0[k6]->iptr[idx][i]; }
    for(int i=0 ; i<np1*2 ; i++){ zdat[i] = lib0[k6]->xptr[idx][i]; }
    idx++;

    if( (law == 1) || (law == 2) || (law == 5) ){
      ctab2 = lib0[k6]->rdata[idx];
      nr2 = ctab2.n1;
      np2 = ctab2.n2;

      /*** copy data in TAB21 */
      for(int i=0 ; i<nr2*2 ; i++){ idat2[i] = lib0[k6]->iptr[idx][i]; }
      idx++;
 
      for(int i0=0 ; i0<np2 ; i0++){
        cont[i0] = lib0[k6]->rdata[idx];
        xdat[i0] = lib0[k6]->xptr[idx]; // copy pointers only
        idx++;
      }
    }

    /*** energy range of MF6 */
    double z0 = zdat[0];
    double z1 = zdat[2*np1-2];

    /*** first, adjust the last point */
    int k1 = np1;
    if(z1 != x1){
      /*** when z1 < x1, duplicate the last point */
      if(z1 < x1){
        zdat[2*k1  ] = x1;
        zdat[2*k1+1] = zdat[2*k1-1];
        cont[k1] = cont[k1-1];
        cont[k1].c2 = x1;
        xdat[k1] = xdat[k1-1];
        k1 ++;
      }
      /*** when z1 > x1, truncate at x1 */
      else{
        for(k1=np1-1 ; k1>0 ; k1--){ if(zdat[k1*2] < x1) break; }
        k1 ++;
        zdat[2*k1  ] = x1;
        cont[k1].c2 = x1;
      }
      np1 = k1;

      message << "ZAP" << setw(8) << zap << " highest energy data changed into (";
      message << setw(13) << setprecision(6) << zdat[2*k1] << ",";
      message << setw(13) << setprecision(6) << zdat[2*k1+1] << ")";
      Notice("DeceMod6:DeceBoundCorrect");
    }

    /*** then adjust the first energy point */
    int k0 = 0;
    if(z0 != x0){
      /*** when z0 < x0, remove points those are less than x0 */
      if(z0 < x0){
        for(k0=1 ; k0<np1 ; k0++){ if(zdat[k0*2] > x0) break; }
        k0 --;
      }
      /*** when z0 > x0, move the first point to x0 */
      zdat[k0*2] = x0;
      cont[k0].c2 = x0;
      np1 -= k0;

      message << "ZAP" << setw(8) << zap << " lowest energy data replaced by (";
      message << setw(13) << setprecision(6) << zdat[k0*2] << ",";
      message << setw(13) << setprecision(6) << zdat[k0*2+1] << ")";
      Notice("DeceMod6:DeceBoundCorrect");
    }


    np2 = np1;

    /*** create new TAB1 */
    ctab1.n2 = np1;
    idat1[(nr1-1)*2] = np1;     // we assume there is only one energy-range, NR = 1
    ENDFPackTAB1(ctab1, idat1, &zdat[k0*2], &lib1);

    if(nr2 > 0){
      /*** create new TAB2 */
      ctab2.n2 = np2;
      idat2[(nr2-1)*2] = np2;
      ENDFPackTAB2(ctab2, &cont[k0], idat2, &xdat[k0], &lib1);
    }
  }

  ENDFLibCopy(&lib1,lib0[k6]);
//ENDFWrite(&lib1);

  delete [] idat1;
  delete [] idat2;
  delete [] zdat;
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


/**********************************************************/
/*      Produce Isotropic Angular Distribution            */
/**********************************************************/
void DeceIsotropicAngularDistribution(ENDFDict *dict, ENDF *lib[], const int mt)
{
  int k3 = dict->getID(3,mt);
  int k6 = dict->getID(6,mt);

  if(k3 < 0) TerminateCode("MT number in MF3 not found",mt);
  if(k6 < 0) TerminateCode("MT number in MF6 not found",mt);

  Record head = lib[k3]->getENDFhead();
  head.l1 = 0;
  head.l2 = 2; // LCT = 2 : CMS
  head.n1 = 1; // NK = 1
  head.n2 = 0;

  lib[k6]->setENDFhead(head);


  int    idat[2];
  double xdat[4];

  idat[0] = 2;
  idat[1] = 2;
  xdat[0] = 1.0;
  xdat[1] = 1.0;
  xdat[2] = 1.0;
  xdat[3] = 1.0;

  Record cont(0.0,0.0,0,3,1,2);
  ENDFPackTAB1(cont,idat,xdat,lib[k6]);

  ENDFWriteHEAD(lib[k6]);
  ENDFWriteTAB1(lib[k6]);
  ENDFWriteSEND(lib[k6]);

}
