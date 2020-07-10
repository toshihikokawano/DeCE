/******************************************************************************/
/**     DeCE TABLE for MF33                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"

static int  DeceTableMF33NC (int, const int, ENDF *);
static int  DeceTableMF33NI (int, const int, ENDF *);
static int  DeceTableMF33NILB0 (int, ENDF *);
static int  DeceTableMF33NILB5 (int, ENDF *);
static int  DeceTableMF33NILB6 (int, ENDF *);
static int  DeceTableMF33NILB8 (int, ENDF *);


/**********************************************************/
/*      Process MF=33                                     */
/**********************************************************/
void DeceTableMF33(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nl   = head.n2;
  int    idx  = 0;

  cout << "# Cross Section Covanriance Matrix" << endl;
  cout << "#           NK" << setw(14) << nl << "  number of subsections" << endl;
  cout << endl;

  for(int n=0 ; n<nl ; n++){
    Record cont = lib->rdata[idx++];
    int mf1 = (int)cont.c1;
    int mt1 = cont.l2;
    int nc = cont.n1;
    int ni = cont.n2;

    cout << "#   Subsection" << setw(14) << n << endl;

    cout << "#          MF1" << setw(14) << mf1 << "  MF number for the second cross section" << endl;
    cout << "#          MT1" << setw(14) << mt1 << "  MT number for the second cross section" << endl;
    cout << "#           NC" << setw(14) << nc << "  derived covariance" << endl;
    cout << "#           NI" << setw(14) << ni << "  relative / absolute covariance"  << endl;

    /*** NC type */
    idx = DeceTableMF33NC(idx,nc,lib);

    /*** NI type */
    idx = DeceTableMF33NI(idx,ni,lib);

    cout << endl;
    cout << endl;
  }
}


/**********************************************************/
/*      NC-Type Covariance                                */
/**********************************************************/
int  DeceTableMF33NC(int idx, const int nc, ENDF *lib)
{
  for(int i=0 ; i<nc ; i++){
    int    lty = lib->rdata[idx++].l2;

    double e1  = lib->rdata[idx].c1;
    double e2  = lib->rdata[idx].c2;
    int    nci = lib->rdata[idx].n2;
    if(lty == 0){
      cout << "# NC-type     " << setw(14) << i << endl;
      cout << "#         Emin" << setw(14) << e1 << endl;
      cout << "#         Emax" << setw(14) << e2 << endl;
      cout << "#          NCI" << setw(14) << nci << endl;
      cout << "#       MT       Constant"<< endl;
      for(int j=0; j<nci ; j++){
        outVal((int)lib->xptr[idx][2*j+1]);
        outVal(     lib->xptr[idx][2*j  ]);
        cout << endl;
      }
    }
    idx++;
  }
  cout << endl;

  return(idx);
}


/**********************************************************/
/*      NI-Type Covariance                                */
/**********************************************************/
int  DeceTableMF33NI(int idx, const int ni, ENDF *lib)
{
  for(int i=0 ; i<ni ; i++){
    int lb = lib->rdata[idx].l2;

    cout << "# NI-type     " << setw(14) << i << endl;
    cout << "#           LB" << setw(14) << lb << endl;

    if( (0 <= lb) && (lb <= 4) ){
      DeceTableMF33NILB0(idx,lib);
    }

    else if(lb == 5){
      DeceTableMF33NILB5(idx,lib);
    }

    else if(lb == 6){
      DeceTableMF33NILB6(idx,lib);
    }

    else if( (lb == 8) || (lb== 9) ){
      DeceTableMF33NILB8(idx,lib);
    }

    idx++;
  }

  return(idx);
}


/**********************************************************/
/*      NI-Type, LB=0, 1, 2, 3, or 4                      */
/**********************************************************/
int DeceTableMF33NILB0(int idx, ENDF *lib)
{
  int lt = lib->rdata[idx].l1;
  int nt = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  cout << "#           NP" << setw(14) << np << "  total number of pairs {Ek,Fk}{El,Fl}" <<  endl;
  cout << "#           NT" << setw(14) << nt << "  NP x 2" << endl;
  cout << "#           LT" << setw(14) << lt << "  number of pairs in the second array {El,Fl}" << endl;

  cout << "# Ek            Fk" << endl;
  for(int i=0; i<np-lt; i++){
    outVal(         lib->xptr[idx][2*i  ]  );
    outVal(sqrt(abs(lib->xptr[idx][2*i+1])));
    cout << endl;
  }

  if(lt > 0){
    cout << "# El            Fl" << endl;
    for(int i=np-lt; i<np; i++){
      outVal(         lib->xptr[idx][2*i  ]  );
      outVal(sqrt(abs(lib->xptr[idx][2*i+1])));
      cout << endl;
    }
  }
  cout << endl;
  
  return(idx);
}


/**********************************************************/
/*      NI-Type, LB=5                                     */
/**********************************************************/
int DeceTableMF33NILB5(int idx, ENDF *lib)
{
  int ls = lib->rdata[idx].l1;
  int nt = lib->rdata[idx].n1;
  int ne = lib->rdata[idx].n2;

  cout << "#           NE" << setw(14) << ne << "  number of energy intervals" <<  endl;
  cout << "#           NT" << setw(14) << nt << "  total number of entries" << endl;
  cout << "#           LS" << setw(14) << ls << "  0: asymmetric / 1: symmetric matrix" << endl;

  double *eptr, *cptr;
  eptr = lib->xptr[idx];
  cptr = &lib->xptr[idx][ne];

  /*** symmetric covariance */
  if(ls == 1){
    for(int i=0; i<ne-1; i++){
      int ki = i+i*(ne-1)-i*(i+1)/2;

      outVal(eptr[i]);
      outVal(sqrt(abs(cptr[ki])));

      for(int j=0; j<=i; j++){
        int kj = j+j*(ne-1)-j*(j+1)/2;
        int k  = (j < i) ?  i+j*(ne-1)-j*(j+1)/2 : j+i*(ne-1)-i*(i+1)/2;
        int c  = 0;
        if(cptr[ki]*cptr[kj] != 0.0){
          c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
          if((i == j) && (c == 999)) c = 1000;
        }
        cout << setw(5) << c;
      }
      cout << endl;
    }
    outVal(eptr[ne-1]); outVal(0.0); cout << endl;
  }

  /*** asymmetric case */
  else{
    int k = 0;
    for(int i=0; i<ne-1; i++){
      outVal(eptr[i]);
      for(int j=0; j<ne-1; j++) outVal(cptr[k++]);
      outVal(0.0);
      cout << endl;
    }
    outVal(eptr[ne-1]);
    outVal(0.0);
    cout << endl;
  }

  cout << endl;

  return(idx);
}


/**********************************************************/
/*      NI-Type, LB=6                                     */
/**********************************************************/
int DeceTableMF33NILB6(int idx, ENDF *lib)
{
  int nt  = lib->rdata[idx].n1;
  int ner = lib->rdata[idx].n2;
  int nec = (nt - 1)/ner;

  cout << "#          NER" << setw(14) << ner << "  number of energies in row" <<  endl;
  cout << "#          NEC" << setw(14) << nec << "  number of energies in column" <<  endl;
  cout << "#           NT" << setw(14) << nt << "  total number of entries" << endl;

  double *eptr1, *eptr2, *cptr;
  eptr1 = &lib->xptr[idx][0];
  eptr2 = &lib->xptr[idx][ner];
  cptr = &lib->xptr[idx][ner+nec];

  cout << "#             ";
  for(int j=0; j<nec; j++) outVal(eptr2[j]);
  cout << endl;

  int k = 0;
  for(int i=0; i<ner-1; i++){
    outVal(eptr1[i]);
    for(int j=0; j<nec-1; j++) outVal(cptr[k++]);
    outVal(0.0);
    cout << endl;
  }
  outVal(eptr1[ner-1]);
  outVal(0.0);
  cout << endl;
  cout << endl;

  return(idx);
}


/**********************************************************/
/*      NI-Type, LB=8 or 9                                */
/**********************************************************/
int DeceTableMF33NILB8(int idx, ENDF *lib)
{
  int nt = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  cout << "#           NP" << setw(14) << np << "  total number of pairs {Ek,Fk}" <<  endl;
  cout << "#           NT" << setw(14) << nt << "  NP x 2" << endl;

  cout << "# Ek            Fk" << endl;
  for(int i=0; i<np; i++){
    outVal(         lib->xptr[idx][2*i  ]  );
    outVal(sqrt(abs(lib->xptr[idx][2*i+1])));
    cout << endl;
  }
  cout << endl;
  
  return(idx);
}


