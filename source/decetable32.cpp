/******************************************************************************/
/**     DeCE TABLE for MF32                                                  **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"
#include "decetable.h"

static int  DeceTableMF32RRR (int, int, ENDF *);
static int  DeceTableMF32LCOMP0 (int, int, ENDF *);
static int  DeceTableMF32LCOMP1LRF3 (int, int, ENDF *);
static int  DeceTableMF32LCOMP1LRF7 (int,      ENDF *);
static int  DeceTableMF32LCOMP2LRF3 (int, int, ENDF *);
static int  DeceTableMF32LCOMP2LRF7 (int, int, ENDF *);
static void DeceTableMF32PrintCorr (int, int, int, double *, double *, double *);
static int  DeceTableMF32URR (int, ENDF *);


/**********************************************************/
/*      Process MF=32                                     */
/**********************************************************/
void DeceTableMF32(ENDF *lib)
{
  int    idx  = 0;
  Record cont = lib->rdata[idx++];
  int    ner  = cont.n1;

  cout << "# Resonance Parameter Covariance Matrix" << endl;
  cout << "#          NER" << setw(14) << ner << "  number of energy ranges" << endl;
  cout << endl;

  for(int i=0 ; i<ner ; i++){
    cont = lib->rdata[idx++];
    double el = cont.c1;
    double eh = cont.c2;
    int lru = cont.l1;
    int lrf = cont.l2;
    int nro = cont.n1;

    cout << "#   Subsection" << setw(14) << i << endl;
    cout << "#         Emin"; outVal(el); cout << endl;
    cout << "#         Emax"; outVal(eh); cout << endl;
    cout << "#          LRU" << setw(14) << lru << endl;
    cout << "#          LRF" << setw(14) << lrf << endl;

    if( !(lrf == 1 || lrf == 2 || lrf == 3 || lrf == 7) ){
      message << "LRF = " << lrf << " not yet implemented";
      WarningMessage();
      return;
    }

    /*** resolved resonance range */
    if(lru == 1){
      /*** when NRO !=0, skip NI-type scattering radius covariance */
      if(nro !=0 ){
        cont = lib->rdata[idx++];
        int ni = cont.n2;
        for(int i=0 ; i<ni ; i++) idx ++;
      }
      idx = DeceTableMF32RRR(idx,lrf,lib);
    }
    /*** unresolved resonance range */
    else{
      idx = DeceTableMF32URR(idx,lib);
    }

    cout << endl;
    cout << endl;
  }
}


/**********************************************************/
/*      Resolved Resonance Parameter Covariance           */
/**********************************************************/
int DeceTableMF32RRR(int idx, int lrf, ENDF *lib)
{
  Record cont = lib->rdata[idx++];
  int lcomp = cont.l2;
  int nls   = cont.n1; int njs = nls;
  int isr   = cont.n2;

  cout << "#        LCOMP" << setw(14) << lcomp << "  compatible flag 0:ENDF/B-V, 1:general case, 2:compact format" << endl;
  cout << "#          ISR" << setw(14) << isr << "  0: no scattering radius uncertainty, 1: data given" << endl;

  /*** scattering radius uncertainty */
  if(isr != 0){
    double ap  = cont.c2;
    if( (lrf == 1) || (lrf == 2) ){
      double dap = lib->rdata[idx++].c2;
      cout << "# AP            DAP" << endl;
      outVal(ap);
      outVal(dap);
      cout << endl;
    }
    else if(lrf == 3){
      int mls = lib->rdata[idx].n1;
      for(int i=0 ; i<=mls ; i++){
        cout << "# AP            DAP" << endl;
        outVal(ap);
        outVal(lib->xptr[idx][i]);
        cout << endl;
      }
      idx++;
    }
    else{
      int njch = lib->rdata[idx].n1;
      cout << "# JS            DAP" << endl;
      for(int j=0 ; j<njs ; j++){
        outVal(j);
        for(int i=0 ; i<njch ; i++) outVal(lib->xptr[idx][i]);
        cout << endl;
      }
      idx++;
    }
  }

  /*** ENDF/B-V Compatible Resolved Resonance Format */
  if(lcomp == 0){
    idx = DeceTableMF32LCOMP0(idx,nls,lib);
  }
  /*** General Resolved Resonance Format */
  else if(lcomp == 1){
    if(lrf == 7) idx = DeceTableMF32LCOMP1LRF7(idx,lib);
    else         idx = DeceTableMF32LCOMP1LRF3(idx,lrf,lib);
  }
  else{
    if(lrf == 7) idx = DeceTableMF32LCOMP2LRF7(idx,njs,lib);
    else         idx = DeceTableMF32LCOMP2LRF3(idx,lrf,lib);
  }
  return(idx);
}


/**********************************************************/
/*      Resolved Resonance Parameter, LCOMP = 0           */
/**********************************************************/
int DeceTableMF32LCOMP0(int idx, int nls, ENDF *lib)
{
  double p[4], e[4], c[10];

  cout << "#          NLS" << setw(14) << nls << "  total number of orbital angular momentum" << endl;


  for(int inls=0 ; inls<nls ; inls++){
    Record cont = lib->rdata[idx];
    int nrs = cont.n2;
    cout << "#      L-value" << setw(14) << inls << endl;
    cout << "#          NRS" << setw(14) << nrs << "  number of resonances" << endl;

    for(int i=0 ; i<nrs ; i++){
      int j = 18*i;

      p[0] = lib->xptr[idx][j  ];
      p[1] = lib->xptr[idx][j+3];
      p[2] = lib->xptr[idx][j+4];
      p[3] = lib->xptr[idx][j+5];

      c[0] = lib->xptr[idx][j+6];
      c[1] = 0.0;
      c[2] = lib->xptr[idx][j+7];
      c[3] = 0.0;
      c[4] = lib->xptr[idx][j+8];
      c[5] = lib->xptr[idx][j+9];
      c[6] = 0.0;
      c[7] = lib->xptr[idx][j+10];
      c[8] = lib->xptr[idx][j+11];
      c[9] = lib->xptr[idx][j+12];

      e[0] = sqrt(c[0]);
      e[1] = sqrt(c[2]);
      e[2] = sqrt(c[5]);
      e[3] = sqrt(c[9]);

      cout << "#    Resonance" << setw(14) << i << endl;
      cout << "# Parameter    Uncertainty" << endl;

      for(int k1=0 ; k1<4 ; k1++){
        outVal(p[k1]);

        if(p[k1] != 0.0) outVal(e[k1]/abs(p[k1]));
        else outVal(0.0);
            
        for(int k2=0 ; k2<=k1 ; k2++){
          int k = 0;
          if( (e[k1]*e[k2] != 0.0) ) k = c[k1*(k1+1)/2 + k2] / e[k1] / e[k2] * 1000.0;
          if( (k1 == k2) && (k == 999) ) k = 1000;
          cout << setw(5) << k;
        }
        cout << endl;
      }
    }
    idx++;
  }

  return(idx);
}


/**********************************************************/
/*      Resolved Resonance  LCOMP = 1, LRF = 1,2,3        */
/**********************************************************/
int DeceTableMF32LCOMP1LRF3(int idx, int lrf, ENDF *lib)
{
  Record cont = lib->rdata[idx++];
  int nsrs = cont.n1;
  int nlrs = cont.n2;

  cout << "#         NSRS" << setw(14) << nsrs << "  short range covariance" << endl;
  cout << "#         NLRS" << setw(14) << nlrs << "  long range covariance" << endl;

  /*** short range covariance */
  for(int insrs=0 ; insrs<nsrs ; insrs++){
    cont = lib->rdata[idx];
    int mpar = cont.l1;
    int nrb  = cont.n2;

    cout << "# ShortRange  " << insrs << endl;
    cout << "#         MPAR" << setw(14) << mpar << "  number of parameters par resonance (ER, GN, GG, ...)" << endl;
    cout << "#          NRB" << setw(14) << nrb << "  number of resonances" << endl;

    int n = nrb * mpar;
    double *p = new double [n];

    /*** copy parameters */
    int k = 0;
    for(int i=0; i<nrb; i++){
      int i0 = 6*i;
      p[k++] = lib->xptr[idx][i0]; if(mpar == 1) continue; // ER
      if(lrf == 3){
        p[k++] = lib->xptr[idx][i0 + 2]; if(mpar == 2) continue; // GN
        p[k++] = lib->xptr[idx][i0 + 3]; if(mpar == 3) continue; // GG
        p[k++] = lib->xptr[idx][i0 + 4]; if(mpar == 4) continue; // GFA
        p[k++] = lib->xptr[idx][i0 + 5];                         // GFB
      }
      else{
        p[k++] = lib->xptr[idx][i0 + 3]; if(mpar == 2) continue; // GN
        p[k++] = lib->xptr[idx][i0 + 4]; if(mpar == 3) continue; // GG
        p[k++] = lib->xptr[idx][i0 + 5]; if(mpar == 4) continue; // GF
        p[k++] = lib->xptr[idx][i0 + 2] - (lib->xptr[idx][i0 + 3] + lib->xptr[idx][i0 + 4] + lib->xptr[idx][i0 + 5]); // GX
      }
    }

    cout << "#     Parameter     Uncertainty   Correlation" << endl;

    double *cptr;
    cptr = &lib->xptr[idx][6*nrb];

    /*** print correlation matrix */
    for(int i0=0; i0<nrb; i0++){
      for(int i1=0; i1<mpar; i1++){
        int i  = i0*mpar + i1;
        int ki = i + i*n - i*(i+1)/2;

        cout << setw(5) <<i;
        outVal(p[i]);
        if(p[i] != 0.0) outVal(sqrt(cptr[ki] / p[i] / p[i]));
        else outVal(0.0);

        for(int j0=0; j0<=i0; j0++){
          int j1max = (j0 == i0) ? i1 : mpar-1;
          for(int j1=0; j1<=j1max; j1++){
            int j  = j0*mpar + j1;
            int kj = j + j*n - j*(j+1)/2;

            int k  = (j < i) ?  i + j*n - j*(j+1)/2 : j + i*n - i*(i+1)/2;

            int c  = 0;
            if(cptr[ki]*cptr[kj] != 0.0){
              c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
              if((i == j) && (c == 999)) c = 1000;
            }
            cout << setw(5) << c;
          }
        }
        cout << endl;
      }
    }
    idx++;

    delete [] p;
  }

  return(idx);
}


/**********************************************************/
/*      Resolved Resonance  LCOMP = 1, LRF = 7            */
/**********************************************************/
int DeceTableMF32LCOMP1LRF7(int idx, ENDF *lib)
{
  Record cont = lib->rdata[idx++];
  int nsrs = cont.n1;

  cout << "#         NSRS" << setw(14) << nsrs << "  short range covariance" << endl;

  /*** short range covariance */
  for(int insrs=0 ; insrs<nsrs ; insrs++){

    cont = lib->rdata[idx++];
    int njsx = cont.l1;

    cout << "# ShortRange  " << insrs << endl;
    cout << "#         NJSX" << setw(14) << njsx << "  number of J-Pi groups"  << endl;

    /*** count total number of parameters */
    int id1 = idx;
    int n = 0;
    for(int j=0; j<njsx ; j++){
      cont = lib->rdata[id1++];
      int nch = cont.l1;
      int nrb = cont.l2;

      cout << "#          NCH" << setw(14) << nch << "  number of channels" << endl;
      cout << "#          NRB" << setw(14) << nrb << "  number of blocks" << endl;

      n += nrb * (nch + 1);
    }

    cout << "#        NPARB" << setw(14) << n << "  total parameters in this group" << endl;

    double *p = new double [n];

    /*** copy parameters */
    int k = 0;
    for(int j=0; j<njsx ; j++){
      cont = lib->rdata[idx];
      int nch = cont.l1;
      int nrb = cont.l2;

      int m0 = (nch+1)/6+1; if((nch+1)%6 == 0) m0--;
      int m1 = m0 * 6;

      for(int i=0 ; i<nrb ; i++){
        int i0 = i * m1;
        p[k++] = lib->xptr[idx][i0];
        for(int c=0 ; c<nch ; c++){
          p[k++] = lib->xptr[idx][i0 + c + 1];
        }
      }
      idx++;
    }

    cout << "#     Parameter     Uncertainty   Correlation" << endl;

    double *cptr;
    cptr = lib->xptr[idx];

    /*** print correlation matrix */
    for(int i=0; i<n; i++){
      int ki = i + i*n - i*(i+1)/2;

      cout << setw(5) <<i;
      outVal(p[i]);
      if(p[i] != 0.0) outVal(sqrt(cptr[ki] / p[i] / p[i]));
      else outVal(0.0);

      for(int j=0; j<=i; j++){
        int kj = j + j*n - j*(j+1)/2;

        int k  = (j < i) ?  i + j*n - j*(j+1)/2 : j + i*n - i*(i+1)/2;

        int c  = 0;
        if(cptr[ki]*cptr[kj] != 0.0){
          c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
          if((i == j) && (c == 999)) c = 1000;
        }
        cout << setw(5) << c;
      }
      cout << endl;
    }

    idx++;

    delete [] p;
  }

  return(idx);
}


/**********************************************************/
/*      Resolved Resonance LCOMP = 2, LRF = 1,2,3         */
/**********************************************************/

static int    row[]  = {0, 0, 18, 13, 11, 9, 8};
static double base[] = {1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6};

int DeceTableMF32LCOMP2LRF3(int idx, int lrf, ENDF *lib)
{
  Record cont = lib->rdata[idx];
  int  nrsa = cont.n2;
  cout << "#         NRSA" << setw(14) << nrsa << "  total number of resonances" <<  endl;

  double *pptr = lib->xptr[idx];
  double *p = new double [nrsa * 4];
  double *e = new double [nrsa * 4];

  idx++;

  cont = lib->rdata[idx];
  int nd = cont.l1;
  int nn = cont.l2;
  int nm = cont.n1;

  int mpar = nn / nrsa;

  cout << "#        NDIGT" << setw(14) << nd << "  number of digits" << endl;
  cout << "#          NNN" << setw(14) << nn << "  total number of resonance parameters" << endl;
  cout << "#           NM" << setw(14) << nm << "  number of INTG data lines" <<  endl;
  cout << "#         MPAR" << setw(14) << mpar << "  number of parameters per resonance (ER, GN, GG, ...)" <<  endl;


  int k = 0;
  for(int i=0; i<nrsa; i++){
    int i0 = 12*i;
    int i1 = i0 + 6;

    p[k  ] = pptr[i0];
    e[k++] = pptr[i1];

    if(lrf == 3){
      p[k  ] = pptr[i0 + 2];  e[k++] = pptr[i1 + 2]; // GN
      p[k  ] = pptr[i0 + 3];  e[k++] = pptr[i1 + 3]; // GG
      if(mpar == 4){
        p[k  ] = pptr[i0 + 4];  e[k++] = pptr[i1 + 4]; // GFA
      }
      else if(mpar == 5){
        p[k  ] = pptr[i0 + 5];  e[k++] = pptr[i1 + 5]; // GFB
      }
    }
    else{
      p[k  ] = pptr[i0 + 3];  e[k++] = pptr[i1 + 3]; // GN
      p[k  ] = pptr[i0 + 4];  e[k++] = pptr[i1 + 4]; // GG
      if(mpar == 4){
        p[k  ] = pptr[i0 + 5];  e[k++] = pptr[i1 + 5]; // GF
      }
    }
  }

  DeceTableMF32PrintCorr(nn,nm,nd,p,e,lib->xptr[idx]);

  delete [] p;
  delete [] e;

  idx++;

  return(idx);
}


/**********************************************************/
/*      Resolved Resonance LCOMP = 2, LRF = 7             */
/**********************************************************/
int DeceTableMF32LCOMP2LRF7(int idx, int njs, ENDF *lib)
{
  idx ++;

  cout << "#          NJS" << setw(14) << njs << "  total number of J-Pi" << endl;

  /*** count total number of parameters */
  int id1 = idx;
  int np = 0;
  for(int j=0; j<njs ; j++){
    Record cont = lib->rdata[id1++];
    int nch = cont.n2;
    double aj = cont.c1;
    double pj = cont.c2;
    cont = lib->rdata[id1++];
    int nrsa = cont.l2;

    int parity = (aj < 0.0) ? -1 : 1;
    if(aj == 0.0) parity = (int)pj;

    cout << "#     JP-Value"; outVal(4,1,abs(aj)); cout << ((parity < 0) ? "(-)" : "(+)") << endl;
    cout << "#          NCH" << setw(14) << nch << "  number of channels" << endl;
    cout << "#         NRSA" << setw(14) << nrsa << "  total number of resonances" <<  endl;

    np += nrsa * (nch + 1);
  }

  double *p = new double [np];
  double *e = new double [np];

  /*** copy parameters */
  int k = 0;
  for(int j=0; j<njs ; j++){
    Record cont = lib->rdata[idx++];
    int nch = cont.n2;
    cont = lib->rdata[idx];
    int nrsa = cont.l2;

    int m0 = (nch+1)/6+1; if((nch+1)%6 == 0) m0--;
    int m1 = m0 * 6;

    for(int i=0 ; i<nrsa ; i++){
      int i0 = 2 * i * m1;
      p[k  ] = lib->xptr[idx][i0];
      e[k++] = lib->xptr[idx][i0 + m1];

      for(int c=0 ; c<nch ; c++){
        p[k  ] = lib->xptr[idx][i0      + c + 1];
        e[k++] = lib->xptr[idx][i0 + m1 + c + 1];
      }
    }
    idx++;
  }

  Record cont = lib->rdata[idx];
  int nd = cont.l1;
  int nn = cont.l2;
  int nm = cont.n1;

  cout << "#        NDIGT" << setw(14) << nd << "  number of digits" << endl;
  cout << "#          NNN" << setw(14) << nn << "  total number of resonance parameters" << endl;
  cout << "#           NM" << setw(14) << nm << "  number of INTG data lines" <<  endl;


  DeceTableMF32PrintCorr(nn,nm,nd,p,e,lib->xptr[idx]);

  delete [] p;
  delete [] e;

  idx++;

  return(idx);
}


/**********************************************************/
/*      Print Correlation Matrix for LCOMP = 2            */
/**********************************************************/
void DeceTableMF32PrintCorr(int nn, int nm, int nd, double *p, double *e, double *cptr)
{
  int *cor = new int [nn * (nn+1) /2];

  for(int i=0 ; i<nn ; i++){
    for(int j=0 ; j<=i ; j++){
      cor[i*(i+1)/2 + j] = (i == j) ? 1000 : 0;
    }
  }

  int i = 0;
  for(int j=0 ; j<nm ; j++){
    int ii = (int)cptr[i++] - 1;
    int jj = (int)cptr[i++] - 1;
    int k0 = ii*(ii+1)/2 + jj;
    int km = ii*(ii+1)/2 + ii;

    for(int k=0 ; k<row[nd] ; k++){
      int c = (int)cptr[i++];

      int ij = k0 + k;
      if(ij < km) cor[ij] = (int)(1000.0 * c / base[nd]);
    }
  }

  cout << "#     Parameter     Uncertainty   Correlation" << endl;

  for(int i=0 ; i<nn ; i++){
    cout << setw(5) <<i;
    outVal(p[i]);
    if(p[i] == 0.0) outVal(0.0);
    else            outVal(abs(e[i] / p[i]));

    for(int j=0 ; j<=i ; j++){
      cout << setw(5) << cor[i*(i+1)/2 + j];
    }
    cout << endl;
  }

  delete [] cor;
}


/**********************************************************/
/*      Unresolved Resonance Parameter Covariance         */
/**********************************************************/
int DeceTableMF32URR(int idx, ENDF *lib)
{
  Record cont = lib->rdata[idx++];
  int nls = cont.n1;

  int idx0 = idx;
  for(int i=0 ; i<nls ; i++) idx0 ++;

  cont = lib->rdata[idx0];
  int mpar = cont.l1;
  int npar = cont.n2;

  double *p = new double [npar];

  int k = 0;
  for(int i=0 ; i<nls ; i++){
    cont = lib->rdata[idx];
    int njs = cont.n2; 

    for(int j=0 ; j<njs ; j++){
      int j0 = 6*j;
      p[k++] = lib->xptr[idx][j0    ]; if(mpar == 1) continue; // ER
      p[k++] = lib->xptr[idx][j0 + 2]; if(mpar == 2) continue; // GN
      p[k++] = lib->xptr[idx][j0 + 3]; if(mpar == 3) continue; // GG
      p[k++] = lib->xptr[idx][j0 + 4]; if(mpar == 4) continue; // GF
      p[k++] = lib->xptr[idx][j0 + 5]; // GX
    }
    idx++;
  }

  cout << "# Parameter    Uncertainty" << endl;

  double *cptr = lib->xptr[idx];

  for(int i=0 ; i<npar ; i++){
    int ki = i + i*npar - i*(i+1)/2;

    outVal(p[i]);
    outVal(sqrt(cptr[ki]));

    for(int j=0 ; j<=i ; j++){
      int kj = j + j*npar - j*(j+1)/2;

      int k  = (j < i) ?  i + j*npar - j*(j+1)/2 : j + i*npar - i*(i+1)/2;
      int c  = 0;
      if(cptr[ki]*cptr[kj] != 0.0){
        c = (int) (cptr[k] / sqrt(cptr[kj]) / sqrt(cptr[ki]) * 1000);
        if((i == j) && (c == 999)) c = 1000;
      }
      cout << setw(5) << c;
    }
    cout << endl;
  }

  delete [] p;

  return(idx);
}
