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

static int DeceTableMF32LCOMP0(int, int, ENDF *);
static int DeceTableMF32LCOMP1(int, int, int, ENDF *);
static int DeceTableMF32LCOMP2(int, int, ENDF *);
static int DeceTableMF32URR   (int, int, ENDF *);


/**********************************************************/
/*      Process MF=32                                     */
/**********************************************************/
void DeceTableMF32(ENDF *lib)
{
  int    idx  = 0;
  Record cont = lib->rdata[idx++];
  int    ner  = cont.n1;

  cout << "#   NER:" << setw(6) << ner << endl;

  for(int i=0 ; i<ner ; i++){
    cont = lib->rdata[idx++];
    double el = cont.c1;
    double eh = cont.c2;
    int lru = cont.l1;
    int lrf = cont.l2;
    int nro = cont.n1;

    cout << "#    EL:"; outVal(el); cout << endl;
    cout << "#    EH:"; outVal(eh); cout << endl;
    cout << "#   LRU:" << setw(6) << lru << endl;
    cout << "#   LRF:" << setw(6) << lrf << endl;

    if( !(lrf == 1 || lrf == 2 || lrf == 3) ){
      message << "LRF = " << lrf << " not yet implemented";
      WarningMessage();
      return;
    }

    if(lru == 1){
      if(nro !=0 ){
        cont = lib->rdata[idx++];
        int ni = cont.n2;
        for(int i=0 ; i<ni ; i++) idx ++;
      }

      cont = lib->rdata[idx++];
      int lcomp = cont.l2;
      int nls   = cont.n1;
      int isr   = cont.n2;

      cout << "# LCOMP:" << setw(6) << lcomp << endl;
      cout << "#   NLS:" << setw(6) << nls << endl;
      cout << "#   ISR:" << setw(6) << isr << endl;

      if(isr != 0) idx++;

      /*** ENDF/B-V Compatible Resolved Resonance Format */
      if(lcomp == 0){
        idx = DeceTableMF32LCOMP0(idx,nls,lib);
      }
      /*** General Resolved Resonance Format */
      else if(lcomp == 1){
        idx = DeceTableMF32LCOMP1(idx,lru,lrf,lib);
      }
      else{
        idx = DeceTableMF32LCOMP2(idx,lrf,lib);
      }
    }
    else{
      cont = lib->rdata[idx++];
      int nls = cont.n1;
      idx = DeceTableMF32URR(idx,nls,lib);
    }
    cout << endl;
    cout << endl;
  }
}


/**********************************************************/
/*      Resolved Resonance Parameter, LCOMP = 0           */
/**********************************************************/
int DeceTableMF32LCOMP0(int idx, int nls, ENDF *lib)
{
  double p[4], e[4], c[10];

  for(int inls=0 ; inls<nls ; inls++){
    Record cont = lib->rdata[idx];
    int nrs = cont.n2;
    cout << "#         L:" << setw(6) << inls << endl;
    cout << "#       NRS:" << setw(6) << nrs << endl;

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

      for(int k1=0 ; k1<4 ; k1++){
        outVal(p[k1]);
        outVal(e[k1]/abs(p[k1]));
            
        for(int k2=0 ; k2<=k1 ; k2++){
          int k = c[k1*(k1+1)/2 + k2] / e[k1] / e[k2] * 1000.0;
          if(k1 == k2) k = 1000;
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
/*      Resolved Resonance Parameter, LCOMP = 1           */
/**********************************************************/
int DeceTableMF32LCOMP1(int idx, int lru, int lrf, ENDF *lib)
{
  Record cont = lib->rdata[idx++];
  int nsrs = cont.n1;
  int nlrs = cont.n2;

  cout << "#  NSRS:" << setw(6) << nsrs << endl;
  cout << "#  NLRS:" << setw(6) << nlrs << endl;

  for(int insrs=0 ; insrs<nsrs ; insrs++){
    cont = lib->rdata[idx];
    int mpar = cont.l1;
    int nrb  = cont.n2;
    int n    = nrb * mpar;

    cout << "#     INSRS:" << setw(6) << insrs << endl;
    cout << "#      MPAR:" << setw(6) << mpar << endl;
    cout << "#       NRB:" << setw(6) << nrb << endl;

    double *p = new double [n];

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

    double *cptr;
    cptr = &lib->xptr[idx][6*nrb];

    for(int i0=0; i0<nrb; i0++){
      for(int i1=0; i1<mpar; i1++){
        int i  = i0*mpar + i1;
        int ki = i + i*n - i*(i+1)/2;

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

  if(lru != 2){
    for(int inlrs=0 ; inlrs<nlrs ; inlrs++) idx++;
  }

  return(idx);
}


/**********************************************************/
/*      Resolved Resonance Parameter, LCOMP = 2           */
/**********************************************************/
int DeceTableMF32LCOMP2(int idx, int lrf, ENDF *lib)
{
  static int    row[]  = {0, 0, 18, 13, 11, 9, 8};
  static double base[] = {1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6};

  Record cont = lib->rdata[idx];
  int nrsa = cont.n2;

  cout << "#  NRSA:" << setw(6) << nrsa << endl;

  double *pptr = lib->xptr[idx];
  double *p = new double [nrsa * 4];
  double *e = new double [nrsa * 4];

  idx++;

  cont = lib->rdata[idx];
  int nd = cont.l1;
  int nn = cont.l2;
  int nm = cont.n1;

  cout << "# NDIGT:" << setw(6) << nd << endl;
  cout << "#    NM:" << setw(6) << nm << endl;

  int mpar = nn / nrsa;

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

  int *cor = new int [nn * (nn+1) /2];
  for(int i=0 ; i<nn ; i++){
    for(int j=0 ; j<=i ; j++){
      cor[i*(i+1)/2 + j] = (i == j) ? 1000 : 0;
    }
  }

  int i = 0;
  for(int j=0 ; j<nm ; j++){
    int ii = lib->iptr[idx][i++] - 1;
    int jj = lib->iptr[idx][i++] - 1;
    int k0 = ii*(ii+1)/2 + jj;
    int km = ii*(ii+1)/2 + ii;

    for(int k=0 ; k<row[nd] ; k++){
      int c = lib->iptr[idx][i++];

      int ij = k0 + k;
      if(ij < km) cor[ij] = (int)(1000.0 * c / base[nd]);
    }
  }

  for(int i=0 ; i<nn ; i++){
    outVal(p[i]);
    if(p[i] == 0.0) outVal(0.0);
    else            outVal(abs(e[i] / p[i]));

    for(int j=0 ; j<=i ; j++){
      cout << setw(5) << cor[i*(i+1)/2 + j];
    }
    cout << endl;
  }

  delete [] p;
  delete [] e;

  idx++;

  return(idx);
}


/**********************************************************/
/*      Unresolved Resonance Parameter Covariance         */
/**********************************************************/
int DeceTableMF32URR(int idx, int nls, ENDF *lib)
{
  int idx0 = idx;
  for(int i=0 ; i<nls ; i++) idx0 ++;

  Record cont = lib->rdata[idx0];
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
