/******************************************************************************/
/**     GFR, Resonance Formulae                                              **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "endflib.h"
#include "gfr.h"

static inline complex<double> breit_wigner_profile  (complex<double>, double, double);
static inline double          arrange_matrixSLBW    (double, double *, BWResonance *);
static inline double          arrange_matrixRM      (double, double *, RMResonance *);
static inline void            invert_matrix         (complex<double> *);

extern double gcLorentzianWidth;
extern Smatrix Smat;

/**********************************************************/
/*      Single Level Breit-Wigner                         */
/**********************************************************/
Pcross gfrSLBreitWigner(const int kmax, const int l, const int j2, const double e,
                        Wfunc *wfn, BWResonance *res)
{
  const int msize = 10;
  Pcross  z;
  double  x[msize];

  double s2 = imag(wfn->phase) * imag(wfn->phase);
  
  for(int k=0 ; k<kmax ; k++){
    if( (res[k].l == l) && (res[k].j2 == j2) ){

      double gt = arrange_matrixSLBW(imag(wfn->d),x,&res[k]);
      double w  = e - res[k].er;
      double d  = w*w + 0.25*gt*gt;
      
      z.elastic += (x[0]*x[0] - 2.0*x[0]*gt*s2 + 2.0*w*x[0]*imag(wfn->phase2)) / d;
      z.fission += x[1]*x[1] / d;
      z.capture += x[3]*x[3] / d;
      z.other   += x[6]*x[6] / d;
    }
  }

  return(z);
}


/**********************************************************/
/*      MLBW defined in ENDF                              */
/**********************************************************/
Pcross gfrMLBreitWignerENDF(const int kmax, const int l, const int s2, const int j2,
                            const double e, Wfunc *wfn, BWResonance *res)
{
  const int msize = 10;
  Pcross  z;
  complex<double> w,y,p,r(0.0,0.0),t(0.0,0.0);
  double  x0[msize],x1[msize];

  for(int k0=0 ; k0<kmax ; k0++){
    if( (res[k0].l == l) && (res[k0].j2 == j2) ){
      double de  = res[k0].gn * (res[k0].s - real(wfn->d)) / (2*res[k0].p);
      double gt0 = arrange_matrixSLBW(imag(wfn->d),x0,&res[k0]);

      p = breit_wigner_profile(complex<double>(e,gcLorentzianWidth),res[k0].er+de,gt0);

      z.total   += (real(p) * real(wfn->phase2) - imag(p) * imag(wfn->phase2))*x0[0]/gt0;
      z.fission +=  real(p)*(x0[1]*x0[1])/gt0/gt0;
      z.capture +=  real(p)* x0[3]*x0[3] /gt0/gt0;

      w = complex<double>(0.0,0.0);
      for(int k1=0 ; k1<kmax ; k1++){
        if( (res[k1].l == l) && (res[k1].j2 == j2) ){
          double gt1 = arrange_matrixSLBW(imag(wfn->d),x1,&res[k1]);
          y = complex<double>(res[k1].er-res[k0].er, -(gt0+gt1)*0.5);
          y = 1.0 /y;
          w += y*x1[0];
        }
      }
      w = complex<double>(-imag(w) + 1.0, -real(w));
      r += x0[0]*w*p/gt0;

      /*** sum i Gn /(Ek + Delta - E -iG/2) */
      t += complex<double>(0.0,x0[0]) / complex<double>(res[k0].er + de - e,-gt0/2.0 - gcLorentzianWidth);
    }
  }

  z.total    = 4*(imag(wfn->phase)*imag(wfn->phase) + z.total);
  z.elastic  = z.total - 4* real(r);
  z.fission *= 4;
  z.capture *= 4;
  z.total    = z.elastic + z.fission + z.capture;  // total re-calculate

  /*** S-matrix element for elastic scattering channel */
  Smat.setElement(l,j2,s2, wfn->phase2 * (1.0 + t));

  return(z);
}


/**********************************************************/
/*      General Multilevel Breit-Wigner                   */
/**********************************************************/
Pcross gfrBreitWigner(const int kmax, const int l, const int j2, const double e,
                      Wfunc *wfn, BWResonance *res)
{
  const int msize = 10;
  Pcross  z;
  complex<double> rmat[msize];
  complex<double> w[msize],y,p;
  double  x0[msize],x1[msize];

  for(int i=0 ; i<msize ; i++)  rmat[i] = complex<double>(0.0,0.0);

  for(int k0=0 ; k0<kmax ; k0++){
    if( (res[k0].l == l) && (res[k0].j2 == j2) ){
      double gt0 = arrange_matrixSLBW(imag(wfn->d),x0,&res[k0]);
      p = breit_wigner_profile(complex<double>(e,gcLorentzianWidth),res[k0].er,gt0);

      z.total += (real(p) * real(wfn->phase2) - imag(p) * imag(wfn->phase2))*x0[0]/gt0;

      for(int i=0 ; i<msize ; i++)  w[i] = complex<double>(0.0,0.0);
      for(int k1=0 ; k1<kmax ; k1++){
        if( (res[k1].l == l) && (res[k1].j2 == j2) ){
          double gt1 = arrange_matrixSLBW(imag(wfn->d),x1,&res[k1]);
          y = complex<double>(res[k1].er-res[k0].er, -(gt0+gt1)*0.5);
          y = 1.0 /y;
          for(int i=0 ; i<msize ; i++) w[i] += y*x1[i];
        }
      }
      for(int i=0 ; i<msize ; i++) w[i] = complex<double>(-imag(w[i]),-real(w[i]));
      w[0] += 1.0;

      for(int i=0 ; i<msize ; i++) rmat[i] += x0[i]*w[i]*p/gt0;
    }
  }

  z.total   = 4*(imag(wfn->phase)*imag(wfn->phase) + z.total);
  z.elastic = z.total - 4 * real(rmat[0]);
  z.fission =         - 4 * real(rmat[1]);
  z.capture =         - 4 * real(rmat[3]);
  z.other   =         - 4 * real(rmat[6]);

  return(z);
}


/**********************************************************/
/*      Reich-Moore R-Matrix                              */
/**********************************************************/
Pcross gfrReichMoore(const int kmax, const int l, const int s2, const int j2,
                     const int tspin2, 
                     const double e, Wfunc *wfn, RMResonance *res)
{
  const int msize = 6;
  complex<double> smat[msize],rmat[msize],wmat[msize];
  double x[msize];

  bool sdep = false; // channel spin dependent
  for(int k=0 ; k<kmax ; k++){ if(res[k].j2 < 0) sdep = true; }

  /*** R-matrix elements */
  for(int i=0 ; i<msize ; i++)  smat[i] = rmat[i] = wmat[i] = complex<double>(0.0,0.0);
  for(int k=0 ; k<kmax ; k++){

    if(res[k].l == l){
      if(sdep){
        int j2r = (int)abs(res[k].j2);
        int sr  = (res[k].j2 < 0) ? -1 : 1;
                      
        if( (j2r == j2) && (sr == s2) ){
          arrange_matrixRM(imag(wfn->d),x,&res[k]);
          complex<double> w(res[k].er-e, -res[k].gg/2.0 - gcLorentzianWidth);
          w = 1.0 / w;
          for(int i=0 ; i<msize ; i++)  rmat[i] += w*x[i];
        }
      }

      else{
        if(res[k].j2 == j2){
          if((l > 0) && (s2 != -1) && (tspin2 != 0)) continue;
          arrange_matrixRM(imag(wfn->d),x,&res[k]);
          complex<double> w(res[k].er-e, -res[k].gg/2.0 - gcLorentzianWidth);
          w = 1.0 / w;
          for(int i=0 ; i<msize ; i++)  rmat[i] += w*x[i];
        }
      }
    }
  }

  for(int i=0 ; i<msize ; i++){
    wmat[i] = complex<double>(imag(rmat[i])/2.0,-real(rmat[i])/2.0);
    if(i==0 || i==2 || i==5) wmat[i] += 1.0;
  }

  invert_matrix(wmat);

  /*** S-matrix elements */
  smat[0] = wfn->phase2 * (2.0*wmat[0]-1.0);
  for(int i=1 ; i<msize ; i++) smat[i] = wfn->phase * 2.0 * wmat[i];

  /*** copy S-matrix element for elastic */
  Smat.setElement(l,j2,s2,smat[0]);

  /*** cross sections from S-matrix */
  Pcross  z;

  z.total    = (1.0-real(smat[0]))*2.0;
  z.elastic  = (1.0-real(smat[0]))*(1.0-real(smat[0])) + imag(smat[0])*imag(smat[0]);
  z.fission  = real(smat[1])*real(smat[1]) + imag(smat[1])*imag(smat[1])
             + real(smat[3])*real(smat[3]) + imag(smat[3])*imag(smat[3]);
  z.capture  = z.total - z.elastic - z.fission;

  return(z);
}


/**********************************************************/
/*      Collision Matrix for MLBW                         */
/**********************************************************/
Pcross gfrBreitWignerUmatrix(const int kmax, const int l, const int s2, const int j2,
                             const double e, Wfunc *wfn, BWResonance *res)
{
  const int msize = 10;
  complex<double> w, tmat[msize], smat[msize];
  double  x[msize];

  for(int i=0 ; i<msize ; i++)  tmat[i] = complex<double>(0.0,0.0);
  for(int k=0 ; k<kmax ; k++){
    if( (res[k].l == l) && (res[k].j2 == j2) ){
      double de  = res[k].gn * (res[k].s - real(wfn->d)) / (2*res[k].p);
      double gt  = arrange_matrixSLBW(imag(wfn->d),x,&res[k]);
      w = complex<double>(res[k].er + de - e, -gt/2.0 - gcLorentzianWidth);
      for(int i=0 ; i<msize ; i++) tmat[i] += complex<double>(0.0,x[i]) / w;
    }
  }

  /*** S-matrix elements */
  for(int i=0 ; i<msize ; i++) if(i==0 || i==2 || i==5 || i==9) tmat[i] += 1.0;

  smat[0] = wfn->phase2 * tmat[0];
  smat[1] = wfn->phase  * tmat[1];
  smat[3] = wfn->phase  * tmat[3];
  smat[6] = wfn->phase  * tmat[6];

  /*** copy S-matrix element for elastic */
  Smat.setElement(l,j2,s2,smat[0]);

  /*** cross sections from S-matrix, not unitary */
  Pcross  z;

  z.total    = (1.0-real(smat[0]))*2.0;
  z.elastic  = (1.0-real(smat[0]))*(1.0-real(smat[0]))  + imag(smat[0])*imag(smat[0]);
  z.fission  = real(smat[1])*real(smat[1]) + imag(smat[1])*imag(smat[1]);
  z.capture  = real(smat[3])*real(smat[3]) + imag(smat[3])*imag(smat[3]);
  z.other    = real(smat[6])*real(smat[6]) + imag(smat[6])*imag(smat[6]);

  return(z);
}


/**********************************************************/
/*      Matrix Element Re-arrangement for SLBW            */
/*        elast  fiss   capt   other                      */
/*        (n,n)  (n,f)  (n,g)  (n,x)                      */
/*      n  x0                                             */
/*      f  x1     x2                                      */
/*      g  x3     x4     x5                               */
/*      x  x6     x7     x8    x9                         */
/**********************************************************/
double arrange_matrixSLBW(double d, double *x, BWResonance *res)
{
  double gn  = res->gn /res->p * d;
  double gg  = res->gg;
  double gf  = abs(res->gf);
  double gx  = res->gx;
  double si  = (res->gf < 0.0) ? -1.0 : 1.0;
  
  x[0] = gn;
  x[1] = sqrt(gn * gf) * si;
  x[2] = gf;
  x[3] = sqrt(gn * gg);
  x[4] = sqrt(gf * gg) * si;
  x[5] = gg;
  x[6] = sqrt(gn * gx);
  x[7] = sqrt(gf * gx) * si;
  x[8] = sqrt(gg * gx);
  x[9] = gx;

  return(gn+gg+gf+gx);
}


/**********************************************************/
/*      Matrix Element Re-arrangement for RM              */
/**********************************************************/
double arrange_matrixRM(double d, double *x, RMResonance *res)
{
  double gn  = res->gn /res->p * d;
  double gg  = res->gg;
  double gfa = abs(res->gf1);
  double gfb = abs(res->gf2);

  x[0] = gn;
  x[1] = sqrt(gn  * gfa) * ((res->gf1 < 0.0) ? -1:1);
  x[2] = gfa;
  x[3] = sqrt(gn  * gfb) * ((res->gf2 < 0.0) ? -1:1);
  x[4] = sqrt(gfa * gfb) * ((res->gf1 < 0.0) ? -1:1) * ((res->gf2 < 0.0) ? -1:1);
  x[5] = gfb;

  return(gn+gg+gfa+gfb);
}


/**********************************************************/
/*      Resonance Shape                                   */
/**********************************************************/
inline complex<double> breit_wigner_profile(complex<double> e, double e0, double g)
{
  complex<double> x  = complex<double>(0.0, g / 2.0);
  complex<double> y  = e - complex<double>(e0,0.0) + x;
  complex<double> bw = x / y;

  return(bw);
}


/**********************************************************/
/*      Invert 3x3 Complex Matrix                         */
/**********************************************************/
void invert_matrix(complex<double> *a)
{
  complex<double> d,unity(1.0,0.0);
  complex<double> b[6];

  if( (a[2] == unity && a[5] == unity) ){
    a[0] = 1.0 / a[0];
    return;
  }

  d = a[0]*a[2]*a[5] + a[1]*a[3]*a[4]*2.0
    -(a[2]*a[3]*a[3] + a[0]*a[4]*a[4] + a[5]*a[1]*a[1]);

  b[0] = a[2]*a[5] - a[4]*a[4];
  b[1] = a[3]*a[4] - a[1]*a[5];
  b[2] = a[0]*a[5] - a[3]*a[3];
  b[3] = a[1]*a[4] - a[2]*a[3];
  b[4] = a[1]*a[3] - a[0]*a[4];
  b[5] = a[0]*a[2] - a[1]*a[1];

  for(int i=0 ; i<6 ; i++) a[i] = b[i]/d;
}
