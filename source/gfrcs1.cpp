/******************************************************************************/
/**     GFR, Calculate SLBW Cross Sections                                   **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "constant.h"

static int gfrLoadBWResonances (int, System *, BWResonance *, ENDF *);
static void gfrResonancePenetrability (const int, System *, const int, BWResonance *, double, double *, double *, ENDF *);
extern Smatrix Smat;

#undef DEBUG_RESONANCE

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection1(const int lrf, const int ner, const double elab, System *sys, ENDF *lib)
{
  ChannelWaveFunc wfn;
  Pcross sig, z;
  BWResonance *res = NULL;

  res = new BWResonance[MAX_RESONANCE];

  /*** load resonance parameters into Resonance class */
  int idx =sys->idx[ner] + 2;
  if(sys->nro[ner] == 1) idx ++;
  int kmax = gfrLoadBWResonances(idx,sys,res,lib);

  double ap_pen = 0.0, ap_phi = 0.0;

  gfrResonancePenetrability(ner,sys,kmax,res,elab,&ap_pen,&ap_phi,lib);

  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01 / ((sys->target_spin2+1)*2);
  double alpha_pen = sys->wave_number * ap_pen;
  double alpha_phi = sys->wave_number * ap_phi;

  /*** for all L partial waves */
  Smat.resetIndex();
  for(int l=0 ; l<sys->nl ; l++){

    /*** calculate phi and P */
    ChannelWaveFunc wpen, wphi;
    gfrPenetrability(l,alpha_phi,&wphi);
    wfn.setPhase(wphi.H);

    gfrPenetrability(l,alpha_pen,&wpen);
    wfn.setData(wpen.a,wpen.H,wpen.D);

    /*** channel spin */
    int smin = abs(sys->target_spin2-1);
    int smax =     sys->target_spin2+1;
    for(int ss=smin ; ss<=smax ; ss+=2){

      /*** total spin */
      int jmin = abs(2*l-ss);
      int jmax =     2*l+ss ;
      for(int jj=jmin ; jj<=jmax ; jj+=2){
        double x3 = (jj+1)*x1;

        if(lrf == 1) z = gfrSLBreitWigner(kmax,l,jj,elab,&wfn,res);
        else         z = gfrMLBreitWignerENDF(kmax,l,ss-smin-1,jj,elab,&wfn,res);

        sig.total    += x3*z.total;
        sig.elastic  += x3*z.elastic;
        sig.capture  += x3*z.capture;
        sig.fission  += x3*z.fission;
        sig.other    += x3*z.other;

        Smat.inclIndex();
      }
    }
  }

  delete [] res;
  return(sig);
}


/**********************************************************/
/*     Copy All Resoance Parameters                       */
/**********************************************************/
int gfrLoadBWResonances(int idx, System *sys, BWResonance *res, ENDF *lib)
{
  int k = 0;

  for(int l=0 ; l<sys->nl ; l++){

    int lx  = lib->rdata[idx].l1;
    int lrx = lib->rdata[idx].l2;
    int nrs = lib->rdata[idx].n2;

    for(int i=0 ; i<nrs ; i++){
      int j = 6*i;

      res[k].l   = lx;
      res[k].er  = lib->xptr[idx][j];
      res[k].j2  = (int)(lib->xptr[idx][j+1] * 2.0);
      res[k].gt  = lib->xptr[idx][j+2];
      res[k].gn  = lib->xptr[idx][j+3];
      res[k].gg  = lib->xptr[idx][j+4];
      res[k].gf  = lib->xptr[idx][j+5];

      if(lrx != 0) res[k].gx = res[k].gt - (res[k].gn + res[k].gg + res[k].gf);

      k++;
      if(k > MAX_RESONANCE) return(0);
    }
    idx++;
  }

#ifdef DEBUG_RESONANCE
  for(int j=0 ; j<k ; j++){
    cout << setprecision(3);
    cout << setw( 3) << res[j].l;
    cout << setw( 3) << res[j].j2; 
    cout << setw(11) << res[j].er;
    cout << setw(11) << res[j].gt;
    cout << setw(11) << res[j].gn;
    cout << setw(11) << res[j].gg;
    cout << setw(11) << res[j].gf << endl;
  }
#endif

  return(k);
}


/**********************************************************/
/*     Penetrability at Each Resoance                     */
/**********************************************************/
void gfrResonancePenetrability(const int ner, System *sys, const int kmax, BWResonance *res, double elab, double *a0, double *a1, ENDF *lib)
{
  int apidx = 0;
  double ap_pen = 0.0, ap_phi = 0.0;

  if(sys->nro[ner] == 1) apidx = sys->idx[ner] + 1;

  /*** NRO = 0: energy independent radius case */
  if(sys->nro[ner] == 0){
    /*** if NAPS = 0, calculate L+iS at 0.123 AWRI**1/3 + 0.08,
         but hard-sphere phase is still at alpha = k x AP */
    ap_pen = (sys->naps[ner] == 0) ? gfrENDFChannelRadius(sys->target_A) : sys->radius;
    ap_phi = sys->radius;
  }
  /*** NRO = 1: energy dependent radius case */
  else{
    ap_phi = ENDFInterpolation(lib,elab,true,apidx) * 10.0;

    if(     sys->naps[ner] == 0) ap_pen = gfrENDFChannelRadius(sys->target_A);
    else if(sys->naps[ner] == 1) ap_pen = ap_phi;
    else                         ap_pen = sys->radius;
  }

  /*** penetrability at each resonance */
  for(int k=0 ; k<kmax ; k++){
    double rho = ap_pen;
    if((sys->nro[ner] == 1) && (sys->naps[ner] == 1)){
      rho = ENDFInterpolation(lib,fabs(res[k].er),true,apidx);
    }
    double alpha = gcKfactor * sqrt(fabs(res[k].er) * 1.0e-06 * sys->reduced_mass) * rho;

    complex<double> q = gfrLfunction(res[k].l,alpha,0.0);
    res[k].s = q.real();
    res[k].p = q.imag();
  }

  *a0 = ap_pen;
  *a1 = ap_phi;
}

