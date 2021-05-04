/******************************************************************************/
/**     GFR, Calculate Reich-Moore Cross Sections                            **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "constant.h"

static Pcross RMMainCalc (const double, System *);
static int    RMLoadResonances (int, System *, ENDF *);
static void   RMScatteringRadius (const int, System *, const double, ENDF *);
static void   RMResonancePenetrability (const int, System *, ENDF *);

extern Smatrix Smat;

static double ap_pen = 0.0, ap_phi = 0.0;
static int kmax = 0;
static RMResonance *res;
static bool dataload = false;

#undef DEBUG_RESONANCE

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection3(const int ner, const double elab, System *sys, ENDF *lib)
{
  /*** release allocated memory */
  if(sys->isLastCall()){
    if(dataload){
      delete [] res;
      kmax = 0;
      dataload = false;
    }
    Pcross s;
    return(s);
  }

  /*** set AP and rho */
  RMScatteringRadius(ner,sys,elab,lib);
  if(sys->isFirstCall()){
    if(dataload){
      delete [] res;
      kmax = 0;
      dataload = false;
    }

    res = new RMResonance [MAX_RESONANCE];

    /*** load resonance parameters into Resonance class */
    int idx =sys->idx[ner] + 2;
    if(sys->nro[ner] == 1) idx ++;
    kmax = RMLoadResonances(idx,sys,lib);
    RMResonancePenetrability(ner,sys,lib);
    dataload = true;
    sys->OnceCalled();
  }

  Pcross sig = RMMainCalc(elab,sys);

  return(sig);
}


/**********************************************************/
/*      Main Calculatioin of Reich Moore Formula          */
/**********************************************************/
Pcross RMMainCalc(const double elab, System *sys)
{
  ChannelWaveFunc wfn;
  Pcross sig, z;

  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01 / ((sys->target_spin2+1)*2);
  double alpha_pen = sys->wave_number * ap_pen;
  double alpha_phi = sys->wave_number * ap_phi;

  /*** for all L partial waves */
  Smat.resetIndex();
  for(int l=0 ; l<sys->nl ; l++){

    /*** this part, it is not so clear which AP should be used */
    if(sys->apl[l] > 0.0){
      alpha_phi = sys->wave_number * sys->apl[l];
//    if(sys->naps[ner] != 0) alpha_pen = alpha_phi;
    }

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

        z = gfrReichMoore(kmax,l,ss-smin-1,jj,sys->target_spin2,elab,&wfn,res);

        sig.total    += x3*z.total;
        sig.elastic  += x3*z.elastic;
        sig.capture  += x3*z.capture;
        sig.fission  += x3*z.fission;

        Smat.inclIndex();
      }
    }
  }

  return(sig);
}


/**********************************************************/
/*     Copy All Resoance Parameters                       */
/**********************************************************/
int RMLoadResonances(int idx, System *sys, ENDF *lib)
{
  int k = 0;

  for(int l=0 ; l<sys->nl ; l++){

    int lx  = lib->rdata[idx].l1;
    int nrs = lib->rdata[idx].n2;

    for(int i=0 ; i<nrs ; i++){
      int j = 6*i;

      res[k].l   = lx;
      res[k].er  = lib->xptr[idx][j];
      res[k].j2  = (int)(lib->xptr[idx][j+1] * 2.0);
      res[k].gn  = lib->xptr[idx][j+2];
      res[k].gg  = lib->xptr[idx][j+3];
      res[k].gf1 = lib->xptr[idx][j+4];
      res[k].gf2 = lib->xptr[idx][j+5];
      
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
    cout << setw(11) << res[j].gn;
    cout << setw(11) << res[j].gg;
    cout << setw(11) << res[j].gf1;
    cout << setw(11) << res[j].gf2 << endl;
  }
#endif

  return(k);
}


/**********************************************************/
/*     Penetrability at Each Resoance                     */
/**********************************************************/
void RMScatteringRadius(const int ner, System *sys, const double elab, ENDF *lib)
{
  int apidx = 0;

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
}


/**********************************************************/
/*     Penetrability at Each Resonance                    */
/**********************************************************/
void RMResonancePenetrability(const int ner, System *sys, ENDF *lib)
{
  int apidx = 0;
  double c2 = 2.0 * AMUNIT / (VLIGHTSQ * HBARSQ);
  double mu = sys->reduced_mass * MNEUTRON;

  if(sys->nro[ner] == 1) apidx = sys->idx[ner] + 1;

  /*** penetrability at each resonance */
  for(int k=0 ; k<kmax ; k++){
    double r = ap_pen;
    if((sys->nro[ner] == 1) && (sys->naps[ner] == 1)){
      r = ENDFInterpolation(lib,fabs(res[k].er),true,apidx);
    }
    double ecm = res[k].er * mu;
    double alpha = sqrt(c2 * mu * fabs(ecm) * 1e-6) * r;

    complex<double> q = gfrLfunction(res[k].l,alpha,0.0);
    res[k].s = q.real();
    res[k].p = q.imag();
  }
}

