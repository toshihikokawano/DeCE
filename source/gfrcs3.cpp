/******************************************************************************/
/**     GFR, Calculate Reich-Moore Cross Sections                            **/
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

static int gfrLoadRMResonances (int, System *, RMResonance *, ENDF *);
extern Smatrix Smat;

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection3(const int ner, const double elab, System *sys, ENDF *lib)
{
  Wfunc  wfn;
  Pcross sig, z;
  RMResonance *res = NULL;

  res = new RMResonance[MAX_RESONANCE];

  /*** load resonance parameters into Resonance class */
  int kmax = gfrLoadRMResonances(sys->idx[ner],sys,res,lib);

  double ap = (sys->naps[ner] == 0) ? gfrENDFChannelRadius(sys->target_A) : sys->radius;
  for(int k=0 ; k<kmax ; k++){
    complex<double> q = gfrLfunction(res[k].l,res[k].er,sys->reduced_mass,ap);
    res[k].s = q.real();
    res[k].p = q.imag();
  }

  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01 / ((sys->target_spin2+1)*2);
  double x2 = sys->wave_number * ap;

  /*** for all L partial waves */
  Smat.resetIndex();
  for(int l=0 ; l<sys->nl ; l++){

    double phase = gfrPenetrability(l,sys->alpha,&wfn);
    wfn.phase  = complex<double>(cos(  phase), -sin(  phase));
    wfn.phase2 = complex<double>(cos(2*phase), -sin(2*phase));

    /*** if NAPS = 0, calculate L+iS at 0.123 AWRI**1/3 + 0.08,
         but hard-sphere phase is still at alpha = ka */
    if(sys->naps[ner] == 0) gfrPenetrability(l,x2,&wfn);

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

  delete [] res;
  return(sig);
}


/**********************************************************/
/*     Copy All Resoance Parameters                       */
/**********************************************************/
int gfrLoadRMResonances(int idx, System *sys, RMResonance *res, ENDF *lib)
{
  /*** first, increment idx to point resonance parameter data */
  idx ++;
  
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

  return(k);
}


