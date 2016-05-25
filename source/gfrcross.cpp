/******************************************************************************/
/**     GFR, Calculate Cross Sections                                        **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <iomanip>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "constant.h"

static void   gfrSubsectionRRR   (const int, const double, System *, ENDF *);
static void   gfrSubsectionURR   (const int, const int, System *, ENDF *);
static int    gfrFindEnergyRange (const int, const double, System *);


extern const double gcLorentzianWidth;
extern Smatrix Smat;

/**********************************************************/
/*      Pointwise Cross Section int Resonance Range       */
/**********************************************************/
Pcross gfrCrossSection(const int urr, const double elab, System *sys, ENDF *lib)
{
  Pcross sig;

  /*** find an energy range for a given incident energy */
  int ner = gfrFindEnergyRange(urr,elab,sys);
  if(ner < 0) return(sig);
 
  int lru = sys->lru[ner];
  int lrf = sys->lrf[ner];

  /*** potential scattering only case, return zero */
  if(lru == 0) return(sig);

  /*** resolved resonance region */
  else if(lru == 1){
    if(lrf == 1){
      gfrSubsectionRRR(sys->idx[ner],elab,sys,lib);
      sig = gfrCrossSection1(1,ner,elab,sys,lib);
    }
    else if(lrf == 2){
      gfrSubsectionRRR(sys->idx[ner],elab,sys,lib);
      sig = gfrCrossSection1(2,ner,elab,sys,lib);
    }
    else if(lrf == 3){
      gfrSubsectionRRR(sys->idx[ner],elab,sys,lib);
      sig = gfrCrossSection3(ner,elab,sys,lib);
    }
    else if(lrf == 7){
      sig = gfrCrossSection7(ner,elab,sys,lib);
    }
  }

  /*** unresolved resonance region */
  else if(lru == 2){
    gfrSubsectionURR(lrf,sys->idx[ner],sys,lib);

    /*** if LSSF = 1, do not calculate cross section */
    if(sys->selfshield_flag != 1){
      sig = gfrCrossSectionURR(ner,elab,sys,lib);
    }
  }

  return(sig);
}


/**********************************************************/
/*     Find Resoance Range for Given Energy               */
/**********************************************************/
int gfrFindEnergyRange(const int urr, const double elab, System *sys)
{
  int ner = -1;
  for(int i=0 ; i<sys->nrange ; i++){
    if( (sys->emin[i] <= elab) &&  (elab <= sys->emax[i]) ){

      if((urr > 0) && (urr != sys->lru[i])) continue;

      ner = i;
      break;
    }
  }
  return(ner);
}


/**********************************************************/
/*     RRR Subsection Specific Parameters                 */
/**********************************************************/
void gfrSubsectionRRR(const int idx, const double elab, System *sys, ENDF *lib)
{
  /*** read CONT at each subsection */
  Record cont = lib->rdata[idx];

  sys->target_spin2 = (int)(2.0*cont.c1);  // 2 x target spin
  sys->radius       = cont.c2 * 10.0;      // scattering radius in fm
  sys->nl           = cont.n1;             // number of orbital angular momentum

  gfrSetEnergy(elab, sys);
}


/**********************************************************/
/*     Set Up Energy Dependent Quantities                 */
/**********************************************************/
void gfrSetEnergy(const double elab, System *sys)
{
  sys->ecms        = elab * sys->reduced_mass;
  sys->wave_number = sqrt(sys->reduced_mass) * gcKfactor*sqrt(sys->ecms * 1.0e-06);
  sys->alpha       = sys->wave_number * sys->radius;
}


/**********************************************************/
/*     URR Subsection Specific Parameters                 */
/**********************************************************/
void gfrSubsectionURR(const int lrf, const int idx, System *sys, ENDF *lib)
{
  /*** read the first CONT at each subsection */
  Record cont = lib->rdata[idx];

  /*** case A
       fission width not given
       all parameters are energy independent */
  if(lrf == 1){
    if(sys->avefission_flag == 0){
      sys->target_spin2    = (int)(2.0*cont.c1);
      sys->radius          = cont.c2 * 10.0;
      sys->selfshield_flag = cont.l1;
      sys->nl              = cont.n1;
    }
  /*** case B
       fission width given
       only fission width is energy-dependent, others energy independent */
    else{
      sys->target_spin2    = (int)(2.0*cont.c1);
      sys->radius          = cont.c2 * 10.0;
      sys->selfshield_flag = cont.l1;
      sys->nfw             = cont.n1;
      sys->nl              = cont.n2;
      sys->fwx = lib->xptr[idx];
    }
  }
  /*** case C
       all energy-dependent */
  else{
    sys->target_spin2    = (int)(2.0*cont.c1);
    sys->radius          = cont.c2 * 10.0;
    sys->selfshield_flag = cont.l1;
    sys->nl              = cont.n1;
  }
}


/**********************************************************/
/*     Hankel Function and Penetrability                  */
/*     ****                                               */
/*        h     = G+iF = O                                */
/*        d     = (G+iF)' = O'                            */
/*        wfn.d = a (G+iF)'/(G+iF) = L                    */
/*        phase = atan(Im O/Re O) = atan F/G              */
/**********************************************************/
double gfrPenetrability(const int l, const double a, Wfunc *wfn)
{
  complex<double> h,d;

  double g0 = cos(a);
  double f0 = sin(a);
  double g1 =  f0 + g0/a;
  double f1 = -g0 + f0/a;

  if(l == 0){
    h = complex<double>( g0,f0);
    d = complex<double>(-f0,g0);
  }
  else if(l==1){
    h = complex<double>(g1,f1);
    d = complex<double>(-f1-g0/(a*a),g1-f0/(a*a));
  }
  else{
    double g2 = 0.0, f2 = 0.0;
    for(int k=2 ; k<=l ; k++){
      g2 = (2*k-1)/a * g1 - g0;  g0 = g1;  g1 = g2;
      f2 = (2*k-1)/a * f1 - f0;  f0 = f1;  f1 = f2;
    }
    h = complex<double>(g2,f2);
    d = complex<double>(g0-l/a*g1, f0-l/a*f1);
  }

  wfn->d = a*d/h;

  double phase = atan( imag(h) / real(h) );

  return(phase);
}


/**********************************************************/
/*     Penetrability at Resoance Energy                   */
/**********************************************************/
complex<double> gfrLfunction(const int l, const double er, const double mu, const double ap)
{
  Wfunc  wfn;
  double alpha = gcKfactor * sqrt(fabs(er) * 1.0e-06 * mu) * ap;

  complex<double> q(0.0, alpha);
  if(l > 0){
    gfrPenetrability(l,alpha,&wfn);
    q = wfn.d;
  }

  return(q);
}


