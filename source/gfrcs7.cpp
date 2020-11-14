/******************************************************************************/
/**     GFR, Calculate R-Matrix Limited Cross Sections                       **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "coulomb.h"
#include "constant.h"
#include "matrix.h"

static int    RMLLoadParticlePairs (int, System *, ParPair *, ENDF *);
static int    RMLLoadRMLParameters (int, System *, RMLParameter *, ENDF *);
static void   RMLMainCalc (const double, System *, ParPair *, RMLParameter *, GFRcross *);
static void   RMLMatrices (const double, const int, ParPair *, RMLParameter *);
static void   RMLCrossSection (const int, RMLChannel *, GFRcross *);

static double RMLIncidentChannel (const double, System *, ParPair *);
static void   RMLStoreChannelParameter (System *, ParPair *, RMLParameter *, RMLChannel *);
static void   RMLStorePenetrability (System *, ParPair *, RMLParameter *, RMLChannel *);
static void   RMLStorePhaseShift (ParPair *, RMLParameter *, RMLChannel *);
static int    RMLArrangeMatrix (const int, ParPair *, RMLParameter *);

extern Smatrix Smat;

static complex<double> *rm, *sm, *wm, *phi0, *phiC;
static double *gm, **pen;
static ChannelWaveFunc *wf;
static int *mtid, *dptr, msize = 0;

#undef DEBUG_INCH
#undef DEBUG_CHANNEL
#undef DEBUG_PHASE
#undef DEBUG_WIDTH
#undef DEBUG_MATRIX
#undef DEBUG_PEN

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection7(const int ner, const double elab, System *sys, ENDF *lib)
{
  GFRcross sig;
  RMLParameter *res;
  ParPair *ppr;

  /*** two particle pair data */
  ppr = new ParPair [MAX_PAIRS];
  RMLLoadParticlePairs(sys->idx[ner],sys,ppr,lib);

  /*** resonance parameters for each spin group */
  res = new RMLParameter [sys->nj];
  RMLLoadRMLParameters(sys->idx[ner],sys,res,lib);

  /*** look for max number of channels for matrix size N(N+1)/2 */
  int chmax = 0;
  for(int j=0 ; j<sys->nj ; j++) if(chmax < res[j].nchannel) chmax = res[j].nchannel;

  /*** memory allocation */
  msize = chmax * (chmax+1) / 2;

  wf = new ChannelWaveFunc [chmax];
  rm = new complex<double> [msize];   // R-matrix
  sm = new complex<double> [msize];   // S-matrix
  wm = new complex<double> [msize];   // W = I +2i X
  gm = new double [msize];            // Gammas
  mtid = new int [chmax];             // MT numbers for each channel
  dptr = new int [chmax];             // pointer to the ENDF data row
  pen = new double * [chmax];         // penetrability [c x Nres]
  for(int i=0 ; i<chmax ; i++){
    pen[i] = new double [MAX_RESONANCE];
  }
  phi0 = new complex<double> [chmax]; // hard-sphare phase factor
  phiC = new complex<double> [chmax]; // Coulomb phase factor

  /*** calculate cross section */
  sig.memalloc(sys->npair);
  sig.clear();
  RMLMainCalc(elab,sys,ppr,res,&sig);

  Pcross s;

  s.energy  = elab;
  s.total   = sig.sum();
  s.elastic = sig.get(2);
  s.capture = sig.get(102);
  s.fission = sig.get(18);

  /*** sum partial proton and alpha cross sections if given */
  s.proton = 0.0;
  if(sig.get(103) > 0.0) s.proton = sig.get(103);
  else if(sig.get(600) > 0.0){
    for(int m = 600 ; m<=649 ; m++) s.proton += sig.get(m);
  }

  s.alpha = 0.0;
  if(sig.get(107) > 0.0) s.alpha = sig.get(107);
  else if(sig.get(800) > 0.0){
    for(int m = 800 ; m<=849 ; m++) s.alpha +=  sig.get(m);
  }

  sig.memfree();

  delete [] ppr;
  delete [] res;

  delete [] wf;
  delete [] rm;
  delete [] sm;
  delete [] wm;
  delete [] gm;
  delete [] mtid;
  delete [] dptr;
  for(int i=0 ; i<chmax ; i++) delete [] pen[i];
  delete [] pen;
  delete [] phi0;
  delete [] phiC;

  return(s);
}


/**********************************************************/
/*      Main Calculatioin of RML Formula                  */
/**********************************************************/
void RMLMainCalc(const double elab, System *sys, ParPair *ppr, RMLParameter *res, GFRcross *sig)
{
  GFRcross z(sys->npair);
  RMLChannel *chn = new RMLChannel [sys->npair];

  /*** set reaction type in GFRcross class */
  for(int p=0 ; p<sys->npair ; p++) sig->type[p] = z.type[p] = ppr[p].mt;

  /*** find incident channel and kinetic parameters */
  double x1 = RMLIncidentChannel(elab,sys,ppr);

  /*** for all spin groups */
  Smat.resetIndex();
  for(int j=0 ; j<sys->nj ; j++){

    /*** channel parameters */
    RMLStoreChannelParameter(sys,ppr,&res[j],chn);

    /*** penetrabilities at resonance energies */
    RMLStorePenetrability(sys,ppr,&res[j],chn);

    /*** hard-sphere phase */
    RMLStorePhaseShift(ppr,&res[j],chn);

    /*** R and S-matrices */
    int nch = res[j].nchannel;
    int mch = nch - 1;  // capture eliminated total channels
    RMLMatrices(elab,mch,ppr,&res[j]);

    /*** cross section */
    RMLCrossSection(mch,chn,&z);

    double x2 = (res[j].j2 + 1.0) * x1; // (2J+1) g Pi/k^2
    for(int i0 = 0 ; i0<nch ; i0++){
      sig->add(mtid[i0], x2 * z.get(mtid[i0]));
    }

    /*** copy S-matrix elements for elastic */
    for(int i=0 ; i<mch ; i++){
      if(mtid[i] == 2){
        int ss = (int)abs(res[j].s2[dptr[i]]) - sys->target_spin2;
        Smat.setElement(res[j].l[dptr[i]], res[j].j2, ss, sm[(i+2)*(i+1)/2-1]);
        Smat.inclIndex();
      }
    }
  }

  delete [] chn;
}


/**********************************************************/
/*      Calculate R and S Matres                          */
/**********************************************************/
void RMLMatrices(const double elab, const int mch, ParPair *ppr, RMLParameter *res)
{
  for(int i=0 ; i<msize ; i++) sm[i] = rm[i] = wm[i] = complex<double>(0.0,0.0);

  /*** for all resonances */
  for(int k=0 ; k<res->nresonance ; k++){

    /*** store width parameters in the matrix */
    RMLArrangeMatrix(k,ppr,res);

    /*** capture width found at the last element */
    int cpt = (mch + 2) * (mch + 1) / 2 - 1; // index of capture width, Gamma_g
    complex<double> w(res->energy[k] - elab, -gm[cpt]/2.0);
    w = 1.0 / w;
    for(int i=0 ; i<mch*(mch+1)/2 ; i++) rm[i] += gm[i] * w; // R-matrix
  }

  /*** W = delta - i/2 R = I - K */
  for(int i=0 ; i<mch*(mch+1)/2 ; i++){
    wm[i] = complex<double>(imag(rm[i])/2.0,-real(rm[i])/2.0);
  }
  for(int i=0 ; i<mch ; i++) wm[(i+2)*(i+1)/2-1] += 1.0;

  /*** (I - K)^{-1} */
  MatrixInverse(mch,wm);

  /*** S-matrix elements = exp^{-(pc + pc')} [2 (I-K)^{-1} - delta(cc')] */
  for(int i0=0 ; i0<mch ; i0++){
    complex<double> p0 = phi0[i0] * phiC[i0];

    for(int i1=0 ; i1<=i0 ; i1++){
      complex<double> p1 = phi0[i1] * phiC[i1];

      int ij = i0*(i0+1)/2 + i1;
      if(i0 == i1) sm[ij] = p0 * p1 * (2.0*wm[ij] - 1.0);
      else         sm[ij] = p0 * p1 *  2.0*wm[ij];
    }
  }

#ifdef DEBUG_MATRIX
  cout << setprecision(12);
  for(int i0=0 ; i0<mch ; i0++){
    for(int i1=0 ; i1<=i0 ; i1++){
      int ij = i0*(i0+1)/2 + i1;
      cout << setw(20) << rm[ij].real() << setw(20) << rm[ij].imag();
//    cout << setw(20) << wm[ij].real() << setw(20) << wm[ij].imag();
//    cout << setw(20) << sm[ij].real() << setw(20) << sm[ij].imag();
    }
    cout << endl;
  }
#endif
}


/**********************************************************/
/*      Calculate Cross Sections from S-Matrix            */
/**********************************************************/
void RMLCrossSection(const int mch, RMLChannel *chn, GFRcross *z)
{
  z->zero();
  double sigtot = 0.0;

  /*** elastic scattering case */
  for(int i0=0 ; i0<mch ; i0++){
    if(mtid[i0] == 2){
      int ij = (i0+2)*(i0+1)/2 - 1;

      /*** add elastic scattering */
      z->add(mtid[i0],norm(phiC[i0] - sm[ij]));

      /*** for neutron, total cross section is defined */
      if(chn[dptr[i0]].coulomb == 0.0) sigtot += (1.0 - real(sm[ij]))*2.0;
    } 
  }

  /*** all other channels */
  for(int i0=0 ; i0<mch ; i0++){    if(mtid[i0] != 2) continue;  // incoming channel = elastic
    for(int i1=0 ; i1<mch ; i1++){  if(mtid[i1] == 2) continue;  // outgoing channel = other

      if(!chn[dptr[i1]].open) continue;

      int ij = i1*(i1+1)/2 + i0;
      double cx = real(sm[ij])*real(sm[ij]) + imag(sm[ij])*imag(sm[ij]);
      z->add(mtid[i1],cx);
    }
  }

  /*** calculate eliminated capture cross section */
  double x = 0.0;
  for(int i0 = 0 ; i0<mch ; i0++){
    if(mtid[i0] != 1) x += z->get(mtid[i0]); // skip total cross section
  }
  z->set(102, sigtot - x);
}


/**********************************************************/
/*      Incident Channal Kinetic Parameters               */
/**********************************************************/
double RMLIncidentChannel(const double elab, System *sys, ParPair *ppr)
{
  /*** elastic channel as the incoming-particle channel */
  for(int p=0 ; p<sys->npair ; p++){
    if(ppr[p].mt == 2){
      /*** incident and target spins */
      sys->incident_spin2 = ppr[p].spin2[0];
      sys->target_spin2   = ppr[p].spin2[1];
      sys->target_parity  = ppr[p].parity[1];

      /*** reduced mass and CMS energy */
      sys->reduced_mass   =  ppr[p].mass[0] * ppr[p].mass[1] 
                          / (ppr[p].mass[0] + ppr[p].mass[1]) * MNEUTRON;
      sys->ecms = elab * ppr[p].mass[1] / (ppr[p].mass[0] + ppr[p].mass[1]);

      sys->wave_number = sqrt(2.0 * AMUNIT * sys->reduced_mass * sys->ecms * 1e-6) / VLIGHT / HBAR;
      break;
    }
  }

  double x = PI / (sys->wave_number*sys->wave_number) * 0.01
                / ((sys->target_spin2 + 1.0) * (sys->incident_spin2 + 1.0));

#ifdef DEBUG_INCH
  cout << "#  2i: " << setw(4) << sys->incident_spin2;
  cout << "   2I: " << setw(4) << sys->target_spin2;
  cout << setprecision(4);
  cout << "   mu: " << setw(11) << sys->reduced_mass;
  cout << " Ecms: " << setw(11) << sys->ecms;
  cout << "    k: " << setw(11) << sys->wave_number << endl;
#endif

  return x;
}


/**********************************************************/
/*      Calculate Channel Parameters                      */
/**********************************************************/
void RMLStoreChannelParameter(System *sys, ParPair *ppr, RMLParameter *res, RMLChannel *chn)
{
  for(int c=0 ; c<res->nchannel ; c++){

    /*** exit channel reduced mass, masses are given as ratios to neutron */
    chn[c].reduced_mass =  ppr[res->pidx[c]].mass[0] * ppr[res->pidx[c]].mass[1] 
                        / (ppr[res->pidx[c]].mass[0] + ppr[res->pidx[c]].mass[1]) * MNEUTRON;

    /*** save mass ratio */
    chn[c].mratio = (ppr[res->pidx[c]].mass[0] + ppr[res->pidx[c]].mass[1]) / ppr[res->pidx[c]].mass[1];

    /*** exit channel energy */
    chn[c].ecms = sys->ecms + ppr[res->pidx[c]].qvalue;
    chn[c].open = (chn[c].ecms > 0.0) ? true : false;

    /*** exit channel wave number, sqrt(c2) is slightly different from 0.21968 */
    double c2 = 2.0 * AMUNIT / (VLIGHTSQ * HBARSQ);
    double k2 = c2 * chn[c].reduced_mass * chn[c].ecms * 1e-6;

    chn[c].wave_number = (chn[c].ecms <= 0.0) ? -sqrt(-k2) : sqrt(k2); 

    /*** Coulomb parameter */
    int zz = ppr[res->pidx[c]].znum[0] * ppr[res->pidx[c]].znum[1];
    chn[c].coulomb = PERMITTIV * COULOMBSQ * zz
                   * sqrt(AMUNIT * chn[c].reduced_mass / (2.0 * abs(chn[c].ecms) * 1e-6)) / VLIGHT / HBAR;
    chn[c].charge = (zz > 0) ? true : false;

    /*** alpha for effective and true channel radii */
    chn[c].alpha_effective = chn[c].wave_number * res->radius_effective[c];
    chn[c].alpha_true      = chn[c].wave_number * res->radius_true[c];

#ifdef DEBUG_CHANNEL
    cout << "# " << setw(3) << c << setw(2) << ppr[res->pidx[c]].znum[0];
    cout << setw(2) << chn[c].open;
    cout << setprecision(3);
    cout <<   " mu: " << setw(10) << chn[c].reduced_mass;
    cout << " Ecms: " << setw(10) << chn[c].ecms;
    cout <<    " k: " << setw(10) << chn[c].wave_number;
    cout <<   " k2: " << setw(10) << k2;
    cout <<  " eta: " << setw(10) << chn[c].coulomb;
    cout << "  Tru: " << setw(10) << res->radius_true[c]      << setw(10) << chn[c].alpha_true;
    cout << "  Eff: " << setw(10) << res->radius_effective[c] << setw(10) << chn[c].alpha_effective << endl;
#endif
  }
}


/**********************************************************/
/*      Calculate Penetrability at Resonance Energy       */
/**********************************************************/
void RMLStorePenetrability(System *sys, ParPair *ppr, RMLParameter *res, RMLChannel *chn)
{
  double c2 = 2.0 * AMUNIT / (VLIGHTSQ * HBARSQ);

  /*** mass ratio, (M + m)/M, for incident channel */
  double mr = 1.0;
  for(int c=0 ; c<res->nchannel ; c++){
    if(ppr[res->pidx[c]].mt == 2){
      mr = chn[c].mratio;
      break;
    }
  }

  for(int c=0 ; c<res->nchannel ; c++){
    int idx = res->pidx[c];

    /*** set default penetration factor */
    for(int k=0 ; k<res->nresonance ; k++) pen[c][k] = 1.0;

    /*** flag for penetrability calculation */
    bool pencalc = false;
    if(ppr[idx].fpen == 1) pencalc = true;
    else if(ppr[idx].fpen == -1) pencalc = false;
    else{
      /*** dont calculate when fission or capture */
      if(ppr[idx].mt == 18 || ppr[idx].mt == 102) pencalc = false;
      else pencalc = true;
    }
    if(!pencalc) continue;

    if(sys->gammaunit_flag == 0){
      for(int k=0 ; k<res->nresonance ; k++){

        /*** resonance energy is LAB */
        double ecm = abs(res->energy[k] / mr + ppr[idx].qvalue);
        double rho = sqrt(c2 * chn[c].reduced_mass * ecm * 1e-6) * res->radius_true[c];
        double eta = 0.0;
        if(chn[c].coulomb != 0.0) eta = chn[c].coulomb * sqrt( abs(chn[c].ecms / ecm) );

        complex<double> q = gfrLfunction(res->l[c],rho,eta);
        pen[c][k] = q.imag();  // P = imag(L)

#ifdef DEBUG_PEN
        cout << setw(3) << c << setw(3) << k << setw(3) << res->l[c];
        cout <<setprecision(3);
        cout << setw(11) << res->energy[k];
        cout << setw(11) << ecm;
        cout << setw(11) << rho;
        cout << setw(11) << eta;
        cout << setw(11) << pen[c][k] << endl;
#endif
      }
    }
  }
}


/**********************************************************/
/*      Calculate Phase Shift                             */
/**********************************************************/
void RMLStorePhaseShift(ParPair *ppr, RMLParameter *res, RMLChannel *chn)
{
  for(int c=0 ; c<res->nchannel ; c++){

    int idx = res->pidx[c];

    if(ppr[idx].mt == 18 || ppr[idx].mt == 102) continue;

    /*** neutron case, wf includes G'+iF' in wf.d */
    if(chn[c].coulomb == 0.0){
      ChannelWaveFunc tmp;
      /*** penetrability calculated with the true radius */
      gfrPenetrability(res->l[c],chn[c].alpha_true,&wf[c]);

      /*** hard-sphare phase by the effective radius */
      gfrPenetrability(res->l[c],chn[c].alpha_effective,&tmp);
      wf[c].setPhase(tmp.H);
      wf[c].setCoulombPhase(0.0);
    }

    /*** charged particle case */
    else{
      /*** C0 = G + iF, C1 = G' + iF' */
      complex<double> C0, C1;
      coulomb(res->l[c],chn[c].alpha_true,chn[c].coulomb,&C0,&C1);
      wf[c].setData(chn[c].alpha_effective,C0,C1);

      /*** hard-sphase phase */
      coulomb(res->l[c],chn[c].alpha_effective,chn[c].coulomb,&C0,&C1);
      wf[c].setPhase(C0);

      /*** Coulomb phase */
      wf[c].setCoulombPhase(coulomb_phaseshift(res->l[c],chn[c].coulomb));
    }

#ifdef DEBUG_PHASE
    cout << setw(3) << c << setw(3) << idx << setw(4) <<  ppr[idx].mt << setw(3) << res->l[c];
    cout << setw(2) << ppr[idx].fpen;
    cout << setprecision(4);
    cout << "   ph0: " << setw(11) << wf[c].p;
    cout << "   eta: " << setw(11) << chn[c].coulomb;
    cout << "   phC: " << setw(11) << wf[c].phaseC.real() << setw(11) << wf[c].phaseC.imag();
    cout << "     P: " << setw(11) << wf[c].P();
    cout << "     S: " << setw(11) << wf[c].S() << endl;
#endif
  }
}


/**********************************************************/
/*      Matrix Element Re-arrangement for RM              */
/**********************************************************/
int RMLArrangeMatrix(const int k, ParPair *ppr, RMLParameter *res)
{
  complex<double> pzero(0.0,0.0);

  /*** diagonal elements */
  /*** elastic should be first */
  int idx = 0;
  for(int c=0 ; c<res->nchannel ; c++){
    if(ppr[res->pidx[c]].mt == 2){
      int i = (idx+2)*(idx+1)/2 - 1;
      gm[i] = res->gamma[c][k] * wf[c].P() /  pen[c][k]; // Gamma x P(E) / P(Eres)
      mtid[idx] = 2;
      dptr[idx] = c;
      phi0[idx] = wf[c].phase;
      phiC[idx] = wf[c].phaseC;
      idx ++;
    }
  }

  /*** other channels */
  for(int c=0 ; c<res->nchannel ; c++){
    if((ppr[res->pidx[c]].mt == 2) || (ppr[res->pidx[c]].mt == 102)) continue;
    int i = (idx+2)*(idx+1)/2 - 1;
    gm[i] = (pen[c][k] == 0.0) ? 0.0 : res->gamma[c][k] * wf[c].P() /  pen[c][k];
    mtid[idx] = ppr[res->pidx[c]].mt;
    dptr[idx] = c;
    phi0[idx] = wf[c].phase;
    phiC[idx] = wf[c].phaseC;
    idx ++;
  }

  /*** capture channel last */
  for(int c=0 ; c<res->nchannel ; c++){
    if(ppr[res->pidx[c]].mt == 102){
      int i = (idx+2)*(idx+1)/2 - 1;
      gm[i] = res->gamma[c][k];
      mtid[idx] = ppr[res->pidx[c]].mt;
      dptr[idx] = c;
      phi0[idx] = pzero;
      phiC[idx] = pzero;
      idx ++;
    }
  }
  int nch = idx;

  /*** off-diagonal elements */
  for(int i=1 ; i<nch ; i++){
    int di = (i+2)*(i+1)/2 - 1; // index for diagonal element
    int si = (gm[di] < 0.0) ? -1 : 1;

    for(int j=0 ; j<=i-1 ; j++){
      int dj = (j+2)*(j+1)/2 - 1;
      int sj = (gm[dj] < 0.0) ? -1 : 1;

      int k = (i+1)*i/2 + j;
      gm[k] = sqrt( abs(gm[di]) * abs(gm[dj]) ) * si * sj;
    }
  }

  /*** make all diagonal elements positive */
  for(int i=0 ; i<nch ; i++){
    int di = (i+2)*(i+1)/2 - 1;
    gm[di] = abs(gm[di]);
  }

#ifdef DEBUG_WIDTH
  double gt = 0.0; // total width
  for(int i=0 ; i<nch ; i++){
    int di = (i+2)*(i+1)/2 - 1;
    cout << setprecision(3);
    cout << setw(4) << mtid[i];
    for(int j=0 ; j<=i ; j++) cout << setw(11) << gm[(i+1)*i/2 + j];
    cout << endl;
    gt += gm[di];
  }
#endif

  return nch;
}


/**********************************************************/
/*     Read Particle Pair Section                         */
/**********************************************************/
int RMLLoadParticlePairs(int idx, System *sys, ParPair *ppr, ENDF *lib)
{
  /*** first line */
  sys->gammaunit_flag = lib->rdata[idx].l1; // IFG: 0, maybe
  sys->format         = lib->rdata[idx].l2; // KRM: 3, Reich-Moore
  sys->nj             = lib->rdata[idx].n1; // NJS: number of J-Pi, Nj
  sys->relativ_flag   = lib->rdata[idx].n2; // KRL: 0, non-relativistic
  idx++;

  sys->npair          = lib->rdata[idx].l1; // NPP: number of pairs
  if(sys->npair > MAX_PAIRS){
    cerr << "too many particle pairs " << sys->npair << endl;
    return(0);
  }

  int k = 0, p[2];
  for(int i=0 ; i<sys->npair ; i++){

    ppr[i].mass[0]  = lib->xptr[idx][k++];
    ppr[i].mass[1]  = lib->xptr[idx][k++];
    ppr[i].znum[0]  = (int)lib->xptr[idx][k++];
    ppr[i].znum[1]  = (int)lib->xptr[idx][k++];
    ppr[i].spin2[0] = (int)(2.0*lib->xptr[idx][k++]);
    ppr[i].spin2[1] = (int)(2.0*lib->xptr[idx][k++]);
    ppr[i].qvalue   = lib->xptr[idx][k++];
    ppr[i].fpen     = (int)lib->xptr[idx][k++];
    ppr[i].fsft     = (int)lib->xptr[idx][k++];
    ppr[i].mt       = (int)lib->xptr[idx][k++];

    /*** determine parity, since parity is given only when I=0 */
    p[0] = (int)lib->xptr[idx][k++];
    p[1] = (int)lib->xptr[idx][k++];

    for(int j=0 ; j<2 ; j++){
      if(ppr[i].spin2[j] == 0) ppr[i].parity[j] = p[j];
      else{
        ppr[i].parity[j] = (ppr[i].spin2[j] < 0) ? -1 : 1;
        ppr[i].spin2[j] = abs(ppr[i].spin2[j]);
      }
    }

    /*** find target spin, look for elastic (MT=2) channel */
    if(ppr[i].mt == 2){
      for(int j=0 ; j<2 ; j++){

        /*** check Mtarg instead of Znum, since sometimes Ztarg is set to zero */
        if(ppr[i].mass[j] > 1.0){ // this should be target
          sys->target_spin2  = ppr[i].spin2[j];
          sys->target_parity = ppr[i].parity[j];
        }
      }
    }
  }

  return sys->npair;
}


/**********************************************************/
/*     Copy All Resoance Parameters                       */
/**********************************************************/
int RMLLoadRMLParameters(int idx, System *sys, RMLParameter *res, ENDF *lib)
{
  idx += 2; // skip first and second CONTs

  int restot = 0;
  for(int j=0 ; j<sys->nj ; j++){

    /*** channel data */
    res[j].j2 = (int)(2.0*lib->rdata[idx].c1);
    int pj    = (int)lib->rdata[idx].c2; //  AJ: J-pi
//  int kbk   = lib->rdata[idx].l1;      // KBK: background R-matrix, maybe zero
//  int kps   = lib->rdata[idx].l2;      // KPS: hard-sphere phase shift specified, maybe zero
    int nch   = lib->rdata[idx].n2;      // NCH: number of channels
    int nres  = lib->rdata[idx+1].l2;    // NRS: number of resonances

    /*** determine resonance parity */
    if(res[j].j2 == 0) res[j].parity = pj;
    else{
      res[j].parity = (res[j].j2 < 0) ? -1 : 1;
      res[j].j2 = abs(res[j].j2);
    }

    res[j].memalloc(nch,nres);
    restot += nres;
    
    /*** for each channel */
    for(int c=0 ; c<nch ; c++){
      int k = 6*c;
      res[j].pidx[c] = lib->xptr[idx][k] - 1;          // PPI: pair index
      res[j].l[c]    = lib->xptr[idx][k+1];            //   L: orbital angular momentum
      res[j].s2[c]   = (int)(2.0*lib->xptr[idx][k+2]); // SCH: channel spin
      res[j].radius_effective[c] = lib->xptr[idx][k+4] * 10.0;  // APE: effective channel radius
      res[j].radius_true[c]      = lib->xptr[idx][k+5] * 10.0;  // APT: true channel radius
    }
    idx ++;

    /*** for each resonance */
    int nline = nch/6 + 1; // required number of lines for each resonance
    for(int i=0 ; i<nres ; i++){
      int k = 6*nline*i;
      res[j].energy[i] = lib->xptr[idx][k];
      for(int c=0 ; c<nch ; c++) res[j].gamma[c][i] = lib->xptr[idx][k+c+1];
    }
    idx ++;
  }

  return restot;
}

