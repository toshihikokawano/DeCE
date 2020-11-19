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

static Pcross RMLCopyCrossSection (const double, GFRcross *);
static int    RMLLoadParticlePairs (int, System *, ENDF *);
static int    RMLLoadResonanceParameters (int, System *, ENDF *);
static void   RMLMainCalc (const double, System *, GFRcross *);
static void   RMLMatrices (const double, const int, RMLParameter *, double **, ChannelWaveFunc *, int *, int *, complex<double> *, complex<double> *, complex<double> *, complex<double> *);
static void   RMLCrossSection (const int, RMLChannel *, GFRcross *, int *, int *, complex<double> *, complex<double> *, complex<double> *);
static int    RMLArrangeMatrix (const int, RMLParameter *, double **, ChannelWaveFunc *, int *, int *, double *, complex<double> *, complex<double> *);
static void   RMLStorePenetrability (System *, RMLParameter *, RMLChannel *, double **);
static double RMLIncidentChannel (const double, System *);
static void   RMLStoreChannelParameter (System *, RMLParameter *, RMLChannel *);
static void   RMLStorePhaseShift (RMLParameter *, RMLChannel *, ChannelWaveFunc *);
static void   RMLAllocateMemory (const int, System *, ENDF *);
static void   RMLFreeMemory (const int);

extern Smatrix Smat;

static RMLParameter *res;
static ParPair *ppr;
static double ***pen;
static int msize = 0, chmax = 0;
static bool dataload = false;

#undef DEBUG_PAIR
#undef DEBUG_RESONANCE
#undef DEBUG_PENETRABILITY
#undef DEBUG_INCIDENT
#undef DEBUG_CHANNEL
#undef DEBUG_PHASE
#undef DEBUG_WIDTH
#undef DEBUG_MATRIX

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection7(const int ner, const double elab, System *sys, ENDF *lib)
{
  /*** when this is the first call, allocate memory, and keep them until the last call */
  if(sys->isFirstCall()){
    if(dataload) RMLFreeMemory(sys->nj);

    RMLAllocateMemory(ner,sys,lib);
    dataload = true;
  }

  /*** calculate cross section */
  GFRcross sig;
  sig.memalloc(sys->npair);
  sig.clear();

  RMLMainCalc(elab,sys,&sig);

  /*** Mapping MT numbers in GFR object to cross section data */
  Pcross s = RMLCopyCrossSection(elab,&sig);

  /*** release allocated memories */
  sig.memfree();
  sys->OnceCalled();
  if(sys->isLastCall()) RMLFreeMemory(sys->nj);

  return(s);
}


/**********************************************************/
/*      Copy Calculated Result to Pcross Object           */
/**********************************************************/
Pcross RMLCopyCrossSection(const double elab, GFRcross *sig)
{
  const double sigcut = 1e-99;
  Pcross s;

  s.clear();
  s.energy  = elab;
  s.total   = sig->sum();     // total cross section
  s.elastic = sig->get(2);    // elastic scattering
  s.capture = sig->get(102);  // capture cross section
  s.fission = sig->get(18);   // fission cross section

  /*** sum partial proton, alpha, and inelastic scattering cross sections if given */
  /*** since we don't know which MT number is assigned to the reaction channel,
       either MT = 4 or MT = 51, 52, ..., scan all MT numbers for inelastic */
  if(sig->get(4) > 0.0) s.inelastic = sig->get(4);
  else if(sig->get(51) > 0.0){
    for(int m = 51 ; m<=91 ; m++) s.inelastic += sig->get(m);
  }

  if(sig->get(103) > 0.0) s.proton = sig->get(103);
  else if(sig->get(600) > 0.0){
    for(int m = 600 ; m<=649 ; m++) s.proton += sig->get(m);
  }

  if(sig->get(107) > 0.0) s.alpha = sig->get(107);
  else if(sig->get(800) > 0.0){
    for(int m = 800 ; m<=849 ; m++) s.alpha += sig->get(m);
  }

  /*** truncate too small cross section */
  if(s.inelastic < sigcut) s.inelastic = 0.0;
  if(s.proton    < sigcut) s.proton = 0.0;
  if(s.alpha     < sigcut) s.alpha = 0.0;

  return s;
}


/**********************************************************/
/*      Main Calculatioin of RML Formula                  */
/**********************************************************/
void RMLMainCalc(const double elab, System *sys, GFRcross *sig)
{
  GFRcross z(sig->getNch());
  RMLChannel      *chn  = new RMLChannel      [chmax]; // channel data
  complex<double> *sm   = new complex<double> [msize]; // S-matrix
  complex<double> *xm   = new complex<double> [msize]; // X-matrix x 2i
  ChannelWaveFunc *wf   = new ChannelWaveFunc [chmax]; // G and F
  complex<double> *phi0 = new complex<double> [chmax]; // hard-sphare phase factor
  complex<double> *phiC = new complex<double> [chmax]; // Coulomb phase factor


  int *mtid = new int [chmax]; // MT numbers for each channel
  int *dptr = new int [chmax]; // pointer to the ENDF data row

  /*** set reaction type in GFRcross class */
  for(int p=0 ; p<sys->npair ; p++) sig->type[p] = z.type[p] = ppr[p].mt;

  /*** find incident channel and kinetic parameters */
  double x1 = RMLIncidentChannel(elab,sys);

  if(sys->isFirstCall()){
    /*** pre-calculate penetrabilities */
    for(int j=0 ; j<sys->nj ; j++){

      /*** channel parameters */
      RMLStoreChannelParameter(sys,&res[j],chn);

      /*** penetrabilities at resonance energies */
      RMLStorePenetrability(sys,&res[j],chn,pen[j]);
    }
  }

  /*** for all spin groups */
  Smat.resetIndex();
  for(int j=0 ; j<sys->nj ; j++){

    /*** channel parameters */
    RMLStoreChannelParameter(sys,&res[j],chn);

    /*** hard-sphere phase */
    RMLStorePhaseShift(&res[j],chn,wf);

    /*** R and S-matrices */
    int nch = res[j].nchannel;
    int mch = nch - 1;  // capture eliminated total channels
    RMLMatrices(elab,mch,&res[j],pen[j],wf,mtid,dptr,sm,xm,phi0,phiC);

    /*** cross section */
    RMLCrossSection(mch,chn,&z,mtid,dptr,sm,xm,phiC);

    double x2 = (res[j].j2 + 1.0) * x1; // (2J+1) g Pi/k^2

    for(int i0 = 0 ; i0<z.getNch() ; i0++){
      int mt = z.type[i0];
      sig->add(mt, x2 * z.get(mt));
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
  delete [] sm;
  delete [] xm;
  delete [] wf;
  delete [] phi0;
  delete [] phiC;
  delete [] mtid;
  delete [] dptr;
}


/**********************************************************/
/*      Calculate R and S Matres                          */
/**********************************************************/
#undef METHOD_A
void RMLMatrices(const double elab, const int mch, RMLParameter *r, double **p, ChannelWaveFunc *wf, int *mtid, int *dptr, complex<double> *sm, complex<double> *xm, complex<double> *phi0, complex<double> *phiC)
{
  complex<double> *wm   = new complex<double> [msize];   // W = I +2i X
  complex<double> *rm   = new complex<double> [msize];   // R-matrix
  double          *gm   = new double          [msize];   // Gammas
  double          *pm   = new double          [msize];   // P

  for(int i=0 ; i<msize ; i++){
    sm[i] = rm[i] = wm[i] = complex<double>(0.0,0.0);
    gm[i] = 0.0;
  }

  /*** for all resonances */
  for(int k=0 ; k<r->nresonance ; k++){

    /*** store width parameters in the matrix */
    RMLArrangeMatrix(k,r,p,wf,mtid,dptr,gm,phi0,phiC);

    /*** capture width found at the last element */
    int cpt = (mch + 2) * (mch + 1) / 2 - 1; // index of capture width, Gamma_g
    complex<double> w(r->energy[k] - elab, -gm[cpt]/2.0);
    w = 1.0 / w;
    for(int i=0 ; i<mch*(mch+1)/2 ; i++) rm[i] += gm[i] * w; // R-matrix
  }

#ifdef METHOD_A
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
#else
  for(int i=0 ; i<mch ; i++) pm[i] = wf[dptr[i]].P();

  /*** W = {L^{-1} - R}^{-1} */
  for(int i=0 ; i<mch ; i++){
    for(int j=0 ; j<=i ; j++){
      int ij = i*(i+1)/2 + j;
      double p = pm[i] * pm[j];

      rm[ij] = (p != 0.0) ? rm[ij] / (2.0 * sqrt(p)) : 0.0;

      if(i == j) wm[ij] = 1.0/complex<double>(0.0,pm[i]) - rm[ij];
      else       wm[ij] = - rm[ij];
    }
  }
  MatrixInverse(mch,wm);

  /*** X = W R */
  for(int i=0 ; i<mch ; i++){
    for(int j=0 ; j<=i ; j++){
      int ij = i*(i+1)/2 + j;
      xm[ij] = complex<double>(0.0,0.0);
      for(int k=0 ; k<mch ; k++){
        int ik = i*(i+1)/2 + k;  if(k > i) ik = k*(k+1)/2 + i;
        int kj = k*(k+1)/2 + j;  if(j > k) kj = j*(j+1)/2 + k;
        xm[ij] += wm[ik] * rm[kj];
      }
    }
  }

  /*** 2i X = 2i P^{1/2}L^{-1} WR P^{1/2} */
  for(int i=0 ; i<mch ; i++){
    double ai = (pm[i] != 0.0) ? 2.0/sqrt(pm[i]) : 0.0;
    for(int j=0 ; j<=i ; j++){
      double aj = (pm[j] != 0.0) ? sqrt(pm[j]) : 0.0;
      int ij = i*(i+1)/2 + j;
      xm[ij] *= ai * aj;
    }
  }

  /*** S-matrix elements = exp^{-(pc + pc')} [delta(cc') + 2iX] */
  for(int i0=0 ; i0<mch ; i0++){
    complex<double> p0 = phi0[i0] * phiC[i0];

    for(int i1=0 ; i1<=i0 ; i1++){
      complex<double> p1 = phi0[i1] * phiC[i1];

      int ij = i0*(i0+1)/2 + i1;
      if(i0 == i1) sm[ij] = p0 * p1 * (xm[ij] + 1.0);
      else         sm[ij] = p0 * p1 *  xm[ij];
    }
  }

  for(int i0=0 ; i0<mch ; i0++){
    for(int i1=0 ; i1<=i0 ; i1++){
      int ij = i0*(i0+1)/2 + i1;
      xm[ij] = xm[ij] / complex<double>(0.0,2.0);
    }
  }

#endif

#ifdef DEBUG_MATRIX
  cout << setprecision(12);
  for(int i0=0 ; i0<mch ; i0++){
    for(int i1=0 ; i1<=i0 ; i1++){
      int ij = i0*(i0+1)/2 + i1;
//    cout << " " << setw(20) << rm[ij].real() << setw(20) << rm[ij].imag();
//    cout << " " << setw(20) << wm[ij].real() << setw(20) << wm[ij].imag();
      cout << " " << setw(20) << sm[ij].real() << setw(20) << sm[ij].imag();
    }
    cout << endl;
  }
#endif

  delete [] wm;
  delete [] rm;
  delete [] gm;
  delete [] pm;
}


/**********************************************************/
/*      Calculate Cross Sections from S-Matrix            */
/**********************************************************/
void RMLCrossSection(const int mch, RMLChannel *chn, GFRcross *z, int *mtid, int *dptr, complex<double> *sm, complex<double> *xm, complex<double> *phiC)
{
  z->zero();
  double sigtot = 0.0;

  /*** elastic scattering case */
  double sigc = 0.0;
  for(int i0=0 ; i0<mch ; i0++){
    if(mtid[i0] == 2){
      int ii = (i0+2)*(i0+1)/2 - 1;

      /*** add elastic scattering */
      z->add(mtid[i0],norm(phiC[i0] - sm[ii]));

      /*** for neutron, total cross section is defined */
      if(chn[dptr[i0]].coulomb == 0.0) sigtot += (1.0 - real(sm[ii]))*2.0;

      double cx = 0.0;
      for(int i1=0 ; i1<mch ; i1++){
        int ij = i1*(i1+1)/2 + i0;
        cx += norm(xm[ij]);
      }
      sigc += 4*(xm[ii].imag() - cx);
    } 
  }

  /*** eliminated capture cross section */
  z->set(102, sigc);

  /*** all other channels */
  for(int i0=0 ; i0<mch ; i0++){    if(mtid[i0] != 2) continue;  // incoming channel = elastic
    for(int i1=0 ; i1<mch ; i1++){  if(mtid[i1] == 2) continue;  // outgoing channel = other
      if(!chn[dptr[i1]].open) continue;
      z->add(mtid[i1],norm(sm[i1*(i1+1)/2 + i0]));
    }
  }
}


/**********************************************************/
/*      Matrix Element Re-arrangement for RM              */
/**********************************************************/
int RMLArrangeMatrix(const int k, RMLParameter *r, double **p, ChannelWaveFunc *wf, int *mtid, int *dptr, double *gm, complex<double> *phi0, complex<double> *phiC)
{
  complex<double> pzero(0.0,0.0);

  /*** diagonal elements */
  /*** elastic should be first */
  int idx = 0;
  for(int c=0 ; c<r->nchannel ; c++){
    if(ppr[r->pidx[c]].mt == 2){
      int i = (idx+2)*(idx+1)/2 - 1;
      gm[i] = r->gamma[c][k] * wf[c].P() /  p[c][k]; // Gamma x P(E) / P(Eres)
      mtid[idx] = 2;
      dptr[idx] = c;
      phi0[idx] = wf[c].phase;
      phiC[idx] = wf[c].phaseC;
      idx ++;
    }
  }

  /*** other channels */
  for(int c=0 ; c<r->nchannel ; c++){
    if((ppr[r->pidx[c]].mt == 2) || (ppr[r->pidx[c]].mt == 102)) continue;
    int i = (idx+2)*(idx+1)/2 - 1;
    gm[i] = (p[c][k] == 0.0) ? 0.0 : r->gamma[c][k] * wf[c].P() /  p[c][k];
    mtid[idx] = ppr[r->pidx[c]].mt;
    dptr[idx] = c;
    phi0[idx] = wf[c].phase;
    phiC[idx] = wf[c].phaseC;
    idx ++;
  }

  /*** capture channel last */
  for(int c=0 ; c<r->nchannel ; c++){
    if(ppr[r->pidx[c]].mt == 102){
      int i = (idx+2)*(idx+1)/2 - 1;
      gm[i] = r->gamma[c][k];
      mtid[idx] = ppr[r->pidx[c]].mt;
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
/*      Calculate Penetrability at Resonance Energy       */
/**********************************************************/
void RMLStorePenetrability(System *sys, RMLParameter *r, RMLChannel *chn, double **p)
{
  double c2 = 2.0 * AMUNIT / (VLIGHTSQ * HBARSQ);

  /*** mass ratio, (M + m)/M, for incident channel */
  double mr = 1.0;
  for(int c=0 ; c<r->nchannel ; c++){
    if(ppr[r->pidx[c]].mt == 2){
      mr = chn[c].mratio;
      break;
    }
  }

  for(int c=0 ; c<r->nchannel ; c++){
    int idx = r->pidx[c];

    /*** set default penetration factor */
    for(int k=0 ; k<r->nresonance ; k++) p[c][k] = 1.0;

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
      for(int k=0 ; k<r->nresonance ; k++){

        /*** resonance energy is LAB */
        double ecm = abs(r->energy[k] / mr + ppr[idx].qvalue);
        double rho = sqrt(c2 * chn[c].reduced_mass * ecm * 1e-6) * r->radius_true[c];
        double eta = 0.0;
        if(chn[c].coulomb != 0.0) eta = chn[c].coulomb * sqrt( abs(chn[c].ecms / ecm) );

        complex<double> q = gfrLfunction(r->l[c],rho,eta);
        p[c][k] = q.imag();  // P = imag(L)

#ifdef DEBUG_PENETRABILITY
        cout << setw(3) << c << setw(3) << k << setw(3) << r->l[c];
        cout <<setprecision(3);
        cout << setw(11) << r->energy[k];
        cout << setw(11) << ecm;
        cout << setw(11) << rho;
        cout << setw(11) << eta;
        cout << setw(11) << p[c][k] << endl;
#endif
      }
    }
  }
}


/**********************************************************/
/*      Incident Channal Kinetic Parameters               */
/**********************************************************/
double RMLIncidentChannel(const double elab, System *sys)
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

#ifdef DEBUG_INCIDENT
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
void RMLStoreChannelParameter(System *sys, RMLParameter *r, RMLChannel *chn)
{
  for(int c=0 ; c<r->nchannel ; c++){

    /*** exit channel reduced mass, masses are given as ratios to neutron */
    chn[c].reduced_mass =  ppr[r->pidx[c]].mass[0] * ppr[r->pidx[c]].mass[1] 
                        / (ppr[r->pidx[c]].mass[0] + ppr[r->pidx[c]].mass[1]) * MNEUTRON;

    /*** save mass ratio */
    chn[c].mratio = (ppr[r->pidx[c]].mass[0] + ppr[r->pidx[c]].mass[1]) / ppr[r->pidx[c]].mass[1];

    /*** exit channel energy */
    chn[c].ecms = sys->ecms + ppr[r->pidx[c]].qvalue;
    chn[c].open = (chn[c].ecms > 0.0) ? true : false;

    /*** exit channel wave number, sqrt(c2) is slightly different from 0.21968 */
    double c2 = 2.0 * AMUNIT / (VLIGHTSQ * HBARSQ);
    double k2 = c2 * chn[c].reduced_mass * chn[c].ecms * 1e-6;

    chn[c].wave_number = (chn[c].ecms <= 0.0) ? -sqrt(-k2) : sqrt(k2); 

    /*** Coulomb parameter */
    int zz = ppr[r->pidx[c]].znum[0] * ppr[r->pidx[c]].znum[1];
    chn[c].coulomb = PERMITTIV * COULOMBSQ * zz
                   * sqrt(AMUNIT * chn[c].reduced_mass / (2.0 * abs(chn[c].ecms) * 1e-6)) / VLIGHT / HBAR;
    chn[c].charge = (zz > 0) ? true : false;

    /*** alpha for effective and true channel radii */
    chn[c].alpha_effective = chn[c].wave_number * r->radius_effective[c];
    chn[c].alpha_true      = chn[c].wave_number * r->radius_true[c];

#ifdef DEBUG_CHANNEL
    cout << "# " << setw(3) << c << " z:" << setw(2) << ppr[r->pidx[c]].znum[0];
    cout << setw(2) << ( (chn[c].open) ? 'o' : 'x' );
    cout << setprecision(3);
    cout <<   " mu: " << setw(10) << chn[c].reduced_mass;
    cout << " Ecms: " << setw(10) << chn[c].ecms;
    cout <<    " k: " << setw(10) << chn[c].wave_number;
    cout <<  " eta: " << setw(10) << chn[c].coulomb;
    cout << "  Tru: " << setw(10) << r->radius_true[c]      << setw(11) << chn[c].alpha_true;
    cout << "  Eff: " << setw(10) << r->radius_effective[c] << setw(11) << chn[c].alpha_effective << endl;
#endif
  }
}


/**********************************************************/
/*      Calculate Phase Shift                             */
/**********************************************************/
void RMLStorePhaseShift(RMLParameter *r, RMLChannel *chn, ChannelWaveFunc *wf)
{
  for(int c=0 ; c<r->nchannel ; c++){

    int idx = r->pidx[c];

    if(ppr[idx].mt == 18 || ppr[idx].mt == 102) continue;

    /*** open neutron channel, wf includes G'+iF' in wf.d */
    if(chn[c].coulomb == 0.0 && chn[c].open){
      ChannelWaveFunc tmp;
      /*** penetrability calculated with the true radius */
      gfrPenetrability(r->l[c],chn[c].alpha_true,&wf[c]);

      /*** hard-sphare phase by the effective radius */
      gfrPenetrability(r->l[c],chn[c].alpha_effective,&tmp);
      wf[c].setPhase(tmp.H);
      wf[c].setCoulombPhase(0.0);
    }

    /*** charged particle or closed neutron channel case */
    else{
      /*** C0 = G + iF, C1 = G' + iF' */
      complex<double> C0, C1;
      coulomb(r->l[c],chn[c].alpha_true,chn[c].coulomb,&C0,&C1);
      wf[c].setData(chn[c].alpha_effective,C0,C1);

      /*** hard-sphase phase */
      coulomb(r->l[c],chn[c].alpha_effective,chn[c].coulomb,&C0,&C1);
      wf[c].setPhase(C0);

      /*** Coulomb phase */
      if(chn[c].coulomb != 0.0) wf[c].setCoulombPhase(coulomb_phaseshift(r->l[c],chn[c].coulomb));
      else wf[c].setCoulombPhase(0.0);
    }

#ifdef DEBUG_PHASE
    cout << setw(3) << c << setw(3) << idx << setw(4) <<  ppr[idx].mt << setw(3) << r->l[c];
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
/*     Read Particle Pair Section                         */
/**********************************************************/
int RMLLoadParticlePairs(int idx, System *sys, ENDF *lib)
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
#ifdef DEBUG_PAIR
  for(int i=0 ; i<sys->npair ; i++){
    cout << setprecision(3);
    cout << "# MT : " << setw(4) << ppr[i].mt << "  Q: " << setw(11) << ppr[i].qvalue << endl;
    for(int j=0 ; j<2 ; j++){
      cout << " M: " << setw(11) << ppr[i].mass[j] << " Z: " << setw(4) << ppr[i].znum[j];
      cout << " I: " << setw(4) << ppr[i].spin2[j] << " p: " << setw(4) << ppr[i].parity[j] << endl;
    }
  }
#endif

  return sys->npair;
}


/**********************************************************/
/*     Copy All Resoance Parameters                       */
/**********************************************************/
int RMLLoadResonanceParameters(int idx, System *sys, ENDF *lib)
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
    restot += res[j].nresonance;
    
    /*** for each channel */
    for(int c=0 ; c<res[j].nchannel ; c++){
      int k = 6*c;
      res[j].pidx[c] = lib->xptr[idx][k] - 1;          // PPI: pair index
      res[j].l[c]    = lib->xptr[idx][k+1];            //   L: orbital angular momentum
      res[j].s2[c]   = (int)(2.0*lib->xptr[idx][k+2]); // SCH: channel spin
      res[j].radius_effective[c] = lib->xptr[idx][k+4] * 10.0;  // APE: effective channel radius
      res[j].radius_true[c]      = lib->xptr[idx][k+5] * 10.0;  // APT: true channel radius
    }
    idx ++;

    /*** for each resonance */
    int nline = res[j].nchannel/6 + 1; // required number of lines for each resonance
    for(int i=0 ; i<res[j].nresonance ; i++){
      int k = 6*nline*i;
      res[j].energy[i] = lib->xptr[idx][k];
      for(int c=0 ; c<res[j].nchannel ; c++) res[j].gamma[c][i] = lib->xptr[idx][k+c+1];
    }
    idx ++;
  }

#ifdef DEBUG_RESONANCE
  for(int j=0 ; j<sys->nj ; j++){
    cout << "# J2 : " << setw(4) << res[j].j2 << setw(6) << res[j].nresonance << endl;
    cout << setprecision(3);
    for(int i=0 ; i<res[j].nresonance ; i++){
      cout << setw(11) << res[j].energy[i];
      for(int c=0 ; c<res[j].nchannel ; c++){
        cout << setw(11) << res[j].gamma[c][i];
      }
      cout << endl;
    }
  }
#endif

  return restot;
}


/**********************************************************/
/*      Allocate Memory and Load Parameters               */
/**********************************************************/
void RMLAllocateMemory(const int ner, System *sys, ENDF *lib)
{
  if(dataload) return;

  int idx =sys->idx[ner] + 1;

  /*** two particle pair data */
  ppr = new ParPair [MAX_PAIRS];
  RMLLoadParticlePairs(idx,sys,lib);

  /*** resonance parameters for each spin group */
  res = new RMLParameter [sys->nj];
  RMLLoadResonanceParameters(idx,sys,lib);

  /*** look for max number of channels for matrix size N(N+1)/2 */
  for(int j=0 ; j<sys->nj ; j++) if(chmax < res[j].nchannel) chmax = res[j].nchannel;

  /*** memory allocation */
  msize = chmax * (chmax+1) / 2;

  pen = new double ** [sys->nj];      // penetrability [Nj x c x Nres]
  for(int j=0 ; j<sys->nj; j++){
    pen[j] = new double * [chmax];
    for(int i=0 ; i<chmax ; i++){
      pen[j][i] = new double [res[j].nresonance];
    }
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void RMLFreeMemory(const int nj)
{
  if(!dataload) return;

  delete [] ppr;
  delete [] res;

  for(int j=0 ; j<nj; j++){
    for(int i=0 ; i<chmax ; i++) delete [] pen[j][i];
    delete [] pen[j];
  }
  delete [] pen;

  dataload = false;
}
