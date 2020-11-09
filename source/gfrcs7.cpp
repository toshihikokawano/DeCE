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

static int gfrLoadParticlePairs (int, System *, ParPair *, ENDF *);
static int gfrLoadRMLParameters (int, System *, RMLParameter *, ENDF *);
static void gfrRMatrixLimited (const double, System *, ParPair *, RMLParameter *);
static void gfrStoreRMLPenetrability (System *, ParPair *, RMLParameter *);
static int arrange_matrixRML (const int, double **, Wfunc *, ParPair *, RMLParameter *);

extern Smatrix Smat;

static complex<double> *rm, *sm, *wm, *phi0, *phiC;
static double *gm, **pen;
static Wfunc *wf;
static int *mtid, *dptr, msize = 0;
static GFRcross sig;

#undef DEBUG

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection7(const int ner, const double elab, System *sys, ENDF *lib)
{
  RMLParameter *res;
  ParPair *ppr;

  /*** two particle pair data */
  ppr = new ParPair [MAX_PAIRS];
  gfrLoadParticlePairs(sys->idx[ner],sys,ppr,lib);

  /*** resonance parameters for each spin group */
  res = new RMLParameter [sys->nj];
  gfrLoadRMLParameters(sys->idx[ner],sys,res,lib);

  /*** look for max number of channels for matrix size N(N+1)/2 */
  int chmax = 0;
  for(int j=0 ; j<sys->nj ; j++) if(chmax < res[j].nch) chmax = res[j].nch;

  /*** memory allocation */
  msize = chmax * (chmax+1) / 2;

  wf = new Wfunc [chmax];
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
  gfrRMatrixLimited(elab,sys,ppr,res);

  Pcross s;

  s.energy  = elab;
  s.total   = sig.total();
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
  for(int i=0 ; i<chmax ; i++){
    delete [] pen[i];
  }
  delete [] pen;
  delete [] phi0;
  delete [] phiC;

  return(s);
}


/**********************************************************/
/*      Main Calculatioin of RML Formula                  */
/**********************************************************/
void gfrRMatrixLimited(const double elab, System *sys, ParPair *ppr, RMLParameter *res)
{
  GFRcross z(sys->npair);
  int nch, mch, cpt;

  /*** elastic channel as the incoming-particle channel */
  for(int p=0 ; p<sys->npair ; p++){
    if(ppr[p].mt == 2){
      /*** incident and target spins */
      sys->incident_spin2 = ppr[p].spin2[0];
      sys->target_spin2   = ppr[p].spin2[1];
      sys->target_parity  = ppr[p].parity[1];

      /*** reduced mass and CMS energy */
      sys->reduced_mass   =  ppr[p].mass[0] * ppr[p].mass[1] 
                          / (ppr[p].mass[0] + ppr[p].mass[1]);
      sys->ecms = elab * sys->reduced_mass;
      sys->wave_number = sqrt(sys->reduced_mass) * gcKfactor * sqrt(sys->ecms * 1.0e-06);
      break;
    }
  }

  /*** set reaction type in GFRcross class */
  for(int p=0 ; p<sys->npair ; p++) sig.type[p] = z.type[p] = ppr[p].mt;

  /*** kinematic factor */
  double x1 = PI / (sys->wave_number*sys->wave_number) * 0.01
                 / ((sys->target_spin2 + 1.0) * (sys->incident_spin2 + 1.0));

  /*** for all spin groups */
  Smat.resetIndex();
  for(int j=0 ; j<sys->nj ; j++){

    /*** (2J+1) g Pi/k^2 */
    double x3 = (res[j].j2 + 1.0) * x1;

    /*** channel radius and hard-sphere phase */
    gfrStoreRMLPenetrability(sys,ppr,&res[j]);

    /*** R-matrix elements */
    mch = 0; nch = 0; cpt = 0;
    for(int i=0 ; i<msize ; i++) sm[i] = rm[i] = wm[i] = complex<double>(0.0,0.0);
    for(int k=0 ; k<res[j].nres ; k++){

      /*** store width parameters in the matrix */
      nch = arrange_matrixRML(k,pen,wf,ppr,&res[j]);
      mch = nch - 1;                 // capture eliminated total channels
      cpt = (nch + 1) * nch / 2 - 1; // index of capture width, Gamma_g

      /*** capture width found at the last element */
      complex<double> w(res[j].er[k] - elab, -gm[cpt]/2.0);
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

    /*** cross section */
    z.clear();
    double sigtot = 0.0;

    /*** elastic scattering case */
    for(int i0=0 ; i0<mch ; i0++){
      if(mtid[i0] == 2){
        int ij = (i0+2)*(i0+1)/2 - 1;
        z.add(mtid[i0],norm(phiC[i0] - sm[ij]));
        /*** for neutron, total cross section is defined */
        if(phiC[i0].imag() == 0.0) sigtot += (1.0 - real(sm[ij]))*2.0;
      } 
    }

    /*** all other channels */
    for(int i0=0 ; i0<mch ; i0++){    if(mtid[i0] != 2) continue;  // incoming channel = elastic
      for(int i1=0 ; i1<mch ; i1++){  if(mtid[i1] == 2) continue;  // outgoing channel = other

        int ij = i1*(i1+1)/2 + i0;
        double cx = real(sm[ij])*real(sm[ij]) + imag(sm[ij])*imag(sm[ij]);
        z.add(mtid[i1],cx);
      }
    }

    /*** calculate eliminated capture cross section */
    double x = 0.0;
    for(int i0 = 0 ; i0<mch ; i0++){
      if(mtid[i0] != 1) x += z.get(mtid[i0]); // skip total cross section
    }
    z.set(102, sigtot - x);

    for(int i0 = 0 ; i0<nch ; i0++){
      sig.add(mtid[i0], x3 * z.get(mtid[i0]));
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
}


/**********************************************************/
/*      Calculate Phase Shift and Penetrability at Res.   */
/**********************************************************/
void gfrStoreRMLPenetrability(System *sys, ParPair *ppr, RMLParameter *res)
{
  /*** set default penetration factor */
  for(int c=0 ; c<res->nch ; c++){
    for(int k=0 ; k<res->nres ; k++) pen[c][k] = 1.0;
  }

  for(int c=0 ; c<res->nch ; c++){

    int    l   = res->l[c];    // orbital angular momentum
    double ape = res->rade[c]; // effective channel radius
    double apt = res->radt[c]; // true channel radius

    double mu  =  ppr[res->pidx[c]].mass[0] * ppr[res->pidx[c]].mass[1] 
               / (ppr[res->pidx[c]].mass[0] + ppr[res->pidx[c]].mass[1]); // exit channel reduced mass
    double eb  = sys->ecms + ppr[res->pidx[c]].qvalue;                    // exit channel energy
    double k2  = gcKfactor*gcKfactor * mu * eb * 1.0e-06;                 // exit channel wave number squared
    double kw  = (eb <= 0.0) ? -sqrt(-k2) : sqrt(k2);                     // exit wave number
    int    zzp = ppr[res->pidx[c]].znum[0] * ppr[res->pidx[c]].znum[1];   // charge product
    double eta = PERMITTIV * COULOMBSQ * zzp * sqrt(AMUNIT * mu / (2.0 * abs(eb) * 1e-6)) / VLIGHT / HBAR;

    /*** flag for penetrability calculation */
    bool pencalc = false;
    if(ppr[res->pidx[c]].fpen == 1) pencalc = true;
    else if(ppr[res->pidx[c]].fpen == -1) pencalc = false;
    else{
      /*** dont calculate when fission or capture */
      if(ppr[res->pidx[c]].mt == 18 || ppr[res->pidx[c]].mt == 102) pencalc = false;
      else pencalc = true;
    }
    if(!pencalc) continue;


    /*** potential scattering sphere phase shift, use effective radius APE */
    double phi = 0.0, wc = 0.0;
    double rho = kw * ape;
    /*** neutron case, wf includes G'+iF' in wf.d */
    if(zzp == 0) phi = gfrPenetrability(l,rho,&wf[c]);
    /*** charged particle case */
    else{
      complex<double> C0, C1;
      coulomb(l,rho,eta,&C0,&C1);  // C0 = G + iF, C1 = G' + iF'
      /*** phi = acos G / sqrt(G^2 + F^2) */
      phi = acos(C0.real() / abs(C0));
      wc  = coulomb_phaseshift(l,eta);
      wf[c].d = rho * C1 / C0;
    }
    wf[c].phase  = complex<double>(cos(-phi), sin(-phi));  // exp^i{-phi}
    wf[c].phaseC = complex<double>(cos( wc ), sin(wc  ));  // exp^{iw}

    /*** penetrability at the resonance energy, calculate L(Er) = S(Er) + iP(Er), use true radius APT */
    if(sys->gammaunit_flag == 0){
      for(int k=0 ; k<res->nres ; k++){
        complex<double> q;
        if(zzp == 0) q = gfrLfunction(l,abs(res->er[k]),mu,apt);
        else         q = gfrLfunctionCoul(l,abs(res->er[k]),mu,apt,eta);
        pen[c][k] = q.imag();
      }
    }
#ifdef DEBUG
    cout << setw(3) << c << setw(3) << res->pidx[c] << setw(4) <<  ppr[res->pidx[c]].mt << setw(3) << l;
    cout << setw(2) << ppr[res->pidx[c]].fpen;
    cout << setprecision(3) << setw(11) << ape << setw(11) << apt;
    cout << "  mu:" << setw(11) << mu;
    cout << "   E:" << setw(11) << eb;
    cout << "   k:" << setw(11) << kw;
    cout << " eta:" << setw(11) << eta;
    cout << " ph0:" << setw(11) << - phi;
    cout << " phC:" << setw(11) << wc;
    cout << "   P:" << setw(11) << wf[c].d.imag() << endl;
#endif
  }
}


/**********************************************************/
/*      Matrix Element Re-arrangement for RM              */
/**********************************************************/
int arrange_matrixRML(const int k, double **pen, Wfunc *wf, ParPair *ppr, RMLParameter *res)
{
  complex<double> pzero(0.0,0.0);

  /*** diagonal elements */
  /*** elastic should be first */
  int idx = 0;
  for(int c=0 ; c<res->nch ; c++){
    if(ppr[res->pidx[c]].mt == 2){
      int i = (idx+2)*(idx+1)/2 - 1;
      gm[i] = res->gam[c][k] * imag(wf[c].d)/ pen[c][k]; // Gamma x P(E) / P(Eres)
      mtid[idx] = 2;
      dptr[idx] = c;
      phi0[idx] = wf[c].phase;
      phiC[idx] = wf[c].phaseC;
      idx ++;
    }
  }

  /*** other channels */
  for(int c=0 ; c<res->nch ; c++){
    if((ppr[res->pidx[c]].mt == 2) || (ppr[res->pidx[c]].mt == 102)) continue;
    int i = (idx+2)*(idx+1)/2 - 1;
    gm[i] = res->gam[c][k];
    mtid[idx] = ppr[res->pidx[c]].mt;
    dptr[idx] = c;
    phi0[idx] = wf[c].phase;
    phiC[idx] = wf[c].phaseC;
    idx ++;
  }

  /*** capture channel last */
  for(int c=0 ; c<res->nch ; c++){
    if(ppr[res->pidx[c]].mt == 102){
      int i = (idx+2)*(idx+1)/2 - 1;
      gm[i] = res->gam[c][k];
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
  double gt = 0.0; // total width
  for(int i=0 ; i<nch ; i++){
    int di = (i+2)*(i+1)/2 - 1;
    gm[di] = abs(gm[di]);
    gt += gm[di];
//    cout << setw(4) << mtid[i] << setw(11) << gt;
//    for(int j=0 ; j<=i ; j++) cout << setw(11) << gm[(i+1)*i/2 + j];
//    cout << endl;
  }

  return(nch);
}


/**********************************************************/
/*     Read Particle Pair Section                         */
/**********************************************************/
int gfrLoadParticlePairs(int idx, System *sys, ParPair *ppr, ENDF *lib)
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

  return(sys->npair);
}



/**********************************************************/
/*     Copy All Resoance Parameters                       */
/**********************************************************/
int gfrLoadRMLParameters(int idx, System *sys, RMLParameter *res, ENDF *lib)
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
      res[j].rade[c] = lib->xptr[idx][k+4] * 10.0;     // APE: effective channel radius
      res[j].radt[c] = lib->xptr[idx][k+5] * 10.0;     // APT: true channel radius
    }
    idx ++;

    /*** for each resonance */
    int nline = nch/6 + 1; // required number of lines for each resonance
    for(int i=0 ; i<nres ; i++){
      int k = 6*nline*i;
      res[j].er[i] = lib->xptr[idx][k];
      for(int c=0 ; c<nch ; c++) res[j].gam[c][i] = lib->xptr[idx][k+c+1];
    }
    idx ++;
  }

  return(restot);
}

