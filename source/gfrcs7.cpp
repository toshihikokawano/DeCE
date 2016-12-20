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
#include "constant.h"
#include "matrix.h"

static int    gfrLoadParticlePairs (int, System *, ParPair *, ENDF *);
static int    gfrLoadRMLParameters (int, System *, RMLParameter *, ENDF *);
static Pcross gfrRMatrixLimited    (const int, double, System *, ParPair *, RMLParameter *);
static int    arrange_matrixRML    (int, int, int *, int *, double **, Wfunc *, double *, ParPair *, RMLParameter *);

extern Smatrix Smat;

/**********************************************************/
/*      Pointwise Cross Section in Resonance Range        */
/**********************************************************/
Pcross gfrCrossSection7(const int ner, const double elab, System *sys, ENDF *lib)
{
  Pcross sig;
  RMLParameter *res;
  ParPair      *ppr;

  gfrSetEnergy(elab, sys);

  /*** two particle pair data */
  ppr = new ParPair [MAX_PAIRS];
  gfrLoadParticlePairs(sys->idx[ner],sys,ppr,lib);

  /*** resonance parameters for each spin group */
  res = new RMLParameter [sys->nj];
  gfrLoadRMLParameters(sys->idx[ner],sys,res,lib);

  //double gk = PI/(sys->wave_number*sys->wave_number) * 0.01 / ((sys->target_spin2+1)*2);


  /*** look for max L and Nch */
  int lmax = 0, cmax = 0;
  for(int j=0 ; j<sys->nj ; j++){

    if(cmax < res[j].nch) cmax = res[j].nch;

    /*** find neutron channel for this spin group */
    for(int c=0 ; c<res[j].nch ; c++){
      if(ppr[res[j].pidx[c]].mt == 2){
        if(lmax < res[j].l[c]) lmax = res[j].l[c];
      }
    }
  }
  sys->nl = lmax+1;

  /*** calculate cross section */
  sig = gfrRMatrixLimited(cmax,elab,sys,ppr,res);


  delete [] ppr;
  delete [] res;

  return(sig);
}



Pcross gfrRMatrixLimited(const int cmax, double e, System *sys, ParPair *ppr, RMLParameter *res)
{
  const int max_elastic_channel = 2;

  Pcross sig,z;
  complex<double> *rmat, *smat, *wmat;
  Wfunc  wf[max_elastic_channel];
  int    nch, mch, *mtid, elidx[max_elastic_channel];
  double *gmat, *pen[max_elastic_channel];

  const int msize = cmax * (cmax+1) / 2;
  rmat = new complex<double> [msize];
  smat = new complex<double> [msize];
  wmat = new complex<double> [msize];
  gmat = new double [msize];
  mtid = new int [msize];
  for(int i=0 ; i<max_elastic_channel ; i++){
    pen[i] = new double [MAX_RESONANCE];
  }

  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01 / ((sys->target_spin2+1)*2);

  /*** for all spin groups */
  Smat.resetIndex();
  for(int j=0 ; j<sys->nj ; j++){

    double x3 = (res[j].j2 + 1.0)*x1;

    /*** look for elastic channels */
    int nel = 0;
    for(int c=0 ; c<res[j].nch ; c++) if( ppr[res[j].pidx[c]].mt == 2 ) elidx[nel++] = c;

    /*** incident channel radius and hard-sphere phase */
    for(int i=0 ; i<nel ; i++){
      int    l   = res[j].l[elidx[i]];
      double apt = res[j].radt[elidx[i]];
      double phi = gfrPenetrability(l,sys->wave_number*apt,&wf[i]);
      wf[i].phase  = complex<double>(cos(  phi), -sin(  phi));

      /*** penetrability at the resonance energy, calculate L(Er) = S(Er) + iP(Er) */
      if(sys->gammaunit_flag == 0){
        for(int k=0 ; k<res[j].nres ; k++){
          complex<double> q = gfrLfunction(l,res[j].er[k],sys->reduced_mass,apt);
          pen[i][k] = q.imag();
        }
      }
    }

    /*** R-matrix elements */
    mch = 0; nch = 0;
    for(int i=0 ; i<msize ; i++) smat[i] = rmat[i] = wmat[i] = complex<double>(0.0,0.0);
    for(int k=0 ; k<res[j].nres ; k++){
      nch = arrange_matrixRML(k,nel,elidx,mtid,pen,wf,gmat,ppr,&res[j]);
      mch = (nch-1)*nch / 2; // capture eliminated

      /*** capture width found at the last element */
      complex<double> w(res[j].er[k]-e, -gmat[nch]/2.0);
      w = 1.0 / w;
      for(int i=0 ; i<mch ; i++) rmat[i] += gmat[i]*w;
    }

    for(int i=0 ; i<mch ; i++){
      wmat[i] = complex<double>(imag(rmat[i])/2.0,-real(rmat[i])/2.0);
    }
    for(int i=0 ; i<nch-1 ; i++) wmat[(i+2)*(i+1)/2-1] += 1.0;

    MatrixInverse(nch-1,wmat);

    /*** S-matrix elements */
    complex<double> ph0(0.0,0.0), ph1(0.0,0.0),pzero(0.0,0.0);
    for(int i0=0 ; i0<mch ; i0++){
      ph0 = (i0 < nel) ? wf[i0].phase : pzero;

      for(int i1=0 ; i1<=i0 ; i1++){
        ph1 = (i1 < nel) ? wf[i1].phase : pzero;

        int ij = i0*(i0+1)/2 + i1;
        if(i0 == i1) smat[ij] = ph0 * ph1 * (2.0*wmat[ij]-1.0);
        else         smat[ij] = ph0 * ph1 * 2.0*wmat[ij];
      }
    }

    /*** cross section */
    z.clear();
    for(int i0=0 ; i0<mch ; i0++){

      if(i0 < nel){
        int ij = (i0+2)*(i0+1)/2 - 1;
        z.total   += (1.0 - real(smat[ij]))*2.0;
        z.elastic += (1.0-real(smat[ij]))*(1.0-real(smat[ij])) + imag(smat[ij])*imag(smat[ij]);
      }

      else{
        for(int i1=0 ; i1<nch ; i1++){
          int ij = i0*(i0+1)/2 + i1;
          double cx = real(smat[ij])*real(smat[ij]) + imag(smat[ij])*imag(smat[ij]);

          if(mtid[i0] == 18) z.fission += cx;
          else               z.other   += cx;
        }
      }
    }

    z.capture  = z.total - z.elastic - z.fission - z.other;

    sig.total    += x3*z.total;
    sig.elastic  += x3*z.elastic;
    sig.capture  += x3*z.capture;
    sig.fission  += x3*z.fission;

    /*** copy S-matrix elements for elastic */
    for(int i=0 ; i<nel ; i++){

      int ss = (int)abs(res[j].s2[elidx[i]]) - sys->target_spin2;

      Smat.setElement(res[j].l[elidx[i]], res[j].j2, ss, smat[(i+2)*(i+1)/2-1]);
      Smat.inclIndex();
    }
  }
  
  delete [] rmat;
  delete [] smat;
  delete [] wmat;
  delete [] gmat;
  delete [] mtid;
  for(int i=0 ; i<max_elastic_channel ; i++){
    delete [] pen[i];
  }

  return(sig);
}



/**********************************************************/
/*      Matrix Element Re-arrangement for RM              */
/**********************************************************/
int arrange_matrixRML(int k, int nel, int *elidx, int *mt, double **pen,
                      Wfunc *wf, double *x, ParPair *ppr, RMLParameter *res)
{
  /*** diagonal elements */
  /*** neutron elastic should be first */
  int idx = 0;
  for(int i=0 ; i<nel ; i++){
    x[(idx+2)*(idx+1)/2 - 1]  = res->gam[k][elidx[i]] * imag(wf[i].d)/ pen[i][k];
    mt[idx] = 2;
    idx ++;
  }

  /*** other channels */
  for(int c=0 ; c<res->nch ; c++){
    if((ppr[res->pidx[c]].mt == 2) || (ppr[res->pidx[c]].mt == 102)) continue;
    x[(idx+2)*(idx+1)/2 - 1] = res->gam[k][c];
    mt[idx] = ppr[res->pidx[c]].mt;
    idx ++;
  }

  /*** capture channel last */
  for(int c=0 ; c<res->nch ; c++){
    if(ppr[res->pidx[c]].mt == 102){
      x[(idx+2)*(idx+1)/2 - 1] = res->gam[k][c];
      mt[idx] = ppr[res->pidx[c]].mt;
      idx ++;
    }
  }
  int nch = idx;

  /*** off-diagonal elements */
  for(int i=1 ; i<nch ; i++){
    int di = (i+2)*(i+1)/2 - 1; // index for diagonal element
    int si = (x[di] < 0.0) ? -1 : 1;

    for(int j=0 ; j<=i-1 ; j++){
      int dj = (j+2)*(j+1)/2 - 1;
      int sj = (x[dj] < 0.0) ? -1 : 1;

      int k = (i+1)*i/2 + j;
      x[k] = sqrt( abs(x[di]) * abs(x[dj]) ) * si * sj;
    }
  }

  /*** diagonal elements */
  double gt = 0.0;
  for(int i=0 ; i<nch ; i++){
    int di = (i+2)*(i+1)/2 - 1;
    x[di] = abs(x[di]);
    gt += x[di];
  }

  return(nch);
}


/**********************************************************/
/*     Read Particle Pair Section                         */
/**********************************************************/
int gfrLoadParticlePairs(int idx, System *sys, ParPair *ppr, ENDF *lib)
{
  /*** first line */
  sys->gammaunit_flag = lib->rdata[idx].l1; // 0, maybe
  sys->format         = lib->rdata[idx].l2; // 3
  sys->nj             = lib->rdata[idx].n1; // Nj
  sys->relativ_flag   = lib->rdata[idx].n2; // 0
  idx++;

  sys->npair          = lib->rdata[idx].l1;
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

    p[0] = (int)lib->xptr[idx][k++];
    p[1] = (int)lib->xptr[idx][k++];

    for(int j=0 ; j<2 ; j++){
      if(ppr[i].spin2[j] == 0) ppr[i].parity[j] = p[j];
      else{
        ppr[i].parity[j] = (ppr[i].spin2[j] < 0) ? -1 : 1;
        ppr[i].spin2[j] = abs(ppr[i].spin2[j]);
      }
    }

    /*** find target spin */
    if(ppr[i].mt == 2){
      for(int j=0 ; j<2 ; j++){
        if(ppr[i].znum[j] == 0){
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
  idx += 2;

  int restot = 0;
  for(int j=0 ; j<sys->nj ; j++){

    /*** channel data */
    res[j].j2 = (int)(2.0*lib->rdata[idx].c1);
    int pj    = (int)lib->rdata[idx].c2;
//  int kbk   = lib->rdata[idx].l1;
//  int kps   = lib->rdata[idx].l2;
    int nch   = lib->rdata[idx].n2;
    int nres  = lib->rdata[idx+1].l2;

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
      res[j].pidx[c] = lib->xptr[idx][k] - 1;
      res[j].l[c]    = lib->xptr[idx][k+1];
      res[j].s2[c]   = (int)(2.0*lib->xptr[idx][k+2]);
      res[j].rade[c] = lib->xptr[idx][k+4] * 10.0;
      res[j].radt[c] = lib->xptr[idx][k+5] * 10.0;
    }
    idx ++;

    /*** for each resonance */
    int nline = nch/6 + 1; // required number of lines for each resonance
    for(int i=0 ; i<nres ; i++){
      int k = 6*nline*i;
      res[j].er[i] = lib->xptr[idx][k];
      for(int c=0 ; c<nch ; c++) res[j].gam[i][c] = lib->xptr[idx][k+c+1];
    }
    idx ++;
  }

  return(restot);
}

