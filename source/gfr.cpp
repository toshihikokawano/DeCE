/******************************************************************************/
/**     GFR, Reconstruct Cross Sections from Resoance Parameters             **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "coupling.h"
#include "terminate.h"

#define printSmat

double gcLorentzianWidth = 0.0, *fact;

Smatrix Smat;

static Pcross  gfrPtCrossFILE    (ifstream *, ENDFDict *, const double);
static Pcross  gfrBackGroundFILE (ifstream *, ENDF *, const double);
static Pcross  gfrBackGround     (ENDFDict *, ENDF **, const double);
static void gfrPrintCrossSection (const double, Pcross);

/**********************************************************/
/*      GFR Interface, Scan File for the Thermal Value    */
/**********************************************************/
void gfrScanThermal(ifstream *fp, ENDFDict *dict, double elab)
{
  const double eth = 0.0253;

  if(elab <= 0.0) elab = eth;

  Pcross crs = gfrPtCrossFILE(fp,dict,elab);
  gfrPrintCrossSection(elab,crs);
}


/**********************************************************/
/*      GFR Interface, Cross Section at Given Energy      */
/**********************************************************/
double gfrGetOnePoint(ifstream *fp, ENDFDict *dict, const int mt, const double elab)
{
  Pcross crs = gfrPtCrossFILE(fp,dict,elab);
  double c = 0.0;
  if(     mt ==   1) c = crs.total;
  else if(mt ==   2) c = crs.elastic;
  else if(mt ==  18) c = crs.fission;
  else if(mt == 102) c = crs.capture;

  return(c);
}


/**********************************************************/
/*      One Point Cross Section from File                 */
/**********************************************************/
Pcross gfrPtCrossFILE(ifstream *fp, ENDFDict *dict, const double elab)
{
  ENDF   librs,libbg;
  System sys;
  Pcross crs, cbg;

  /*** Read beginning of the file */
  ENDFSeekHead(fp,&librs,1,451);

  /*** fissile and resonance flags */
  int lrp  = dict->getLRP();   //  resonance parameter flag
  if(lrp < 0) return(crs);

  ENDFReadMF2(fp,&librs);

  gfrReadHEADData(&sys,&librs);

  /*** resonance contribution */
  crs = gfrCrossSection(0,elab,&sys,&librs);

  /*** background contribution */
  cbg = gfrBackGroundFILE(fp,&libbg,elab);
  crs = crs + cbg;

  return(crs);
}


/**********************************************************/
/*      GFR Interface to DeCE, Tabulate Cross Sections    */
/**********************************************************/
void gfrPtCross(ENDFDict *dict, ENDF *lib[], double emin, double emax, double de)
{
  System sys;
  Pcross crs, cbg;
  double *elab;
  int np = 0;

  if(dict->emaxRe == 0.0){
    message << "no resonance region given";
    WarningMessage();
    return;
  }

  elab = new double [MAX_DBLDATA/2];

  int kres = dict->getID(2,151);
  gfrReadHEADData(&sys,lib[kres]);

  gfrPrintCrossSection(-1.0,crs);

  /*** determine energy grid automatically */
  if((emin == 0.0) && (emax == 0.0)){

    np = gfrAutoEnergyRRR(&sys,lib[kres],elab,dict->emaxRR);

    for(int i=0 ; i<np ; i++){

      /*** resonance cross sections */
      crs = gfrCrossSection(1,elab[i],&sys,lib[kres]);
      /*** background cross section in MF3 */
      cbg = gfrBackGround(dict,lib,elab[i]);

      /*** add background */
      crs = crs + cbg;

      /*** print cross section */
      gfrPrintCrossSection(elab[i],crs);
    }

    np = gfrAutoEnergyURR(elab,dict->emaxRR,dict->emaxUR);

    for(int i=0 ; i<np ; i++){

      /*** unresolved resonance cross sections */
      crs = gfrCrossSection(2,elab[i],&sys,lib[kres]);

      /*** print cross section */
      gfrPrintCrossSection(elab[i],crs);
    }
  }
  /*** equidistant energy grid case */
  else{
    np = gfrFixedEnergyRRR(emin,emax,de,elab,dict->emaxRR,dict->emaxUR);

    for(int i=0 ; i<np ; i++){

      /*** resonance cross sections */
      crs = gfrCrossSection(0,elab[i],&sys,lib[kres]);
      if(elab[i] < dict->emaxRR){
        /*** background cross section in MF3 */
        cbg = gfrBackGround(dict,lib,elab[i]);
        crs = crs + cbg;
      }

      /*** print cross section */
      gfrPrintCrossSection(elab[i],crs);
    }
  }

  delete [] elab;
}


/**********************************************************/
/*      GFR Generate Legendre Coefficient from Resonances */
/**********************************************************/
void gfrAngDist(ENDFDict *dict, ENDF *lib[], double emin, double emax, double de)
{
  const int ndiv = 100;
  System sys;
  Pcross c;
  double *elab, *pleg;

  elab = new double [MAX_DBLDATA/2];
  pleg = new double [LMAX*2];

  Smat.memalloc(2*(LMAX+1)*(LMAX+1)-1);

  /*** factorial */
  factorial_allocate();

  int kres = dict->getID(2,151);
  gfrReadHEADData(&sys,lib[kres]);

  cout <<"# Energy[eV]   ";
  for(int l=1 ; l<LMAX ; l++)  cout <<"P"<< setw(1) << l <<"/P0        ";
  cout << endl;

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(5);

  if( (emin == 0.0) && (emax == 0.0) && (de == 0.0) ){
    emax = dict->emaxRR;
    de   = emax / ndiv;
    emin = de;
  }

  /*** always equi-distant energy grid */
  int np = gfrFixedEnergyRRR(emin,emax,de,elab,dict->emaxRR,dict->emaxUR);

  for(int i=0 ; i<np ; i++){

    /*** cross section calculation provides S-matrix elements */
    c = gfrCrossSection(0,elab[i],&sys,lib[kres]);

    /*** Legendre coefficients calculated from the S-matrix elements */
    gfrLegendreCoefficient(&sys,pleg);
    if(gcLorentzianWidth > 0.0)  pleg[0] = gfrCompoundReaction(&sys);

    cout << setw(13) << elab[i];
    for(int l=1 ; l<LMAX ; l++){
      if(pleg[0] == 0.0) cout << setw(13) << 0.0;
      else               cout << setw(13) << pleg[l]/pleg[0];
    }
    cout << endl;
  }

  Smat.memfree();
  delete [] elab;
  delete [] pleg;
  factorial_delete();
}


/**********************************************************/
/*      GFR Smoothed Legendre Coefficient from Resonances */
/**********************************************************/
void gfrAngDistSmooth(ENDFDict *dict, ENDF *lib[], double width)
{
  const int npoint = 100;
  const int ndiv   = 100;
  const int wf     =   1;
  System sys;
  Pcross c;
  double **pave, *wsum, *pleg;

  pave = new double * [npoint];
  wsum = new double [npoint];
  for(int i=0 ; i<npoint ; i++){
    pave[i] = new double [LMAX*2];
  }
  pleg = new double [LMAX*2];

  Smat.memalloc(2*(LMAX+1)*(LMAX+1)-1);

  /*** factorial */
  factorial_allocate();

  int kres = dict->getID(2,151);
  gfrReadHEADData(&sys,lib[kres]);

  cout <<"# Energy[eV]   ";
  for(int l=1 ; l<LMAX ; l++)  cout <<"P"<< setw(1) << l <<"/P0        ";
  cout << endl;

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(5);

  /*** divide resolved resonance region by Npoint x Ndiv */
  int    np   = npoint * ndiv;
  double emin = 1e-05;
  double emax = dict->emaxRR;
  double de0  = emax / np;
  double de1  = emax / npoint;

  if(width == 0.0) width = de1;

  cout << setw(13) << emin;
  for(int l=1 ; l<LMAX ; l++) cout << setw(13) << 0.0;
  cout << endl;

  for(int n=0 ; n<npoint ; n++){
    wsum[n] = 0.0;
    for(int l=1 ; l<LMAX*2 ; l++) pave[n][l] = 0.0;
  }
  

  for(int i=1 ; i<=np ; i++){
    double elab = i * de0;
    c = gfrCrossSection(0,elab,&sys,lib[kres]);
    gfrLegendreCoefficient(&sys,pleg);
  
    for(int n=0 ; n<npoint ; n++){
      double e0 = de1 * (n+1);

      /*** weighted by Gaussian (1), Lorentzian (2), or constant (3) */
      double weight = 0.0;
      double x = e0 - elab;
      switch(wf){
      case   1: weight = exp(-x*x/(width*width)); break;
      case   2: weight = 1.0/(x*x + width*width); break;
      case   3: weight = 1.0; break;
      default : weight = 1.0; break;
      }
      wsum[n] += weight;

      for(int l=0 ; l<LMAX*2 ; l++) pave[n][l] += weight * pleg[l];
    }
  }
  for(int n=0 ; n<npoint ; n++){
    for(int l=0 ; l<LMAX*2 ; l++) pave[n][l] /= wsum[n];
  }

  for(int n=0 ; n<npoint ; n++){
    double e0 = de1 * (n+1);
    cout << setw(13) << e0;
    for(int l=1 ; l<LMAX ; l++){
      if(pave[n][0] == 0.0) cout << setw(13) << 0.0;
      else                  cout << setw(13) << pave[n][l]/pave[n][0];
    }
    cout << endl;
  }

  Smat.memfree();
  for(int i=0 ; i<npoint ; i++){
    delete [] pave[i];
  }
  delete [] pave;
  delete [] wsum;
  delete [] pleg;
  factorial_delete();
}


/**********************************************************/
/*      GFR Print Scattering Matrix of Resonances         */
/**********************************************************/
void gfrSmatrixElement(ENDFDict *dict, ENDF *lib[])
{
  const int ndiv = 10000;
  System sys;
  Wfunc  wfn;
  double *elab;

  double emax = dict->emaxRR;
  double de   = emax / ndiv;
  double emin = de;

  elab = new double [MAX_DBLDATA/2];
  int np = gfrFixedEnergyRRR(emin,emax,de,elab,dict->emaxRR,dict->emaxUR);

  Smat.memalloc(2*(LMAX+1)*(LMAX+1)-1);

  int kres = dict->getID(2,151);
  gfrReadHEADData(&sys,lib[kres]);

  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(4);
  int l = 0, j2 = 0, s2 = 0;

  gfrCrossSection(0,elab[0],&sys,lib[kres]);

  cout <<"# Energy[eV] ";
  for(int j=0 ; j<Smat.getIndex() ; j++){
    Smat.getElement(j,&l,&j2,&s2);
    cout << setw(3) << l;
    cout << setw(3) << j2 << "/2                ";
  }
  cout <<" Hard-Sphere ";
  cout << endl;

  for(int i=0 ; i<np ; i++){

    gfrCrossSection(0,elab[i],&sys,lib[kres]);

    cout << setw(13) << elab[i];
    for(int j=0 ; j<Smat.getIndex() ; j++){
      cout << setw(12) << Smat.getElement(j).real();
      cout << setw(12) << Smat.getElement(j).imag();
    }

    for(l=0 ; l<sys.nl ; l++){
      double phase = gfrPenetrability(l,sys.alpha,&wfn);
      cout << setw(12) <<  cos(2*phase);
      cout << setw(12) << -sin(2*phase);
    }

    cout << endl;
  }

  Smat.memfree();
  delete [] elab;
}


/**********************************************************/
/*      Read HEAD section: GFR Interface to DeCE          */
/**********************************************************/
void gfrReadHEADData(System *sys, ENDF *lib)
{
  /*** HEAD section in MF2 MT151 */
  Record head = lib->getENDFhead();
  double za   = head.c1;
  double awr  = head.c2;

  sys->target_Z = (int)za/1000;
  sys->target_A = (int)za - sys->target_Z*1000;
  sys->reduced_mass = awr/(awr+1.0);

  /*** CONT0 */
  Record cont = lib->rdata[0];
  sys->avefission_flag  = cont.l2;
  sys->nrange = cont.n1;

  /*** for each resonance range */
  int idx = 1;
  for(int i=0 ; i<sys->nrange ; i++){

    /*** CONT in each NER range */
    cont = lib->rdata[idx];
    sys->emin[i]     = cont.c1;
    sys->emax[i]     = cont.c2;
    sys->lru[i]      = cont.l1; // 0: scat.rad only, 1: resolved, 2: unresolved
    sys->lrf[i]      = cont.l2; // 1: SLBW, 2: MLBW, 3: RM, 7: RML
    sys->nro[i]      = cont.n1;
    sys->naps[i]     = cont.n2;

    if(sys->nro[i] == 1) TerminateCode("NRO=1, not inplemented yet");
    
    idx ++;
    sys->idx[i]      = idx;

    /*** skip all data section */
    /*** Resolved Resonance Region */
    if(sys->lru[i] == 1){

      if(sys->lrf[i] <= 3){
        cont = lib->rdata[idx++];
        int nls = cont.n1;
        if(sys->nro[i] !=0) idx++;
        idx += nls;
      }

      else if(sys->lrf[i] == 7){
        cont = lib->rdata[idx];
        int njs = cont.n1;
        idx ++;
        idx += njs;
      }
    }

    /*** Unresolved Resonance Region */
    else if(sys->lru[i] == 2){
      if(sys->lrf[i] == 2){
        cont = lib->rdata[idx++];
        int nls = cont.n1;
        for(int inls=0 ; inls<nls ; inls++){
          cont = lib->rdata[idx++];
          int njs = cont.n1;
          idx += njs;
        }
      }
      else TerminateCode("LRU=2, LRF=1, not implemented yet");
    }
  }
}


/**********************************************************/
/*      Background Cross Section in MF3                   */
/**********************************************************/
Pcross  gfrBackGroundFILE(ifstream *fp, ENDF *lib, const double elab)
{
  Pcross cbg;

  ENDFReadMF3(fp,lib,  1); cbg.total   = ENDFInterpolation(lib,elab,true,0);
  ENDFReadMF3(fp,lib,  2); cbg.elastic = ENDFInterpolation(lib,elab,true,0);
  ENDFReadMF3(fp,lib, 18); cbg.fission = ENDFInterpolation(lib,elab,true,0);
  ENDFReadMF3(fp,lib,102); cbg.capture = ENDFInterpolation(lib,elab,true,0);

  return(cbg);
}


Pcross  gfrBackGround(ENDFDict *dict, ENDF **lib, const double elab)
{
  Pcross cbg;

  if(dict->getID(3,  1) >= 0) cbg.total   = ENDFInterpolation(lib[dict->getID(3,  1)],elab,true,0);
  if(dict->getID(3,  2) >= 0) cbg.elastic = ENDFInterpolation(lib[dict->getID(3,  2)],elab,true,0);
  if(dict->getID(3, 18) >= 0) cbg.fission = ENDFInterpolation(lib[dict->getID(3, 18)],elab,true,0);
  if(dict->getID(3,102) >= 0) cbg.capture = ENDFInterpolation(lib[dict->getID(3,102)],elab,true,0);

  return(cbg);
}


/**********************************************************/
/*      Print Crooss Section                              */
/**********************************************************/
static void gfrPrintCrossSection (const double elab, Pcross crs)
{
  /*** prinr heading */
  if(elab < 0.0){
    cout.setf(ios::scientific, ios::floatfield);
    cout <<"# Energy[eV]    Total[b]      Elastic[b]    Capture[b]    Fission[b]    Proton[b]     Alpha[b]      Other[b]";
    cout << endl;
  }
  else{
    cout << setw(14) << setprecision(7) << elab;
    cout << setw(14) << setprecision(6) << crs.total;
    cout << setw(14) << crs.elastic;
    cout << setw(14) << crs.capture;
    cout << setw(14) << crs.fission;
    cout << setw(14) << crs.proton;
    cout << setw(14) << crs.alpha;
    cout << setw(14) << crs.other;
    cout << endl;
  }
}
