/******************************************************************************/
/**     GFR, Calculate Unresolved Resonance Cross Sections                   **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <iomanip>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "constant.h"
#include "terminate.h"

static int gfrURetrieveParameter (const int, int, ENDF *, URResonance *);
static int gfrUFindRange         (double, URResonance *, bool *);
static Pcross gfrCrossSectionURR (int, int, int, System *, URResonance *);
static Pcross gfrBreitWignerURR  (int, double, double, Wfunc *, URResonance *);
static double gfrMoldauer        (int, double *, double *, double *);
static Pcross gfrUInterpolation  (double, Pcross, Pcross);

static const int MAX_URESONANCE =  20;
static const int MAX_EPOINTS    = 100;

static int MAX_GAUSS20 = 10;
static double
gauss20_x[]={
  0.9931285991850949,  0.9639719272779138,  0.9122344282513259,  0.8391169718222188,
  0.7463319064601508,  0.6360536807265150,  0.5108670019508271,  0.3737060887154196,
  0.2277858511416451,  0.0765265211334973},
gauss20_a[]={
  0.0176140071391521,  0.0406014298003869,  0.0626720483341091,  0.0832767415767047,
  0.1019301198172404,  0.1181945319615184,  0.1316886384491766,  0.1420961093183821,
  0.1491729864726037,  0.1527533871307258};



/**********************************************************/
/*      Pointwise Cross Section in Unresolved Resonance   */
/**********************************************************/
Pcross gfrCrossSectionURR(const int ner, const double elab, System *sys, ENDF *lib)
{
  Pcross z0, z1, z2;
  URResonance *res;

  res = new URResonance [MAX_URESONANCE];
  for(int i=0 ; i<MAX_URESONANCE ; i++) res[i].memalloc(MAX_EPOINTS);

  int lrf = sys->lrf[ner];

  /*** case A: energy independent parameter, no fission */
  if(lrf == 1){
    if(sys->avefission_flag == 0){
    }
  /*** case B: energy independent, but for fission */
    else{
    }

    WarningMessage("LRF=0,1 not yet supported");
  }
  /*** case C: all parameters are energy dependent */
  else{
    bool itp = false;
    int pm = gfrURetrieveParameter(sys->nl,sys->idx[ner]+1,lib,res);
    int ke = gfrUFindRange(elab,res,&itp);

    /*** need interpolation */
    if(itp){
      z1 = gfrCrossSectionURR(pm,ke  ,ner,sys,res);
      z2 = gfrCrossSectionURR(pm,ke+1,ner,sys,res);
      z0 = gfrUInterpolation(elab,z1,z2);
    }
    /*** no interpolation case */
    else{
      z0 = gfrCrossSectionURR(pm,ke,ner,sys,res);
    }
  }

  delete [] res;

  return(z0);
}


/**********************************************************/
/*      Unresolved Resonance Parameters at Given Energy   */
/**********************************************************/
int gfrURetrieveParameter(const int nl, int idx, ENDF *lib, URResonance *res)
{
  int p = 0; // index for a given (L,J) pair

  /*** for all L partial waves */
  for(int l=0 ; l<nl ; l++){
    /*** CONT */
    int njs = lib->rdata[idx].n1;
    idx ++;

    for(int j=0 ; j<njs ; j++){
      /*** LIST */
      res[p].l   = l;
      res[p].j2  = (int)(2.0*lib->rdata[idx].c1);
      res[p].itp = lib->rdata[idx].l1;
      res[p].ne  = lib->rdata[idx].n2;;

      res[p].dfx = lib->xptr[idx][2];
      res[p].dfn = lib->xptr[idx][3];
      res[p].dfg = lib->xptr[idx][4];
      res[p].dff = lib->xptr[idx][5];

      for(int k=0 ; k<lib->rdata[idx].n2 ; k++){
        int q = (k+1)*6;
        res[p].bw[k].er = lib->xptr[idx][q];     // energy point
        res[p].bw[k].d  = lib->xptr[idx][q+1];   // average spacing
        res[p].bw[k].gx = lib->xptr[idx][q+2];   // competitive width
        res[p].bw[k].gn = lib->xptr[idx][q+3];   // reduced neutron width
        res[p].bw[k].gg = lib->xptr[idx][q+4];   // gamma width
        res[p].bw[k].gf = lib->xptr[idx][q+5];   // fission width
      }
      p ++;
      idx ++;

      if(p >= MAX_URESONANCE){
        cerr << "too many L,J groups in URR" << endl;
        return(0);
      }
    }
  }
  return(p);
}


/**********************************************************/
/*      Find Energy Grids for Given Elab                  */
/**********************************************************/
int gfrUFindRange(double elab, URResonance *res, bool *itp)
{
  /*** look at the first L,J section (k=0) to determine the energy grid */
  const int k = 0;

  /*** check if the highest point */
  if(res[k].bw[res[k].ne-1].er == elab){
    *itp = false;
    return(res[k].ne-1);
  }

  /*** search for adjacent grids that contain Elab */
  int kp = 0;
  *itp = true;
  for(int i=0 ; i<res[k].ne-1 ; i++){

    if(res[k].bw[i].er == elab){
      kp = i;
      *itp = false;
      break;
    }
    if((res[k].bw[i].er < elab) && (elab < res[k].bw[i+1].er) ){
      kp =  i;
      break;
    }
  }

  return(kp);
}


/**********************************************************/
/*      Calculate Cross Section                           */
/**********************************************************/
Pcross gfrCrossSectionURR(int km, int ke, int ner, System *sys, URResonance *res)
{
  Wfunc  wfn;
  Pcross sig, z;

  sig.energy = res[0].bw[ke].er; // we hope energy grids are always the same.
  gfrSetEnergy(sig.energy,sys);

  double ap = (sys->naps[ner] == 0) ? gfrENDFChannelRadius(sys->target_A) : sys->radius;
  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01;
  double x2 = sys->wave_number * ap;
  double p2 = 2.0*PI;

  for(int l=0 ; l<sys->nl ; l++){

    double phase = gfrPenetrability(l,sys->alpha,&wfn);
    wfn.phase  = complex<double>(cos(  phase), -sin(  phase));
    wfn.phase2 = complex<double>(cos(2*phase), -sin(2*phase));

    if(sys->naps[ner] == 0) gfrPenetrability(l,x2,&wfn);

    int smin = abs(sys->target_spin2-1);
    int smax =     sys->target_spin2+1;

    for(int ss=smin ; ss<=smax ; ss+=2){

      int jmin = abs(2*l-ss);
      int jmax =     2*l+ss ;

      for(int jj=jmin ; jj<=jmax ; jj+=2){
        double gj = (jj+1.0) / ((sys->target_spin2+1)*2.0);

        int kp = 0;
        for(int k=0 ; k<km ; k++){
          if((l == res[k].l) && (jj == res[k].j2)){ kp = k; break; }
        }

        z = gfrBreitWignerURR(ke,gj,x2,&wfn,&res[kp]);

        sig.elastic  += p2 * x1 * z.elastic;
        sig.capture  += p2 * x1 * z.capture;
        sig.fission  += p2 * x1 * z.fission;
        sig.other    += p2 * x1 * z.other;
      }
    }

    sig.elastic += 4.0*x1*(2.0*l+1.0) * imag(wfn.phase)*imag(wfn.phase);
  }

  sig.setTotal();

  return(sig);
}


/**********************************************************/
/*      Breit-Wigner form in Unresolved Range             */
/**********************************************************/
Pcross gfrBreitWignerURR(int ke, double gj, double x, Wfunc *wfn, URResonance *res)
{
  Pcross  z;
  double tr[4],df[4],wf[4];

  double gf = imag(wfn->d) / x * sqrt(res->bw[ke].er);
  double gn = res->bw[ke].gn * gf;
  double gt = gn + res->bw[ke].gg + res->bw[ke].gf + res->bw[ke].gx;
  double s2 = imag(wfn->phase) * imag(wfn->phase);

  double r  = gj / res->bw[ke].d;

  tr[0] = gn;
  tr[1] = res->bw[ke].gf;
  tr[2] = res->bw[ke].gg;
  tr[3] = res->bw[ke].gx;

  df[0] = res->dfn;  if(df[0] == 0.0) df[0] =  1.0; // ad hot
  df[1] = res->dff;  if(df[1] == 0.0) df[2] =  4.0; // ad hot
  df[2] = res->dfg;  if(df[2] == 0.0) df[2] = 20.0; // ad hot
  df[3] = res->dfx;
 
 
  gfrMoldauer(4,tr,df,wf);
  
  z.elastic = r *(wf[0] * gn * gn             / gt - 2.0*gn*s2);
  z.fission = r * wf[1] * gn * res->bw[ke].gf / gt;
  z.capture = r * wf[2] * gn * res->bw[ke].gg / gt;
  z.other   = r * wf[3] * gn * res->bw[ke].gx / gt;

  return(z);
}


/**********************************************************/
/*      Moldauer's Width Fluctuation Factor               */
/**********************************************************/
double gfrMoldauer(int n, double *tc, double *nu, double *wfc)
{
  double ts = 0.0;
  double w[2*MAX_GAUSS20];

  for(int i=0 ; i<2*MAX_GAUSS20 ; i++) w[i] = 1.0;

  /*** sum all transmission coefficients */
  for(int i=0 ; i<n ; i++) ts += tc[i];

  /*** calculate integrant at each Gauss-Legendre point */
  for(int i=0 ; i<n ; i++){
    double nu2 = nu[i]/2.0;
    double x   = 1.0-tc[i]/(ts*nu2);

    for(int ig=0 ; ig<MAX_GAUSS20 ; ig++){
      double x1 = 0.5*(1+gauss20_x[ig]);
      double x2 = 0.5*(1-gauss20_x[ig]);
      /*** W at Gauss-Legendre inegration */
      w[MAX_GAUSS20+ig  ] *= pow( ((1-x*x1)/x2),-nu2 );
      w[MAX_GAUSS20-ig-1] *= pow( ((1-x*x2)/x1),-nu2 );
    }
  }

  /*** calculate width fluctuation correction for each channel */
  double df0 = nu[0];

  for(int i=0 ; i<n ; i++){
    if(tc[i] == 0.0) wfc[i] = 1.0;
    else{
      double df1 = nu[i];
      double x3  = 2.0/(df0*ts);
      double x4  = 2.0/(df1*ts);

      double sum = 0.0;
      for(int ig=0 ; ig<MAX_GAUSS20 ; ig++){
        double x1 = 0.5*(1+gauss20_x[ig]);
        double x2 = 0.5*(1-gauss20_x[ig]);
        double g  = 1.0 - tc[0]*x3;
        double f  = 1.0 - tc[i]*x4;
        sum += gauss20_a[ig]* ( w[MAX_GAUSS20+ig  ]/(1-f*x1)/(1-g*x1)
                               +w[MAX_GAUSS20-ig-1]/(1-f*x2)/(1-g*x2));
      }
      wfc[i] = sum/2.0;
    }
  }

  wfc[0] *= 1.0 + 2.0/df0;

  return(ts);
}


/**********************************************************/
/*      Interpolate Two Cross Sections                    */
/**********************************************************/
Pcross gfrUInterpolation(double elab, Pcross z1, Pcross z2)
{
  Pcross z0;

  /*** linear interpolation used */
  double dx = (elab - z1.energy) / (z2.energy - z1.energy);

  z0.total   = (z2.total   - z1.total  ) * dx + z1.total;
  z0.elastic = (z2.elastic - z1.elastic) * dx + z1.elastic;
  z0.capture = (z2.capture - z1.capture) * dx + z1.capture;
  z0.fission = (z2.fission - z1.fission) * dx + z1.fission;
  z0.other   = (z2.other   - z1.other  ) * dx + z1.other;

  z0.energy  = elab;

  return(z0);
}
