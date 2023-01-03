/******************************************************************************/
/**     DeCE TABLE for MF6                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "decemisc.h"
#include "global.h"
#include "constant.h"
#include "kalbach.h"

static int DeceTableMF6Law1(ENDF *, int);
static void DeceTableMF6Law1Lang1  (ENDF *, int);
static void DeceTableMF6Law1Lang2  (ENDF *, int);
static void DeceTableMF6Law1Lang11 (ENDF *, int, const int);
static int DeceTableMF6Law2(ENDF *, int);
static int DeceTableMF6Law5(ENDF *, int);
static int DeceTableMF6Law6(ENDF *, int);
static int DeceTableMF6Law7(ENDF *, int);


/**********************************************************/
/*      Process MF=6                                      */
/**********************************************************/
void DeceTableMF6(ENDF *lib3, ENDF *lib6)
{
  Record head = lib6->getENDFhead();
  int    nk   = head.n1;
  int    lct  = head.l2;
  int    idx  = 0;

  cout << "# Energy and angle distribution" << endl;
  cout << "#           NK" << setw(14) << nk << "  number of subsections" << endl;
  cout << "#          LCT" << setw(14) << lct << " 1:LAB, 2:CMS, 3:CMS for particle and LAB for recoil, 4:CMS for initial and LAB for subsequent breakup" << endl;
  cout << endl;

  /*** for each sub block, make lib for (E,yield) */
  for(int ik=0 ; ik<nk ; ik++){
    Record cont = lib6->rdata[idx];
    int    zap  = (int)cont.c1;
    int    lip  = cont.l1;
    int    law  = cont.l2;
    int    np   = cont.n2;
 
    cout << "#   Subsection" << setw(14) << ik << endl;
    cout << "#          ZAP" << setw(14) << zap << endl;
    cout << "#          LIP" << setw(14) << lip << "  product modifier" << endl;
    cout << "#          LAW" << setw(14) << law << "  distribution law" << endl;
    cout << "#           NP" << setw(14) << np  << "  yield energy points" << endl;
    cout << "#   Energy[eV]  CrossSec [b]  Multiplicity  CrossSec x M"<< endl;
    for(int i=0 ; i<np ; i++){
      double x  = lib6->xptr[idx][i*2  ];
      double y6 = lib6->xptr[idx][i*2+1];
      double y3 = ENDFInterpolation(lib3,x,false,0);
      outVal(x); outVal(y3); outVal(y6); outVal(y3*y6);
      cout << endl;
    }
    cout << endl;
    cout << endl;

    /*** increment index */
    idx++;
    if(     law == 1) idx = DeceTableMF6Law1(lib6,idx);
    else if(law == 2) idx = DeceTableMF6Law2(lib6,idx);
    else if(law == 5) idx = DeceTableMF6Law5(lib6,idx);
    else if(law == 6) idx = DeceTableMF6Law6(lib6,idx);
    else if(law == 7) idx = DeceTableMF6Law7(lib6,idx);
  }
}


/**********************************************************/
/*      Continuum Energy-Angle Distribution               */
/**********************************************************/
int DeceTableMF6Law1(ENDF *lib6, int idx)
{
  int    lang = lib6->rdata[idx].l1;
  int    lep  = lib6->rdata[idx].l2;
  int    nr   = lib6->rdata[idx].n1;
  int    ne   = lib6->rdata[idx].n2; idx++;

  cout << "#           NR" << setw(14) << nr << endl;
  cout << "#           NE" << setw(14) << ne << endl;
  cout << "#         LANG" << setw(14) << lang << "  1: Legendre, 2:KM systematics, 11-15: tabulated" << endl;
  cout << "#          LEP" << setw(14) << lep << "  secondary energy interpolation" << endl;

  /*** when DDX is calculated from KM systematics, setup parameters */
  if(lang == 2 && opt.AngleStep > 0.0){
    /*** neutron incident assumed */
    double za  = (lib6->getENDFhead()).c1;
    double zap = lib6->rdata[idx-2].c1;
    ddxKalbachSetParm(1.0, za, zap);
  }
  
  for(int i0=0 ; i0<ne ; i0++){
    double e1  = lib6->rdata[idx].c2;
    int    nd  = lib6->rdata[idx].l1;
    int    na  = lib6->rdata[idx].l2;
    int    nep = lib6->rdata[idx].n2;

    cout << "#           E1"; outVal(e1); cout << "  incident energy" << endl;
    cout << "#          NEP" << setw(14) << nep << "  number of secondary energy points" << endl;
    cout << "#           ND" << setw(14) << nd << "  number of discrete lines" << endl;
    cout << "#           NA" << setw(14) << na << "  number of angular parameters" << endl;

    /*** discrete lines if given */
    if(nd > 0){
      cout << "# Energy       Production" << endl;
      for(int i1=0 ; i1<nd ; i1++){
        outVal(lib6->xptr[idx][2*i1  ]);
        outVal(lib6->xptr[idx][2*i1+1]);
        cout << endl;
      }
      cout << endl;
      cout << endl;
    }

    if(     lang == 1) DeceTableMF6Law1Lang1(lib6,idx);
    else if(lang == 2) DeceTableMF6Law1Lang2(lib6,idx);
    else if( (11 <= lang) && (lang <= 15) ) DeceTableMF6Law1Lang11(lib6,idx,lang);
    else{
      cout <<"LANG = " << lang << " not yet implemented";
    }

    cout << endl;
    cout << endl;
    idx++;
  }
  cout << endl;

  return(idx);
}


/*** continuum spectra */
void DeceTableMF6Law1Lang1(ENDF *lib6, int idx)
{
  int    nd  = lib6->rdata[idx].l1;
  int    na  = lib6->rdata[idx].l2;
  int    nep = lib6->rdata[idx].n2;
  double da  = opt.AngleStep;

  /*** if da is given then calculate actual double-differential cross sections */
  if(da > 0.0 && na >= 1){
    int np = 180.0/da;
    if( np*da == 180.0 ) np++;
    cout << "# Angle         Energy        Spectrum" << endl;

    double t = 0.0;
    while(t <= 180.0){
      for(int i1=nd ; i1<nep ; i1++){
        outVal(t);
        outVal(lib6->xptr[idx][(na+2)*i1]);
        double f = 0.0;
        for(int i2=0 ; i2<=na ; i2++){
          f += (i2+0.5) * lib6->xptr[idx][(na+2)*i1+i2+1] * legendre(i2,t) * 2.0;
        }
        outVal(f);
        cout << endl;
      }
      t += da;
      cout << endl;
    }
  }
  /*** print coefficients */
  else{
    /*** Legendre coefficients */
    cout << "# Energy      ";
    if(na == 0) cout << " Spectrum";
    else{
      for(int i0=0 ; i0 <= na ; i0++) cout << " PL(" << setw(2) << i0 << ")       ";
    }
    cout << endl;

    for(int i1=nd ; i1<nep ; i1++){
      outVal(lib6->xptr[idx][(na+2)*i1]);
      for(int i2=1 ; i2<=na+1 ; i2++) outVal(lib6->xptr[idx][(na+2)*i1+i2]);
      cout << endl;
    }
  }
}


/*** K-M systematics */
void DeceTableMF6Law1Lang2(ENDF *lib6, int idx)
{
  double e1  = lib6->rdata[idx].c2;
  int    nd  = lib6->rdata[idx].l1;
  int    na  = lib6->rdata[idx].l2;
  int    nep = lib6->rdata[idx].n2;
  double da  = opt.AngleStep;

  /*** if da is given then calculate actual double-differential cross sections */
  if(da > 0.0){

    int np = 180.0/da;
    if( np*da == 180.0 ) np++;

    double *ang, **ddx, *e2;
    ang = new double [np];
    ddx = new double * [nep];
    e2  = new double [nep];
    for(int i1=0 ; i1<nep ; i1++) ddx[i1] = new double [np];
    for(int i0=0 ; i0<np ; i0++) ang[i0] = (double)i0 * da;

    for(int i1=nd ; i1<nep ; i1++){
      e2[i1]    = lib6->xptr[idx][(na+2)*i1];
      double f  = lib6->xptr[idx][(na+2)*i1+1];
      double r  = lib6->xptr[idx][(na+2)*i1+2];
      ddxKalbach(np,e1,e2[i1],r,f,ang,ddx[i1]);
    }

    cout << "# Angle         Energy        Spectrum" << endl;
    for(int i0=0 ; i0<np ; i0++){
      for(int i1=nd ; i1<nep ; i1++){
        outVal(ang[i0]); outVal(e2[i1]); outVal(ddx[i1][i0]);
        cout << endl;
      }
      cout << endl;
    }

    delete [] ang;
    delete [] e2;
    for(int i1=0 ; i1<nep ; i1++) delete [] ddx[i1];
    delete [] ddx;
  }

  /*** print coefficients */
  else{
    cout << "# Energy       Spectrum      PEfraction" << endl;
    for(int i1=nd ; i1<nep ; i1++){
      outVal(lib6->xptr[idx][(na+2)*i1]);
      outVal(lib6->xptr[idx][(na+2)*i1+1]);
      outVal(lib6->xptr[idx][(na+2)*i1+2]);
      cout << endl;
    }
  }
}


/*** tabluated with given interpolation scheme */
void DeceTableMF6Law1Lang11(ENDF *lib6, int idx, const int lang)
{
  int    nd  = lib6->rdata[idx].l1;
  int    na  = lib6->rdata[idx].l2;
  int    nep = lib6->rdata[idx].n2;

  /*** print probabilities */

  cout << "#     Interpol" << setw(14) << lang - 10 << endl;
  cout << "#  cos(th)      Probability   Angle         Energy        Spectrum" << endl;

  for(int i1=nd ; i1<nep ; i1++){
    int i0 = (na+2)*i1;
    double e1 = lib6->xptr[idx][i0];
    double f0 = lib6->xptr[idx][i0 + 1];
    for(int i2=2 ; i2<=na ; i2+=2){

      double mu = lib6->xptr[idx][i0 + i2];
      double f1 = lib6->xptr[idx][i0 + i2 + 1];
      double s1 = f1;
      if(i2 != 1) s1 *= 2.0 * f0;

      outVal(mu);
      outVal(f1);

      outVal(acos(mu) * 180.0 / PI);
      outVal(e1);
      outVal(s1);
      cout << endl;
    }
    cout << endl;
  }
}


/**********************************************************/
/*      Two-Body Reaction Angular Distribution            */
/**********************************************************/
int DeceTableMF6Law2(ENDF *lib6, int idx)
{
  int    nr   = lib6->rdata[idx].n1;
  int    ne   = lib6->rdata[idx].n2; idx++;
  double da   = opt.AngleStep;

  cout << "#           NR" << setw(14) << nr << endl;
  cout << "#           NE" << setw(14) << ne << endl;

  for(int i0=0 ; i0<ne ; i0++){
    double e1   = lib6->rdata[idx].c2;
    int    lang = lib6->rdata[idx].l1;
    int    nw   = lib6->rdata[idx].n1;
    int    nl   = lib6->rdata[idx].n2;

    cout << "#           E1"; outVal(e1); cout << "  incident energy" << endl;
    cout << "#         LANG" << setw(14) << lang << "  0: Legendre, 12,14: tabulated probability" << endl;
    cout << "#           NW" << setw(14) << nw << "  number of parameters in LIST" << endl;
    cout << "#           NL" << setw(14) << nl << "  the highest Legendre order, or number of cosines tabulated" << endl;

    if(lang == 0){
      if(da > 0){
        int np = 180/da;
        if( np*da == 180.0 ) np++;
        cout << "#           NP" << setw(14) << np << endl;
        cout << "# Energy        Angle         Probability" << endl;
        double t = 0.0;
        while(t <= 180.0){
          double f = 0.5;
          for(int i1=0 ; i1<nl ; i1++) f += (i1+1.5)*lib6->xptr[idx][i1]*legendre(i1+1,t);
          outVal(e1);
          outVal(t);
          outVal(f);
          cout << endl;
          t += da;
        }
      }
      else{
        cout << "# Energy       Leg. Order    Coefficient" << endl;
        outVal(e1); outVal(0); outVal(1.0); cout << endl;
        for(int i1=0 ; i1<nl ; i1++){
          outVal(e1);
          outVal(i1+1);
          outVal(lib6->xptr[idx][i1]);
          cout << endl;
        }
      }
    }
    else{
      for(int i1=0 ; i1<nl ; i1++){
        double p = lib6->xptr[idx][2*i1+1];
        if(lang == 14) p = exp(p);
        outVal(lib6->xptr[idx][2*i1]);
        outVal(p);
        cout << endl;
      }
    }
    cout << endl;
    cout << endl;
    idx++;
  }
  cout << endl;

  return(idx);
}


/**********************************************************/
/*      Charged Particle Elastic Scattering               */
/**********************************************************/
int DeceTableMF6Law5(ENDF *lib6, int idx)
{
  double spi  = lib6->rdata[idx].c1;
  int    lidp = lib6->rdata[idx].l1;
  int    nr   = lib6->rdata[idx].n1;
  int    ne   = lib6->rdata[idx].n2; idx++;

  cout << "#           NR" << setw(14) << nr << endl;
  cout << "#           NE" << setw(14) << ne << endl;
  cout << "#         LIDP" << setw(14) << lidp << "  0: different, 1: identical particles" << endl;
  cout << "#          SPI" << setw(14) << spi << "  spin of the particle" << endl;

  for(int i0=0 ; i0<ne ; i0++){
    double e1  = lib6->rdata[idx].c2;
    int    ltp = lib6->rdata[idx].l1;
    int    nw  = lib6->rdata[idx].n1;
    int    nl  = lib6->rdata[idx].n2;

    cout << "#           E1"; outVal(e1); cout << "  incident energy" << endl;
    cout << "#          LTP" << setw(14) << ltp << " representation flag" << endl;
    cout << "#           NW" << setw(14) << nw << "  number of parameters in LIST" << endl;
    cout << "#           NL" << setw(14) << nl << "  the highest Legendre order, or number of cosines tabulated" << endl;

    for(int i1=0 ; i1<nw ; i1++){
      outVal(i1);
      outVal(lib6->xptr[idx][i1]);
      cout << endl;
    }
    cout << endl;
    cout << endl;
    idx++;
  }
  cout << endl;

  return(idx);
}


/**********************************************************/
/*      N-Body Phase Space Distribution                   */
/**********************************************************/
int DeceTableMF6Law6(ENDF *lib6, int idx)
{
  double apsx = lib6->rdata[idx].c1;
  int    npsx = lib6->rdata[idx].n2; idx++;

  cout << "#         APSX"; outVal(apsx); cout << "  total mass of n-particles" << endl;
  cout << "#         NPSX" << setw(14) << npsx << "  number of particles" << endl;
  cout << endl;
  cout << endl;

  return(idx);
}


/**********************************************************/
/*      Laboratory Angle-Energy Law                       */
/**********************************************************/
int DeceTableMF6Law7(ENDF *lib6, int idx)
{
  int    nr   = lib6->rdata[idx].n1;
  int    ne   = lib6->rdata[idx].n2; idx++;

  cout << "#           NR" << setw(14) << nr << endl;
  cout << "#           NE" << setw(14) << ne << endl;

  for(int i0=0 ; i0<ne ; i0++){
    double e1  = lib6->rdata[idx].c2;
    int    nrm = lib6->rdata[idx].n1;
    int    nmu = lib6->rdata[idx].n2; idx++;

    cout << "#           E1"; outVal(e1); cout << "  incident energy" << endl;
    cout << "#          NRM" << setw(14) << nrm << "  NR for cosine" << endl;
    cout << "#          NMU" << setw(14) << nmu << "  number of cosines" << endl;

    for(int i1=0 ; i1<nmu ; i1++){

      double mu  = lib6->rdata[idx].c2;
      int    nrp = lib6->rdata[idx].n1;
      int    nep = lib6->rdata[idx].n2;

      cout << "#          NRP" << setw(14) << nrp << "  NR for secondary energy" << endl;
      cout << "#          NEP" << setw(14) << nep << "  number of secondary energy points" << endl;

      for(int i2=0 ; i2<nep ; i2++){
        outVal(mu);
        outVal(lib6->xptr[idx][2*i2  ]);
        outVal(lib6->xptr[idx][2*i2+1]);
        cout << endl;
      }
      cout << endl;
      idx++;
    }
    cout << endl;
    cout << endl;
  }
  cout << endl;

  return(idx);
}
