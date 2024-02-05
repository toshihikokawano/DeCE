/******************************************************************************/
/**     DeCE ANGDIST                                                         **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <ostream>
#include <cmath>
#include <cstring>

using namespace std;

#include "constant.h"
#include "dece.h"
#include "global.h"
#include "terminate.h"
#include "polysq.h"
#include "masstable.h"

static const int ANGLE_POINTS = 180;
static const int MAX_ENERGY   = 500;
static const int MAX_LEGCOEF  =  60;

static int  readADdata (char *, int, int, double *, double **, double *);
static int  geneADdata (int, int, double, double *, double *, double **, double  **, Record *);
static void storeMF4   (int, int, double **, Record *, ENDF *);
static void storeMF6   (int, int, double **, Record *, ENDF *, char *);
static int  discretegammaspec (const int, const int, const int, double *, double *, double *, char *);

static bool   addata = true;
static double za = 0.0, awr = 0.0;
static int    mat = 0;

/**********************************************************/
/*      Read in External Data from a File                 */
/**********************************************************/
void DeceAngdist(ENDFDict *dict, ENDF *lib[], const int mf, const int mt, char *datafile, char *gammafile, int ofset)
{
  double   *cx=NULL, **cy=NULL, *en=NULL, **lg=NULL;
  Record   *cont;

  if( (mf != 4) && (mf !=6 ) ) return;
  if( !((mt == 2) || (51 <= mt && mt <= 90) || (mt >= 600)) ) return;

  int k3 = dict->getID(3,mt);
  if(k3 < 0){
    message << "process skipped since no MT" << mt << " found in MF3";
    WarningMessage();
    DeceDelete(dict,mf,mt);
    return;
  }

  int k4 = dict->getID(4,mt);
  if(mf == 6) k4 = dict->getID(6,mt);

  /*** allocate data array and open data file */
  en = new double   [MAX_ENERGY];
  cx = new double   [ANGLE_POINTS];
  cy = new double * [MAX_ENERGY];
  lg = new double * [MAX_ENERGY];
  cont = new Record [MAX_ENERGY];

  for(int i=0 ; i<MAX_ENERGY ; i++){
    cy[i] = new double [ANGLE_POINTS];
    lg[i] = new double [MAX_LEGCOEF];
  }

  /*** read angular distribution data */
  int ne = readADdata(datafile,ofset,mt,cx,cy,en);
  if(ne == 0){
    message << "no data to be added from " << datafile << " for MT = " << mt;
    WarningMessage();
  }

  if(ne > 0){
    /*** threshold energy */
    double eth = (mt == 2) ? 1e-05 : lib[k3]->xdata[0];

    /*** generate floating point data */
    ne = geneADdata(mt,ne,eth,en,cx,cy,lg,cont);

    /*** ZA, AWR, and MAT number from Dictionary */
    mat  = dict->getMAT();
    za   = dict->getZA();
    awr  = dict->getAWR();

    if(mf == 4) storeMF4(mt,ne,lg,cont,lib[k4]);
    else        storeMF6(mt,ne,lg,cont,lib[k4],gammafile);
  }
  else{
    DeceDelete(dict,mf,mt);
  }

  /*** Clean all */
  for(int i=0 ; i<MAX_ENERGY ; i++){
    delete [] cy[i];
    delete [] lg[i];
  }
  delete [] en;
  delete [] cx;
  delete [] cy;
  delete [] lg;
  delete [] cont;

  return;
}


/**********************************************************/
/*      Read in Angular Distribution Data                 */
/**********************************************************/
int readADdata(char *file, int ofset, int mt, double *x, double **y, double *en)
{
  ifstream fp;
  string   line,legflag;

  if(ofset == 0){
    if(mt == 2) ofset = 1;
    else if( ( 51 <= mt) && (mt <=  91) ) ofset = mt -  50 + 1;
    else if( (600 <= mt) && (mt <= 640) ) ofset = mt - 600 + 1;
    else if( (650 <= mt) && (mt <= 690) ) ofset = mt - 650 + 1;
    else if( (700 <= mt) && (mt <= 740) ) ofset = mt - 700 + 1;
    else if( (750 <= mt) && (mt <= 790) ) ofset = mt - 750 + 1;
    else if( (800 <= mt) && (mt <= 840) ) ofset = mt - 800 + 1;
  }

  if(ofset == 0){
    message << "MF4:MT" << mt << " cannot read from " << file << " because ofset is not given";
    Notice("DeceAngdist:readADdata");
    return 0;
  }

  fp.open(file);
  if(!fp){ message << "cannot open data file " << file; TerminateCode("readADdata"); }

  int ne=0;
  while(getline(fp,line)){

    istringstream s1(&line[1]);  // skip comment #
    s1 >> en[ne] >> legflag;

    /*** convert energy unit into eV */
    en[ne] *= opt.ReadXdataConversion;

    bool nonzero = false;
    if(legflag == "Legcoef"){
      for(int i=0 ; i<MAX_LEGCOEF ; i++){
        getline(fp,line);
        istringstream s2(line);
        s2 >> x[i];
        for(int j=0 ; j<ofset ; j++) s2 >> y[ne][i];
      }
      if(y[ne][0] != 0.0) nonzero = true;
      addata = false;
    }else{
      for(int i=0 ; i<ANGLE_POINTS ; i++){
        getline(fp,line);
        istringstream s2(line);
        s2 >> x[i];
        for(int j=0 ; j<ofset ; j++) s2 >> y[ne][i];
      }
      if(y[ne][0] != 0.0) nonzero = true;
      addata = true;
    }

    /*** skip data if range is set by options */
    if(DeceCheckReadRange(x[ne])) continue;

    if(nonzero) ne++;

    if(ne >= MAX_ENERGY){ message << "too many energy points for angular distributions, " << ne; TerminateCode("readADdata"); }
  }

  fp.close();

  message << "MT" << mt << " angular distributions at " << ne << " energy points imported from " << file;
  Notice("DeceAngdist:readADdata");

  return ne;
}


/**********************************************************/
/*      Generate 2-Dim Array of Legendre Coefficients     */
/**********************************************************/
int geneADdata(int mt, int ne, double eth, double *en, double *x, double **y, double **lg, Record *cont)
{
  const double eps1 = 1.0e-03, eps2 = 1.0e-12;
  double a[MAX_LEGCOEF];

  for(int i=0 ; i<MAX_ENERGY ; i++){
    for(int l=0 ; l<MAX_LEGCOEF ; l++) lg[i][l] = 0.0;
  }

  int idx = 0;

  /*** insert the first point (zero or threshold energy) */
  cont[idx].setRecord(0.0,eth,0,0,2,0);
  lg[idx][0] = 0.0;
  lg[idx][1] = 0.0;  idx++;

  int mmax = 0;
  for(int i=0 ; i<ne ; i++){

    /*** copy data if larger than Eth */
    if(en[i] <= eth) continue;

    a[0] = 1.0;
    for(int j=1 ; j<MAX_LEGCOEF ; j++) a[j] = 0.0;

    int m = 0;
    /*** data file given by tabulated distribution */
    if(addata) m = polysq(ANGLE_POINTS,MAX_LEGCOEF,x,y[i],a);
    /*** data given in Legendre coefficients */
    else{
      /*** make the number of points even */
      for(int j=0 ; j<MAX_LEGCOEF ; j+=2){
        a[j  ] = y[i][j  ];
        a[j+1] = y[i][j+1];
        if(a[j] == 0.0){ m = j-1; break; }
      }
    }

    if(m <= 1){
      cont[idx].setRecord(0.0,en[i],0,0,2,0);
      lg[idx][0] = 0.0;
      lg[idx][1] = 0.0;
    }
    else{
      double p = a[1]/(3.0*a[0]);

      /*** eliminate numerical noise for elastic at low energies */
      if( (x[i] <= 10000.0) && (fabs(p) < eps1) && (mt == 2) ){
        cont[idx].setRecord(0.0,en[i],0,0,2,0);
        lg[idx][0] = 0.0;
        lg[idx][1] = 0.0;
      }
      else{
        if(m >= mmax){
          cont[idx].setRecord(0.0,en[i],0,0,m-1,0);
          for(int i=1 ; i<=m ; i++){
            p = a[i]/((2*i+1)*a[0]);
            lg[idx][i-1] = (fabs(p) > eps2) ? p : 0.0;
          }
        }
        else{
          cont[idx].setRecord(0.0,en[i],0,0,mmax-1,0);
          for(int i=1 ; i<=mmax ; i++){
            if(i <= m){
              p = a[i]/((2*i+1)*a[0]);
              lg[idx][i-1] = (fabs(p) > eps2) ? p : 0.0;
            }
            else{
              lg[idx][i-1] = 0.0;
            }
          }
        }
      }
    }
    /*** save the largest L */
    if(m >= mmax) mmax = m;

    idx++;
  }
/*
  for(int i=0 ; i<idx ; i++){
    cout << cont[i].c2 <<" " << cont[i].n1 << "  "<< lg[i][0] << " " << lg[i][1] << endl;
  }
*/
  return(idx);
}


/**********************************************************/
/*      Store Legendre Coefficients in ENDF lib (MF4)     */
/**********************************************************/
void storeMF4(const int mt, const int ne, double **xdat, Record *xcont, ENDF *lib)
{
  Record cont;
  int    idat[2];

  /*** Make HEAD and CONT */
  int    lvt  = 0; // transformation matrix not given
  int    ltt  = 1; // Legendre parameters given
  int    li   = 0; // not isotropic
  int    lct  = 2; // center-of-mass system

  lib->setENDFmat(mat);
  lib->setENDFmf(4);
  lib->setENDFmt(mt);
  lib->setENDFhead(za, awr, lvt, ltt, 0, 0);

  /*** extra card for MF4 */
  cont.setRecord(0.0, awr, li, lct, 0, 0);
  ENDFPackCONT(cont,lib);

  /*** Make TAB2 (LIST) */
  cont.setRecord(0.0, 0.0, 0, 0, 1, ne);
  idat[0] = ne;     // there are NE incident energies
  idat[1] = 2;      // lin-lin interpolation
  ENDFPackTAB2(cont, xcont, idat, xdat, lib);
}


/**********************************************************/
/*      Store Legendre Coefficients in ENDF lib (MF6)     */
/**********************************************************/
void storeMF6(const int mt, const int ne, double **xdat, Record *xcont, ENDF *lib, char *gfile)
{
  double  zap  = 1.0; // product identifier (assume neutron)
  double  awp  = 1.0; // product mass

  int pid = 0, mt0 = mt;
  if(      ( 50 <= mt) && (mt <=  91)) { zap =    1.0; awp = 1.0;                   pid = 1; mt0 =  50; }
  else if( (600 <= mt) && (mt <= 649)) { zap = 1001.0; awp = MPROTON   / MNEUTRON;  pid = 2; mt0 = 600; }
  else if( (650 <= mt) && (mt <= 699)) { zap = 1002.0; awp = MDEUTERON / MNEUTRON;  pid = 4; mt0 = 650; }
  else if( (700 <= mt) && (mt <= 749)) { zap = 1003.0; awp = MTRITON   / MNEUTRON;  pid = 5; mt0 = 700; }
  else if( (750 <= mt) && (mt <= 799)) { zap = 2003.0; awp = MHELIUM3  / MNEUTRON;  pid = 6; mt0 = 750; }
  else if( (800 <= mt) && (mt <= 849)) { zap = 2004.0; awp = MALPHA    / MNEUTRON;  pid = 3; mt0 = 800; }

  bool gammaspec = false;
  if( (strlen(gfile) > 0) && (mt != mt0) ) gammaspec = true;

  Record cont;
  int    idat[2];
  double zdat[4];

  int     lct  = 2;   // center-of-mass system
  int     nk   = (gammaspec) ? 3 : 2;   // number of subsection, incl. recoil and gamma lines
  int     lip  = 0;   // isomer flag
  int     law  = 2;   // two-body scattering

  /*** Make HEAD CONT */
  lib->setENDFmat(mat);
  lib->setENDFmf(6);
  lib->setENDFmt(mt);
  lib->setENDFhead(za, awr, 0, lct, nk, 0);


  /*** TAB1 for fraction = 1.0 in the [emin,emax] range */
  cont.setRecord(zap, awp, lip, law, 1, 2);
  idat[0] = 2;              // two points
  idat[1] = 2;              // lin-lin interpolation
  zdat[0] = xcont[0].c2;    // min incident energy
  zdat[1] = 1.0;            // fraction = 1.0
  zdat[2] = xcont[ne-1].c2; // max incident energy
  zdat[3] = 1.0;            // fraction = 1.0
  ENDFPackTAB1(cont, idat, zdat, lib);


  /*** Make TAB2 (LIST) */
  cont.setRecord(0.0, 0.0, 0, 0, 1, ne);
  idat[0] = ne;             // NE incident energies
  idat[1] = 2;              // lin-lin interpolation


  /*** change NW in the case of MF6 */
  for(int i=0 ; i<ne ; i++){
    Record r =   xcont[i];
    r.n2 = r.n1;
    xcont[i] = r;
  }
  ENDFPackTAB2(cont, xcont, idat, xdat, lib);


  /*** Make TAB1 for recoil */
  zap = za + 1.0 - zap;     // neutron incident assumed
  int z = zap/1000;
  int a = zap - z*1000;
  double mass = mass_excess(z,a);
  awp = (mass / AMUNIT + a) / MNEUTRON;  // recoil mass
  law = 4;                  // discrete two body recoil
  cont.setRecord(zap, awp, lip, law, 1, 2);
  idat[0] = 2;              // two points
  idat[1] = 2;              // lin-lin interpolation
  ENDFPackTAB1(cont, idat, zdat, lib);

  /*** add discrete gamma-ray spectrum */
  if(gammaspec){
    int nlev = mt - mt0;

    /*** possible total number of transitions */
    int ntot = nlev * (nlev + 1)/2; 

    double *x = new double [ntot];
    double *y = new double [ntot];
    double q = 0.0; // multiplicity

    /*** gamma-ray spectrum data */
    int ngamma = discretegammaspec(mt,mt0,pid,x,y,&q,gfile);

    zap = 0.0;
    awp = 0.0;
    lip = 0;
    law = 1;

    zdat[1] = zdat[3] = q;

    /*** first CONT and yield */
    cont.setRecord(zap, awp, lip, law, 1, 2);
    ENDFPackTAB1(cont, idat, zdat, lib);

    /*** TAB2 for two points */
    int lang = 1;
    int lep  = 2;
    cont.setRecord(0.0, 0.0, lang, lep, 1, 2);

    int nd  = ngamma;
    int na  = 0;
    int nep = ngamma;
    int nw  = nep * (na + 2);
    Record cdat[2];
    double *xtab[2];
    for(int i=0 ; i<2 ; i++) xtab[i] = new double [2*ngamma];

    cdat[0].setRecord(0.0, zdat[0], nd, na, nw, nep);
    cdat[1].setRecord(0.0, zdat[2], nd, na, nw, nep);
    for(int i=0 ; i<ngamma ; i++){
      xtab[0][2*i  ] = xtab[1][2*i  ] = x[i];
      xtab[0][2*i+1] = xtab[1][2*i+1] = y[i];
    }
    ENDFPackTAB2(cont, cdat, idat, xtab, lib);

    delete [] x;
    delete [] y;
    for(int i=0 ; i<2 ; i++) delete [] xtab[i];
  }
}



/**********************************************************/
/*      Produce Discrete Transision Spectrum              */
/**********************************************************/
int discretegammaspec(const int mt, const int mt0, const int pid, double *x, double *y, double *q, char *datafile)
{
  const double pcut = 1e-11; // cut-off for gamma-ray production probability 
  const int mlev = 50; // max number of levels allowed by ENDF

  int nlev = mt - mt0;
  int ntot = nlev * (nlev + 1)/2;

  double **br = new double * [mlev];
  double **gp = new double * [mlev];
  int    **fs = new int * [mlev];
  int     *ng = new int [mlev];
  double  *ex = new double [mlev];

  for(int i=0 ; i<mlev ; i++){
    br[i] = new double [mlev]; // branching ratio
    gp[i] = new double [mlev]; // gamma-ray probability
    fs[i] = new int [mlev];    // decay final state
    for(int j=0 ; j< mlev ; j++){
      br[i][j] = gp[i][j] = 0.0;
      fs[i][j] = 0;
    }
    ex[i] = 0.0;
    ng[i] = 0;
  }

  /*** read gamma-ray branching ratio data */
  ifstream fpin;
  fpin.open(datafile);
  if(!fpin){
    message << "data file " << datafile << " cannot open";
    TerminateCode("DeCEMod6:DeceAddDiscrete:discretegammaspec");
  }

  bool   withconv = false;
  string dat1, dat2, dat3, line;
  while(1){
    getline(fpin,line);
    if(fpin.eof() != 0) break;
    if(line == "#gammaray"){
      withconv = false; break;
    }
    else if(line == "#gammaraywithconversion"){
      withconv = true; break;
    }
  }

  while(1){
    getline(fpin,line);
    if(fpin.eof() != 0) break;
    if(line.length() == 0) break;

    dat1 = line.substr( 5, 4);  int n = atoi(&dat1[0]); // number of levels
    dat2 = line.substr( 0, 5);  int p = atoi(&dat2[0]); // particle ID

    /*** for each discrete level */
    for(int k=0 ; k<n ; k++){
      getline(fpin,line);

      /*** when found, copy branching ratio data to array */
      if( (p == pid) && (k < mlev) ){

        dat1 = line.substr( 9, 4);  // number of gamma-rays
        dat2 = line.substr(13,13);  // level energy

        ng[k] = atoi(&dat1[0]);
        ex[k] = atof(&dat2[0]) * 1e+6;

        for(int g=0 ; g<ng[k] ; g++){
          if(withconv){
            dat1 = line.substr(26 + 30*g,  4); // final state
            dat2 = line.substr(30 + 30*g, 13); // branching ratio
            dat3 = line.substr(43 + 30*g, 13); // gamma-ray emission probability
          }
          else{
            dat1 = line.substr(26 + 17*g,  4); // final state
            dat2 = line.substr(30 + 17*g, 13); // branching ratio
          }
          fs[k][g] = atoi(&dat1[0]);
          br[k][g] = atof(&dat2[0]);
          gp[k][g] = (withconv) ? atof(&dat3[0]) : 1.0;
        }
      }
    }
  }
  fpin.close();

  /*** level population for gamma-ray cascade calculation */
  double *lpop = new double [nlev + 1];
  for(int k=0; k<nlev ; k++) lpop[k] = 0.0;
  lpop[nlev] = 1.0;

  /*** generate discrete gamma spectrum from given level */
  int n = 0;
  for(int k=nlev ; k>=0 ; k--){
    for(int g=0 ; g<ng[k] ; g++){
      double eg = ex[k] - ex[ fs[k][g] ];  // gamma-ray energy
      double dp = br[k][g] * gp[k][g];     // emission probability
      lpop[ fs[k][g] ] += dp;              // increment population of the final state

      /*** discrete gamma spectrum */
      double p = dp * lpop[k];
      if(p > pcut){
        x[n] = eg;
        y[n] = dp * lpop[k];
        n++;
        if(n >= ntot) break;
      }
    }
    if(n >= ntot) break;
  }
  delete [] lpop;


  /*** sort the gamma lines by energies */
  double z;
  for(int j=0 ; j<n ; j++){
    int k = j;
    for(int i=j ; i<n ; i++){
      if(x[i] > x[k]) k = i;
    }
    z = x[j];  x[j] = x[k];  x[k] = z;
    z = y[j];  y[j] = y[k];  y[k] = z;
  }

  /*** normalize spectrum and calculate multiplicity */
  double psum = 0.0, eave = 0.0;
  for(int i=0 ; i<n ; i++){
    psum += y[i];
    eave += x[i] * y[i];
  }
  if(psum > 0.0){
    eave /= psum;
    for(int i=0 ; i<n ; i++) y[i] /= psum;
    *q = ex[nlev] / eave;
  }
  else *q = 0.0;

  for(int i=0 ; i< mlev ; i++){
    delete [] br[i];
    delete [] gp[i];
    delete [] fs[i];
  }
  delete [] br;
  delete [] gp;
  delete [] fs;

  delete [] ng;
  delete [] ex;

  message << "number of " << n << " gammas included in " << mt << " produced by " << datafile << " from ";
  message << setw(13) << setprecision(6) << x[n-1] << " to ";
  message << setw(13) << setprecision(6) << x[0];
  Notice("DeceMod6:discretegammaspec");

  return n;
}
