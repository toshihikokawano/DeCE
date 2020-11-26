/******************************************************************************/
/**     DeCE ANGDIST                                                         **/
/******************************************************************************/

#include <iostream>
#include <ostream>
#include <cmath>

using namespace std;

#include "dece.h"
#include "global.h"
#include "terminate.h"
#include "polysq.h"

static const int ANGLE_POINTS = 180;
static const int MAX_ENERGY   = 100;
static const int MAX_LEGCOEF  =  60;

static int  readADdata (char *, int, int, double *, double **, double *);
static int  geneADdata (int, int, double, double *, double *, double **, double  **, Record *);
static void storeMF4   (int, int, double **, Record *, ENDF *);
static void storeMF6   (int, int, double **, Record *, ENDF *);

static bool   addata = true;
static double za = 0.0, awr = 0.0;
static int    mat = 0;

/**********************************************************/
/*      Read in External Data from a File                 */
/**********************************************************/
void DeceAngdist(ENDFDict *dict, ENDF *lib[], const int mf, const int mt, char *datafile, int ofset)
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
    else        storeMF6(mt,ne,lg,cont,lib[k4]);
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

  fp.open(file);
  if(!fp) TerminateCode("cannot open data file",(string)file);

  if(ofset == 0){
    if(mt == 2) ofset = 1;
    else if( ( 51 <= mt) && (mt <=  91) ) ofset = mt -  50 + 1;
    else if( (600 <= mt) && (mt <= 640) ) ofset = mt - 600 + 1;
    else if( (650 <= mt) && (mt <= 690) ) ofset = mt - 650 + 1;
    else if( (700 <= mt) && (mt <= 740) ) ofset = mt - 700 + 1;
    else if( (750 <= mt) && (mt <= 790) ) ofset = mt - 750 + 1;
    else if( (800 <= mt) && (mt <= 840) ) ofset = mt - 800 + 1;
  }

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
    if(ne >= MAX_ENERGY) TerminateCode("too many energy points for angular distributions");
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
void storeMF4(int mt, int ne, double **xdat, Record *xcont, ENDF *lib)
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
void storeMF6(int mt, int ne, double **xdat, Record *xcont, ENDF *lib)
{
  Record cont;
  int    idat[2];
  double zdat[4];

  int     lct  = 2;   // center-of-mass system
  int     nk   = 2;   // number of subsection, incl. recoil
  int     lip  = 0;   // isomer flag
  int     law  = 2;   // two-body scattering
  double  zap  = 1.0; // product identifier (assume neutron)
  double  awp  = 1.0; // product mass

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
  zap = za;
  awp = awr;
  law = 4;                  // discrete two body recoil
  cont.setRecord(zap, awp, lip, law, 1, 2);
  idat[0] = 2;              // two points
  idat[1] = 2;              // lin-lin interpolation
  ENDFPackTAB1(cont, idat, zdat, lib);
}

