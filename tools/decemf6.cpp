/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF6 from CoH ECLIPSE output                    **/
/**                                                            Oct. 2010     **/
/**                                                            T. Kawano     **/
/**                                       Los Alamos National Laboratory     **/
/**                                                                          **/
/**                                                                          **/
/**     usage: decemf6 -t MTnumber -e ECLIPSE.dat -f ENDF.dat                **/
/**                                                                          **/
/**     notes: this code is a part of the DeCE package                       **/
/**            ENDF.dat must have corresponding subsection in MF3 section    **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <unistd.h>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);

static const int NEIN   =   100; // max number of incident energies
static const int NDAT   =  5000; // max number of floating point data
static const int NPAR   =     7; // max number of particles

static const int WFIELD =    13; // data field width

static int mat = 9999;           // MAT number
static int mt  =    0;           // MT number

int    dataread   (ifstream *, string);
int    datarecord (ifstream *, int, int, double);

void   processMF6 (int, ENDF *);
void   mf6spec    (int, int, double, ENDF *);
void   mf6yield   (int, int, int, double, double, ENDF *);
string reactionMT (void);

inline Record mf6head    (ENDF *);
inline double coltodbl   (string, int);
inline string coltostr   (string, int);

static double *elab, *gyield, ***spc;
static int    *ng, **ns, **nl;


/**********************************************************/
/*      DeCE MF6 Main Program                             */
/**********************************************************/
int main(int argc, char *argv[])
{
  ifstream fpin;
  string   libname = "", eclname = "", reacid = "";
  ENDF     lib;

  /*** command line options */
  int p = 0;
  while((p=getopt(argc,argv,"e:f:t:"))!=-1){
    switch(p){
    case 'e':
      eclname = optarg;
      break;
    case 'f':
      libname = optarg;
      break;
    case 't':
      mt = atoi(optarg);
      break;
    default:
      break;
    }
  }

  if(mt <= 0 || mt > 849){
    cerr << "MT number " << mt << " out of range" << endl;
    exit(-1);
  }
  reacid = reactionMT();

  /*** Read ENDF MF3 */
  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;
    exit(-1);
  }
  int icond = ENDFReadMF3(&fpin,&lib,mt);
  fpin.close();
  if(icond < 0){
    cerr << "MT number not given in MF3 for " << mt << endl;
    exit(-1);
  }

  mat = lib.getENDFmat();

  /*** open ECLIPSE file */
  fpin.open(eclname.c_str());
  if(!fpin){
    cerr << "ECLIPISE file cannot open" << endl;
    exit(-1);
  }

  /*** data arrays */
  elab   = new double [NEIN];
  gyield = new double [NEIN];

  ng     = new int [NEIN];
  ns     = new int * [NPAR];
  nl     = new int * [NPAR];
  spc    = new double ** [NPAR];
  for(int p=0 ; p<NPAR ; p++){
    ns[p]  = new int [NEIN];
    nl[p]  = new int [NEIN];
    spc[p] = new double * [NEIN];
    for(int n=0 ; n<NEIN ; n++) spc[p][n] = new double [NDAT];
  }

  for(int n=0 ; n<NEIN ; n++){
    ng[n] = 0;
    elab[n] = gyield [n] = 0.0;
  }
  for(int p=0 ; p<NPAR ; p++){
    for(int n=0 ; n<NEIN ; n++) ns[p][n] = nl[p][n] = 0;
  }

  int ne = dataread(&fpin,reacid);
  fpin.close();

  if(ne > 0) processMF6(ne,&lib);

  delete [] elab;
  delete [] gyield;

  for(int p=0 ; p<NPAR ; p++){
    for(int n=0 ; n<NEIN ; n++) delete [] spc[p][n];
    delete [] ns[p];
    delete [] nl[p];
    delete [] spc[p];
  }
  delete [] ng;
  delete [] ns;
  delete [] nl;
  delete [] spc;

  return(0);
}


/**********************************************************/
/*      Generage Subsection in MF6                        */
/**********************************************************/
void processMF6(int ne, ENDF *lib3)
{
  double emin,emax;
  ENDF   lib;

  lib.setENDFmat(mat);
  lib.setENDFmf(6);
  lib.setENDFmt(mt);
  lib.setENDFhead( mf6head(lib3) );


  /*** neutron capture case */
  if( mt == 102 ){
    emin = elab[0]; // at zero energy
    emax = elab[ne-1];
  }

  /*** particle emission case */
  else{
    emin = lib3->xdata[0];
    emax = elab[ne-1];

    string part = reactionMT();

    for(int pid=1 ; pid<=NPAR-1 ; pid++){
      int n = (int)part[2*pid] - '0';
      if(n > 0){
        mf6yield(ne,pid,n,emin,emax,&lib);
        mf6spec(ne,pid,emin,&lib);
      }
    }
  }

  /*** finally, gamma-rays */
  mf6yield(ne,0,0,emin,emax,&lib);
  mf6spec(ne,0,emin,&lib);

  ENDFWriteMF6(&lib);
}


/**********************************************************/
/*      MF6 Particle Yield Section                        */
/**********************************************************/
void mf6yield(int nelab, int pid, int nyield, double emin, double emax, ENDF *lib)
{
  int    idat[2], lip = 0, law = 1, nr = 1, np = 2;
  double *xdat, zap = 0.0, awr = 0.0, yield;

  xdat = new double [NEIN*2];

  if(     pid == 0){ zap =    0.0;   awr = 0.000000; }
  else if(pid == 1){ zap =    1.0;   awr = 1.000000; }
  else if(pid == 2){ zap = 1001.0;   awr = 0.999167; }
  else if(pid == 3){ zap = 2004.0;   awr = 3.968220; }
  else if(pid == 4){ zap = 1002.0;   awr = 1.996800; }
  else if(pid == 5){ zap = 1003.0;   awr = 2.990140; }
  else if(pid == 6){ zap = 2003.0;   awr = 3.000034; }

  int ki = 0;
  idat[ki++] = 2;
  idat[ki++] = 2;

  int kx = 0;
  if(pid == 0){
    if(mt == 102){
      np  = nelab;
      for(int i=0 ; i<nelab ; i++){
        xdat[kx++] = elab[i];
        xdat[kx++] = gyield[i];
      }
    }
    else{
      yield = gyield[0];
      for(int i=0 ; i<nelab ; i++){
        if(elab[i] > emin){
          yield = gyield[i];
          break;
        }
      }
      xdat[kx++] = emin;
      xdat[kx++] = yield;

      for(int i=0 ; i<nelab ; i++){
        if(elab[i] > emin){
          xdat[kx++] = elab[i];
          xdat[kx++] = gyield[i];
        }
      }
      np = kx/2;
    }
    idat[0] = np;
  }
  else{
    np  = 2;
    yield = (double)nyield;
    xdat[kx++] = emin;
    xdat[kx++] = yield;
    xdat[kx++] = emax;
    xdat[kx++] = yield;
  }

  Record cont(zap,awr,lip,law,nr,np);
  ENDFPackTAB1(cont,idat,xdat,lib);
//ENDFPrint1Dim(lib,0);

  delete [] xdat;
}


/**********************************************************/
/*      MF6 Particle Energy Spectra Section               */
/**********************************************************/
void mf6spec(int nelab, int pid, double emin, ENDF *lib)
{
  double **xdat;
  Record  *xcon;
  int      idat[2], lang = 1, lep = 2, nr = 1, ne = 1;
  int      na, nd, nw, np;

  /*** count number of block for a real NE */
  ne = 1;
  for(int i=0 ; i<nelab ; i++){
    if(elab[i] > emin) ne++;
  }

  xcon = new Record [ne+1];
  xdat = new double * [ne+1];
  for(int i=0 ; i<ne+1 ; i++) xdat[i] = new double [NDAT];

  /*** at threshold, dummy spectrum */
  int kx = 0, idx = 0;
  if(mt != 102){
    xdat[idx][kx++] = 0.0;
    xdat[idx][kx++] = 0.0;
    xdat[idx][kx++] = 0.5;
    xdat[idx][kx++] = 2.0;
    xdat[idx][kx++] = 1.0;
    xdat[idx][kx++] = 0.0;
    nw = kx;
    np = kx/2;
    nd = 0;
    na = 0;
    xcon[idx++].setRecord(0.0,emin,nd,na,nw,np);
  }

  /*** for each incident energy */
  for(int i=0 ; i<nelab ; i++){

    if(elab[i] >= emin){
      nd = (pid == 0) ? ng[i] : 0;
      na = nl[pid][i] - 1;
      np = ns[pid][i] + nd;
      nw = (na+2)*np;

      /*** skip this energy if zero-sum for a1
           except for the last energy point */
      kx = 0;
      double sum = 0.0;
      for(int k=0 ; k<np ; k++){
        for(int l=0 ; l<=nl[pid][i] ; l++){
          if(l == 1){
            sum += spc[pid][i][kx];
          }
          kx++;
        }
      }
      if((sum == 0.0) && (i != nelab-1)){
        ne --;
        continue;
      }

      /*** when no spectrum given, copy a dummy */
      kx = 0;
      if(np == 0){
        xdat[idx][kx++] = 0.0;
        xdat[idx][kx++] = 0.0;
        xdat[idx][kx++] = 0.5;
        xdat[idx][kx++] = 2.0;
        xdat[idx][kx++] = 1.0;
        xdat[idx][kx++] = 0.0;
        nw = kx;
        np = kx/2;
        nd = 0;
        na = 0;
      }
      else{
        for(int k=0 ; k<np ; k++){
          for(int l=0 ; l<=nl[pid][i] ; l++){
            xdat[idx][kx] = spc[pid][i][kx];
            kx++;
          }
        }
      }
      xcon[idx++].setRecord(0.0,elab[i],nd,na,nw,np);
    }
  }
 
  idat[0] = ne;
  idat[1] = 22;

  Record cont(0.0,0.0,lang,lep,nr,ne);
  ENDFPackTAB2(cont,xcon,idat,xdat,lib);

  for(int i=0 ; i<ne+1 ; i++) delete [] xdat[i];
  delete [] xdat;
  delete [] xcon;
}


/**********************************************************/
/*      Analyze ECLIPSE Data                              */
/**********************************************************/
int dataread(ifstream *fp, string rid)
{
  string line;
  int    m = -1;  // number of energy points in the file

  bool rec = false;
  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;
    if( line.length() ==0 ) continue;

    string key = line.substr(0,WFIELD);

    /*** new incident energy section starts */
    if(key == "# ECLIPSE    "){
      m++;
      if(m >= NEIN){ m = -1; break; }

      elab[m] = coltodbl(line,1) * 1e+06;
      rec = false;
      continue;
    }

    /*** sub-section in the energy loop, particle emission identifier */
    else if(key == "# Nucleus    "){
      rec  = (coltostr(line,1) == rid) ? true : false;
      continue;
    }

    /*** if the MT is what we are looking for, get them */
    else if(key == "# Particle   "){
      if(rec){
        int    pid  = (int)coltodbl(line,1);
        int    nsec = (int)coltodbl(line,2);
        double yld  =      coltodbl(line,3);

        if(nsec > 0){
          if( (datarecord(fp, m, pid, yld)) < 0 ){
            cerr << "too many data point" << endl;  return(-1);
          }
        }
      }
    }
  }

  return(m+1);
}


/**********************************************************/
/*      Read ECLIPSE Output                               */
/**********************************************************/
int datarecord(ifstream *fp, int m, int pid, double yld)
{
  string line;

  /*** gamma-ray */
  if(pid == 0){
    gyield[m] = yld;

    /*** read in discrete gamma-rays */
    getline(*fp,line);
    ng[m] = (int)coltodbl(line,1);

    int k=0;
    for(int i=0 ; i<ng[m] ; i++){
      if(k >= NDAT) return(-1);
      getline(*fp,line);
      spc[pid][m][k++] = coltodbl(line,0) * 1e+06;
      spc[pid][m][k++] = coltodbl(line,1);
    }

    getline(*fp,line);
    getline(*fp,line);

    /*** read in continuum gamma-rays */
    getline(*fp,line);
    ns[pid][m] = (int)coltodbl(line,1);
    nl[pid][m] = (int)coltodbl(line,2);

    for(int i=0 ; i<ns[pid][m] ; i++){
      if(k >= NDAT) return(-1);
      getline(*fp,line);
      spc[pid][m][k++] = coltodbl(line,0) * 1e+06;
      spc[pid][m][k++] = coltodbl(line,1) * 1e-06;
    }
  }

  /*** other particles */
  else if(pid < NPAR){

    getline(*fp,line);
    ns[pid][m] = (int)coltodbl(line,1);
    nl[pid][m] = (int)coltodbl(line,2);

    int k=0;
    for(int i=0 ; i<ns[pid][m] ; i++){
      if(k >= NDAT) return(-1);
      getline(*fp,line);
      spc[pid][m][k++] = coltodbl(line,0) * 1e+06;
      for(int j=1 ; j<=nl[pid][m] ; j++){
        if(k >= NDAT) return(-1);
        spc[pid][m][k++] = coltodbl(line,j) * 1e-06;
      }
    }
  }
  return(0);
}


/**********************************************************/
/*      Determine MT Number from Number of Particles      */
/**********************************************************/
string reactionMT(void)
{
  string p = "";

  switch(mt){
  case  16: p = " n2p0a0d0t0h0"; break;
  case  17: p = " n3p0a0d0t0h0"; break;
  case  22: p = " n1p0a1d0t0h0"; break;
  case  23: p = " n1p0a3d0t0h0"; break;
  case  24: p = " n2p0a1d0t0h0"; break;
  case  25: p = " n3p0a1d0t0h0"; break;
  case  28: p = " n1p1a0d0t0h0"; break;
  case  29: p = " n1p0a2d0t0h0"; break;
  case  30: p = " n2p0a2d0t0h0"; break;
  case  32: p = " n1p0a0d1t0h0"; break;
  case  33: p = " n1p0a0d0t1h0"; break;
  case  34: p = " n1p0a0d0t0h1"; break;
  case  35: p = " n1p0a2d1t0h0"; break;
  case  36: p = " n1p0a2d0t1h0"; break;
  case  37: p = " n4p0a0d0t0h0"; break;
  case  41: p = " n2p1a0d0t0h0"; break;
  case  42: p = " n3p1a0d0t0h0"; break;
  case  44: p = " n1p2a0d0t0h0"; break;
  case  45: p = " n1p1a1d0t0h0"; break;
  case  91: p = " n1p0a0d0t0h0"; break;
  case 102: p = " n0p0a0d0t0h0"; break;
  case 103: p = " n0p1a0d0t0h0"; break;
  case 104: p = " n0p0a0d1t0h0"; break;
  case 105: p = " n0p0a0d0t1h0"; break;
  case 106: p = " n0p0a0d0t0h1"; break;
  case 107: p = " n0p0a1d0t0h0"; break;
  case 108: p = " n0p0a2d0t0h0"; break;
  case 109: p = " n0p0a3d0t0h0"; break;
  case 111: p = " n0p2a0d0t0h0"; break;
  case 112: p = " n0p1a1d0t0h0"; break;
  case 113: p = " n0p0a2d0t1h0"; break;
  case 114: p = " n0p0a2d1t0h0"; break;
  case 115: p = " n0p1a0d1t0h0"; break;
  case 116: p = " n0p1a0d0t1h0"; break;
  case 117: p = " n0p0a1d1t0h0"; break;
  case 649: p = " n0p1a0d0t0h0"; break;
  case 699: p = " n0p0a0d1t0h0"; break;
  case 749: p = " n0p0a0d0t1h0"; break;
  case 799: p = " n0p0a0d0t0h1"; break;
  case 849: p = " n0p0a1d0t0h0"; break;
  default:  p = ""; break;
  }

  return (p);
}


/**********************************************************/
/*      MF6 HEAD Record                                   */
/**********************************************************/
inline Record mf6head(ENDF *lib)
{
  Record head = lib->getENDFhead();
  double za  = head.c1;
  double awr = head.c2;
  int    lct = 2;   // CMS
  int    nk  = 1;   // number of subsections

  string part = reactionMT();
  nk = 1;
  for(int i=1 ; i<=6 ; i++){
    int n = (int)part[2*i] - '0';
    if(n > 0) nk++;
  }

  head.setRecord(za, awr, 0, lct, nk, 0);
  return(head);
}


/**********************************************************/
/*      Convert Data Field into Numerics or Text          */
/**********************************************************/
inline double coltodbl(string c, int k)
{
  char d[WFIELD+1];
  for(int j=0 ; j<WFIELD ; j++) d[j] = c[k*WFIELD + j];
  d[WFIELD] = '\0';
  return (atof(d));
}


inline string coltostr(string c, int k)
{
  char d[WFIELD+1];
  for(int j=0 ; j<WFIELD ; j++) d[j] = c[k*WFIELD + j];
  d[WFIELD] = '\0';
  return ((string)d);
}

