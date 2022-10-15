/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF12 and MF15 for Fission Gamma from Data File **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);

static const int NEIN   =   100; // max number of incident energies
static const int NDAT   = 10000; // max number of floating point data
static const int WFIELD =    14; // data field width

/*
  In this program, the spectrum data in a file are assumed 
  in the following format, with the constant data field width of WFIELD.

#               0.000000e+00  7.916936e+00     // first energy and multiplicity
  0.000000e+00  1.711458e-01
  5.000000e-03  6.622625e-01
  1.500000e-02  6.802960e-01
  ...           ...
  3.945000e+01  7.237446e-27
  3.955000e+01  0.000000e+00                   // data end by a blank line

#               1.000000e+00  7.946835e+00     // second energy
  0.000000e+00  1.669042e-01
  5.000000e-03  6.758213e-01
  ...           ...
 */

// if not zero, the highest blocks will be duplicated to make a dummy data block
static const double duplicatepoint = 0.0;


void   processMF12 (int, ENDF *);
void   processMF14 (     ENDF *);
void   processMF15 (int, ENDF *);
int    dataread    (ifstream *);
int    datadummy   (int);
inline double coltodbl (string, int);

static double **xtab;
static Record  *ctab;
static double  *ydat;

int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: decemf15 data_file ENDF_file" << endl;  exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     lib12, lib14, lib15;

  datname = argv[1];
  libname = argv[2];

  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;  exit(-1);
  }
  ENDFSeekHead(&fpin,&lib12,1,451);
  fpin.close();

  /*** copy MAT number */
  lib14.setENDFmat(lib12.getENDFmat());
  lib15.setENDFmat(lib12.getENDFmat());

  /*** data arrays */
  ctab = new Record [NEIN];
  ydat = new double [NEIN * 2];
  xtab = new double * [NEIN];
  for(int i=0 ; i<NEIN ; i++) xtab[i] = new double [NDAT];

  for(int i=0 ; i<NEIN ; i++){
    ctab[i].setRecord(0.0, 0.0, 0, 0, 0,0);
    for(int j=0 ; j<NDAT ; j++) xtab[i][j] = 0.0;
  }

  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }
  int ne = dataread(&fpin);

  if(ne > 0){
    ne = datadummy(ne);
    processMF12(ne,&lib12);
    processMF14(   &lib14);
    processMF15(ne,&lib15);
  }

  delete [] ctab;
  for(int i=0 ; i<NEIN ; i++) delete [] xtab[i];
  delete [] xtab;
  delete [] ydat;

  return(0);
}


void processMF12(int ne, ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nk   = 1; // number of subsections
  int    lo   = 1; // multiplicities

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,lo,0,nk,0);
  lib->setENDFhead(head);
  lib->setENDFmf(12);
  lib->setENDFmt(18);

  Record cont(0.0,0.0,0,1,1,ne); // LF = 1, spectrum given in MF15

  int idat[2];
  idat[0] = ne; // there will be two points
  idat[1] =  2; // lin-lin interpolation

  /*** pack the data into TAB1 format */
  ENDFPackTAB1(cont,idat,ydat,lib);
  ENDFWriteMF12(lib);
}


void processMF14(ENDF *lib)
{
  Record head = lib->getENDFhead();
  int li = 1; // isotropic angular distributoin
  int nk = 0;

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,li,0,nk,0);
  lib->setENDFhead(head);
  lib->setENDFmf(14);
  lib->setENDFmt(18);
  ENDFWriteMF14(lib);
}


void processMF15(int ne, ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    nc   = 1; // number of subsections
  int    lf   = 1; // arbitrary tabulated function

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,0,0,nc,0);
  lib->setENDFhead(head);
  lib->setENDFmf(15);
  lib->setENDFmt(18);

  int    idat[2], **itab;
  double xdat[4];

  itab = new int * [NEIN];
  for(int i=0 ; i<NEIN ; i++) itab[i] = new int [2];

  /*** TAB1 for spectrum fraction.
       fraction = 1.0 in the [emin,emax] range */
  Record cont(0.0, 0.0, 0, lf, 1, 2);
  idat[0] = 2;             // there will be two points
  idat[1] = 2;             // lin-lin interpolation
  xdat[0] = ctab[0].c2;    // min incident energy
  xdat[1] = 1.0;           // fraction
  xdat[2] = ctab[ne-1].c2; // max incident energy
  xdat[3] = 1.0;           // fraction
  ENDFPackTAB1(cont,idat,xdat,lib);  // make a TAB1

  /*** TAB2 for the spectrun data */
  /*** TAB1 for ourter loop */
  cont.setRecord(0.0, 0.0, 0, 0, 1, ne);
  idat[0] = ne;            // there are NE incident energies
  idat[1] = 2;             // lin-lin interpolation

  /*** for all incident energies */
  for(int i=0 ; i<ne ; i++){
    itab[i][0] = ctab[i].n2;  // number of outgoing energies
    itab[i][1] = 1;           // histogram
  }
  ENDFPackTAB21(cont,idat,ctab,itab,xtab,lib);  // make a TAB2

  /*** output */
  ENDFWriteMF15(lib);

  for(int i=0 ; i<NEIN ; i++) delete [] itab[i];
  delete [] itab;
}


int dataread(ifstream *fp)
{
  int ne = -1, k = 0;
  string line;

  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;

    /*** data start with #-line, with the incident energy
         in the second column. */
    if(line.substr(0,1) == "#"){
      ne++;
      if(ne >= NEIN) return(0);
      double e = coltodbl(line,1) * 1e+6;
      double y = coltodbl(line,2);
      
      if(e == 0.0) e= 1e-11;

      ctab[ne].c2  = e;
      ydat[2*ne  ] = e;
      ydat[2*ne+1] = y;
      continue;
    }
    /*** data end at the first blank line,
         set number of outgoing energy points */
    else if( line.length() == 0 ){
      if(k != 0){
        ctab[ne].n1 = 1; // NR
        ctab[ne].n2 = k; // NP
        k = 0;
      }
      continue;
    }
    /*** secondary energy and spectrum data */
    xtab[ne][2*k  ] = coltodbl(line,0) * 1e+6;
    xtab[ne][2*k+1] = coltodbl(line,1) * 1e-6;
    k++;
  }

  return(ne+1);
}


int datadummy(int ne)
{
  if(duplicatepoint == 0.0) return ne;
  if(ne == NEIN) return ne;

  ctab[ne] = ctab[ne-1];
  ctab[ne].c2 = duplicatepoint;

  ydat[2*ne  ] = duplicatepoint;
  ydat[2*ne+1] = ydat[2*(ne-1)+1];

  for(int k=0 ; k<ctab[ne].n2 ; k++){
    xtab[ne][2*k]   = xtab[ne-1][2*k];
    xtab[ne][2*k+1] = xtab[ne-1][2*k+1];
  }

  return(ne+1);
}


inline double coltodbl(string c, int k)
{
  char d[WFIELD+1];
  for(int j=0 ; j<WFIELD ; j++) d[j] = c[k*WFIELD + j];
  d[WFIELD] = '\0';
  return (atof(d));
}
