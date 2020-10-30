/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF5 for Fission Spectrum from Tabulated Data   **/
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
static const int NDAT   =  2000; // max number of floating point data
static const int WFIELD =    14; // data field width
/*
  In this program, the spectrum data in a file are assumed 
  in the following format, with the constant data field width of WFIELD.

#               1.000000e-05   // start with #-line, energy in 
  0.000000e+00  0.000000e+00   // secondary energy, spectrum data
  1.000000e+01  1.964834e-09
  1.100000e+01  2.060735e-09
  1.200000e+01  2.152367e-09

  ...           ...

  2.980000e+07  2.763210e-17
  3.000000e+07  2.309210e-17
  3.020000e+07  0.000000e+00   // data end by a blank line



#               5.000000e+05   // start the next energy data
  0.000000e+00  0.000000e+00
  1.000000e+01  1.951657e-09
  1.100000e+01  2.046915e-09
 */

int    dataread   (ifstream *);
void   processMF5 (int, ENDF *);
inline double coltodbl   (string, int);

static double **xtab;
static Record  *ctab;

int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: decemf5 data_file ENDF_file" << endl;  exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     lib;

  datname = argv[1];
  libname = argv[2];

  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;  exit(-1);
  }
  ENDFSeekHead(&fpin,&lib,1,451);
  fpin.close();

  /*** data arrays */
  ctab = new Record [NEIN];
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

  if(ne > 0) processMF5(ne,&lib);

  delete [] ctab;
  for(int i=0 ; i<NEIN ; i++) delete [] xtab[i];
  delete [] xtab;

  return(0);
}


void processMF5(int ne, ENDF *lib)
{
  Record cont = lib->getENDFhead();
  int    nk   = 1; // number of subsections
  int    lf   = 1; // arbitrary tabulated function

  /*** reset index, set HEAD record, MF and MT */
  cont.setRecord(cont.c1,cont.c2,0,0,nk,0);
  lib->setENDFhead(cont);
  lib->setENDFmf( 5);
  lib->setENDFmt(18);

  int    idat[2], **itab;
  double xdat[4];

  itab = new int * [NEIN];
  for(int i=0 ; i<NEIN ; i++) itab[i] = new int [2];

  /*** TAB1 for spectrum fraction.
       fraction = 1.0 in the [emin,emax] range */
  cont.setRecord(0.0, 0.0, 0, lf, 1, 2);
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
    itab[i][1] = 2;           // lin-lin interpolation
  }
  ENDFPackTAB21(cont,idat,ctab,itab,xtab,lib);  // make a TAB2

  /*** output */
  ENDFWriteMF5(lib);

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
      ctab[ne].c2 = coltodbl(line,1);
      continue;
    }
    /*** data end at the first blank line,
         set number of outgoing energy points */
    else if( line.length() ==0 ){
      if(k != 0){
        ctab[ne].n1 = 1; // NR
        ctab[ne].n2 = k; // NF
        k = 0;
      }
      continue;
    }
    /*** secondary energy and spectrum data */
    xtab[ne][2*k  ] = coltodbl(line,0);
    xtab[ne][2*k+1] = coltodbl(line,1);
    k++;
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
