/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF4 for Resonance Reconstructed P(L)           **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);

static const int NEIN   = 10000; // number of incodent energies
static const int LMAX   =     4; // max order of Legendre coefficients
static const int WFIELD =    13; // data field width
/*
  In this program, the Pl data are givne in the following format,
  with the constant data field width of WFIELD.

# Energy[eV]   P1/P0        P2/P0        P3/P0        P4/P0        
  1.00000e-05  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
  5.00000e+03 -6.33430e-03  1.16188e-03 -1.20863e-08  8.32435e-12
  1.00000e+04  2.45790e-02  1.22970e-04  2.63004e-07  1.21758e-10
  1.50000e+04  1.78008e-03  3.94202e-04  1.88554e-07  1.88837e-10
  2.00000e+04 -4.36081e-06  1.88401e-05 -2.80346e-08  4.26228e-10
  ...
 */

int    dataread   (ifstream *);
void   processMF4 (int, ENDF *);
inline double coltodbl   (string, int);

static double **xtab;
static Record  *ctab;

int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: decemf4 data_file ENDF_file" << endl;  exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     lib(L);

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
  for(int i=0 ; i<NEIN ; i++) xtab[i] = new double [LMAX+1];

  for(int i=0 ; i<NEIN ; i++){
    ctab[i].setRecord(0.0, 0.0, 0, 0, 0,0);
    for(int j=0 ; j<LMAX+1 ; j++) xtab[i][j] = 0.0;
  }

  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }
  int ne = dataread(&fpin);

  if(ne > 0) processMF4(ne,&lib);

  delete [] ctab;
  for(int i=0 ; i<NEIN ; i++) delete [] xtab[i];
  delete [] xtab;

  return(0);
}


void processMF4(int ne, ENDF *lib)
{
  Record head = lib->getENDFhead();
  int    idat[2];

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,0,1,0,0);
  lib->setENDFhead(head);
  lib->setENDFmf(4);
  lib->setENDFmt(2);

  Record cont;
  cont.setRecord(0.0,head.c2,0,2,0,0);
  ENDFPackCONT(cont,lib);

  cont.setRecord(0.0,0.0,0,0,1,ne);
  idat[0] = ne;
  idat[1] = 2;

  ENDFPackTAB2(cont,ctab,idat,xtab,lib);

  ENDFWriteMF4(lib);
}


int dataread(ifstream *fp)
{
  int k = 0;
  string line;

  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;

    /*** skip comment lines */
    if(line.substr(0,1) == "#") continue;

    for(int l=0 ; l<LMAX ; l++) xtab[k][l] = coltodbl(line,l+1);
    int lm = LMAX;
    bool zero = true;
    for(int l=LMAX-1 ; l>=0 ; l--){
      if(xtab[k][l] != 0.0){
        lm = l+1;
        zero = false;
        break;
      }
    }
    if(zero) lm = 2;
    ctab[k].setRecord(0.0,coltodbl(line,0),0,0,lm,0);
    k++;
  }

  return(k);
}


inline double coltodbl(string c, int k)
{
  char d[WFIELD+1];
  for(int j=0 ; j<WFIELD ; j++) d[j] = c[k*WFIELD + j];
  d[WFIELD] = '\0';
  return (atof(d));
}
