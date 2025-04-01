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

/*
  The spectrum data in a file are assumed in the following format.

# 0.000000e+00  7.916936e+00                   // first energy and multiplicity
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


  For the multi-group delayed photon, the data include 7 columns,
  first the neutron energy, and spectra for 6-groups follow.
  The first row include the fractions.

#   1.00000e-11  8.531008e-04  2.346744e-03  1.692697e-03  3.844134e-03  1.438060e-03  2.744023e-04
   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
   1.000000e-02  1.003570e-07  4.929373e-01  2.794212e-02  6.112030e-01  6.074969e-01  1.103672e+00
   2.000000e-02  7.095672e-07  1.122969e-01  1.243022e-02  2.426452e+00  3.939898e-01  5.359372e-01
   3.000000e-02  2.285918e-06  8.838770e-02  1.934914e+01  2.933705e+00  3.124705e-01  3.865900e-01
 */


static const int NEIN   =   500; // max number of incident energies
static const int NDAT   = 10000; // max number of floating point data
static const int NNF    =     6; // number of delayed neutron precursor groups

static int mt = 18; // default MT number is fission

int    main(int, char *[]);
int    dataread (ifstream *, const int, int **, double ***, double **, double **);
int    datadummy   (int);
void   processMF1  (const int, const int, ENDF *, double *);
void   processMF12 (const int, const int, const int, ENDF *, double **);
void   processMF14 (const int, const int, ENDF *);
void   processMF15 (const int, const int, const int, ENDF *, int **, double **, double ***);


int main(int argc, char *argv[])
{
  if(argc <= 2){
    cerr << "usage: decemf15 data_file ENDF_file" << endl;  exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     lib1, lib12, lib14, lib15;

  datname = argv[1];
  libname = argv[2];

  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;  exit(-1);
  }
  ENDFSeekHead(&fpin,&lib1,1,451);
  fpin.close();

  /*** copy MAT number */
  lib12.setENDFmat(lib1.getENDFmat());
  lib14.setENDFmat(lib1.getENDFmat());
  lib15.setENDFmat(lib1.getENDFmat());

  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }
  else{
    /*** determine if MF = 18 (prompt) or 460 (delayed photon) */
    string dummy, line;
    getline(fpin,line);
    istringstream ss(line);
    double d0, d1, d2;
    ss >> dummy >> d0 >> d1 >> d2; 
    if(d2 > 0.0) mt = 460;
    else         mt =  18;

    fpin.seekg(0,ios_base::beg);
  }

  /*** determine number of subsections */
  int nk = (mt == 460) ? NNF : 1;

  int     **ndat = new int * [nk];
  double  **pdat = new double  * [nk];
  double  **fdat = new double  * [nk];
  double ***xdat = new double ** [nk];
  for(int k=0 ; k<nk ; k++){
    ndat[k] = new int [NEIN];
    xdat[k] = new double * [NEIN];
    pdat[k] = new double   [NEIN * 2];
    fdat[k] = new double   [NEIN * 2];
    for(int i=0 ; i<NEIN ; i++) xdat[k][i] = new double [NDAT];
  }

  /*** read data and process them */
  int ne = dataread(&fpin,nk,ndat,xdat,pdat,fdat);
  double lambda[NNF];
  if(mt == 460){
    for(int k=0 ; k<nk ; k++){
      /*** because energy-dependent lambda is not allowed, we take the zero energy data */
      lambda[k] = fdat[k][0];
      for(int i=0 ; i<ne ; i++){
        fdat[k][2*i] = pdat[k][2*i]; // replace by incident energies
      }
    }
  }


  /*** for each of subsections */
  for(int k=0 ; k<nk ; k++){

    /*** if the first point is not zero, insert it */
    for(int i=0 ; i<ne ; i++){
      if(xdat[k][i][0] != 0.0){
        /*** shift data */
        for(int j=ndat[k][i]-1 ; j>=0 ; j--){
          xdat[k][i][j*2+2] = xdat[k][i][j*2  ];
          xdat[k][i][j*2+3] = xdat[k][i][j*2+1];
        }
        xdat[k][i][0] = 0.0;
        xdat[k][i][1] = 0.0;
        ndat[k][i] ++;
      }
    }

    /*** if the lastt point is not zero, set it zero */
    for(int i=0 ; i<ne ; i++){
      if(xdat[k][i][(ndat[k][i]-1)*2+1] != 0.0) xdat[k][i][(ndat[k][i]-1)*2+1] = 0.0;
    }

    /*** eliminate zeros in the spectrum */
    for(int i=0 ; i<ne ; i++){

      int nj = ndat[k][i];
      for(int j=ndat[k][i]-2 ; j>=0 ; j--){
        if(xdat[k][i][j*2+1] != 0.0){
          nj = j + 2;
          break;
        }
      }
      ndat[k][i] = nj;
    }

    /*** re-normalize spectrum */
    for(int i=0 ; i<ne ; i++){

      double sum = 0.0;
      for(int j=1 ; j<ndat[k][i] ; j++){
        sum += (xdat[k][i][j*2  ] - xdat[k][i][j*2-2]) * (xdat[k][i][j*2+1] + xdat[k][i][j*2-1]) * 0.5;
      }

      if(sum != 0.0){
        sum = 1.0 / sum;
        for(int j=0 ; j<ndat[k][i] ; j++) xdat[k][i][2*j+1] *= sum;
      }
    }
  }

  if(ne > 0){
    if(mt == 460){
      processMF1(mt,nk,&lib1,lambda);
      processMF12(mt,nk,ne,&lib12,pdat);
      processMF14(mt,nk,&lib14);
      processMF15(mt,nk,ne,&lib15,ndat,fdat,xdat);

      ENDFWriteMF1(&lib1);
      ENDFWriteMF12(&lib12);
      ENDFWriteMF14(&lib14);
      ENDFWriteMF15(&lib15);
    }
    else{
      processMF12(mt,nk,ne,&lib12,pdat);
      processMF14(mt,nk,&lib14);
      processMF15(mt,nk,ne,&lib15,ndat,pdat,xdat);

      ENDFWriteMF12(&lib12);
      ENDFWriteMF14(&lib14);
      ENDFWriteMF15(&lib15);
    }
  }

  for(int k=0 ; k<nk ; k++){
    for(int i=0 ; i<NEIN ; i++) delete [] xdat[k][i];
    delete [] xdat[k];
    delete [] fdat[k];
    delete [] pdat[k];
    delete [] ndat[k];
  }
  delete [] pdat;
  delete [] fdat;
  delete [] xdat;
  delete [] ndat;

  return(0);
}


void processMF1(const int mt, const int nk, ENDF *lib, double *lambda)
{
  Record head = lib->getENDFhead();
  int lo = 2; // continuous representation

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,lo,0,0,0);
  lib->setENDFhead(head);
  lib->setENDFmf(1);
  lib->setENDFmt(mt);

  Record cont(0.0,0.0,0,0,nk,0);
  ENDFPackLIST(cont,lambda,lib);
}



void processMF12(const int mt, const int nk, const int ne, ENDF *lib, double **pdat)
{
  Record head = lib->getENDFhead();
  int    lo   = 1; // multiplicities
  int    lf   = 1; // normalized spectrum given in MF15

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,lo,0,nk,0);
  lib->setENDFhead(head);
  lib->setENDFmf(12);
  lib->setENDFmt(mt);

  Record cont;
  int idat[2];

  if(nk > 1){
    double *x = new double [2*ne];
    for(int i=0 ; i<ne ; i++){
      x[2*i  ] = pdat[0][2*i];
      x[2*i+1] = 0.0;
      for(int k=0 ; k<nk ; k++) x[2*i+1] += pdat[k][2*i+1]; // total yield
    }
    cont.setRecord(0.0, 0.0, 0, 0, 1, ne);
    idat[0] = ne; // there are NE energy points
    idat[1] = 2;  // lin-lin interpolation
    ENDFPackTAB1(cont,idat,x,lib);  // make a TAB1 for fractions
  }

  /*** for each of subsections */
  for(int k=0 ; k<nk ; k++){

    /*** TAB1 for spectrum fraction */
    cont.setRecord(0.0, 0.0, 0, lf, 1, ne);
    idat[0] = ne;  // there are NE energy points
    idat[1] = 2;   // lin-lin interpolation
    ENDFPackTAB1(cont,idat,pdat[k],lib);  // make a TAB1 for fractions
  }
}


void processMF14(const int mt, const int nk, ENDF *lib)
{
  Record head = lib->getENDFhead();
  int li = 1; // isotropic angular distributoin

  /*** reset index, set HEAD record, MF and MT */
  head.setRecord(head.c1,head.c2,li,0,nk,0);
  lib->setENDFhead(head);
  lib->setENDFmf(14);
  lib->setENDFmt(mt);
}


void processMF15(const int mt, const int nk, const int ne, ENDF *lib, int **ndat, double **pdat, double ***xdat)
{
  Record cont = lib->getENDFhead();
  int    lf   = 1; // arbitrary tabulated function

  /*** reset index, set HEAD record, MF and MT */
  cont.setRecord(cont.c1,cont.c2,0,0,nk,0);
  lib->setENDFhead(cont);
  lib->setENDFmf(15);
  lib->setENDFmt(mt);

  int idat[2];
  int **itab = new int * [ne];
  Record *ctab = new Record [ne];
  for(int i=0 ; i<ne ; i++) itab[i] = new int [2];

  /*** for each of subsections */
  for(int k=0 ; k<nk ; k++){

    if(nk == 1){
      double emin = pdat[0][0];
      double emax = pdat[0][2*ne-2];

      double xdum[4];

      /*** TAB1 for spectrum fraction.
           fraction = 1.0 in the [emin,emax] range */
      cont.setRecord(0.0, 0.0, 0, lf, 1, 2);
      idat[0] = 2;             // there will be two points
      idat[1] = 2;             // lin-lin interpolation
      xdum[0] = emin;          // min incident energy
      xdum[1] = 1.0;           // fraction = 1.0
      xdum[2] = emax;          // max incident energy
      xdum[3] = 1.0;           // fraction = 1.0
      ENDFPackTAB1(cont,idat,xdum,lib);  // make a TAB1 for fractions
    }
    else{
      /*** TAB1 for spectrum fraction */
      cont.setRecord(0.0, 0.0, 0, lf, 1, ne);
      idat[0] = ne;            // there are NE energy points
      idat[1] = 2;             // lin-lin interpolation
      ENDFPackTAB1(cont,idat,pdat[k],lib);  // make a TAB1 for fractions
    }

    /*** TAB2 for the spectrun data */
    cont.setRecord(0.0, 0.0, 0, 0, 1, ne);

    for(int i=0 ; i<ne ; i++){
      itab[i][0] = ndat[k][i];            // number of outgoing energies
      itab[i][1] = 2;                     // lin-lin interpolation
    }
    for(int i=0 ; i<ne ; i++) ctab[i].setRecord(0.0, pdat[k][i*2], 0, 0, 1, ndat[k][i]);
    ENDFPackTAB21(cont,ctab,idat,itab,xdat[k],lib);  // make a TAB2 for spectrum
  }

  for(int i=0 ; i<ne ; i++) delete [] itab[i];
  delete [] itab;
  delete [] ctab;
}



int dataread(ifstream *fp, const int ncol, int **ndat, double ***xdat, double **pdat, double **fdat)
{
/*  xdat[ group 1..ncol][ inc energy ][ spec data in 1dim ]
    pdat[ group 1..ncol][ (inc energy, fraction) ] */

  int i = 0, j = 0;
  string line, dummy;
  double x, y[3*ncol];

  while(getline(*fp,line)){
    if(line[0] == '#'){
      istringstream ss(line);
      ss >> dummy >> x;
      for(int k=0 ; k<3*ncol ; k++) ss >> y[k];

      if(ncol == 1){
        pdat[0][2*i  ] = x * 1.0e+6;  // incident energy, MeV assumed
        pdat[0][2*i+1] = y[0];        // multiplicity
        fdat[0][2*i  ] = x * 1.0e+6;  // incident energy, MeV assumed
        fdat[0][2*i+1] = 1.0;         // fraction
      }
      else{
        for(int k=0 ; k<ncol ; k++){
          pdat[k][2*i  ] = x * 1.0e+6;      // incident energy, MeV assumed
          pdat[k][2*i+1] = y[k];            // multiplicities
          fdat[k][2*i  ] = y[ncol + 2*k];   // lambda
          fdat[k][2*i+1] = y[ncol + 2*k+1]; // fraction
        }
      }
      continue;
    }
    else if(line.length() == 0){
      if(j != 0){
        for(int k=0 ; k<ncol ; k++) ndat[k][i] = j;
        i ++;
        j = 0;

        if(i >= NEIN){
          cerr << "too many energies" << endl;
          exit(-1);
        }
      }
      continue;
    }

    /*** secondary energy and spectrum data */
    istringstream ss(line);
    ss >> x;
    for(int k=0 ; k<ncol ; k++) ss >> y[k];

    for(int k=0 ; k<ncol ; k++){
      xdat[k][i][2*j  ] = x    * 1.0e+6;
      xdat[k][i][2*j+1] = y[k] * 1.0e-6;
    }
    j++;

    if(j >= NDAT/2){
      cerr << "too many data points" << endl;
      exit(-1);
    }
  }

  int ne = i;
/*
  for(int k=0 ; k<ncol ; k++){
    for(int i=0 ; i<ne ; i++){
      cout << "### " << k << "  " << ndat[k][i] << " " <<  pdat[k][2*i] << " " << pdat[k][2*i+1] << endl;

      for(int j=0 ; j<ndat[k][i] ; j++){
        cout << " " <<  xdat[k][i][2*j] << " " << xdat[k][i][2*j+1] << endl;;
      }
      cout << endl;
    }
    cout << endl;
  }
*/
  return(ne);
}



