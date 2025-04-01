/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF5 for Prompt and Delayed Fission Spectrum    **/
/**                  from Tabulated Data                                     **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

/*-----------------------------------------------------------
  In this program, the prompt / delayed spectrum data in a file are assumed 
  in the following format.

# 1.000000e-05   // start with #-line, energy in 
  0.000000e+00  0.000000e+00   // secondary energy, spectrum data
  1.000000e+01  1.964834e-09
  ...           ...
  3.000000e+07  2.309210e-17
  3.020000e+07  0.000000e+00   // data end by a blank line

#               5.000000e+05   // start the next energy data
  0.000000e+00  0.000000e+00
  1.000000e+01  1.951657e-09

  For the multi-group delayed neutron, the data include 7 columns,
  first the neutron energy, and spectra for 6-groups follow.
  The first row include the fractions.

# 1.0000e-05  8.3283e-02  2.7600e-01  2.0660e-01  3.0814e-01  1.0511e-01  2.0871e-02
  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00
  1.0000e-02  1.9485e-05  3.5025e-05  2.5216e-05  4.1858e-05  2.8676e-05  6.0682e-06
  2.0000e-02  1.9009e-05  5.5210e-05  2.8938e-05  4.6476e-05  3.2495e-05  7.7854e-06
  ...
  1.2240e+01  0.0000e+00  0.0000e+00  ...    // data end by a blank line

# 5.0000e+05  8.2560e-02  2.7625e-01  ...    // start the next energy data
  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00
  1.0000e-02  1.9500e-05  3.5733e-05  2.5372e-05  4.1917e-05  2.8781e-05  6.3048e-06 
  ...
 */

static const int NEIN   =   500; // max number of incident energies
static const int NDAT   =  5000; // max number of floating point data
static const int NNF    =     6; // number of delayed neutron precursor groups

static int mt = 18; // default MT number is fission

int    main(int, char *[]);
int    dataread (ifstream *, const int, int **, double ***, double **);
void   processMF5 (const int, const int, ENDF *, int **, double ***, double **);


int main(int argc, char *argv[])
{
  if(argc <= 2){
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

  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }
  else{
    /*** determine if MF = 18 (prompt) or 455 (delayed) */
    string dummy, line;
    getline(fpin,line);
    istringstream ss(line);
    double d0, d1;
    ss >> dummy >> d0 >> d1; 
    if(d1 > 0.0) mt = 455;
    else         mt =  18;

    fpin.seekg(0,ios_base::beg);
  }

  /*** determine number of subsections */
  int nk = (mt == 455) ? NNF : 1;
    
  int     **ndat = new int * [nk];
  double  **pdat = new double  * [nk];
  double ***xdat = new double ** [nk];
  for(int k=0 ; k<nk ; k++){
    ndat[k] = new int [NEIN];
    xdat[k] = new double * [NEIN];
    pdat[k] = new double   [NEIN * 2];
    for(int i=0 ; i<NEIN ; i++) xdat[k][i] = new double [NDAT];
  }

  /*** read data and process them */
  int ne = dataread(&fpin,nk,ndat,xdat,pdat);

  if(ne > 0) processMF5(nk,ne,&lib,ndat,xdat,pdat);

  /*** output */
  ENDFWriteMF5(&lib);

  for(int k=0 ; k<nk ; k++){
    for(int i=0 ; i<NEIN ; i++) delete [] xdat[k][i];
    delete [] xdat[k];
    delete [] pdat[k];
    delete [] ndat[k];
  }
  delete [] pdat;
  delete [] xdat;
  delete [] ndat;

  return(0);
}


void processMF5(const int nk, const int ne, ENDF *lib, int **ndat, double ***xdat, double **pdat)
{
  Record cont = lib->getENDFhead();
  int    lf   = 1;   // arbitrary tabulated function

  /*** reset index, set HEAD record, MF and MT */
  cont.setRecord(cont.c1,cont.c2,0,0,nk,0);
  lib->setENDFhead(cont);
  lib->setENDFmf( 5);
  lib->setENDFmt(mt);

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

    Record *ctab = new Record [ne];
    int idat[2], **itab;
    itab = new int * [ne];
    for(int i=0 ; i<ne ; i++) itab[i] = new int [2];

    if(nk == 1){
      double emin = pdat[0][0];
      double emax = pdat[0][2*ne-2];

      /*** TAB1 for spectrum fraction.
           fraction = 1.0 in the [emin,emax] range */
      cont.setRecord(0.0, 0.0, 0, lf, 1, 2);
      idat[0] = 2;             // there will be two points
      idat[1] = 2;             // lin-lin interpolation
      pdat[0][0] = emin;       // min incident energy
      pdat[0][1] = 1.0;        // fraction = 1.0
      pdat[0][2] = emax;       // max incident energy
      pdat[0][3] = 1.0;        // fraction = 1.0
    }
    else{
      /*** TAB1 for spectrum fraction */
      cont.setRecord(0.0, 0.0, 0, lf, 1, ne);
      idat[0] = ne;            // there are NE energy points
      idat[1] = 2;             // lin-lin interpolation
    }
    ENDFPackTAB1(cont,idat,pdat[k],lib);  // make a TAB1 for fractions

    /*** TAB2 for the spectrun data */
    cont.setRecord(0.0, 0.0, 0, 0, 1, ne);

    for(int i=0 ; i<ne ; i++){
      itab[i][0] = ndat[k][i];            // number of outgoing energies
      itab[i][1] = 2;                     // lin-lin interpolation
    }
    for(int i=0 ; i<ne ; i++) ctab[i].setRecord(0.0, pdat[k][i*2], 0, 0, 1, ndat[k][i]);
    ENDFPackTAB21(cont,ctab,idat,itab,xdat[k],lib);  // make a TAB2 for spectrum

    for(int i=0 ; i<ne ; i++) delete [] itab[i];
    delete [] itab;
    delete [] ctab;
  }
}


int dataread(ifstream *fp, const int ncol, int **ndat, double ***xdat, double **pdat)
{
/*  xdat[ group 1..ncol][ inc energy ][ spec data in 1dim ]
    pdat[ group 1..ncol][ (inc energy, fraction) ] */

  int i = 0, j = 0;
  string line, dummy;
  double x, y[ncol];

  while(getline(*fp,line)){
    if(line[0] == '#'){
      istringstream ss(line);
      ss >> dummy >> x;
      if(ncol > 1){
        for(int k=0 ; k<ncol ; k++) ss >> y[k];
      }
      else y[0] = 1.0;

      for(int k=0 ; k<ncol ; k++){
        pdat[k][2*i  ] = x * 1.0e+6;  // incident energy, MeV assumed
        pdat[k][2*i+1] = y[k];        // fractions
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

