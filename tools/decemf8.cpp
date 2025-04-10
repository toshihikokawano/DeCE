/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF8 Fission Product Yield Data                 **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"


/*-----------------------------------------------------------
  This program assumes that the FPY data are in the following format.

#   1.00000e-11   // start with #-line followed by the incident energy
     31  74   0   3.824614e-08   9.908241e-06  // Z, A, meta, indep, cumul yields
     31  74   1   8.836146e-09   8.836146e-09
     31  75   0   4.224868e-06   2.063774e-05
...
     64 159   0   0.000000e+00   1.519093e-06
     64 160   0   0.000000e+00   6.408843e-07
     64 161   0   0.000000e+00   0.000000e+00   // data end by a blank line

#   5.00000e-01  // start the next energy data
     31  74   0   5.547339e-08   1.166481e-05
     31  74   1   1.292866e-08   1.292866e-08

  The number of FPYs at each of inciden energies should be the same.
*/

static const int NEIN   =   100; // max number of incident energies
static const int NDAT   =  5000; // max number of floating point data


int    main (int, char *[]);
int    dataread (ifstream *, const int, unsigned int *, double *, double *, double *);
void   processMF8 (const int, const double, ENDF *, unsigned int *, double *);
double uncertainty (const double);


int main(int argc, char *argv[])
{
  if(argc <= 2){
    cerr << "usage: decemf8 data_file ENDF_file" << endl;  exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     lib,libi,libc;

  datname = argv[1];
  libname = argv[2];

  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;  exit(-1);
  }

  ENDFSeekHead(&fpin,&lib,1,451);
  Record head = lib.getENDFhead();
  fpin.close();

  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }

  /*** data arrays */
  unsigned int *id = new unsigned int [NDAT];  // nuclide ID
  double *en = new double [NEIN];              // incident energy array
  double *xi = new double [NDAT];              // independent yield data, [energy] x [yield]
  double *xc = new double [NDAT];              // cumulative yield data
  for(int i=0 ; i<NEIN ; i++){ en[i] = 0.0; }
  for(int j=0 ; j<NDAT ; j++){ xi[j] = xc[j] = 0.0; }

  libi.setENDFmat(lib.getENDFmat());
  libc.setENDFmat(lib.getENDFmat());

  libi.setENDFmf(8);
  libc.setENDFmf(8);
  
  libi.setENDFmt(454);
  libc.setENDFmt(459);

  Record cont(head.c1,head.c2,0,0,0,0);
  libi.setENDFhead(cont);
  libc.setENDFhead(cont);

  /*** read FPY data */
  int ne = 0;
  for(int n=0 ; n<NEIN ; n++){
    int nk = dataread(&fpin,n,id,en,xi,xc);
    if(nk < 0) break;
    processMF8(nk,en[n],&libi,id,xi);
    processMF8(nk,en[n],&libc,id,xc);
    ne ++;
  }
  fpin.close();

  int le = ne - 1;

  /*** change the number of energy points in HEAD */
  cont.setRecord(head.c1,head.c2,le+1,0,0,0);
  libi.setENDFhead(cont);
  libc.setENDFhead(cont);

  /*** change the first CONT */
  cont = libi.rdata[0];
  libi.rdata[0].setRecord(cont.c1,cont.c2,le,cont.l2,cont.n1,cont.n2);

  cont = libc.rdata[0];
  libc.rdata[0].setRecord(cont.c1,cont.c2,le,cont.l2,cont.n1,cont.n2);

  ENDFWrite(&libi);
  ENDFWrite(&libc);

  /*** post process */
  delete [] xi;
  delete [] xc;
  delete [] id;
  delete [] en;

  return(0);
}


void processMF8(const int nfp, const double en, ENDF *lib, unsigned int *id, double *fpy)
{
  int    intp = 22;
  int    nn = 4 * nfp;
  double *xdat = new double [nn];

  for(int i=0 ; i<nfp ; i++){
    int zafp = id[i]/10;
    int fps  = id[i] - zafp*10;

    xdat[i*4  ] = (double)zafp;
    xdat[i*4+1] = (double)fps;
    xdat[i*4+2] = fpy[i];
    xdat[i*4+3] = uncertainty(fpy[i]);
  }

  Record cont(en, 0.0, 0, intp, nn, nfp);
  ENDFPackLIST(cont,xdat,lib);

  delete [] xdat;
}


double uncertainty(const double x)
{
  double z = 0.0;
  if(     x >  0.04){ z = 0.10; }
  else if(x >  0.01){ z = 0.20; }
  else if(x > 0.005){ z = 0.50; }
  else              { z = 0.90; }

  return(z * x);
}


int dataread(ifstream *fp, const int n, unsigned int *id, double *en, double *xi, double *xc)
{
  int i = 0, k = -1;
  string line, dummy;
  unsigned int z = 0, a = 0, m = 0;

  while(getline(*fp,line)){
    if(line[0] == '#'){
      istringstream ss(line);
      ss >> dummy >> en[n];
      en[n] *= 1.0e+6;
      continue;
    }
    else if( line.length() == 0 ){
      if(i != 0){
        k = i;
        break;
      }
      continue;
    }
    istringstream ss(line);
    ss >> z >> a >> m >> xi[i] >> xc[i];

    id[i] = (unsigned int)(10000 * z + 10 * a + m);
    i++;
  }

  // for(i = 0 ; i < k ; i++){
  //   std::cout << n << " " << en[n] << "  " << i << " " << id[i] << " " << xi[i] << " " << xc[i] << std::endl;
  // }

  return k;
}
