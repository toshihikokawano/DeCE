/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Calculate photo reaction cross section from ENDF Data   **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);
static void DecePhotoProduction(const int, const int, ENDF *);
static inline void printdata(const string);

static int np0, counter = 0;
static double *x0, *y0, *y1, *y2;

int main(int argc, char *argv[])
{
  if(argc <=  1){
    cerr << "usage: decephoto ENDF_file" << endl; exit(-1);
  }

  ifstream fpin;    // file pointer to input library
  ENDF     lib3(M); // allocate cross section data block
  ENDF     lib6(L); // allocate residual production data block

  x0 = new double [MAX_DBLDATA]; // common grid from photo-absorption cross section
  y0 = new double [MAX_DBLDATA]; // photo-abs cross section
  y1 = new double [MAX_DBLDATA]; // multiplicity
  y2 = new double [MAX_DBLDATA]; // cumulative multiplicities

  /*** read in MF6 MT5 */
  fpin.open(argv[1]);
  if(!fpin){
    cerr << "ENDF file cannot open" << endl; exit(-1);
  }

  int id3 = ENDFReadMF3(&fpin,&lib3,3);
  int id5 = ENDFReadMF3(&fpin,&lib3,5);

  if(     id3 < 0) ENDFReadMF3(&fpin,&lib3,5);
  else if(id5 < 0) ENDFReadMF3(&fpin,&lib3,3);
  else             ENDFReadMF3(&fpin,&lib3,3);

  ENDFReadMF6(&fpin,&lib6,5);
  fpin.close();

  /*** ZA */
  Record head = lib3.getENDFhead();
  int zat =(int) head.c1;

  for(int i=0 ; i<MAX_DBLDATA ; i++) x0[i] = y0[i] = 0.0;


  cout.setf(ios::scientific, ios::floatfield);
  cout << "# " << setw(3) << counter << " photo absorption" << endl;
  np0 = lib3.rdata[0].n2;
  for(int i=0 ; i<np0 ; i++){
    x0[i] = lib3.xdata[2*i  ];
    y0[i] = lib3.xdata[2*i+1];
    cout << setprecision(6) << setw(14) << x0[i] * 1e-6;
    cout << setprecision(6) << setw(14) << y0[i] * 1e+3 << endl;
  }
  cout << endl;
  cout << endl;
  counter ++;

  /*** neutron production, ZAP = 1 */
  DecePhotoProduction(1,np0,&lib6);
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("neutron production");

  /*** (g,n) */
  DecePhotoProduction(zat - 1,np0,&lib6);    // (g,n)
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("(g,n)");

  /*** (g,1nx) */
  DecePhotoProduction(zat - 1002,np0,&lib6); // (g,n+p)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 2005,np0,&lib6); // + (g,n+A)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  printdata("(g,1nx)");

  /*** (g,2n) */
  DecePhotoProduction(zat - 2,np0,&lib6);    // (g,2n)
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("(g,2n)");

  /*** (g,2nx) */
  DecePhotoProduction(zat - 1003,np0,&lib6); // + (g,2n+p)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 2006,np0,&lib6); // + (g,2n+A)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  printdata("(g,2nx)");

  /*** (g,3n) */
  DecePhotoProduction(zat - 3,np0,&lib6);    // (g,3n)
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("(g,3n)");

  /*** (g,sn) */
  DecePhotoProduction(zat - 1,np0,&lib6);    // (g,n)
  for(int i=0 ; i<np0 ; i++) y2[i]  = y1[i];
  DecePhotoProduction(zat - 2,np0,&lib6);    // + (g,2n)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 3,np0,&lib6);    // + (g,3n)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 1002,np0,&lib6); // + (g,n+p)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 2005,np0,&lib6); // + (g,n+A)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 1003,np0,&lib6); // + (g,2n+p)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  DecePhotoProduction(zat - 2006,np0,&lib6); // + (g,2n+A)
  for(int i=0 ; i<np0 ; i++) y2[i] += y1[i];
  printdata("(g,sn)");

  /*** (g,p) */
  DecePhotoProduction(zat - 1001,np0,&lib6); // (g,p)
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("(g,p)");

  /*** (g,alpha) */
  DecePhotoProduction(zat - 2004,np0,&lib6); // (g,A)
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("(g,a)");

  /*** (g,4n) */
  DecePhotoProduction(zat - 4,np0,&lib6);    // (g,4n)
  for(int i=0 ; i<np0 ; i++) y2[i] = y1[i];
  printdata("(g,4n)");


  delete [] x0;
  delete [] y0;
  delete [] y1;
  delete [] y2;

  return 0;
}


static inline void printdata(string title)
{
  bool out = false;
  for(int i=0 ; i<np0 ; i++){
    if(y2[i] > 0.0){ out = true; break; }
  }

  if(out){
    cout << "# " << setw(3) << counter << " " << title << endl;
    for(int i=0 ; i<np0 ; i++){
      if(y2[i] == 0.0) continue;
      cout << setprecision(6) << setw(14) << x0[i] * 1e-6;
//    cout << setprecision(6) << setw(14) << y2[i];
      cout << setprecision(6) << setw(14) << y0[i] * y2[i] * 1e+3<< endl;
    }
    cout << endl;
    cout << endl;
    counter ++;
  }
}


/**********************************************************/
/*      Process MF=6                                      */
/**********************************************************/
void DecePhotoProduction(const int k, const int np0, ENDF *lib6)
{
  for(int i=0 ; i<np0 ; i++) y1[i] = 0.0;

  Record head = lib6->getENDFhead();
  int    nk   = head.n1;
  int    idx  = 0;

  for(int ik=0 ; ik<nk ; ik++){
    Record cont = lib6->rdata[idx];
    int    zap  = (int)cont.c1;
//  int    lip  = cont.l1;
    int    law  = cont.l2;
    int    np1  = cont.n2;

    if(k == zap){
      for(int i0=0 ; i0<np0 ; i0++){

        if(x0[i0] == lib6->xptr[idx][(np1-1)*2]) y1[i0] += lib6->xptr[idx][(np1-1)*2+1];
        else{
          for(int i1=0 ; i1<np1-1 ; i1++){
            double xa = lib6->xptr[idx][i1*2  ];
            double ya = lib6->xptr[idx][i1*2+1];
            double xb = lib6->xptr[idx][i1*2+2];
            double yb = lib6->xptr[idx][i1*2+3];
            if((xa == x0[i0]) && (xa == xb)){
              y1[i0] += yb;
              break;
            }
            if((xa <= x0[i0]) && (x0[i0] < xb)){
              if(xb > xa) y1[i0] += (yb - ya)/(xb - xa) * (x0[i0] - xa) + ya;
              break;
            }
          }
        }
      }
    }

    /*** increment index */
    idx++;
    if(law == 1){
      int ne = lib6->rdata[idx].n2;
      idx++;
      for(int i0=0 ; i0<ne ; i0++) idx++;
    }
  }
}


