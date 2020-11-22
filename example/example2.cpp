// Example for use of endflib.cpp and endfio.cpp

#include <iostream>

using namespace std;

#include "../source/endflib.h"

int main(void);

static const int NELAB = 13;
static double x[NELAB] 
  = {8.083816e+6, 9.000000e+6, 1.000000e+7, 1.100000e+7,  1.200000e+7,
     1.300000e+7, 1.400000e+7, 1.500000e+7, 1.600000e+7,  1.700000e+7,
     1.800000e+7, 1.900000e+7, 2.000000e+7};

static double y[NELAB] 
  = {0.000000e+0, 1.978970e-1, 5.880670e-1, 8.896460e-1, 1.060970e+0,
     1.149730e+0, 1.183550e+0, 1.188810e+0, 1.178160e+0, 1.153660e+0,
     1.105680e+0, 1.053960e+0, 9.783820e-1};

int main(void)
{
  ENDF   lib, cpy;  // objects

  /*** data to be stored */
  double za  = 3.307400e+4;
  double awr = 7.328890e+1;
  double q   = -7.975000e+6;
  int    nr  =  1;
  int    np  = NELAB;

  /*** head and cont */
  Record head(za,awr,0,0,0,0);
  Record cont(q,q,0,0,nr,np);

  /*** int and double data arrays */
  int    idat[2];
  double xdat[2*NELAB];

  /*** make an interpolation scheme data */
  idat[0] = NELAB;
  idat[1] = 2;

  /*** store floating point data */
  for(int i=0 ; i<NELAB ; i++){
    xdat[2*i  ] = x[i];
    xdat[2*i+1] = y[i];
  }

  /*** set MAT, MF, MT, and HEAD */
  lib.setENDFmat(3322);
  lib.setENDFmf(3);
  lib.setENDFmt(16);
  lib.setENDFhead(head);

  /*** pack the data into TAB1 format */
  ENDFPackTAB1(cont,idat,xdat,&lib);

  /*** print it in ENDF format */
  ENDFWriteMF3(&lib);
  ENDFLibPeek(&lib);

  /*** copy the whole data into another */
  ENDFLibCopy(&lib,&cpy);

  /*** change MT into 17 */
  cpy.setENDFmt(17);

  /*** print it in ENDF format */
  ENDFWriteMF3(&cpy);
  ENDFLibPeek(&cpy);

  return(0);
}
