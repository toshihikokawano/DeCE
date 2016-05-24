// Example for use of endflib.cpp

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "endflib.h"

int main(int, char *[]);

int main(int argc, char *argv[])
{
  ifstream fpin;
  ENDF     lib('M');

  if(argc<=1){
    cerr << "ENDF file not supplied" << endl;
    exit(-1);
  }

  fpin.open(argv[1]);
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;
    exit(-1);
  }
  /*** read in MF3 MT1 */
  ENDFReadMF3(&fpin,&lib,1);
  fpin.close();

  /*** print it in table format */
  ENDFPrint1Dim(&lib,0);

  /*** access to the data */
  int n = lib.rdata[0].n2;
  for(int i=0 ; i<n ; i++){
    cout <<  setw(14) << lib.xdata[2*i  ] / 1e+06;
    cout <<  setw(14) << lib.xdata[2*i+1] * 1000.0;
    cout << endl;
  }

  return(0);
}
