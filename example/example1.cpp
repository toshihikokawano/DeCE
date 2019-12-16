// Example for use of endflib.cpp and endfio.cpp

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);

int main(int argc, char *argv[])
{
  ifstream fpin;   // file pointer to input library
  ENDF     lib(M); // allocate medium size data block

  if(argc <= 1){
    cerr << "ENDF file not supplied" << endl;
    exit(-1);
  }

  /*** read in MF3 MT1 */
  fpin.open(argv[1]);
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;
    exit(-1);
  }
  ENDFReadMF3(&fpin,&lib,1);
  fpin.close();

  /*** print it in table format */
  ENDFPrint1Dim(&lib,0);

  /*** print it in ENDF format */
  ENDFWriteMF3(&lib);

  /*** access to the data */
  int n = lib.rdata[0].n2;
  for(int i=0 ; i<n ; i++){
    cout <<  setw(14) << lib.xdata[2*i  ] / 1e+06;
    cout <<  setw(14) << lib.xdata[2*i+1] * 1000.0;
    cout << endl;
  }

  return(0);
}
