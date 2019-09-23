// Example for use of endflib.cpp and endfio.cpp

/*
   This example first scans the entire ENDF file 
   and store the index information in the ENDFDict object.
   Then, for each MF/MT, ENDFExtract copies from the original file to the standard output.

 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);

int main(int argc, char *argv[])
{
  ENDFDict dict;   // allocate ENDF Dictionary

  if(argc <= 1){
    cerr << "ENDF file not supplied" << endl;
    exit(-1);
  }

  /*** scan ENDF file and store data in dict */
  if(ENDFScanLibrary((string)argv[1],&dict) < 0){
    cerr << "ENDF file cannot open " << argv[1] << endl;
    exit(-1);
  }

  /*** read original ENDF file section by section, and print each of them */
  ifstream fpin;
  fpin.open(argv[1]);

  /*** MT numbers in each MF */
  int *mt = new int [1000];

  /*** for all MF numbers */
  for(int mf=1 ; mf <= 40 ; mf++){

    /*** check if MT subsections are given in this MF */
    int k = 0;
    for(int i=0 ; i<dict.getSEC() ; i++) if(dict.mf[i] == mf) mt[k++] = dict.mt[i];

    /*** if given, print it */
    if(k > 0){
      sort(mt,mt+k); // sort MT, just in case
      for(int i=0 ; i<k ; i++) ENDFExtract(&fpin,mf,mt[i]);
    }
  }
  fpin.close();

  return(0);
}
