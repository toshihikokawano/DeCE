// Example for use of endflib.cpp and endfio.cpp

/*
   This example first scans the entire ENDF file 
   and store the index information in the ENDFDict object.
   Then, for each MF/MT, ENDFExtract copies from the original file to the standard output.
 */

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

int main(int, char *[]);

int main(int argc, char *argv[])
{
  if(argc <= 1){
    cerr << "ENDF file not supplied" << endl;
    exit(-1);
  }

  /*** scan ENDF file and store data in dict */
  ENDFDict dict;   // allocate ENDF Dictionary

  if(ENDFScanLibrary((string)argv[1],&dict) < 0){
    cerr << "ENDF file cannot open " << argv[1] << endl;
    exit(-1);
  }

  /*** read original ENDF file section by section, and print each of them */
  ifstream fpin;
  fpin.open(argv[1]);

  /*** write Tape ID */
  ENDFWriteTPID(&dict);
  
  /*** for all MF numbers */
  for(int mf=1 ; mf <= 40 ; mf++){

    /*** check if MT subsections are given in this MF */
    for(int i=0 ; i<dict.getSEC() ; i++){
      if(dict.mf[i] == mf){
        /*** copy the section */
        ENDFExtract(&fpin,mf,dict.mt[i]);
      }
    }
  }

  /*** write closing items */
  ENDFWriteFEND(0);
  ENDFWriteFEND(-1);

  fpin.close();

  return(0);
}
