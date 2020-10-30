// Example for use of endflib.cpp and endfio.cpp

/*
   This example first scans the entire ENDF file 
   and store the index information in the ENDFDict object.
   Then, for each MF/MT, ENDFExtract copies from the original file to the standard output.
   Alternatively, each section is once copied into an ENDF object, and print.
   This is activated when ON_MEMORY is defined.
 */

#include <iostream>
#include <fstream>

using namespace std;

#include "../source/endflib.h"

#define ON_MEMORY

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
  
#ifdef ON_MEMORY
  /*** first, extract MF1 MT451, because DeCE does not read the comment section */
  ENDFExtract(&fpin,1,451);
#endif

  /*** for all MF numbers */
  for(int mf=1 ; mf <= 40 ; mf++){

    bool mfexist = false;

    /*** check if MT subsections are given in this MF */
    for(int i=0 ; i<dict.getSEC() ; i++){

      if(dict.mf[i] == mf){
#ifdef ON_MEMORY
        ENDF lib; // allocate a new object

        /*** read data in the source file, and write them */
        ENDFRead(&fpin,&lib,mf,dict.mt[i]);
        ENDFWrite(&lib);
#else
        /*** copy the section directly from source file */
        ENDFExtract(&fpin,mf,dict.mt[i]);
#endif
        mfexist = true;
      }
    }
    /*** write file end if this MF exists in the source file */
    if(mfexist) ENDFWriteFEND(dict.getMAT());
  }

  /*** write closing items */
  ENDFWriteFEND(0);    // MEND
  ENDFWriteFEND(-1);   // TEND

  fpin.close();

  return(0);
}
