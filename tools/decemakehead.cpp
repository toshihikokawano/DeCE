/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate Header Part                                    **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <unistd.h>

using namespace std;

class MassExcess{
 public:
  unsigned int za;    // Z*1000 + A
  float        mass;  // mass excess
};

#include "../source/endflib.h"
#include "../source/constant.h"
#include "../source/masstable_audi2012_frdm2012.h"

int    main(int, char *[]);

int main(int argc, char *argv[])
{
  int anum = 0, znum = 0, mat = 0;
  string elem = "";

  if(argc <= 4){
    cerr << "usage: decemakehead Element Znumber Anumber MATnumber" << endl;
  }

  /*** command line options */
  elem = argv[1];
  znum = atoi(argv[2]);
  anum = atoi(argv[3]);
  mat  = atoi(argv[4]);

  if(znum < 2 || znum > 100){ cerr << "Z number " << znum << " out of range" << endl; exit(-1); }
  if(anum < 2 || anum > 300){ cerr << "A number " << anum << " out of range" << endl; exit(-1); }

  unsigned int zanum = znum * 1000 + anum;
  double mass  = 0.0;

  bool found = false;
  for(int i=0 ; i<nMassTable ; i++){
    if(MassTable[i].za == zanum){
      found = true;
      mass = MassTable[i].mass;
      break;
    }
  }
  if(!found){ cerr << " mass for Z:" << znum << " A:" << anum << " not found" << endl; exit(-1); }


  double za   = (double)zanum;
  double awr  = (mass / AMUNIT + anum) / MNEUTRON;
  int    lrp  = 1;
  int    lfi  = 0;
  int    nlib = 0;
  int    nmod = 1;
  Record head(za, awr, lrp, lfi, nlib, nmod);

  ENDFPrintLineNumber(false);
  cout << setw(70) << 1 << setw(2) << 0 << setw(3) << 0 << endl;
  ENDF lib;
  lib.setENDFmat(mat);
  lib.setENDFmf(1);
  lib.setENDFmt(451);
  lib.setENDFhead(head);

  ENDFWriteHEAD(&lib);

  double elis = 0.0;
  double sta  = 0.0;
  int lis  = 0;
  int liso = 0;
  int nfor = 6;
  Record cont(elis, sta, lis, liso, 0, nfor);
  ENDFPackCONT(cont,&lib);
  ENDFWriteCONT(&lib);

  double awi  = 1.0;
  double emax = 2.0e+7;
  int    lrel = 0;
  int    nsub = 10;
  int    nver = 1;
  cont.setRecord(awi, emax, lrel, 0, nsub, nver);
  ENDFPackCONT(cont,&lib);
  ENDFWriteCONT(&lib);

  double temp = 0.0;
  int    ldrv = 0;
  int    nwd  = 0;
  int    nxc  = 0;
  cont.setRecord(temp, 0.0, ldrv, 0, nwd, nxc);
  ENDFPackCONT(cont,&lib);
  ENDFWriteCONT(&lib);

  cout << setw(3) << znum << '-' << setw(2) << elem << '-' << setw(3) << anum << endl;
  cout << endl;
  cout << "----" << endl;
  cout << "-----INCIDENT NEUTRON DATA" << endl;
  cout << "------ENDF-6 FORMAT" << endl;

  ENDFWriteSEND(&lib);

  return(0);
}

