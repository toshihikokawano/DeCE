/******************************************************************************/
/**     DeCE HEADER                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"
#include "masstable.h"
#include "constant.h"
#include "endftext.h"

static ENDFText zsymam(11,0,0), alab(11,11,0), edate(11,22,0), auth(33,33,0);
static ENDFText refer(22,0,1), ddate(11,22,1), rdate(11,33,1), endate(11,55,1);
static ENDFText libname(18,4,2);
static ENDFText sublib(TEXT_WIDTH-5,5,3);
static ENDFText format(TEXT_WIDTH-6,6,4);



/**********************************************************/
/*      Print Header Data                                 */
/**********************************************************/
void DeceShowHeaders(ENDFDict *dict)
{
  cout << "header: ZA    " << setw(13) << dict->getZA() << endl;
  cout << "header: AWR   " << setw(13) << dict->getAWR() << endl;
  cout << "header: AWI   " << setw(13) << dict->getAWI() << endl;
  cout << "header: NLIB  " << setw(13) << dict->getNLIB() << "  library identifier" << endl;
  cout << "header: LREL  " << setw(13) << dict->getLREL() << "  release number" << endl;
  cout << "header: NVER  " << setw(13) << dict->getNVER() << "  version number" << endl;
  cout << "header: NSUB  " << setw(13) << dict->getNSUB() << "  sublibrary number" << endl;
  cout << "header: NMOD  " << setw(13) << dict->getNMOD() << "  modification number" << endl;
  cout << "header: LDRV  " << setw(13) << dict->getLDRV() << "  0: primary evaluation, 1: derived" << endl;

  cout << "header: NFOR  " << setw(13) << dict->getNFOR() << "  library format" << endl;
  cout << "header: LRP   " << setw(13) << dict->getLRP() << "  -1: no MF2, 0: no resonances, 1: resolved and/or unresolved, 2: PENDF" << endl;
  cout << "header: LFI   " << setw(13) << dict->getLFI() << "  0: no fission, 1: fission" << endl;

  cout << "header: EMAX  " << setw(13) << dict->getEMAX() << "  upper limit energy" << endl;
  cout << "header: TEMP  " << setw(13) << dict->getTEMP() << "  temperature" << endl;
  cout << "header: ELIS  " << setw(13) << dict->getELIS() << "  excitation energy of target" << endl;
  cout << "header: STA   " << setw(13) << dict->getSTA() << "  0: stable, 1: unstable" << endl;
  cout << "header: LIS   " << setw(13) << dict->getLIS() << "  state number of target" << endl;
  cout << "header: LISO  " << setw(13) << dict->getLISO() << "  isomeric state number" << endl;
}


/**********************************************************/
/*      Edit Header Data                                  */
/**********************************************************/
void DeceEditHeader(ENDFDict *dict, string parameter, const double x)
{
  int c = (int)x;

  if(     parameter == "NLIB") dict->setNLIB(c);
  else if(parameter == "LREL") dict->setLREL(c);
  else if(parameter == "NVER") dict->setNVER(c);
  else if(parameter == "NSUB") dict->setNSUB(c);
  else if(parameter == "NMOD") dict->setNMOD(c);
  else if(parameter == "LDRV") dict->setLDRV(c);
  else if(parameter == "NFOR") dict->setNFOR(c);
  else if(parameter == "LRP" ) dict->setLRP(c);
  else if(parameter == "LFI" ) dict->setLFI(c);
  else if(parameter == "EMAX") dict->setEMAX(x);
  else if(parameter == "TEMP") dict->setTEMP(x);
  else if(parameter == "ELIS") dict->setELIS(x);
  else if(parameter == "STA" ) dict->setSTA(c);
  else if(parameter == "LIS" ) dict->setLIS(c);
  else if(parameter == "LISO") dict->setLISO(c);
  else if(parameter == "ZA" || parameter == "AWR" || parameter == "AWI"){
    message << "header parameter [ " << parameter << " ] not allowed to change";
    WarningMessage();
  }
  else{
    message << "header parameter [ " << parameter << " ] not defined";
    WarningMessage();
  }
}


/**********************************************************/
/*      Fix AWR in DICT                                   */
/**********************************************************/
void DeceFixAWR(ENDFDict *dict)
{
  /*** ZA and AWR from Dictionary */
  double za   = dict->getZA();    //  1000*Z + A number

  int z = za/1000;
  int a = za - z*1000;

  double mass = mass_excess(z,a);
  double awr  = (mass / AMUNIT + a) / MNEUTRON;

  dict->setAWR(awr);
}


