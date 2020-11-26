/******************************************************************************/
/**     DeCE HEADER                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"
#include "terminate.h"
#include "masstable.h"
#include "constant.h"
#include "endftext.h"


static void DeceHeaderCopyData (char *[]);
static void DeceHeaderReplaceData (ENDFDict *);

static ENDFText zsymam(11,0,0), alab(11,11,0), edate(11,22,0), auth(33,33,0);
static ENDFText refer(22,0,1), ddate(11,22,1), rdate(11,33,1), endate(11,55,1);
static ENDFText libname(18,4,2);
static ENDFText sublib(TEXT_WIDTH-5,5,3);
static ENDFText format(TEXT_WIDTH-6,6,4);


static bool datacopied = false;


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

  cout << "[0---+----1----+----2----+----3----+----4----+----5----+----6----+-]" << endl;
  if(dict->getSTDHeader()){
    for(int i=0 ; i<5 ; i++){
      cout << "[" << dict->text[i] << "]" << endl;
    }
  }
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
/*      Show Header Text                                  */
/**********************************************************/
void DeceShowHeaderText(ENDFDict *dict)
{
  if(!dict->getSTDHeader()){
    message << "head text part is not in a standard form";
    WarningMessage();
    return;
  }

  DeceHeaderCopyData(dict->text);

  cout << "  ZSYMAM:"; zsymam.print();  cout << endl;
  cout << "    ALAB:"; alab.print();    cout << endl;
  cout << "    AUTH:"; auth.print();    cout << endl;

  cout << "   REFER:"; refer.print();   cout << endl;
  cout << "   EDATE:"; edate.print();   cout << endl;
  cout << "   DDATE:"; ddate.print();   cout << endl;
  cout << "   RDATE:"; rdate.print();   cout << endl;
  cout << "  ENDATE:"; endate.print();  cout << endl;

  cout << " LIBNAME:"; libname.print(); cout << endl;
  cout << "  SUBLIB:"; sublib.print();  cout << endl;
  cout << "  FORMAT:"; format.print();  cout << endl;
}


/**********************************************************/
/*      Edit Header Text                                  */
/**********************************************************/
void DeceEditHeaderText(ENDFDict *dict, string field, char *text)
{
  if(!dict->getSTDHeader()){
    message << "cannot edit head text part because is not in a standard form";
    WarningMessage();
    return;
  }

  if(!datacopied) DeceHeaderCopyData(dict->text);

  if(     field == "ZSYMAM" ){ zsymam.read( text); }
  else if(field == "ALAB"   ){ alab.read(   text); }
  else if(field == "AUTH"   ){ auth.read(   text); }

  else if(field == "REFER"  ){ refer.read(  text); }
  else if(field == "EDATE"  ){ edate.read(  text); }
  else if(field == "DDATE"  ){ ddate.read(  text); }
  else if(field == "RDATE"  ){ rdate.read(  text); }
  else if(field == "ENDATE" ){ endate.read( text); }

  else if(field == "LIBNAME"){ libname.read(text); }
  else if(field == "SUBLIB" ){ sublib.read( text); }
  else if(field == "FORMAT" ){ format.read( text); }
  else{
    message << "header text field [ " << field << " ] not defined";
    WarningMessage();
  }
 
  DeceHeaderReplaceData(dict);
}


/**********************************************************/
/*      Copy TEXT Filed Data into Objects                 */
/**********************************************************/
void DeceHeaderCopyData(char *line[])
{
  zsymam.copy( line[ zsymam.ypos]);
  alab.copy(   line[   alab.ypos]);
  edate.copy(  line[  edate.ypos]);
  auth.copy(   line[   auth.ypos]);

  refer.copy(  line[  refer.ypos]);
  ddate.copy(  line[  ddate.ypos]);
  rdate.copy(  line[  rdate.ypos]);
  endate.copy( line[ endate.ypos]);

  libname.copy(line[libname.ypos]);
  sublib.copy( line[ sublib.ypos]);
  format.copy( line[ format.ypos]);

  datacopied = true;
}


/**********************************************************/
/*      Replace TEXT Data in DICT Object by Given Data    */
/**********************************************************/
void DeceHeaderReplaceData(ENDFDict *dict)
{
  zsymam.paste( TEXT_WIDTH,dict->text[zsymam.ypos]);
  alab.paste(   TEXT_WIDTH,dict->text[  alab.ypos]);
  edate.paste(  TEXT_WIDTH,dict->text[ edate.ypos]);
  auth.paste(   TEXT_WIDTH,dict->text[  auth.ypos]);

  refer.paste(  TEXT_WIDTH,dict->text[ refer.ypos]);
  ddate.paste(  TEXT_WIDTH,dict->text[ ddate.ypos]);
  rdate.paste(  TEXT_WIDTH,dict->text[ rdate.ypos]);
  endate.paste( TEXT_WIDTH,dict->text[endate.ypos]);

  for(int i=44 ; i<55 ; i++) dict->text[1][i] = ' ';

  libname.paste(TEXT_WIDTH,dict->text[libname.ypos]);
  sublib.paste( TEXT_WIDTH,dict->text[ sublib.ypos]);
  format.paste( TEXT_WIDTH,dict->text[ format.ypos]);

  for(int i=0 ; i<4 ; i++) dict->text[2][i] = '-';
  for(int i=0 ; i<5 ; i++) dict->text[3][i] = '-';
  for(int i=0 ; i<6 ; i++) dict->text[4][i] = '-';

  sprintf(&dict->text[2][22],"MATERIAL "); dict->text[2][31] = ' ';
  sprintf(&dict->text[2][31],"%4d",dict->getMAT());
  for(int i=35 ; i<66 ; i++) dict->text[2][i] = ' ';

//  for(int i=0 ; i<5 ; i++){
//    dict->text[i][66] = '\0';
//    cout << dict->text[i] << "<" << endl;
//  }
}

