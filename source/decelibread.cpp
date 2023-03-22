/******************************************************************************/
/**     DeCE LIBREAD                                                         **/
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "dece.h"
#include "terminate.h"


/**********************************************************/
/*      Copy One Section from a File                      */
/**********************************************************/
void DeceLibRead(ENDFDict *dict, ENDF *lib, char *file)
{
  ifstream fpin;
  int mf = lib->getENDFmf();
  int mt = lib->getENDFmt();

  fpin.open(file);
  if(!fpin){ message << "cannot open data file " << file; TerminateCode("DeceLibRead"); }

  int c = ENDFRead(&fpin,lib,lib->getENDFmf(),lib->getENDFmt());
  fpin.close();

  if(c < 0){
    message << "no section to be imported from " << file << " for MF/MT = " << mf << "/" << mt;
    WarningMessage();
  }
  else{
    /*** replace ZA and AWR from Dictionary */
    Record head = lib->getENDFhead();
    head.c1 = dict->getZA();
    head.c2 = dict->getAWR();

    lib->setENDFhead(head);

    message << "MF " << lib->getENDFmf() << " MT " << lib->getENDFmt() << " imported from " << file;
    Notice("DeceLibRead");
  }
}


/**********************************************************/
/*      Check If Section is Given in a File               */
/**********************************************************/
bool DeceLibScan(const int mf, const int mt, char *file)
{
  ifstream fpin;
  ENDF lib;

  fpin.open(file);
  if(!fpin){ message << "cannot open data file " << file; TerminateCode("DeceLibScan"); }
  int c = ENDFRead(&fpin,&lib,mf,mt);
  fpin.close();

  return ((c < 0) ? false : true);
}


