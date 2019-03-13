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

  fpin.open(file);
  if(!fpin) TerminateCode("cannot open data file",file);

  ENDFRead(&fpin,lib,lib->getENDFmf(),lib->getENDFmt());
  fpin.close();

  /*** replace ZA and AWR from Dictionary */
  Record head = dict->head;
  double za   = head.c1;
  double awr  = head.c2;

  head = lib->getENDFhead();
  head.c1 = za;
  head.c2 = awr;

  lib->setENDFhead(head);
}

