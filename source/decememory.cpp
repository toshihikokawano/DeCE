/******************************************************************************/
/**     DeCE MEMORY USAGE                                                    **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "dece.h"

static void memoryPrint(ENDF *);
static void memoryTotal(void);

static bool firstcall = true;
static unsigned long memsize = 0;
static unsigned long memused = 0;

/**********************************************************/
/*      Print Allocated Memory and Usage                  */
/**********************************************************/
void DeceMemoryUsage(ENDFDict *dict, ENDF *lib[])
{
  firstcall = true;

  /*** for all MF */
  for(int mf=1 ; mf<=40 ; mf++){
    /*** for all MT */
    for(int i=0 ; i<dict->getSEC() ; i++){
      if(dict->mf[i] == mf){
        int id = dict->getID(mf,dict->mt[i]);
        if(id < 0) continue;
        if(lib[id]->isalloc()) memoryPrint(lib[id]);
      }
    }
  }

  /*** print total memory */
  if(!firstcall) memoryTotal();
}


/**********************************************************/
/*      Print Memory Usage Table                          */
/**********************************************************/
void memoryPrint(ENDF *lib)
{
  if(firstcall){
    memsize = 0L;
    memused = 0L;

    cout << "  MF  MT   Record                  Integer                 Double                  Total" << endl;
    cout << "           Allocated        Used   Allocated        Used   Allocated        Used   Allocated        Used" << endl;
    cout << "--------------------------------------------------------------------------------------------------------" << endl;
    firstcall = false;
  }

  /*** allocated buffer size */
  unsigned long rsize = (unsigned long)lib->getRSIZE()   * sizeof(Record);
  unsigned long isize = (unsigned long)lib->getISIZE()   * sizeof(int);
  unsigned long xsize = (unsigned long)lib->getXSIZE()   * sizeof(double);
  unsigned long tsize = rsize + isize + xsize;

  /*** actual memory size */
  unsigned long rused = (unsigned long)(lib->getPOS()+1) * sizeof(Record);
  unsigned long iused = (unsigned long)lib->getNI()      * sizeof(int);
  unsigned long xused = (unsigned long)lib->getNX()      * sizeof(double);
  unsigned long tused = rused + iused + xused;

  /*** print table */
  cout << setw(4) << lib->getENDFmf();
  cout << setw(4) << lib->getENDFmt();

  cout << setw(12) << rsize << setw(12) << setw(12) << rused;
  cout << setw(12) << isize << setw(12) << setw(12) << iused;
  cout << setw(12) << xsize << setw(12) << setw(12) << xused;

  cout << setw(12) << tsize;
  cout << setw(12) << tused << endl;

  /*** add to total */
  memsize += tsize;
  memused += tused;
}


/**********************************************************/
/*      Print Total Memory Usage                          */
/**********************************************************/
void memoryTotal()
{
  cout << "--------------------------------------------------------------------------------------------------------" << endl;
  cout << "   Total                                                                        ";
  cout << setw(12) << memsize;
  cout << setw(12) << memused << endl;
}
