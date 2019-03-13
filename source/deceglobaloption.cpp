/******************************************************************************/
/**     DeCE GLOBAL OPTION                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "dece.h"
#include "global.h"
#include "terminate.h"

static void printOption(void);

/**********************************************************/
/*      set / unset global options                        */
/**********************************************************/
void DeceGlobalOption(string ope, string option, const double x)
{
  if(ope == "showoptions") printOption();
  else{
    if(option == "LineNumber"){
      opt.LineNumber = (ope == "set") ? true : false;
      ENDFPrintLineNumber(opt.LineNumber);
    }
    else if(option == "EnergyConversion"){
      if(ope == "set") opt.EnergyConversion = x;
      else WarningMessage("cannot unset option",option);
    }
    else if(option == "CrossSectionConversion"){
      if(ope == "set") opt.CrossSectionConversion = x;
      else WarningMessage("cannot unset option",option);
    }
    else WarningMessage("option not found",option);
  }
}


/**********************************************************/
/*      print global options                              */
/**********************************************************/
void printOption()
{
  cout << "option: LineNumber             " << ((opt.LineNumber) ? " ON" : "OFF") << endl;
  cout << "option: EnergyConversion       " << setw(13) << opt.EnergyConversion << endl;
  cout << "option: CrossSectionConversion " << setw(13) << opt.CrossSectionConversion << endl;
}
