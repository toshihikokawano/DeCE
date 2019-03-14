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

    if(ope == "set"){
      if(option == "LineNumber"){
        opt.LineNumber = true;
        ENDFPrintLineNumber(opt.LineNumber);
      }
      else if(option == "EnergyConversion"){
        opt.EnergyConversion = x;
      }
      else if(option == "CrossSectionConversion"){
        opt.CrossSectionConversion = x;
      }
      else if(option == "AngleStep"){
        opt.AngleStep = x;
      }
      else WarningMessage("option not found ",option);
    }
    else if(ope == "unset"){
      if(option == "LineNumber"){
        opt.LineNumber = false;
        ENDFPrintLineNumber(opt.LineNumber);
      }
      else WarningMessage("option not found ",option);
    }
  }
}


/**********************************************************/
/*      print global options                              */
/**********************************************************/
void printOption()
{
  cout << "option: LineNumber             " << setw(13) << ((opt.LineNumber) ? " ON" : "OFF") << endl;
  cout << "option: AngleStep              " << setw(13) << opt.AngleStep << endl;
  cout << "option: EnergyConversion       " << setw(13) << opt.EnergyConversion << endl;
  cout << "option: CrossSectionConversion " << setw(13) << opt.CrossSectionConversion << endl;
}
