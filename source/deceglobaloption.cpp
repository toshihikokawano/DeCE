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

static void optionPrint (void);
static void optionToggle (string, string);
static void optionSet (string, string, const double);
static void optionSet (string, string, string);

/**********************************************************/
/*      set / unset global options                        */
/**********************************************************/
void DeceGlobalOption(string ope, string option, string value, const double x)
{
  /*** print current setting */
  if(ope == "showoptions") optionPrint();

  /*** toggle options */
  else if(option == "LineNumber") optionToggle(ope,option);

  /*** set text */
  else if(option == "OutPut") optionSet(ope,option,value);

  /*** set a value */
  else optionSet(ope,option,x);
}


/**********************************************************/
/*      Check Data Reading Range                          */
/**********************************************************/
bool DeceCheckReadRange(const double x)
{
  bool outofrange = false;

  /*** skip data if range is set by options */
  if(opt.ReadRangeMin > 0.0){
    if(x < opt.ReadRangeMin) outofrange = true;
  }
  if(opt.ReadRangeMax > 0.0){
    if(x > opt.ReadRangeMax) outofrange = true;
  }
  return outofrange;
}


/**********************************************************/
/*      Check Data Editing Range                          */
/**********************************************************/
bool DeceCheckEditRange(const double x)
{
  bool outofrange = false;

  /*** skip data if range is set by options */
  if(opt.EditRangeMin > 0.0){
    if(x < opt.EditRangeMin) outofrange = true;
  }
  if(opt.EditRangeMax > 0.0){
    if(x > opt.EditRangeMax) outofrange = true;
  }
  return outofrange;
}


/**********************************************************/
/*      global options, set or unset                      */
/**********************************************************/
void optionToggle(string ope, string option)
{
  if(ope == "set"){
    if(option == "LineNumber"){
      opt.LineNumber = true;
      ENDFPrintLineNumber(opt.LineNumber);
    }
  }
  else if(ope == "unset"){
    if(option == "LineNumber"){
      opt.LineNumber = false;
      ENDFPrintLineNumber(opt.LineNumber);
    }
  }

  if(option == "LineNumber"){
    message << "option" << option << " line number is ";
    if(opt.LineNumber) message << "on";
    else message << "off";
    Notice("DeceGlobalOption:optionToggle");
  }
}


/**********************************************************/
/*      global options, value setting                     */
/**********************************************************/
void optionSet(string ope, string option, string value)
{
  if(ope == "set"){
    if(option == "OutPut"){
      opt.OutPut = value;
      message << "option " << option << " set to " << value;
    }
  }
  else if(ope == "unset"){
    if(option == "OutPut"){
      opt.OutPut = "";
      message << "option " << option << " unset";
    }
  }

  Notice("optionSet");
}

void optionSet(string ope, string option, const double x)
{
  if(ope == "unset"){
    message << "option [ " << option << " ] cannnot unset: ";
    WarningMessage();
    return;
  }

  if(     option == "ReadXdataConversion")  opt.ReadXdataConversion  = x;
  else if(option == "ReadYdataConversion")  opt.ReadYdataConversion  = x;
  else if(option == "WriteXdataConversion") opt.WriteXdataConversion = x;
  else if(option == "WriteYdataConversion") opt.WriteYdataConversion = x;
  else if(option == "ReadRangeMin")         opt.ReadRangeMin         = x;
  else if(option == "ReadRangeMax")         opt.ReadRangeMax         = x;
  else if(option == "EditRangeMin")         opt.EditRangeMin         = x;
  else if(option == "EditRangeMax")         opt.EditRangeMax         = x;
  else if(option == "AngleStep"){
    if((0.0 <= x) && (x < 180.0) ) opt.AngleStep = x;
    else{
      message << "option AngleStep = " << x << " out of range" ;
      WarningMessage();
    }
  }
  else{
    message << "option [ " << option << " ] not found" ;
    WarningMessage();
    return;
  }

  message << "option " << option << " set to " << x;
  Notice("optionSet");
}


/**********************************************************/
/*      print global options                              */
/**********************************************************/
void optionPrint()
{
  cout << "option: LineNumber           " << setw(13) << ((opt.LineNumber) ? " ON" : "OFF") << endl;
  cout << "option: AngleStep            " << setw(13) << opt.AngleStep << endl;
  cout << "option: ReadXdataConversion  " << setw(13) << opt.ReadXdataConversion << endl;
  cout << "option: ReadYdataConversion  " << setw(13) << opt.ReadYdataConversion << endl;
  cout << "option: WriteXdataConversion " << setw(13) << opt.WriteXdataConversion << endl;
  cout << "option: WriteYdataConversion " << setw(13) << opt.WriteYdataConversion << endl;
  cout << "option: ReadRangeMin         " << setw(13) << opt.ReadRangeMin << endl;
  cout << "option: ReadRangeMax         " << setw(13) << opt.ReadRangeMax << endl;
  cout << "option: EditRangeMin         " << setw(13) << opt.EditRangeMin << endl;
  cout << "option: EditRangeMax         " << setw(13) << opt.EditRangeMax << endl;

  cout << "option: OutPut               ";
  if(opt.OutPut.length() == 0) cout << setw(13) << "-none-" << endl;
  else                         cout << setw(13) << opt.OutPut << endl;
}


