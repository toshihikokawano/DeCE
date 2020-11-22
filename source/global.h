/******************************************************************************/
/**                                                                          **/
/**     Code Options                                                         **/
/**             define global variables, which can be changed                **/
/**             by the "set" command                                         **/
/**                                                                          **/
/******************************************************************************/

/**************************************/
/*      Class GlobalOption            */
/**************************************/
class GlobalOption{
 public:
  bool   LineNumber            ;    // print line numbers
  double AngleStep             ;    // calculate angular distribution at this step
  double ReadXdataConversion   ;    // conversion factor of energy when importing data
  double ReadYdataConversion   ;    // conversion factor of cross section
  double WriteXdataConversion  ;    // conversion factor when tabulating energy
  double WriteYdataConversion  ;    // conversion factor when tabulating cross section
  double ReadRangeMin          ;    // data reading range, low-side
  double ReadRangeMax          ;    // data reading range, high-side
  double EditRangeMin          ;    // data modification range, low-side
  double EditRangeMax          ;    // data modification range, high-side
  string Output                ;    // file name table/extract output will be written

  GlobalOption(){
    LineNumber           =  false;
    AngleStep            =  0.0;    // when zero, Legendre coefficients are printed
    ReadXdataConversion  =  1.0e+6; // default in MeV
    ReadYdataConversion  =  0.001;  // default in milli-barn
    WriteXdataConversion =  1.0;    // default in eV
    WriteYdataConversion =  1.0;    // default in barn
    ReadRangeMin         =  0.0;    // default, all data will be read
    ReadRangeMax         =  0.0;    // default, all data will be read
    EditRangeMin         =  0.0;    // default, all data will be modified
    EditRangeMax         =  0.0;    // default, all data will be modified
    Output               =  "";     // default output STDOUT
  }
};


#ifndef DECE_TOPLEVEL
extern GlobalOption opt;
#endif
