/*
   global.h : 
        global variables, can be changed by the "set" command
 */

/**************************************/
/*      Code Options                  */
/**************************************/
class GlobalOption{
 public:
  bool   LineNumber              ;    // print line numbers
  double EnergyConversion        ;    // conversion factor of energy when imported
  double CrossSectionConversion  ;    // conversion factor of cross section

  GlobalOption(){
    LineNumber             =  false;
    EnergyConversion       =  1.0e+6; // default in MeV
    CrossSectionConversion =  0.001;  // default in mb
  }
};


#ifndef DECE_TOPLEVEL
extern GlobalOption opt;
#endif
