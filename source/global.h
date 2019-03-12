/*
   global.h : 
        global variables, can be changed by the "set" command
 */

/**************************************/
/*      Code Options                  */
/**************************************/
class GlobalOption{
 public:
  bool linenumber         ;    // print line numbers

  GlobalOption(){
    linenumber          = false;
  }
};


#ifndef DECE_TOPLEVEL
extern GlobalOption opt;
#endif
