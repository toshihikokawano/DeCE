/*
   terminate.h : 
        prototype of function to print out warning messages, 
        and to terminate code execution.
        the functions are defined in main.cpp
 */

#include <ostream>

#ifndef DECE_TOPLEVEL
extern ostringstream message;
#endif

/**************************************/
/*      dece.cpp                      */
/**************************************/
void    WarningMessage     ();
void    Notice             (std::string);

int     TerminateCode      (std::string);
int     TerminateCode      (std::string, int);
int     TerminateCode      (std::string, double);
int     TerminateCode      (std::string, std::string);
