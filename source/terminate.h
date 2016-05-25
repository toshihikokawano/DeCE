/*
   terminate.h : 
        prototype of function to print out warning messages, 
        and to terminate code execution.
        the functions are defined in main.cpp
 */

/**************************************/
/*      main.cpp                      */
/**************************************/
void    WarningMessage     (std::string);
void    WarningMessage     (std::string, int);
void    WarningMessage     (std::string, double);
void    WarningMessage     (std::string, std::string);

int     TerminateCode      (std::string);
int     TerminateCode      (std::string, int);
int     TerminateCode      (std::string, double);
int     TerminateCode      (std::string, std::string);
