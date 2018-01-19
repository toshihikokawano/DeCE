/*
   dece.h : 
        prototype definition for DeCE functions
 */

#include <fstream>

#ifndef __ENDFLIB_H__
#define __ENDFLIB_H__
#include "endflib.h"
#endif

//#define ENERGY_UNIT_EV
//#define CROSS_SECTION_UNIT_B

#define ENERGY_UNIT_MEV
#define CROSS_SECTION_UNIT_MB

/*** decetable.cpp */
void   DeceTable        (ENDFDict *, ENDF *[], ifstream *, const int, const int, const int);
void   DeceFileToTable  (ifstream *, const int, const int, const int);
void   DeceDataPoint    (ifstream *, const int, const int, const double);

/*** deceangdist */
void   DeceAngdist      (ENDFDict *, ENDF **, const int, const int, char *, int);

/*** dececalc.cpp */
void   DeceCalc         (ENDFDict *, ENDF **, const int, const int, const int, const char);
void   DeceCalc452      (ENDFDict *dict, ENDF *lib[]);
void   DeceDuplicate    (ENDFDict *, ENDF **, const int, const int, const int);

/*** dececopy.cpp */
void   DeceExtract      (ENDFDict *, ENDF *[], ifstream *, const int, const int);
void   DeceLibRead      (ENDFDict *, ENDF *, char *);
void   DeceDelete       (ENDFDict *, const int, const int);
void   DeceScan         (ENDFDict *);

/*** deceheader.cpp */
void   DeceFixAWR       (ENDFDict *);

/*** deceoutput.cpp */
void   DeceOutput       (ifstream *, ENDFDict *, ENDF **);
void   DeceRenumber     (string, string, ENDFDict *);

/*** decepoint.cpp */
void   DecePoint        (ENDFDict *, ENDF **, const int, const int, double, double, string);
void   DeceFactor       (ENDFDict *, ENDF **, const int, const int, double, double, double, double);
void   DeceApplyFunc    (ENDFDict *, ENDF **, const int, const int, int, double, double, double);
void   DeceChangeInt    (ENDFDict *, ENDF **, const int, int, int, int);

/*** deceqvalue.coo */
void   DeceChangeQvalue (ENDFDict *, ENDF **, const int, double, double);
void   DeceCheckEnergy  (ENDFDict *, ENDF **, const bool);

/*** deceread.cpp */
void   DeceRead         (ENDFDict *, ENDF *, const int, const int, char *, int, bool);

/*** decemod6 */
void   DeceBoundCorrect (ENDFDict *, ENDF **, const int);
void   DeceDuplicatePoint (ENDFDict *, ENDF **, const int, double);
void   DeceGenProdCS    (ENDFDict *, ENDF **, const int, const int);

/*** gfr.cpp */
void   gfrScanThermal   (ifstream *, ENDFDict *);
double gfrGetOnePoint   (ifstream *, ENDFDict *, const double);
void   gfrPtCross       (ENDFDict *, ENDF **, double, double, double);
void   gfrAngDist       (ENDFDict *, ENDF **, double, double, double);
void   gfrAngDistSmooth (ENDFDict *, ENDF **, double);

/*** deceprocpoint.cpp */
void   DeceGeneratePointwise (ENDFDict *, ENDF **);

