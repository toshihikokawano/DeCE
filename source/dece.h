/******************************************************************************/
/**                                                                          **/
/**     DeCE Main Header                                                     **/
/**             define DeCE functions                                        **/
/**                                                                          **/
/******************************************************************************/

#include <fstream>

#ifndef __ENDFLIB_H__
#define __ENDFLIB_H__
#include "endflib.h"
#endif

/*** dece.cpp */
void   DeceCreateLib    (ENDFDict *, int, int);
void   DeceCheckMT      (int);

/*** deceoperation.cpp */
void   DeceOperation    (ENDFDict *, ENDF **, ifstream *);

/*** dececalc.cpp */
void   DeceCalc         (ENDFDict *, ENDF **, const int, const int, const int, const char);
void   DeceCalc452      (ENDFDict *, ENDF **);

/*** decedelete.cpp */
void   DeceDelete       (ENDFDict *, const int, const int);

/*** deceread.cpp */
void   DeceRead         (ENDFDict *, ENDF *, const int, const int, char *, const int, const int);

/*** deceangdist.cpp */
void   DeceAngdist      (ENDFDict *, ENDF **, const int, const int, char *, int);

/*** decelibread.cpp */
void   DeceLibRead      (ENDFDict *, ENDF *, char *);

/*** decetable.cpp */
void   DeceTable        (ENDFDict *, ENDF *[], ifstream *, const int, const int);
void   DeceFileToTable  (ifstream *, const int, const int);
void   DeceDataPoint    (ifstream *, const int, const int, const double);

/*** deceextract.cpp */
void   DeceExtract      (ENDFDict *, ENDF *[], ifstream *, const int, const int);

/*** decepoint.cpp */
void   DecePoint        (ENDFDict *, ENDF **, const int, const int, double, double, string);

/*** decefactor.cpp */
void   DeceFactor       (ENDFDict *, ENDF **, const int, const int, double, double, double, double);

/*** deceapplyfunc.cpp */
void   DeceApplyFunc    (ENDFDict *, ENDF **, const int, const int, int, double, double, double);

/*** decereadjust.cpp */
void   DeceReadjust     (ENDFDict *, ENDF **, const int, const int);

/*** decechangeint.cpp */
void   DeceChangeInt    (ENDFDict *, ENDF **, const int, int, int, int);

/*** deceqvalue.cpp */
void   DeceChangeQvalue (ENDFDict *, ENDF **, const int, double, double);
void   DeceCheckThreshold (ENDFDict *, ENDF **, const bool);

/*** deceheader.cpp */
void   DeceShowHeaders  (ENDFDict *);
void   DeceEditHeader   (ENDFDict *, string, const double);
void   DeceShowHeaderText (ENDFDict *);
void   DeceEditHeaderText (ENDFDict *, string, char *);
void   DeceFixAWR       (ENDFDict *);

/*** decemod6 */
void   DeceBoundCorrect (ENDFDict *, ENDF **, const int);
void   DeceDuplicatePoint (ENDFDict *, ENDF **, const int, double);
void   DeceGenProdCS    (ENDFDict *, ENDF **, const int, const int);
void   DeceIsotropicAngularDistribution (ENDFDict *, ENDF **, const int);

/*** gfr.cpp */
void   gfrScanThermal   (ifstream *, ENDFDict *, double);
double gfrGetOnePoint   (ifstream *, ENDFDict *, const double);
void   gfrPtCross       (ENDFDict *, ENDF **, double, double, double);
void   gfrAngDist       (ENDFDict *, ENDF **, double, double, double, double);
void   gfrAngDistSmooth (ENDFDict *, ENDF **, double);
void   gfrSmatrixElement(ENDFDict *, ENDF **, double, double, double);

/*** deceprocpointwise.cpp */
void   DeceGeneratePointwise (ENDFDict *, ENDF **);

/*** deceprocgroup.cpp */
void   DeceGenerateGroup (ENDFDict *, ENDF **, const int,  const int, string);

/*** decescanindex.cpp */
void   DeceScanIndex    (ENDFDict *);

/*** deceglobaloption.cpp */
void   DeceGlobalOption (string, string, string, const double);
bool   DeceCheckReadRange (const double);
bool   DeceCheckEditRange (const double);

/*** decememory.cpp */
void   DeceMemoryUsage  (ENDFDict *, ENDF **);

/*** deceoutput.cpp */
void   DeceOutput       (ifstream *, ENDFDict *, ENDF **);
void   DeceRenumber     (string, string, ENDFDict *);
