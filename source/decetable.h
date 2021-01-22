/******************************************************************************/
/**                                                                          **/
/**     DeCE TABLE                                                           **/
/**             functions for tabulating ENDF contents                       **/
/**                                                                          **/
/******************************************************************************/

/**************************************/
/*      decetable.cpp                 */
/**************************************/

void   DeceLibToTable    (ENDF *, ENDF *);
void   DeceTableMF1      (ENDF *);
void   DeceTableMF2      (ENDF *);
void   DeceTableMF3      (ENDF *);
void   DeceTableMF4      (ENDF *);
void   DeceTableMF5      (ENDF *);
void   DeceTableMF6      (ENDF *, ENDF *);
void   DeceTableMF7      (ENDF *);
void   DeceTableMF8      (ENDF *);
void   DeceTableMF9      (ENDF *);
void   DeceTableMF10     (ENDF *);
void   DeceTableMF12     (ENDF *);
void   DeceTableMF13     (ENDF *);
void   DeceTableMF14     (ENDF *);
void   DeceTableMF15     (ENDF *);
void   DeceTableMF32     (ENDF *);
void   DeceTableMF33     (ENDF *);
void   DeceTableMF34     (ENDF *);
void   DeceTableMF35     (ENDF *);


/**************************************/
/*      general output format         */
/**************************************/

static inline void outVal(int x)
{ cout << setw(10) << x << "    "; }

static inline void outVal(double x)
{ cout.setf(ios::scientific, ios::floatfield);
  if(x >= 0.0) cout << setprecision(7) << setw(14) << x;
  else         cout << setprecision(6) << setw(14) << x; }

static inline void outVal(int w, int p, double x)
{ cout.setf(ios::fixed, ios::floatfield);
  cout << setw(w) << setprecision(p) << x; }

