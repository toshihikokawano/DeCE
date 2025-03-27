/******************************************************************************/
/**                                                                          **/
/**     DeCE READ                                                            **/
/**             functions for reading data files                             **/
/**                                                                          **/
/******************************************************************************/

struct Qval
{
  double qm; // Q-value for mass difference
  double qi; // reaction Q-value
  double et; // threshold energy
};

struct Prod
{
  int      col; // data column ofset
  int       za; // product 1000*Z + A
  int      pid; // out-going particles
  int    nx[3]; // meta state state index
  double ex[3]; // excitation energy of first meta
  double th[3]; // half-life
};


/**************************************/
/*      decereadfile.cpp              */
/**************************************/
void   readSetCharged (void);
int    readCSdata (char *, int, const int, double *, double *);
int    readISdata (char *, int, const int, double *, double *, double *);
int    readNUDdata (char *, int, double *, double *, double **, const int);
int    readNUPdata (char *, int, double *, double *);
struct Prod readMShead (char *, const int, int);
int    readMSdata (char *, const int, const int, double *, double *);
int    geneCSdata  (int, double *, double *, double, double, double *);
int    mergeCSdata (int, double *, double *, double, double *, const int, double *);
