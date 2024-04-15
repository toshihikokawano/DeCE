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


/**************************************/
/*      decereadfile.cpp              */
/**************************************/
void   readSetCharged (void);
int    readCSdata  (char *, int, const int, double *, double *);
int    readISdata  (char *, int, const int, double *, double *, double *);
int    readNUdata  (char *, int,            double *, double *);
int    readMSdata  (char *, int,            double *, double *, int *, int *, double *);
int    geneCSdata  (int, double *, double *, double, double, double *);
int    mergeCSdata (int, double *, double *, double, double *, double *);
