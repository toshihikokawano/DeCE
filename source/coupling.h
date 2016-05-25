/*
   coupling.h : 
        prototype of functions for angular momentum coupling
 */

const int    MAX_FACTORIAL = 200;        // maximal factorial n!  (2 x Lmax)
const double ARRAY_OVER    = 1.0e+300;   // force overflow

/**************************************/
/*      coupling.c                    */
/**************************************/
void    factorial             (int);

double  triangle              (int, int, int);
double  clebsh_gordan         (int, int, int, int, int);
double  wigner_3j             (int, int, int, int, int, int);
double  wigner_6j             (int, int, int, int, int, int);
double  wigner_9j             (int, int, int, int, int, int, int, int, int);
double  racah                 (int, int, int, int, int, int);
double  z_coefficient         (int, int, int, int, int, int);
double  zbar_coefficient      (int, int, int, int, int, int);
double  reduced_matrix_element(int, int, int, int, int, int, int);
