/*
   coupling.h : 
        prototype of functions for angular momentum coupling
 */

const int    MAX_FACTORIAL = 200;        // maximal factorial n!  (2 x Lmax)
const double ARRAY_OVER    = 1.0e+300;   // force overflow

/**************************************/
/*      coupling.c                    */
/**************************************/
void    factorial_allocate    (void);
void    factorial_delete      (void);

double  triangle              (const int, const int, const int);
double  clebsh_gordan         (const int, const int, const int, const int, const int);
double  wigner_3j             (const int, const int, const int, const int, const int, const int);
double  wigner_6j             (const int, const int, const int, const int, const int, const int);
double  wigner_9j             (const int, const int, const int, const int, const int, const int, const int, const int, const int);
double  racah                 (const int, const int, const int, const int, const int, const int);
double  z_coefficient         (const int, const int, const int, const int, const int, const int);
double  zbar_coefficient      (const int, const int, const int, const int, const int, const int);
double  reduced_matrix_element(const int, const int, const int, const int, const int, const int, const int);
