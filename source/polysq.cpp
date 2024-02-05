/******************************************************************************/
/**     POLYSQ, Least-Squares Fitting by Polynomial or Legendre              **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "constant.h"
#include "terminate.h"
#include "polysq.h"
#include "decemisc.h"

static void   polyDesignMatrix (const int, const int, double *, double *);

static const double eps = 1.0e-03;

/**********************************************************/
/*      Least-Squares Fitting to Data                     */
/**********************************************************/
int polysq(const int n, const int m, double *xdata, double *ydata, double *a)
{
  double  *f, *v, *x;
  int     order = 0;
  
  v = new double [n*(n+1)/2];
  x = new double [m*(m+1)/2];
  f = new double [n*m];

  polyDesignMatrix(n,m,xdata,f);

  /*** optimize the highest order */
  for(int mopt=3 ; mopt<=m ; mopt+=2){

    /*** design matrix */
    polyDesignMatrix(n,mopt,xdata,f);

    /*** covariance matrix, just diagonal */
    for(int i=0 ; i<n ; i++){
      for(int j=0 ; j<=i ; j++)  v[i*(i+1)/2+j] = (i==j) ? 1.0 : 0.0;
    }

    if( (least_sq(n,mopt,ydata,a,v,x,f)) < 0.0 ){
      message << "least-squares equation not solved"; TerminateCode("polysq");
    }

    double chi2 = 0.0;
    for(int i=0 ; i<n ; i++){
      double xx = 0.0;
      for(int j=0 ; j<mopt ; j++) xx += f[i*mopt+j]*a[j];
      xx = ydata[i]-xx;
      chi2 += xx*xx;
    }

    if( chi2 <= eps ){
      order = mopt;
      break;
    }
  }

/*
  for(int i=0 ; i<180 ; i++){
    double z = a[0];
    for(int j=1 ; j<m ; j++){
      z += a[j]*legendre(j,(i+1.0));
    }
    cout << setw(12) << i+1 << setw(12) << z <<endl;
  }
*/
  delete [] v;
  delete [] x;
  delete [] f;

  return(order);
}


/**********************************************************/
/*      Least-Squares Fitting to Data                     */
/**********************************************************/
void polyDesignMatrix(const int n, const int m, double *x, double *f)
{
  for(int i=0 ; i<n ; i++){
    f[i*m] = 1.0;
    for(int j=1 ; j<m ; j++) f[i*m+j] = legendre(j,x[i]);
  }
}

