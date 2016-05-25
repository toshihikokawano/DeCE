/******************************************************************************/
/**     Matrix Calculations                                                  **/
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <complex>

using namespace std;

#include "matrix.h"

static int MatrixInverseCholeski(const int, complex<double> *);
static int MatrixInverseCalc1   (const int, complex<double> *);
static int MatrixLUDecomposition(const int, complex<double> **);
static int MatrixInverseCalc2   (const int, complex<double> **, complex<double> *);


void MatrixPrint(const int m, complex<double> **u)
{
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      cout <<setw(12) << real(u[i][j]) << setw(12) << imag(u[i][j]);
    }
    cout <<endl;
  }
  cout <<endl;
}


void MatrixPrint(const int m, complex<double> *u)
{
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<=i ; j++){
      int k = i*(i+1)/2+j;
      cout <<setw(12) << real(u[k]) << setw(12) << imag(u[k]);
    }
    cout <<endl;
  }
  cout <<endl;
}


int MatrixInverse(const int m, complex<double> *a)
{
  int c=  MatrixInverseCholeski(m,a);
  if(c < 0){
    cerr << "pivot zero" << endl;
  }
  else{
    c = MatrixInverseCalc1(m,a);
    if(c == -1) cerr << "matrix not positiv" << endl;
  }
  return(c);
}


int MatrixInverse1(const int m, complex<double> **a, complex<double> **b)
{
  complex<double> *d;
  int c= 0;

  try{
    d = new complex<double> [m*(m+1)/2];
  }
  catch(bad_alloc){
    cerr << "memory allocation error";
    return(-1);
  }

  int k = 0;
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<=i ; j++){
      d[k++] = a[i][j];
//      cout << setw(12) << real(a[i][j]) << endl;
    }
  }


  c = MatrixInverseCholeski(m,d);

  if(c < 0){
    cerr << "pivot zero" << endl;
  }
  else{
    c = MatrixInverseCalc1(m,d);
    if(c == -1) cerr << "matrix not positiv" << endl;
  }

  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      if(j<=i) k = i*(i+1)/2+j;
      else     k = j*(j+1)/2+i;
      b[i][j] = d[k];
    }
  }

  delete [] d;
  return(c);
}


int MatrixInverse2(const int m, complex<double> **a, complex<double> **b)
{
  complex <double> *d;
  try{
    d = new complex<double> [m];
  }
  catch(bad_alloc){
    cerr << "memory allocation error";
    return(-1);
  }

  MatrixLUDecomposition(m,a);

  for(int j=0 ; j<m ; j++){
    for(int i=0 ; i<m ; i++){
      d[i] = (i == j) ? complex<double>(1.0,0.0) : complex<double>(0.0,0.0);
    }
    MatrixInverseCalc2(m,a,d);
    for(int i=0 ; i<m ; i++) b[i][j] = d[i];
  }

  delete [] d;
  return(0);
}



int MatrixInverseCholeski(const int n, complex<double> *a)
{
  const complex<double> eps(1e-72,0.0);
  complex<double> x1,x2,x3,x4;
  int c = 0;

  if(abs(a[0]) == 0.0) return(-2);
  a[0] = complex<double>(1.0,0.0)/a[0];
  if(n == 1) return(c);

  x1 = a[0]*a[1];
  x2 = a[1]*x1;
  a[1] = x1;
  x4 = a[2];
  x3 = x4-x2;
  if( abs(x3) < abs(x4*eps) ) return(-2);
  a[2] = complex<double>(1.0,0.0)/x3;
  if(n == 2) return(c);

  int n2 = n-2;
  int l=4;
  for(int i=1 ; i<=n2 ; i++){
    int jk=2;
    int ik=l;
    for(int j=1 ; j<=i ; j++){
      complex<double> sum1(0.0,0.0);
      int lj = l+j-1;
      for(int ij=l ; ij<=lj ; ij++){
        sum1 += a[ij-1]*a[jk-1]; 
        jk++;
      }
      jk++;
      ik++;
      a[ik-1] -= sum1;
    }

    int i1 = i+2;
    int jj = 1;
    int ij = l;

    complex<double> sum2(0.0,0.0);
    for(int j=2 ; j<=i1 ; j++){
      x1 = a[ij-1]*a[jj-1];
      sum2 += a[ij-1]*x1;
      a[ij-1] = x1;
      jj += j;
      ij++;
    }

    x4 = a[jj-1];
    x3 = x4-sum2;
    if( abs(x3) < abs(x4*eps) ) return(-2);
    a[jj-1] = complex<double>(1.0,0.0)/x3;
    l=jj+1;
  }

  return(c);
}


int MatrixInverseCalc1(const int n, complex<double> *a)
{
  int c = 0;

  if(n == 1) return(c);

  a[1] = -a[1];
  int n1 = n-1;
  if(n!=2){
    int ije= 2;
    for(int i=2 ; i<=n1 ; i++){
      int ij=ije+2;
      ije += i+1;
      int kjs=0;
      for(int j=2 ; j<=i ; j++){
        complex<double> sum1 = a[ij-1];
        ij++;
        int kjsp = j;
        kjs += kjsp;
        int kj = kjs;
        for(int ik=ij ; ik<=ije ; ik++){
          sum1 += a[ik-1]*a[kj-1];
          kj+=kjsp;
          kjsp++;
        }
        a[kj-1] = -sum1;
      }
      a[ije-1] = -a[ije-1];
    }
  }

  int ii= 1;
  int ki= 2;
  complex<double> sum2 = a[0];

  for(int k=2 ; k<=n ; k++){
    ii += k;
    complex<double> x1 = a[ki-1]*a[ii-1];
    sum2 += a[ki-1]*x1;
    a[ki-1] = x1;
    ki += k;
  }
  a[0] = sum2;
  if(n == 2) return(c);

  int ij=1;
  for(int i=2 ; i<=n1 ; i++){
    int i1 = i+1;
    for(int j=2 ; j<=i ; j++){
      ij++;
      int kj=ij;
      complex<double> sum3 = a[kj-1];
      int ijsp=i1-j;
      for(int k=i ; k<=n1 ; k++){
        kj += k;
        ki = kj+ijsp;
        sum3 += a[kj-1]*a[ki-1];
      }
      a[ij-1] = sum3;
    }
    ij++;

    int ki = ij+i;
    int ii = ij;
    complex<double> sum4 = a[ii-1];

    for(int k=i1 ; k<=n ; k++){
      ii += k;
      complex<double> x1 = a[ki-1]*a[ii-1];
      sum4 += a[ki-1]*x1;
      a[ki-1] = x1;
      ki += k;
    }
    a[ij-1] = sum4;
   }

   return(c);
}


int MatrixLUDecomposition(const int n, complex<double> **a)
{
  for(int j=0 ; j<n ; j++){

    for(int i=0 ; i<j ; i++){
      complex<double> s = a[i][j];
      for(int k=0 ; k<i ; k++) s -= a[i][k]*a[k][j];
      a[i][j] = s;
    }

    for(int i=j ; i<n ; i++){
      complex<double> s = a[i][j];
      for(int k=0 ; k<j ; k++) s -= a[i][k]*a[k][j];
      a[i][j] = s;
    }

    if(abs(a[j][j]) == 0.0) a[j][j] = complex<double>(1.0e-32,0.0);
    if(j != n-1){
      complex<double> d = complex<double>(1.0,0.0)/a[j][j];
      for(int i=j+1 ; i<n ; i++) a[i][j] *= d;
    }
  }
  return(0);
}


int MatrixInverseCalc2(const int n, complex<double> **a,  complex<double> *b)
{
  for(int i=0 ; i<n ; i++){
    complex<double>sum = b[i];
    for(int j=0 ; j<i ; j++) sum -= a[i][j]*b[j];
    b[i] = sum;
  }
  for(int i=n-1 ; i>=0 ; i--){
    complex<double>sum = b[i];
    for(int j=i+1 ; j<n ; j++) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
  return(0);
}

