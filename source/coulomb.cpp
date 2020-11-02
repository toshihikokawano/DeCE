/******************************************************************************/
/**                                                                          **/
/**     COULOMB :   Coulomb Function, G+iF, G'+iF'                           **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

using namespace std;

#include "constant.h"
#include "coulomb.h"

static const int    MAX_L              =     60;
static const double CRITERIA_ITERATION =  1e-24;
static const int    MAX_ITERATION      =  10000;
static const double TINY_NUMBER        =  1e-64;

static void  omCoulombBarnett (const int, const double, const double, complex<double> *, complex<double> *);
static void  omCoulombFRatio  (const int, const double, const double, double *);
static void  omCoulombFrecur  (const int, const double, const double, double *, complex<double> *, complex<double> *);
static complex<double> omCoulombWRatio (const int, const double, const double);
static double omCoulombPowerSeries (const double, const double);


/**********************************************************/
/*      Coulomb function                                  */
/*      for given Rho and Eta, up to Lmax                 */
/**********************************************************/
int coulomb(const int lmax, const double rho, const double eta, complex<double> *Coul0, complex<double> *Coul1)
{
  int lmax1 = MAX_L;
  if(lmax > lmax1) lmax1 = lmax;

  complex<double> *C0 = new complex<double>[lmax1+1]; // G + iF
  complex<double> *C1 = new complex<double>[lmax1+1]; // G'+ iF'

  /*** Barnett-Aoki algorithm */
  omCoulombBarnett(lmax1, rho, eta, C0, C1);

  /*** copy calculated results */
  for(int l=0 ; l<=lmax ; l++){
    Coul0[l] = C0[l];
    Coul1[l] = C1[l];
  }

  delete [] C0;
  delete [] C1;

  return lmax1;
}


/**********************************************************/
/*      A.R. Barnett                                      */
/*      Comp. Phys. Comm. 21 (1981), 297 - 314            */
/*      modified by Y. Aoki (Tsukuba U)                   */
/**********************************************************/

/*** R_lambda in Eq.(23), but squared */
inline static double Rfunc(int k, double eta)
{ return(1.0 + eta*eta/((double)(k*k))); }

/*** S in Eq.(23) */
inline static double Sfunc(int k, double eta, double rho)
{ return((double)k / rho + eta / (double)k); }

/*** T in Eq. (28) */
inline static double Tfunc(int k, double eta, double rho)
{ return (2.0*k+1.0) * (1.0/rho + eta/(k*(k+1.0))); }



void omCoulombBarnett(int lmax, double rho, double eta, complex<double> *c0, complex<double> *c1)
{
  double fw[3];

  /*** F(lambda+1) / F(lambda) */
  omCoulombFRatio(lmax, eta, rho, fw);

  /*** Recurrance Formula for F(lambda)
       fw[0] = F1, fw[1] = F0, fw[2] = F'(0)/F(0) */
  omCoulombFrecur(lmax, eta, rho, fw, c0, c1);

  /*** ratio (G'+iF')/(G+iF) = p + iq for l=0 */
  complex<double> w = omCoulombWRatio(0, eta, rho);

  /*** F for l = 0 */
  /*** terning point, rho(TP) = eta + sqrt{eta^2 + lambda(lambda+1)} = 2eta */ 
  double rtp = 2*eta;

  /*** normalization constant, sqrt{q / [(f - q)^2 + q^2]} */
  double f0 = 0.0;
  if(rho >= rtp){
    f0 = sqrt(w.imag() /((fw[2] - w.real())*(fw[2] - w.real()) + w.imag()*w.imag()));
    if(fw[1] <  0.0) f0 *= -1.0;
  }
  else{
    f0 = omCoulombPowerSeries(eta, rho);
  }

  double f1 = f0 * fw[2];                                      // F(0)' = F0 * (F'/F)
  double g0 = (1 - w.imag() * f0 * f0) / (f1 - w.real() * f0); // G = (F' - pF)/q (avoid q=0 case)
  double g1 = w.real() * g0 - w.imag() * f0;                   // G' = pG - qF

  c0[0] = complex<double>(g0,f0);
  c1[0] = complex<double>(g1,f1);

  /*** normalization factor of F, Fnorm = F / F(calc) */
  double fnorm = f0 / fw[1];

  /*** upward recursion relation, Eq.(26) */
  for(int l=1 ; l<=lmax ; l++){
    double r = c0[l].real();
    double s = c1[l].real();

    g0 = (s * c0[l-1].real() - c1[l-1].real()) / r;
    g1 = r * c0[l-1].real() - g0 * s;
    f0 = fnorm * c0[l].imag();   // because F are not normalized
    f1 = r * c0[l-1].imag() - f0 * s;

    c0[l] = complex<double>(g0,f0);
    c1[l] = complex<double>(g1,f1);
  }

  return;
}


/**********************************************************/
/*      Ratio of F(lambda+1)/F(lambda) in Eq.(28)         */
/*      Return un-normalized F(0) and F(1) in fw[]        */
/**********************************************************/
void omCoulombFRatio(const int lambda, const double eta, const double rho, double *fw)
{
  int l = lambda + 1;
  double d0 = 1.0 / Tfunc(l, eta, rho);
  double h0 = Rfunc(l, eta) * d0;
  double f0 = h0;

  bool conv = false, neg = false;
  int itr = 0;
  do{
    l++;
    double t1 = Tfunc(l, eta, rho);
    double d1 = 1.0 / (t1 - d0 * Rfunc(l, eta));
    double h1 = (t1 * d1 - 1.0) * h0;
    double f1 = f0 + h1;

    if(abs(f0 / f1 - 1.0) < CRITERIA_ITERATION){
      /*** save sign */
      if(d1 < 0.0) neg = true;
      conv = true;
      break;
    }

    d0 = d1;
    h0 = h1;
    f0 = f1;

  }while(++itr < MAX_ITERATION);

  if(!conv){
    cerr << "continued fraction calculation for F(l+1)/F(l) didn't converge" << endl;
    exit(-1);
  }

  /*** In Eq.(28), lambda+1 term is not squared */
  f0 = f0 / sqrt(Rfunc(lambda+1, eta));

  /*** set F0 to be very small number, then F1 is from the ratio f0 */
  fw[0] = TINY_NUMBER * ((neg) ? -1.0 : 1.0 );
  fw[1] = fw[0] * f0;
}


/**********************************************************/
/*      Recurrance Formula for F(lambda) in Eq.(23)       */
/**********************************************************/
void omCoulombFrecur(const int lambda, const double eta, const double rho, double *fw, complex<double> *c0, complex<double> *c1)
{
  double f0 = fw[1];  // F(lambda+1)
  double f1 = fw[0];  // F(lambda)

  int l = lambda + 1;

  double r0 = sqrt(Rfunc(l,eta)); // R(lambda + 1)
  double s0 = Sfunc(l, eta, rho); // S(lambda + 1)

  for(int l = lambda ; l>=1 ; l--){

    double r1 = sqrt(Rfunc(l,eta));  // R(lambda)
    double s1 = Sfunc(l, eta, rho);  // S(lambda)
    
    /*** Eq.(23) is modified as
         F(2) = [F(1) (S(0) + S(1)) - F(0) R(0)] / R(1) */
    double f2 = (f1*(s0 + s1) - f0 * r0) / r1;

    /*** save F(lambda), R(lambda), and S(lambda) */
    c0[l] = complex<double>(r1,f1);
    c1[l] = complex<double>(s1,0.0);

    f0 = f1;  f1 = f2;
    r0 = r1;  s0 = s1;
  }

  fw[0] = f0;
  fw[1] = f1;
  /*** store F'(0)/F(0) = S(1) - R(1) F(1)/F(0) */
  fw[2] = s0 - r0 * f0 / f1;
}


/**********************************************************/
/*      Ratio (G'+iF')/(G+iF) = p + iq in Eq.(29)         */
/**********************************************************/
complex<double> omCoulombWRatio(const int lambda, const double eta, const double rho)
{
  complex<double>x0(-lambda,eta);     // (-lambda, eta)
  complex<double>x1(lambda+1.0,eta);  // (lmabda+1,eta)
  complex<double>x2(2*(rho-eta),2.0); // (2(rho-eta), 2)

  complex<double>d0 = 1.0/x2;
  complex<double>h0 = x0*x1/x2;
  complex<double>z0 = h0;

  bool conv = false;
  int itr = 0;
  do{
    x0 += complex<double>(1.0,0.0);
    x1 += complex<double>(1.0,0.0);
    x2 += complex<double>(0.0,2.0);

    complex<double>d1 = 1.0 / (x2 + x0*x1*d0);
    complex<double>h1 = (x2*d1 - 1.0)*h0;
    complex<double>z1 = z0 + h1;

    double c = abs(norm(z0)/norm(z1) - 1.0);
    if(c < CRITERIA_ITERATION){
      conv = true;
      break;
    }

    d0 = d1;
    h0 = h1;
    z0 = z1;

  }while(++itr < MAX_ITERATION);

  if(!conv){
    cerr << "continued fraction calculation for G'+iF'/G+iF didn't converge" << endl;
    exit(-1);
  }

  /*** multiply by i/rho and add the first term of (1 - eta/rho)i */
  z0 = complex<double>(0.0,1.0-eta/rho) + z0 * complex<double>(0.0,1.0/rho);

  return z0;
}


/**********************************************************/
/*      Power Series Expansion of Coulomb Function        */
/*      Y. Aoki note                                      */
/**********************************************************/
double omCoulombPowerSeries(const double eta, const double rho)
{
  double b1 = 2 * PI * eta;
  double b0 = ((b1 <= 40.0) ? sqrt(b1 / (exp(b1) - 1.0)) : sqrt(b1) * exp(-b1/2)) * rho;

  b1 = eta * rho * b0;

  double f0 = b0 + b1;

  int l = 2;
  bool conv = false;
  int itr = 0;
  do{
    double b2 = (2*eta*rho*b1 - rho*rho*b0)/(l*(l+1.0));
    double f1 = f0 + b2;

    double c = abs(f1/f0 - 1.0);
    if(c < CRITERIA_ITERATION){
      conv = true;
      break;
    }
    b0 = b1;
    b1 = b2;
    f0 = f1;
    l ++;

  }while(++itr < MAX_ITERATION);

  if(!conv){
    cerr << "power seriese expansion didn't converge" << endl;
    exit(-1);
  }

  return f0;
}
