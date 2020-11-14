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
static const int    MAX_ITERATION      =  50000;
static const double TINY_NUMBER        =  1e-64;
static const double VALUE_CUTOFF       =  1e-24;

static void  omCoulombBarnett (const int, const double, const double, complex<double> *, complex<double> *);
static void  omCoulombFRatio  (const int, const double, const double, double *);
static void  omCoulombFrecur  (const int, const double, const double, double *, complex<double> *, complex<double> *);
static complex<double> omCoulombWRatio (const int, const double, const double);
static double omCoulombPowerSeries (const double, const double);
static void  omCoulombClosed  (const int, const double, const double, complex<double> *, complex<double> *);
static double omAsymptoticClosed (const double, const double);


/**********************************************************/
/*      Coulomb function                                  */
/*      for given Rho and Eta, up to Lmax                 */
/**********************************************************/
void coulomb(const int l, const double rho, const double eta, complex<double> *Coul0, complex<double> *Coul1)
{
  int lmax1 = MAX_L;
  if(l > lmax1) lmax1 = l;

  complex<double> *C0 = new complex<double>[lmax1+1]; // G + iF
  complex<double> *C1 = new complex<double>[lmax1+1]; // G'+ iF'

  /*** closed channel */
  if(rho < 0.0){
    omCoulombClosed(l, rho, eta, C0, C1);
  }
  /*** Barnett-Aoki algorithm */
  else{
    omCoulombBarnett(lmax1, rho, eta, C0, C1);
  }

  *Coul0 = C0[l];
  *Coul1 = C1[l];

  delete [] C0;
  delete [] C1;
}


void coulomb_function(const int lmax, const double rho, const double eta, complex<double> *Coul0, complex<double> *Coul1)
{
  int lmax1 = MAX_L;

  complex<double> *C0 = new complex<double>[lmax1+1]; // G + iF
  complex<double> *C1 = new complex<double>[lmax1+1]; // G'+ iF'

  /*** closed channel */
  if(rho < 0.0){
    omCoulombClosed(lmax1, rho, eta, C0, C1);
  }
  /*** Barnett-Aoki algorithm */
  else{
    omCoulombBarnett(lmax1, rho, eta, C0, C1);
  }

  /*** copy all calculated results */
  for(int l=0 ; l<=lmax ; l++){
    Coul0[l] = C0[l];
    Coul1[l] = C1[l];
  }
  for(int l=lmax+1 ; l<=lmax1 ; l++){
    Coul0[l] = complex<double>(0.0,0.0);
    Coul1[l] = complex<double>(0.0,0.0);
  }

  delete [] C0;
  delete [] C1;
}


/**********************************************************/
/*      Coulomb Phase Shift                               */
/**********************************************************/
double coulomb_phaseshift(const int l, const double eta)
{
  double w = 0.0;
  if(l > 0){
    for(int n=1 ; n<=l ; n++) w += atan(eta / (double)n);
  }
/*
  double y1 = eta;
  double y2 = y1*y1;
  double y3 = 16.0 + y2;
  double y4 = y3*y3;
  w = -y1+y1*log(y3)/2.+3.5*atan(y1/4.)-(atan(y1)+atan(y1/2.)+atan(y1/3.))
          -y1*(1.+(y2-48.)/(30.*y4)+(y2*y2-160.*y2+1280.)/(105.*y4*y4))
          /(12.*y3);
*/
  return w;
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

  /*** when we cannot calculate F0, default */
  if(f0 == 0.0){
    for(int l=0 ; l<=lmax ; l++){
      c0[l] = complex<double>( 1.0,0.0);
      c1[l] = complex<double>(-1.0,0.0);
    }
  }
  else{
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
    cerr << "continued fraction calculation for F(l+1)/F(l) didnt converge for eta = " << eta << " rho = " << rho << endl;
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
    cerr << "continued fraction calculation for Gp+iFp / G+iF didnt converge for eta = " << eta << " rho = " << rho << endl;
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

  if(f0 < VALUE_CUTOFF) return 0.0;

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
    cerr << "power seriese expansion didnt converge for eta = " << eta << " rho = " << rho << endl;
    exit(-1);
  }

  return f0;
}


/**********************************************************/
/*      Coulomb Function for Negative Energy              */
/**********************************************************/
void omCoulombClosed(const int lmax, const double rhox, const double eta, complex<double> *C0, complex<double> *C1)
{
  double p1, p2, p3;
  double rho = abs(rhox);

  if(eta == 0.0){

    /*** recurrence form, f[n+1] = [2n+1]/x f[n] + f[n-1] */
    p1 = 1.0/rho;
    p2 = p1+1.0/(rho*rho);
    C0[0] = complex<double>(p1,0.0);
    C0[1] = complex<double>(p2,0.0);

    for(int l=1 ; l<lmax ; l++){
      p3 = (2.0*l+1.0)/rho*p2 + p1;
      C0[l+1] = complex<double>(p3,0.0);
      p1 = p2;  p2 = p3;
    }

    /*** convert into k[n], by multiplying exp(-x) */
    double e = exp(-rho);
    for(int l=0 ; l<=lmax ; l++){
      C0[l] *= e;
      C0[l].imag(0.0);
    }

    /*** derivative, (xk[n])' = -n k[n] -x k[n-1] */
    C1[0].real(-exp(-rho));
    C1[0].imag(0.0);

    for(int l=1 ; l<=lmax ; l++){
      double wr = -l*C0[l].real() - rho*C0[l-1].real();
      C1[l] = complex<double>(wr,0.0);
    }
    /*** covert k[n] into x k[n] */
    for(int l=0 ; l<=lmax ; l++){
      double wr = rho*C0[l].real();
      C0[l].real(wr);
    }
  }
  else{
    /*** closed Coulomb function, F' at rho */
    double f1 = omAsymptoticClosed(rho, eta);

    int lp = lmax + 24 + (int)(5.0*eta);
    double x = 1.0;
    for(int l=lp ; l>=0 ; l--){
      double a = eta/(l+1.0);
      double b = a + (l+1.0)/rho;
      x  = -(a*a - 1.0)/(b + x) + b;
      if(l <= lmax) C1[l] = complex<double>(x, 0.0);
    }

    for(int l=0 ; l<=lmax ; l++){
      double g1 = C1[l].real();
      double d  = 1.0/sqrt(abs(g1 - f1));
      C0[l] = complex<double>(d,d);
      C1[l] = complex<double>(g1*d,f1*d);

      double a = eta/(l+1.0);
      double b = a + (l+1.0)/rho;
      f1 = (a*a - 1.0)/(b - f1) - b;
    }
  }
}


/***********************************************************/
/*      J. Raynal                                          */
/*      Closed Channel Coulomb Function, F'                */
/***********************************************************/
double omAsymptoticClosed(const double rm, const double coulomb)
{
  const int max_itr = 1000;
  const double eps = 1e-16;

  double f1 = 0.0;

  if((coulomb + 1.0)*rm > 8.0){

    /*** long range integration */
    if(rm < 10.0*(coulomb+1.0)){
      double s[7];
      double h = 0.25*rm;  if(h > 0.001953125) h = 0.001953125;
      double h12 = h*h/12.0;
      
      int n = 1 + 10.0/h;
      double v1 = 0.0;
      double v2 = h12*(1.0 + 2.0*coulomb/(rm + h* n   ));
      double v3 = h12*(1.0 + 2.0*coulomb/(rm + h*(n-1)));

      for(int j=0 ; j<=4 ; j++) s[j] = 0.0;
      s[5] = exp(-h);
      s[6] = 1.0;

      for(int i=n-2 ; i>=-3 ; i--){
        for(int j=0 ; j<6 ; j++) s[j] = s[j+1]/s[6];
        v1 = v2;
        v2 = v3;
        v3 = h12*(1.0 + 2.0*coulomb/(rm + h*i));
        s[6] = (s[5]*(2.0 + 10.0*v2) - s[4]*(1.0 - v1))/(1.0-v3);
      }
      f1 = ((s[0] - s[6])/60.0 + 0.15*(s[5] - s[1]) + 0.75*(s[2] - s[4]))/(h*s[3]);
    }

    /*** asymptotic expansion */
    else{
      double a = 1.0, b = 1.0, c = 0.0;
      for(int m=1 ; m<=26 ; m++){
        a = -a*0.5*(coulomb + m - 1.0)*(coulomb + m)/(rm*m);
        b += a;
        c -= a*m/rm;
      }
      f1 = c/b - 1.0 - coulomb/rm;
    }
  }

  /*** series expansion */
  else{
    double u0 = 0.0;
    double u1 = rm;
    double v0 = 1.0;
    double v1 = 0.0;

    double u  = u0 + u1;
    double v  = v0 + v1;
    double up = 1.0;
    double vp = 0.0;

    for(int n=2 ; n <max_itr ; n++){
      double cn = n*(n-1.0);
      double u2 = rm*(2.0*coulomb*u1 + rm*u0)/cn;
      double v2 = rm*(2.0*coulomb*v1 + rm*v0)/cn - 2.0*coulomb*(2.0*n-1.0)*u2/cn;
      u  += u2;
      v  += v2;
      up += u2*n/rm;
      vp += v2*n/rm;
      if(abs(u2/u) > eps){
        u0 = u1; u1 = u2;
        v0 = v1; v1 = v2;
      }
      else{
        if(abs(v2/v) < eps) break;
      }
    }

    double psr = 0.0;
    if(coulomb >= 1.0e-8){
      int    k  = (coulomb <= 7.5) ? (int)(8.5 - coulomb) : 0;
      double x  = 1.0 + coulomb + k;
      double uu = 1.0/(x*x);

      psr = log(x) - 0.5/x
          - uu/12.0 + uu*uu/120.0 - pow(uu,3.0)/252.0 + pow(uu,4.0)/240.0
          - pow(uu,5.0)/132.0 + pow(uu,6.0) * 691.0/32760.0;
      if(k != 0){
        for(int i=1 ; i<=k ; i++) psr -= 1.0/(x-(double)i);
      }
      psr += -0.5/coulomb + 2.0*EULER - 1.0;
    }
    else psr = EULER - 1.0;

    double ce = 2.0*coulomb*(psr + log(2.0*rm));
    f1 = (vp + up*ce + 2.0*coulomb*u/rm)/(v + u*ce);
  }

  return(f1);
}

