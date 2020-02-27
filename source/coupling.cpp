#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

/*
 *  coupling.cpp - calculate coupling coefficients of angular momenta
 *                 this source is a part of the CoH code
 *  note that arguments of those functions must be doubled, namely 1/2 is 1, etc.
 *
 *  Functions:
 *  factorial()
 *
 *    This function must be called first, in the main routine, to set up
 *    the array of factorial (logarithm). The maximum size of the array is
 *    defined in coupling.h (currently 160). The array "fact" must be declared as 
 *    an external variable in the main.c. Here is an exaple:
 *
 *    double  *fact;
 *    int main(){
 *       factorial(MAX_FACTORIAL);
 *    }
 *
 *
 *  wigner_3j(j1,j2,j3,j4,j5,j6)
 *    Wigner's 3J symbol (similar to Clebsh-Gordan)
 *               = / j1 j2 j3 \
 *                 \ j4 j5 j6 /
 *
 *  wigner_6j(j1,j2,j3,j4,j5,j6)
 *    Wigner's 6J symbol (similar to Racah)
 *               = { j1 j2 j3 }
 *                 { j4 j5 j6 }
 *
 *  wigner_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
 *    Wigner's 9J symbol
 *                 / j1 j2 j3 \
 *               = | j4 j5 j6 |
 *                 \ j7 j8 j9 /
 *
 *  racah(j1,j2,j3,j4,j5,j6)
 *    Racah coefficient
 *               = W(j1 j2 j3 j4 ; j5 j6)
 *
 *  clebsh_gordan(j1,j2,m1,m2,j3)
 *    Clebsh-Gordan coefficient 
 *               = <j1,j2,m1,m2|j3,m1+m3>
 *
 *  z_coefficient(l1,j1,l2,j2,S,L)
 *    Biedenharn's Z-coefficientn coefficient
 *               =  Z(l1  j1  l2  j2 | S L )
 *
 *  reduced_matrix_element(L,S,J,l0,j0,l1,j1)
 *    Reduced Matrix Element for Tensor Operator
 *               = < l1j1 || T(YL,sigma_S)J || l0j0 >
 *
 */


#include "coupling.h"

extern  double  *fact;
static inline int parity(int x){
  return( (((x)/2)%2 == 0) ? 1 : -1 );
}

static inline bool halfint(int x){
  return( ((x)%2 != 0) ? true : false );
}

static inline int max3(int a, int b, int c){
  if(a < b) a = b;
  if(a < c) a = c;
  return(a);
}

static inline int max4(int a, int b, int c, int d){
  if(a < b) a = b;
  if(a < c) a = c;
  if(a < d) a = d;
  return(a);
}

static inline int min3(int a, int b, int c){
  if(a > b) a = b;
  if(a > c) a = c;
  return(a);
}

static inline double w6j0  (const int, int *);
static inline double w6j1  (           int *);
static inline double cg1   (const int, const int, const int);
static inline double cg2   (const int, const int, const int, const int, const int, const int, const int, const int);
static inline double cg3   (const int, const int, const int, const int, const int, const int);


/***********************************************************/
/*      Factorial Calc. and Store in fact[]                */
/***********************************************************/
void factorial_allocate()
{
  fact = new double [MAX_FACTORIAL];

  fact[0] = 0.0;
  for(int i=1 ; i<MAX_FACTORIAL ; i++){ fact[i] = fact[i-1] + log((double)i); }
}

void factorial_delete()
{
  delete [] fact;
}


/***********************************************************/
/*      Wigner's 3j Symbol                                 */
/***********************************************************/
double wigner_3j(const int j1, const int j2, const int j3, const int j4, const int j5, const int j6)
{
  double cg = 0.0;

  if( (j4+j5+j6) != 0 ) return(0.0);
  if((cg=clebsh_gordan(j1, j2, j4, j5, j3)) == 0.0) return(0.0);
  /*** Brink, page 136 */
  return( (((j1-j2-j6)%4==0) ?  1.0 : -1.0)*cg/sqrt(j3+1.0) );
}


/***********************************************************/
/*      Wigner's 6j Symbol                                 */
/***********************************************************/
double wigner_6j(const int j1,const int j2,const int j3,const int j4,const int j5,const int j6)
{
  int x[6];

  x[0]=j1; x[1]=j2; x[2]=j3; x[3]=j4; x[4]=j5; x[5]=j6;

  /*** When it has zero, use Brink P.142,
       W(abcd,0f) = (-)^{a+c-f} delta(a,b) delta(c,d)/sqrt(2a+1)/sqrt(2c+1) */
  for(int i=0 ; i<6 ; i++) if(x[i] == 0) return(w6j0(i,x));

  /*** general case */
  return(w6j1(x));
}


inline double w6j0(const int i, int *x)
{
  switch(i){
  case 0: if((x[1] != x[2]) || (x[4] != x[5])) return(0.0);
    x[5]=x[3]; x[0]=x[1]; x[3]=x[4];
    break;
  case 1: if((x[0] != x[2]) || (x[3] != x[5])) return(0.0);
    x[5]=x[4];
    break;
  case 2: if((x[0] != x[1]) || (x[3] != x[4])) return(0.0);
    break;
  case 3: if((x[1] != x[5]) || (x[2] != x[4])) return(0.0);
    x[5]=x[0]; x[0]=x[4]; x[3]=x[1];
    break;
  case 4: if((x[0] != x[5]) || (x[2] != x[3])) return(0.0);
    x[5]=x[1];
    break;
  case 5: if((x[0] != x[4]) || (x[1] != x[3])) return(0.0);
    x[5]=x[2];
    break;
  }

  if( x[5] > (x[0]+x[3]) || x[5] < abs(x[0]-x[3]) ) return(0.0);
  if( x[0] > MAX_FACTORIAL || x[3] > MAX_FACTORIAL){
    cerr << "factorial n! too large" << endl;
    return(ARRAY_OVER);
  }

  return(1.0/sqrt((x[0]+1)*(x[3]+1))*( ((x[0]+x[3]+x[5])/2)%2!=0 ? -1 : 1));
}


inline double w6j1(int *x)
{
  int k1,k2,l1,l2,l3,l4,n1,n2,n3,m1,m2,m3,x1,x2,x3,y[4];
  static int a[3][4]={{0,0,3,3},
                      {1,4,1,4},
                      {2,5,5,2}};

  double w6j = 0.0;

  for(int k=0 ; k<4 ; k++){
    x1 = x[(a[0][k])];
    x2 = x[(a[1][k])];
    x3 = x[(a[2][k])];

    int n = (x1+x2+x3)/2;
    if(n > MAX_FACTORIAL){
      cerr << "factorial n! too large" << endl;
      return(ARRAY_OVER);
    }
    else if(n < 0) return(0.0);

    if((n1=n-x3) < 0) return(0.0);
    if((n2=n-x2) < 0) return(0.0);
    if((n3=n-x1) < 0) return(0.0);

    y[k] = n+2;
    w6j += fact[n1]+fact[n2]+fact[n3]-fact[n+1];
  }

  n1 = (x[0]+x[1]+x[3]+x[4])/2;
  n2 = (x[0]+x[2]+x[3]+x[5])/2;
  n3 = (x[1]+x[2]+x[4]+x[5])/2;

  k1 = max4(y[0],y[1],y[2],y[3])-1;
  k2 = min3(n1,n2,n3)+1;

  l1 = k1-y[0]+1;  m1 = n1-k1+1;
  l2 = k1-y[1]+1;  m2 = n2-k1+1;
  l3 = k1-y[2]+1;  m3 = n3-k1+1;
  l4 = k1-y[3]+1;

  w6j = exp(0.5*w6j+fact[k1]-fact[l1]-fact[l2]-fact[l3]-fact[l4]
                   -fact[m1]-fact[m2]-fact[m3]) * ((k1%2)==0 ? -1:1);

  if(k1 != k2){
    double w = w6j;
    int k = k2-k1;
    m1 -= k-1; m2 -= k-1; m3 -= k-1;
    l1 += k  ; l2 += k  ; l3 += k  ; l4 += k;

    for(int i=0 ; i<k ; i++)
           w6j = w - w6j*((k2-i)*(m1+i)*(m2+i)*(m3+i))
                        /((l1-i)*(l2-i)*(l3-i)*(l4-i));
  }
  return(w6j);
}


/***********************************************************/
/*      Wigner's 9j Symbol                                 */
/***********************************************************/
double wigner_9j(const int j1,const int j2,const int j3,
                 const int j4,const int j5,const int j6,
                 const int j7,const int j8,const int j9)
{
  int i0 = max3(abs(j1-j9),abs(j2-j6),abs(j4-j8));
  int i1 = min3(   (j1+j9),   (j2+j6),   (j4+j8));

  double rac = 0.0;
  for(int i=i0 ; i<=i1 ; i+=2){ rac += racah(j1,j4,j9,j8,j7, i)
                                      *racah(j2,j5, i,j4,j8,j6)
                                      *racah(j9, i,j3,j2,j1,j6)*(i+1);
  }

  return( (( (int)((j1+j3+j5+j8)/2+j2+j4+j9)%2==0) ?  1.0 : -1.0)*rac );
}

/***********************************************************/
/*      Racah Coefficient                                  */
/***********************************************************/
double racah(const int a, const int b, const int c, const int d, const int e, const int f)
{
  if( (a+b+c+d+e+f) >= 2*MAX_FACTORIAL){
    cerr << "factorial n! too large" << endl;
    return(ARRAY_OVER);
  }

  int i1 = (a+b+e)/2;
  int i2 = (c+d+e)/2;
  int i3 = (a+c+f)/2;
  int i4 = (b+d+f)/2;
  int i5 = (a+b+c+d)/2;
  int i6 = (a+d+e+f)/2;
  int i7 = (b+c+e+f)/2;
  double tri = triangle(a,b,e) * triangle(c,d,e) * triangle(a,c,f) * triangle(b,d,f);
  if(tri ==0.0 ) return(0.0);

  int k0 = i1; if(i2 > k0) k0 = i2; if(i3 > k0) k0 = i3; if(i4 > k0) k0 = i4;
  int k1 = i5; if(i6 < k1) k1 = i6; if(i7 < k1) k1 = i7;

  double sig = 0.0;
  for(int k=k0 ; k<=k1 ; k++){
    sig+=( ((k+i5)%2==0) ?  1.0 : -1.0 )
             *exp( fact[k+1]-(fact[k-i1]+fact[k-i2]+fact[k-i3]+fact[k-i4]
                                        +fact[i5-k]+fact[i6-k]+fact[i7-k]) );
  }

  return( tri*sig );
/* Wigner-6j: */
//return( tri*sig * ((((a+b+c+d)/2)%2==0) ? 1: -1) );
}


double triangle(const int a, const int b, const int c)
{
  int j1,j2,j3,j4;

  if(c < abs(a-b) || c > (a+b)) return(0.0);
  if((j1 = ( a+b-c)/2) < 0) return(0.0);
  if((j2 = ( a-b+c)/2) < 0) return(0.0);
  if((j3 = (-a+b+c)/2) < 0) return(0.0);
      j4 = ( a+b+c)/2+1;

  return( exp(0.5*(fact[j1]+fact[j2]+fact[j3]-fact[j4])) );
}


/***********************************************************/
/*      Clebsch-Gordan Coefficients                        */
/*     ( J1  J2  M1  M2 | J3  M1+M2 )                      */
/***********************************************************/
double clebsh_gordan(const int j1, const int j2, const int m1, const int m2, const int j3)
{
  int x1,x2,x3,y1,y2,y3;
  double cg=0.0;

  if(j1 < 0 || j2 < 0 || j3 < 0) return(0.0);
  if(j1+j2+j3 > 2*MAX_FACTORIAL){
    cerr << "factorial n! too large" << endl;
    return(ARRAY_OVER);
  }

  /* check int/half-int of J and M */
  if(halfint(j1) != halfint(m1)) return(0.0);
  if(halfint(j2) != halfint(m2)) return(0.0);
  if(halfint(j1+j2) != halfint(j3)) return(0.0);

  int m3 = m1+m2;

  if((x1=(j1+m1)/2+1) <= 0) return(0.0);
  if((x2=(j2+m2)/2+1) <= 0) return(0.0);
  if((x3=(j3-m3)/2+1) <= 0) return(0.0);

  if((y1=x1-m1) <= 0) return(0.0);
  if((y2=x2-m2) <= 0) return(0.0);
  if((y3=x3+m3) <= 0) return(0.0);

  if(j3 == 0){
    if(j1 == j2)      cg=(1.0/sqrt((double)j1+1.0)*((y1%2==0 ) ? -1:1));
  }
  else if( (j1 == 0 || j2 == 0) ){
    if((j1+j2) == j3) cg=1.0;
  }
  else{
    if(     m3 == 0 && abs(m1) <= 1){
      if(m1 == 0)   cg=cg1(x1,x2,x3);
      else          cg=cg2(x1+y1-y2,x3-1,x1+x2-2,x1-y2,j1,j2,j3, m2);
    }
    else if(m2 == 0 && abs(m1) <= 1){
                    cg=cg2(x1-y2+y3,x2-1,x1+x3-2,x3-y1,j1,j3,j3, m1);
    }
    else if(m1 == 0 && abs(m3) <= 1){
                    cg=cg2(x1      ,x1-1,x2+x3-2,x2-y3,j2,j3,j3,-m3);
    }
    else            cg=cg3(x1,x2,x3,y1,y2,y3);
  }

  return(cg);
}


inline double cg1(const int x1, const int x2, const int x3)
{
  int p1 = x1+x2+x3-1; if((p1%2) != 0) return(0.0);
  int p2 = x1+x2-x3;
  int p3 =-x1+x2+x3;
  int p4 = x1-x2+x3;
  if(p2 <= 0 || p3 <= 0 || p4 <= 0) return(0.0);
  if(p1 >= MAX_FACTORIAL){
    cerr << "factorial n! too large" << endl;
    return(ARRAY_OVER);
  }

  int q1 = (p1+1)/2-1; p1--;
  int q2 = (p2+1)/2-1; p2--;
  int q3 = (p3+1)/2-1; p3--;
  int q4 = (p4+1)/2-1; p4--;

  double a = fact[q1]-( fact[q2]+fact[q3]+fact[q4])
                 +0.5*( fact[2*x3-1]-fact[2*x3-2]
                       +fact[p2]+fact[p3]+fact[p4]-fact[p1]);

  return( (((q1+x1-x2)%2==0) ? 1.0 : -1.0) * exp(a) );
}


inline double cg2(const int k, const int q0, const int z1, const int z2, const int w1, const int w2, const int w3, const int mm)
{
  int p1 = z1 + q0 +2;
  int p2 = z1 - q0 +1;
  int p3 = z2 + q0 +1;
  int p4 =-z2 + q0 +1;
  if(p2 <= 0 || p3 <= 0 || p4 <= 0) return(0.0);
  if(p1 >= MAX_FACTORIAL){
    cerr << "factorial n! too large" << endl;
    return(ARRAY_OVER);
  }

  int q1 = (p1+1)/2-1; p1--;
  int q2 = (p2+1)/2-1; p2--;
  int q3 = (p3+1)/2-1; p3--;
  int q4 = (p4+1)/2-1; p4--;

  double a = fact[q1]-(fact[q2]+fact[q3]+fact[q4])
                +0.5*( fact[w3+1]-fact[w3  ]
                      +fact[w1  ]-fact[w1+1]
                      +fact[w2  ]-fact[w2+1]
                      +fact[p2]+fact[p3]+fact[p4]-fact[p1] );

  return( (((q4+k+(mm>0)*(p1+2))%2==0) ? -1.0 : 1.0) * 2.0 * exp(a) );
}


inline double cg3(const int x1, const int x2, const int x3, const int y1, const int y2, const int y3)
{
  int z1,z2,z3;
  double a,cg;

  int nx = x1+x2+x3-1;
  if( (z1=nx-x1-y1) < 0 ) return(0.0);
  if( (z2=nx-x2-y2) < 0 ) return(0.0);
  if( (z3=nx-x3-y3) < 0 ) return(0.0);

  int k1 = x2-y3;
  int k2 = y1-x3;

  int q1 = max3(k1,k2, 0);
  int q2 = min3(y1,x2,z3+1)-1;
  int q3 = q1-k1;
  int q4 = q1-k2;

  int p1 = y1-q1-1;
  int p2 = x2-q1-1;
  int p3 = z3-q1;

  a = cg = exp(0.5*( fact[x3+y3-1]-fact[x3+y3-2]-fact[nx-1]
                    +fact[z1  ]+fact[z2  ]+fact[z3  ]
                    +fact[x1-1]+fact[x2-1]+fact[x3-1]
                    +fact[y1-1]+fact[y2-1]+fact[y3-1])
                    -fact[p1  ]-fact[p2  ]-fact[p3  ]
                    -fact[q1  ]-fact[q3  ]-fact[q4  ])*(((q1%2)==0) ? 1 : -1);

  if(q1 != q2){
    q3 = q2-k1;
    q4 = q2-k2;
    p1 = y1-q2;
    p2 = x2-q2;
    p3 = z3-q2+1;
    for(int i=0; i<(q2-q1); i++)
      cg=a-cg*((p1+i)*(p2+i)*(p3+i))/((q2-i)*(q3-i)*(q4-i));
  }

  return(cg);
}


/***********************************************************/
/*      Biedenharn's Z-coefficientn Coefficient            */
/*      Z( L1  J1  L2  J2 | S L )                          */
/*      Blatt, Biedenharn                                  */
/*      Rev. Mod. Phys. 24, 258 (1952), Eq.(4.3)           */
/***********************************************************/
double z_coefficient(const int l1, const int j1, const int l2, const int j2, const int s, const int ll)
{
  int p = (-l1+l2+ll)/2;
  if(p%2 != 0) return(0.0);

  double z = ((p%4 == 0) ? 1.0 : -1.0)
           * sqrt(l1+1.0) * sqrt(l2+1.0) * sqrt(j1+1.0) * sqrt(j2+1.0)
           * clebsh_gordan(l1,l2,0,0,ll) * racah(l1,j1,l2,j2,s,ll);

  return(z);
}


/***********************************************************/
/*      Z-coefficientn Coefficient in Lane-Thomas          */
/***********************************************************/
double zbar_coefficient(const int l1, const int j1, const int l2, const int j2, const int s, const int ll)
{
  int p = (-l1+l2+ll)/2;
  if(p%2 != 0) return(0.0);

  double z = sqrt(l1+1.0) * sqrt(l2+1.0) * sqrt(j1+1.0) * sqrt(j2+1.0)
           * clebsh_gordan(l1,l2,0,0,ll) * racah(l1,j1,l2,j2,s,ll);

  return(z);
}


/***********************************************************/
/*      Reduced Matrix Element for Tensor Operator         */
/*      < l1j1 || T(YL,sigma_S)J || l0j0 >                 */
/*      M.B.Johnson, L.W.Owen, G.R.Satchler                */
/*      Phys. Rev. 142, 748 (1966)                         */
/*              Note: definition differs from JOS by       */
/*              the factor sqrt(2j1+1)                     */
/***********************************************************/
double reduced_matrix_element(const int lt, const int st, const int jt,
                              const int l0, const int j0, const int l1, const int j1)
{
  double x1,x2,x3,reduced_mat;
  const double FOUR_PI = 12.56637061435917295384;

  if(parity(lt) != parity(l0)*parity(l1)) return(0.0);

  if( abs( l0-l1)    > lt ||  (l0+l1)    < lt ) return(0.0);
  if( abs((j0-j1)/2) > jt || ((j0+j1)/2) < jt ) return(0.0);

  int llt = 2*lt;
  int jjt = 2*jt;
  int sst = 2*st;

  reduced_mat = 1.0/sqrt(FOUR_PI)
               *clebsh_gordan(j1,j0,1,-1,jjt)/sqrt(jjt+1.0)
               *sqrt( (j0+1.0) * (j1+1.0) * (llt+1.0) )
               *parity((j1-j0)/2)*parity((-l0+l1+lt)/2)*parity((j0-1)/2);

  if(sst == 2){
    x1 = (l0-j0/2.0)*(j0+1.0);
    x2 = (l1-j1/2.0)*(j1+1.0);
    if(jjt == llt){
      x3 = (lt==0) ? 0 :    (     x1-x2)/sqrt( lt             *(lt+1.0) );
    }
    else if(jjt == (llt-sst)){
      x3 = (lt == 0) ? 0 : -(lt  +x1+x2)/sqrt( lt*(2.0*lt+1.0)          );
    }
    else if(jjt == (llt+sst)){
      x3 =                  (lt+1-x1-x2)/sqrt(    (2.0*lt+1.0)*(lt+1.0) );
    }
    else{
      x3 = 1.0;
    }
  }
  else x3=1.0;
  reduced_mat *= x3;

  return(reduced_mat);
}
