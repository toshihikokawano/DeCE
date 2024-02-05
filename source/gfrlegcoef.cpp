/******************************************************************************/
/**     GFR, Calculate Legendre Coefficients from Resolved Resonances        **/
/******************************************************************************/

#include <complex>
#include <iostream>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "coupling.h"
#include "constant.h"

extern const double gcLorentzianWidth;
extern Smatrix Smat;


/**********************************************************/
/*      Legendre Coefficient from Resonance Parameters    */
/*      pl[] includes 4pi factor, so that pl[0] is        */
/*      the elastic scattering cross section              */
/**********************************************************/
void gfrLegendreCoefficient(System *sys, double *pl)
{
  int lmax = sys->nl;
  if(lmax == 0) return;
  for(int l=0 ; l<2*LMAX ; l++) pl[l] = 0.0;

  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01 / ((sys->target_spin2+1)*2);

  int smin = abs(sys->target_spin2-1);
  int smax =     sys->target_spin2+1;

  for(int l=0 ; l<2*lmax ; l++){
    double c = x1/(2.0*l + 1.0);

    for(int l0=0 ; l0<lmax ; l0++){
      for(int s0=smin ; s0<=smax ; s0+=2){

        int ss0    = s0 - smin - 1;
        int jj0min = abs(2*l0-s0);
        int jj0max =     2*l0+s0 ;
        for(int jj0=jj0min ; jj0<=jj0max ; jj0+=2){

          int idx0 = Smat.findIndex(l0,jj0,ss0);
          if(idx0 < 0) continue;

          for(int l1=0 ; l1<lmax ; l1++){
            if((l0+l1+l)%2 != 0) continue;
            for(int s1=smin ; s1<=smax ; s1+=2){
              if(s0 != s1) continue;

              int ss1    = s1 - smin - 1;
              int jj1min = abs(2*l1-s1);
              int jj1max =     2*l1+s1 ;
              for(int jj1=jj1min ; jj1<=jj1max ; jj1+=2){

                int idx1 = Smat.findIndex(l1,jj1,ss1);
                if(idx1 < 0) continue;

                /*** Re{(1-S0)(1-S1)^*} */
                double z = z_coefficient(2*l0,jj0,2*l1,jj1,s0,2*l);
                double x2 =  (1.0-real(Smat.getElement(idx0)))*(1.0-real(Smat.getElement(idx1)))
                            +     imag(Smat.getElement(idx0)) *     imag(Smat.getElement(idx1));
                pl[l] += c*z*z*x2;
              }
            }
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Compound Cross Section from Average S-matrix      */
/**********************************************************/
double gfrCompoundReaction(System *sys)
{
  int lmax = sys->nl;
  if(lmax == 0) return(0.0);

  double x1 = PI/(sys->wave_number*sys->wave_number) * 0.01;

  /*** Reaction Cross Section */
  double sigma = 0.0;
  for(int l=0 ; l<lmax ; l++){
    for(int s=-1 ; s<=1 ; s+=2){
      int j = l*2+s; if(j<0) continue;
      int idx = 2*l+(s+1)/2;
      double tj = 1 - (  real(Smat.getElement(idx))*real(Smat.getElement(idx))
                       + imag(Smat.getElement(idx))*imag(Smat.getElement(idx))  );
      if(tj < 0.0) tj = 0.0;

      sigma += ((s < 0) ? l : l+1) * tj;
    }
  }
  sigma *= x1;

  return(sigma);
}


