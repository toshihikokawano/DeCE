/******************************************************************************/
/**     GFR, Determine Energy Grid                                           **/
/******************************************************************************/

#include <complex>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstring>

using namespace std;

#include "endflib.h"
#include "gfr.h"
#include "constant.h"

static void gfrSortResonanceEnergies     (double *, const int);
static int  gfrIncludeResonanceEnergies  (ENDF *, double *);

static const int nstep = 18;
static double estep[nstep] = {
  1.00,  1.10,  1.20,  1.30,  1.50,  1.75,
  2.00,  2.25,  2.50,  2.75,  3.00,  3.50,
  4.00,  4.50,  5.00,  6.00,  7.00,  8.50};

static double energy_step_accuracy  = 0.01;
static double energy_step_expansion = 1.2;
static double energy_step_shrink    = 0.8;

/**********************************************************/
/*      Energy Point in the Resolved Resonance Region     */
/**********************************************************/
int gfrAutoEnergyRRR(System *sys, ENDF *lib2, double *elab, const double ebr)
{
  int k = gfrIncludeResonanceEnergies(lib2,elab);

  if(k <= 0) return(0);
  int nr = k;
  int np = nr;

  gfrSortResonanceEnergies(elab,np);

  /*** half way to the first resonance */
  double e0 = elab[0] * 0.5;

  /*** add thermal point */
  elab[k++] = 0.0253;

  /*** up to the first resonance */
  double emag = 1.0e-5;
  double e1   = emag;
  while(emag < e0){
    for(int i=0 ; i<nstep ; i++){
      e1 = estep[i] * emag;
      if(e1 > e0){ emag = e0; break; }
      elab[k++] = e1;
    }
    emag *= 10.0;
  }
  np = k;
  if(np <= 3) return(np);

  /*** successive four points */
  k = k - 2;

  Pcross cs[4];
  double ex[4], ds[2], de;
  ex[0] = elab[k];
  ex[1] = elab[k+1];
  de    = ex[1] - ex[0];
  ex[2] = ex[1] + de;
  ex[3] = ex[2] + de;
  for(int i=0 ; i<4 ; i++) cs[i] = gfrCrossSection(1,ex[i],sys,lib2);

  int nstep = 1;
  do{
    double dm = ex[0] * 1e-05; // minimum interval

    ds[0] = ((ex[1]-ex[0])*cs[2].elastic + (ex[2]-ex[1])*cs[0].elastic)/(ex[2]-ex[0]);
    ds[1] = ((ex[2]-ex[1])*cs[3].elastic + (ex[3]-ex[2])*cs[1].elastic)/(ex[3]-ex[1]);

    ds[0] = abs(ds[0]/cs[1].elastic - 1.0);
    ds[1] = abs(ds[1]/cs[2].elastic - 1.0);

//  cout << "   Delta = " << ds[0] <<" " << ds[1] << "   " << dm <<  " "  << de << endl;

    /*** if linear approximation is reasonably good for the next point */
    if( (ds[0] < energy_step_accuracy) || (de == dm) ){

      /*** if all 4 points are on a straight line, expand interval */
      if(ds[1] < energy_step_accuracy){
        de *= energy_step_expansion;
      }

      /*** proceed to the next point */
      for(int i=1 ; i<4 ; i++){
        ex[i-1] = ex[i];
        cs[i-1] = cs[i];
      }

      /*** check if resonance exists between ex2 and ex2+dE */
      if(de > dm){
        int kr = -1;
        for(int j=0 ; j<nr ; j++){
          if((ex[2] < elab[j]) && (elab[j] <= ex[2]+de)){ kr = j ; break; }
        }
        if(kr >= 0){
          de = (elab[kr]-ex[2]) * 0.8;
          if(de <= dm) de = dm;
        }
      }

      ex[3] = ex[2] + de;
      cs[3] = gfrCrossSection(1,ex[3],sys,lib2);

    }

    /*** if not linear, shorten the interval */
    else{
      de *= energy_step_shrink;
      if(de <= dm) de = dm;
      for(int i=1 ; i<4 ; i++) ex[i] = ex[i-1] + de;
      for(int i=1 ; i<4 ; i++) cs[i] = gfrCrossSection(1,ex[i],sys,lib2);
    }

    elab[k++] = ex[1];
//  cout << k <<" " << ex[1] << " " << cs[1].elastic << endl;

    nstep ++;
    if(nstep >= MAX_DBLDATA/2){
      cerr << "too many energy points in GFR" << endl;
      break;
    }
  }while(ex[2] < ebr);

  np = k;
  
  gfrSortResonanceEnergies(elab,np);

//   /*** remove duplicated points */
//   char stri[15], strj[15];
//   int j = 0;
//   do{
//     sprintf(strj,"% 13.6e",elab[j]);
//     int i = j+1;
//     do{
//       sprintf(stri,"% 13.6e",elab[i]);
//       if(!strncmp(stri,strj,14)){
//         /*** shift array */
//         for(int m=i+1; m<np ; m++) elab[m-1] = elab[m];
//         np--;
//       }
//       i++;
//     }while(i < np);
//     j++;
//   }while(j < np);


/*
  char strj[15];
  for(int i=0 ; i<np ; i++){
    sprintf(strj,"% 13.6e",elab[i]);
    cout << setw(5) << i<< setw(13) << strj << endl;
  }
*/
  return(np);
}


/**********************************************************/
/*       Energy Points in Unresolved Resonance Range      */
/**********************************************************/
int gfrAutoEnergyURR(double *elab, const double ebr, const double ebu)
{
  if(ebu == 0.0) return(0);

  int k = 0;
  double e0 = ebr;
  double e2 = ebu;
  double e1 = 0.0;

  /*** when e0 is zero, no resolved resonance is given */
  if(e0 == 0.0){
    e0 = 1e-5;
    elab[k++] = e0;
    elab[k++] = 0.0253;

    do{
      e1 = e0 * 2.0; 
      elab[k] = e1;  if(k >= MAX_DBLDATA/2-2 || e1 > e2) break;
      k++;

      e1 = e0 * 5.0; 
      elab[k] = e1;  if(k >= MAX_DBLDATA/2-2 || e1 > e2) break;
      k++;

      e0 *= 10.0;
    }while(e0 < e2);

    elab[k++] = e2;

    gfrSortResonanceEnergies(elab,k);

  }
  else{
    elab[k++] = e0;
  
    /*** find the order of the first point */
    int q0 = (int)(log(e0)/log(10.0));
    int r0 = (int)(e0 /pow(10.0,(double)q0));

    /*** find the order of the last point */
    int q2 = (int)(log(e2)/log(10.0));

    /*** determine energy interval */
    double de = 1.0;
    if((q0 == q2) && (q0 >0)) de = pow(10.0,(double)(q0-1));
    else de = pow(10.0,(double)q0);

    e1 = r0 * pow(10,(double)q0);
    if(e1 <= e0) e1 = (r0+1) * pow(10,(double)q0);

    do{
      elab[k++] = e1;
      if(k >= MAX_DBLDATA/2-2) break;
    e1 += de;
    }while(e1 < e2);

    elab[k++] = e2;
  }

  return(k);
}


/**********************************************************/
/*      Equi-Distant Energy Points, User Input            */
/**********************************************************/
int gfrFixedEnergyRRR(double emin, double emax, double de, double *elab, const double ebr, const double ebu)
{
  if(ebu > 0.0){
    if(emax > ebu) emax = ebu;
  }
  else{
    if(emax > ebr) emax = ebr;
  }

  if(de <= 0.0){
    emax = emin;
    de   = 1.0;
  }

  int i = 0;
  double e = emin;
  while(e <= emax){
    elab[i] = e;
    e += de;
    i++;
  }

  return(i);
}


/**********************************************************/
/*      Include Resonance Energies in the Array           */
/**********************************************************/
int gfrIncludeResonanceEnergies(ENDF *lib2, double *elab)
{
  int k = 0, idx = 0;
  Record cont = lib2->rdata[idx++];
  int nrange = cont.n1;

  for(int i=0 ; i<nrange ; i++){
    cont = lib2->rdata[idx++];
    int lru   = cont.l1;
    int lrf   = cont.l2;
    double e0 = cont.c1;
    double e1 = cont.c2;

    if(lru == 1){

      elab[k++] = e1; // include energy boundary

      if(lrf <= 3){
        cont = lib2->rdata[idx++];
        int nls = cont.n1;

        for(int inls=0 ; inls<nls ; inls++){
          cont = lib2->rdata[idx];
          int nrs = cont.n2;

          for(int irs=0 ; irs<nrs ; irs++){
            double e2 = lib2->xptr[idx][6*irs];
            if((e0 <= e2) && (e2 <= e1)) elab[k++] = e2;
          }
          idx ++;
        }
      }
      else if(lrf == 7){
        cont = lib2->rdata[idx++];
        int njs = cont.n1;
        idx ++;
        for(int ijs=0 ; ijs<njs ; ijs++){

          cont = lib2->rdata[idx++];
          int nch = cont.n2;
          int m = (nch+1)/6+1; if((nch+1)%6 == 0) m--;

          cont = lib2->rdata[idx];
          int nrs = cont.l2;

          for(int irs=0 ; irs<nrs ; irs++){
            double e2 = lib2->xptr[idx][6*m*irs];
            if((e0 <= e2) && (e2 <= e1)) elab[k++] = e2;
          }
          idx ++;
        }
      }
    }
  }

  return(k);
}


/**********************************************************/
/*      Sort Energies                                     */
/**********************************************************/
void gfrSortResonanceEnergies(double *e, const int np)
{
  /*** quick sort */
  for(int j=0 ; j<np ; j++){
    int m  = j;
    for(int i=j ; i<np ; i++){
      if(e[i] <= e[m]) m = i;
    }
    double x = e[j];
    e[j] = e[m];
    e[m] = x;
  }
}
