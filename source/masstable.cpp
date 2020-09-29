/******************************************************************************/
/**     Mass Table                                                           **/
/******************************************************************************/

#include <string>

using namespace std;

#include "constant.h"
#include "terminate.h"
#include "masstable.h"

//#include "masstable_aw95.h"              // Audi Wapstra 1995 table
//#include "masstable_ripl2.h"             // AW95 + FRDM95 from RIPL2
//#include "masstable_ripl3.h"             // AW03 + FRDM95 from RIPL3
//#include "masstable_audi2011.h"          // AW11 + FRDM95 from RIPL3
#include "masstable_audi2012_frdm2012.h"   // AW12 + FRDM2012

double mass_excess(int z, int a)
{
  double    mx  = 0.0;
  unsigned int za = z*1000+a;

  bool found = false;
  for(int i=0 ; i<nMassTable ; i++){
    if(MassTable[i].za == za){
      found = true;
      mx = MassTable[i].mass;
      break;
    }
  }

  if(!found) TerminateCode("mass data for not found", (int)za);

  return(mx);
}

#include <iostream>

/**********************************************************/
/*      Q-value Calculation                               */
/**********************************************************/
double qvalue(const int proj, const int targ, int mt)
{
  double q = 0.0, e1 = 0.0;
  int z = targ/1000;
  int a = targ - z*1000;

  int zp = proj/1000;
  int ap = proj - zp*1000;

  if( (mt == 1) || ((51 <= mt) && (mt <= 91))) return(q);

  if(      (600 <= mt) && (mt <= 649) ) mt = 600;
  else if( (650 <= mt) && (mt <= 699) ) mt = 650;
  else if( (700 <= mt) && (mt <= 749) ) mt = 700;
  else if( (750 <= mt) && (mt <= 799) ) mt = 750;
  else if( (800 <= mt) && (mt <= 849) ) mt = 800;

  double ep = 0.0;
  switch(proj){
  case    0:  ep = 0.0;       break;
  case    1:  ep = ENEUTRON;  break;
  case 1001:  ep = EPROTON;   break;
  case 1002:  ep = EDEUTERON; break;
  case 2004:  ep = EALPHA;    break;
  default  :                  break;
  }

  double e0 = mass_excess(z,a) + ep ;
  z += zp;
  a += ap;

  switch(mt){
  case  16:  e1 = mass_excess(z  ,a- 2) + ENEUTRON*2;                      break; // (x,2n)
  case  17:  e1 = mass_excess(z  ,a- 3) + ENEUTRON*3;                      break; // (x,3n)
  case  22:  e1 = mass_excess(z-2,a- 5) + ENEUTRON + EALPHA;               break; // (x,nA)
  case  23:  e1 = mass_excess(z-6,a-13) + ENEUTRON + EALPHA*3;             break; // (x,n3A)
  case  24:  e1 = mass_excess(z-2,a- 6) + ENEUTRON*2 + EALPHA;             break; // (x,2nA)
  case  25:  e1 = mass_excess(z-2,a- 7) + ENEUTRON*3 + EALPHA;             break; // (x,3nA)
  case  28:  e1 = mass_excess(z-1,a- 2) + ENEUTRON + EPROTON;              break; // (x,np)
  case  29:  e1 = mass_excess(z-4,a- 9) + ENEUTRON + EALPHA*2;             break; // (x,n2A)
  case  30:  e1 = mass_excess(z-4,a-10) + ENEUTRON*2 + EALPHA*2;           break; // (x,2n2A)
  case  32:  e1 = mass_excess(z-1,a- 3) + ENEUTRON + EDEUTERON;            break; // (x,nd)
  case  33:  e1 = mass_excess(z-1,a- 4) + ENEUTRON + ETRITON;              break; // (x,nt)
  case  34:  e1 = mass_excess(z-2,a- 4) + ENEUTRON + EHELIUM3;             break; // (x,nh)
  case  35:  e1 = mass_excess(z-5,a-11) + ENEUTRON + EDEUTERON + EALPHA*2; break; // (x,nd2A)
  case  36:  e1 = mass_excess(z-5,a-12) + ENEUTRON + ETRITON + EALPHA*2;   break; // (x,nt2A)
  case  37:  e1 = mass_excess(z  ,a- 4) + ENEUTRON*4;                      break; // (x,4n)
  case  41:  e1 = mass_excess(z-1,a- 3) + ENEUTRON*2 + EPROTON;            break; // (x,2np)
  case  42:  e1 = mass_excess(z-1,a- 4) + ENEUTRON*3 + EPROTON;            break; // (x,3np)
  case  44:  e1 = mass_excess(z-2,a- 3) + ENEUTRON + EPROTON*2;            break; // (x,n2p)
  case  45:  e1 = mass_excess(z-3,a- 6) + ENEUTRON + EPROTON + EALPHA;     break; // (x,npA)
  case  50:  e1 = mass_excess(z  ,a- 1) + ENEUTRON;                        break; // (x,n0)
  case 102:  e1 = mass_excess(z  ,a   );                                   break; // (x,g)
  case 600:
  case 103:  e1 = mass_excess(z-1,a- 1) + EPROTON;                         break; // (x,p)
  case 650:
  case 104:  e1 = mass_excess(z-1,a- 2) + EDEUTERON;                       break; // (x,d)
  case 700:
  case 105:  e1 = mass_excess(z-1,a- 3) + ETRITON;                         break; // (x,t)
  case 750:
  case 106:  e1 = mass_excess(z-2,a- 3) + EHELIUM3;                        break; // (x,h)
  case 800:
  case 107:  e1 = mass_excess(z-2,a- 4) + EALPHA ;                         break; // (x,A)
  case 108:  e1 = mass_excess(z-4,a- 8) + EALPHA*2;                        break; // (x,2A)
  case 109:  e1 = mass_excess(z-6,a-12) + EALPHA*3;                        break; // (x,3A)
  case 111:  e1 = mass_excess(z-2,a- 2) + EPROTON*2;                       break; // (x,2p)
  case 112:  e1 = mass_excess(z-3,a- 5) + EPROTON + EALPHA;                break; // (x,pA)
  case 113:  e1 = mass_excess(z-5,a-11) + ETRITON + EALPHA*2;              break; // (x,t2A)
  case 114:  e1 = mass_excess(z-5,a-10) + EDEUTERON + EALPHA*2;            break; // (x,d2A)
  case 115:  e1 = mass_excess(z-2,a- 3) + EPROTON + EDEUTERON;             break; // (x,pd)
  case 116:  e1 = mass_excess(z-2,a- 3) + EPROTON + ETRITON;               break; // (x,pt)
  case 117:  e1 = mass_excess(z-3,a- 6) + EDEUTERON + EALPHA;              break; // (x,dA)
  default :  e1 = e0; break;
  }

  q = e0 - e1;

  return(q*1e+6);
}


/**********************************************************/
/*      Threshold Energy                                  */
/**********************************************************/
double threshold(const int za, double q)
{
  int z = za/1000;
  int a = za - z*1000;
  double mt = a + mass_excess(z,a)/AMUNIT;
  double et = -(q * (mt+MNEUTRON)/mt);

  return(et);
}


