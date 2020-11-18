/******************************************************************************/
/**     DeCE TABLE for MF2                                                   **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>

using namespace std;

#include "dece.h"
#include "decetable.h"
#include "terminate.h"

static int DeceTableMF2RRR(ENDF *, int, int);
static int DeceTableMF2RR7(ENDF *, int);
static int DeceTableMF2URA(ENDF *, int);
static int DeceTableMF2URB(ENDF *, int);
static int DeceTableMF2URC(ENDF *, int);
static int DeceTableMF2AP(ENDF *, int);


/**********************************************************/
/*      Process MF=2 Resonance Parameters                 */
/**********************************************************/
void DeceTableMF2(ENDF *lib)
{
  int idx = 0;
  Record cont = lib->rdata[idx];
  int lfw = cont.l2;
  int ner = cont.n1;
  idx++;

  cout << "# Resonance parameters" << endl;
  cout << "#          NER" << setw(14) << ner << "  number of resonance energy ranges" << endl;
  cout << "#          LFW" << setw(14) << lfw << "  average fission width given or not" << endl;
  cout << endl;

  for(int iner=0 ; iner<ner ; iner++){
    cont = lib->rdata[idx++];
    double el   = cont.c1;
    double eh   = cont.c2;
    int    lru  = cont.l1;
    int    lrf  = cont.l2;
    int    nro  = cont.n1;
    int    naps = cont.n2;
    cout << "#           EL"; outVal(el); cout << "  lower limit of this energy range" << endl;
    cout << "#           EH"; outVal(eh); cout << "  upper limit" << endl;
    cout << "#          LRU" << setw(14) << lru  << "  0:scattering radius only, 1:RRR, 2:URR" << endl;

    if(lru == 1){
      cout << "#          LRF" << setw(14) << lrf  << "  1:SLBW, 2:MLBW, 3:RM, 4:AA, 7:R-Matrix" << endl;
    }else{
      cout << "#          LRF" << setw(14) << lrf  << "  1:energy independent widths exept for fission, 2: energy dependent widths" << endl;
    }

    cout << "#          NRO" << setw(14) << nro  << "  energy dependent scattering radius" << endl;
    cout << "#         NAPS" << setw(14) << naps << "  scattering radius control" << endl;

    if(nro == 1) idx = DeceTableMF2AP(lib,idx);

    /*** Resolved Resonance Region */
    if(lru == 1){
      if(lrf <= 3){
        cont = lib->rdata[idx++];
        double spi  = cont.c1;
        double ap   = cont.c2;
        int    nls  = cont.n1;

        cout << "#          SPI"; outVal(spi); cout << "  target spin" << endl;
        cout << "#           AP"; outVal(ap); cout << "  scattering radius" << endl;
        cout << "#          NLS" << setw(14) << nls << "  number of L-values" << endl;

        for(int inls=0 ; inls<nls ; inls++){
          cout << "# L       " << setw(4) << inls << endl;
          idx = DeceTableMF2RRR(lib,lrf,idx);
        }
      }
      else if(lrf == 7){
        idx = DeceTableMF2RR7(lib,idx);
      }
    }

    /*** Unresolved Resonance Region */
    else if(lru == 2){
      cont = lib->rdata[idx];
      double spi  = cont.c1;
      double ap   = cont.c2;
      int    lssf = cont.l1;

      cout << "#          SPI"; outVal(spi); cout << "  target spin" << endl;
      cout << "#           AP"; outVal(ap); cout << "  scattering radius" << endl;
      cout << "#         LSSF" << setw(14) << lssf << "  0: cross section calculated from URR parameters, 1:FILE2 used for self-shielding only" << endl;

      if(lrf == 1){
        // Case A, not tested
        if(lfw == 0) idx = DeceTableMF2URA(lib,idx);
        // Case B
        else         idx = DeceTableMF2URB(lib,idx);
      }
      // Case C
      else if(lrf == 2) idx = DeceTableMF2URC(lib,idx);
    }
    else{
      message << "this resonance parameters cannot be processed";
      WarningMessage();
    }
  }
}


/**********************************************************/
/*      Resolved Resonance Range                          */
/**********************************************************/
int DeceTableMF2RRR(ENDF *lib, int lrf, int idx)
{
  Record cont = lib->rdata[idx];
  int    nrs  = cont.n2;

  cout << "#          NRS" << setw(14) << nrs << "  number of resonances"<< endl;
  cout << "# E             J    ";
  if((lrf == 1) || (lrf == 2))
    cout << "G(total)      G(neutron)    G(gamma)      G(fission)" << endl;
  else if(lrf == 3)
    cout << "G(neutron)    G(gamma)      G(fissionA)   G(fissionB)" << endl;
  for(int i=0 ; i<nrs ; i++){
    int j = 6*i;
    outVal(lib->xptr[idx][j  ]);
    outVal(5,1,lib->xptr[idx][j+1]);
    outVal(lib->xptr[idx][j+2]);
    outVal(lib->xptr[idx][j+3]);
    outVal(lib->xptr[idx][j+4]);
    outVal(lib->xptr[idx][j+5]);
    cout << endl;
  }
  idx++;
  cout << endl;
  cout << endl;

  return(idx);
}


int DeceTableMF2RR7(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];
  int    ifg  = cont.l1;
  int    krm  = cont.l2;
  int    njs  = cont.n1;
  int    krl  = cont.n2;

  cout << endl;
  cout << "#          IFG" << setw(14) << ifg << "  0: GAM in eV, 1: GAM in sqrt(eV)" << endl;
  cout << "#          KRM" << setw(14) << krm << "  1:SLBW, 2:MLBW, 3:RM, 4:R-Matrix" << endl;
  cout << "#          NJS" << setw(14) << njs << "  number of J-pi values" << endl;
  cout << "#          KRL" << setw(14) << krl << "  0: non-relativistic, 1: relativistic" << endl;
  idx ++;

  cont = lib->rdata[idx];
  int    npp  = cont.l1;
  cout << "#          NPP" << setw(14) << npp << "  number of pairs"<< endl;
  cout << endl;

  char **channel;
  channel = new char * [npp];
  for(int i=0 ; i<npp ; i++) channel[i] = new char [11];

  /*** for each pair */
  for(int ipp=0 ; ipp<npp ; ipp++){
    double ma  = lib->xptr[idx][ipp*12];
    double mb  = lib->xptr[idx][ipp*12+1];
    double za  = lib->xptr[idx][ipp*12+2];
    double zb  = lib->xptr[idx][ipp*12+3];
    double ia  = lib->xptr[idx][ipp*12+4];
    double ib  = lib->xptr[idx][ipp*12+5];
    double q   = lib->xptr[idx][ipp*12+6];
    double pnt = lib->xptr[idx][ipp*12+7];
    double shf = lib->xptr[idx][ipp*12+8];
    double mt  = lib->xptr[idx][ipp*12+9];
    double pa  = lib->xptr[idx][ipp*12+10];
    double pb  = lib->xptr[idx][ipp*12+11];

    if(ia != 0.0) pa = (ia < 0.0) ? -1.0 : 1.0;
    if(ib != 0.0) pb = (ib < 0.0) ? -1.0 : 1.0;

    cout << "# Pair " << setw(3) << ipp+1 << "     MT  = " << setw(11) << (int)mt << endl;
    cout << "#              PNT = " << setw(11) << (int)pnt << "  1: calculate penetrability, -1: dont, 0: depends on MT" << endl;
    cout << "#              SHF = " << setw(11) << (int)shf << "  1: calculate shift factor, 0: dont" << endl;
    cout << "#          Q-value = "; outVal(11,2,q);
    cout << endl;

    cout << "#                    Mass     Charge spin  par" << endl;
    cout << "#             "; outVal(11,4,ma); outVal(11,4,za); outVal(5,1,abs(ia)); outVal(5,1,pa);
    cout << endl;
    cout << "#             "; outVal(11,4,mb); outVal(11,4,zb); outVal(5,1,abs(ib)); outVal(5,1,pb);
    cout << endl;

    int A = (int)(ma + 0.1);
    int Z = (int)za;

    if(      (A == 0) && (Z == 0) ) strncpy(channel[ipp],"Photon    ",11);
    else if( (A == 1) && (Z == 0) ) strncpy(channel[ipp],"Neutron   ",11);
    else if( (A == 1) && (Z == 1) ) strncpy(channel[ipp],"Proton    ",11);
    else if( (A == 2) && (Z == 1) ) strncpy(channel[ipp],"Deuteron  ",11);
    else if( (A == 3) && (Z == 1) ) strncpy(channel[ipp],"Triton    ",11);
    else if( (A == 3) && (Z == 2) ) strncpy(channel[ipp],"He-3      ",11);
    else if( (A == 4) && (Z == 2) ) strncpy(channel[ipp],"Alpha     ",11);
    else                            strncpy(channel[ipp],"unknown   ",11);
  }
  idx++;

  /*** J-Pi loop */
  for(int ijs=0 ; ijs<njs ; ijs++){
    cont = lib->rdata[idx];
    double aj   = cont.c1;
    double pj   = cont.c2;
    int    kbk  = cont.l1;
    int    kps  = cont.l2;
    int    nch  = cont.n2;

    int parity = (aj < 0.0) ? -1 : 1;
    if(aj == 0.0) parity = (int)pj;

    cout << "# "; outVal(4,1,abs(aj)); cout << ((parity < 0) ? "(-)" : "(+)") << endl;
    cout << "#          NBK" << setw(14) << kbk << "  non-zero if background R-matrix exist" << endl;
    cout << "#          KPS" << setw(14) << kps << "  non-zero if non-hard-sphere specified" << endl;
    cout << "#          NCH" << setw(14) << nch << "  number of channels"<< endl;

    int m = (nch+1)/6+1; if((nch+1)%6 == 0) m--;

    int ppi[7], l[7];
    for(int ich=0 ; ich<nch ; ich++){
      ppi[ich]    = (int)lib->xptr[idx][ich*6];
      l[ich]      = (int)lib->xptr[idx][ich*6+1];
      double sch  = lib->xptr[idx][ich*6+2];
      double bnd  = lib->xptr[idx][ich*6+3];
      double ape  = lib->xptr[idx][ich*6+4];
      double apt  = lib->xptr[idx][ich*6+5];

      cout << "#          PPI" << setw(14) << ppi[ich] << "  pair index" << endl;
      cout << "#            L" << setw(14) << l[ich]   << "  L-value" << endl;
      cout << "#          SCH"; outVal(sch); cout << "  channel spin" << endl;
      cout << "#          BND"; outVal(bnd); cout << "  boundary condition" << endl;
      cout << "#          APE"; outVal(ape); cout << "  effective channel radius" << endl;
      cout << "#          APT"; outVal(apt); cout << "  true channel radius" << endl;
    }
    idx++;

    cont = lib->rdata[idx];
    int    nrs  = cont.l2;
    cout << "#          NRS" << setw(14) << nrs << "  number of resonances" << endl;

    cout << "#  Energy[eV]  ";
    for(int ich=0 ; ich<nch ; ich++) cout << setw(10) << channel[ppi[ich]-1] << "    ";
    cout << endl;
    cout << "#         L = ";
    for(int ich=0 ; ich<nch ; ich++) cout << setw(14) << l[ich];
    cout << endl;

    for(int irs=0 ; irs<nrs ; irs++){
      int i0 = irs*m*6;
      outVal(lib->xptr[idx][i0]);
      for(int ich=0 ; ich<nch ; ich++) outVal(lib->xptr[idx][i0+ich+1]);
      cout << endl;
    }

    cout << endl;
    cout << endl;

    idx++;
  }

  for(int i=0 ; i<npp ; i++) delete [] channel[i];
  delete [] channel;

  return(idx);
}


/**********************************************************/
/*      Unresolved Resonance Range                        */
/**********************************************************/
int DeceTableMF2URA(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];

  int nls = cont.n1;
  cout << "#          NLS" << setw(14) << nls << "  number of L-values" << endl;
  idx++; 

  for(int inls=0 ; inls<nls ; inls++){
    cout << "# L       " << setw(4) << inls << endl;

    cont = lib->rdata[idx];
    int njs = cont.n2;

    cout << "#          NJS" << setw(14) << njs << "  number of J-values"<< endl;

    for(int injs=0 ; injs<njs ; injs++){

      cout << "# J            D             Deg.Free(n)   G(neutron)    G(gamma)" << endl;
      outVal(lib->xptr[idx][injs*6 + 1]);
      outVal(lib->xptr[idx][injs*6    ]);
      outVal(lib->xptr[idx][injs*6 + 2]);
      outVal(lib->xptr[idx][injs*6 + 3]);
      outVal(lib->xptr[idx][injs*6 + 4]);
      cout << endl;

      cout << endl;
      cout << endl;
    }
    idx++;
  }

  return(idx);
}


int DeceTableMF2URB(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx];

  int ne  = cont.n1;
  int nls = cont.n2;
  cout << "#          NLS" << setw(14) << nls << "  number of L-values" << endl;
  cout << "#           NE" << setw(14) << ne  << "  number of energy points" << endl;

  int m = idx; // index for the energy array

  idx++; 

  for(int inls=0 ; inls<nls ; inls++){
    cout << "# L       " << setw(4) << inls << endl;

    cont = lib->rdata[idx++];
    int njs = cont.n1;

    cout << "#          NJS" << setw(14) << njs << "  number of J-values"<< endl;

    for(int injs=0 ; injs<njs ; injs++){
      cont = lib->rdata[idx];
      double muf = (double)cont.l2;

      cout << "# J            D             Deg.Free(n)   Deg.Free(f)   G(neutron)    G(gamma)" << endl;
      outVal(lib->xptr[idx][1]);
      outVal(lib->xptr[idx][0]);
      outVal(lib->xptr[idx][2]);
      outVal(muf);
      outVal(lib->xptr[idx][3]);
      outVal(lib->xptr[idx][4]);
      cout << endl;

      cout << "# E            G(fission)" << endl;
      for(int i=0 ; i<ne ; i++){
        outVal(lib->xptr[m][i]);
        outVal(lib->xptr[idx][i+6]);
        cout << endl;
      }

      cout << endl;
      cout << endl;

      idx++;
    }
  }

  return(idx);
}


int DeceTableMF2URC(ENDF *lib, int idx)
{
  Record cont = lib->rdata[idx++];
  int nls = cont.n1;
  cout << "#          NLS" << setw(14) << nls << "  number of L-values" << endl;

  for(int inls=0 ; inls<nls ; inls++){
    cout << "# L       " << setw(4) << inls << endl;

    cont = lib->rdata[idx++];
    int njs = cont.n1;

    cout << "#          NJS" << setw(14) << njs << "  number of J-values"<< endl;

    for(int injs=0 ; injs<njs ; injs++){
      cont = lib->rdata[idx];
      double aj   = cont.c1;
      int    ne   = cont.n2;

      cout << "# J       "; outVal(4,1,aj); cout << endl;
      cout << "#           NE" << setw(14) << ne-1 << "  nubmber of energy points" << endl;
      cout << "# Deg. Freedom                other         neutron       gamma"<<endl;
      outVal(lib->xptr[idx][0]);
      outVal(lib->xptr[idx][1]);
      outVal(lib->xptr[idx][2]);
      outVal(lib->xptr[idx][3]);
      outVal(lib->xptr[idx][4]);
      outVal(lib->xptr[idx][5]);
      cout << endl;

      cout << "# E             D             G(other)      G(neutron)";
      cout << "    G(gamma)      G(fission)" << endl;

      for(int i=1 ; i<=ne ; i++){
        int j = 6*i;
        outVal(lib->xptr[idx][j  ]);
        outVal(lib->xptr[idx][j+1]);
        outVal(lib->xptr[idx][j+2]);
        outVal(lib->xptr[idx][j+3]);
        outVal(lib->xptr[idx][j+4]);
        outVal(lib->xptr[idx][j+5]);
        cout << endl;
      }
      idx ++;

      cout << endl;
      cout << endl;
    }
  }
  return(idx);
}


/**********************************************************/
/*      Energy-Dependent Scattering Radius, AP            */
/**********************************************************/
int DeceTableMF2AP(ENDF *lib, int idx)
{
  cout << endl;
  cout << "# Energy-dependent AP" << endl;
  ENDFPrint1Dim(lib,idx,"Energy","Radius");
  idx++;

  return(idx);
}


