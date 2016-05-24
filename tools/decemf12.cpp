/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Generate MF12 and MF14 for Discrete Gamma-Rays          **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "../source/endflib.h"

#define NEUTRON_INELASTIC_ONLY

const int NSUB = 200;  // maximum 200 levels
const int NSEC = 600;  // 200 levels for n, p, and alpha

int    main       (int, char *[]);
int    processMF12(ifstream *, ENDF **);

static int    mat = 0, mf = 12;
static double za = 0.0, awr = 0.0;

int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: decemf12 CoHoutput.dat ENDF_file" << endl;
    cerr << "       CoHoutput generated with -p512 option" << endl;
    exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     wrk(S),*lib[NSEC];

  datname = argv[1];
  libname = argv[2];

  /*** read comment section for ZA, AWR, MAT */
  fpin.open(libname.c_str());
  if(!fpin){
    cerr << "ENDF file cannot open" << endl;  exit(-1);
  }
  ENDFSeekHead(&fpin,&wrk,1,451);
  fpin.close();

  Record head = wrk.getENDFhead();
  za  = head.c1;
  awr = head.c2;
  mat = wrk.getENDFmat();

  /*** read ECLIPSE data */
  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }

  for(int i=0 ; i<NSEC ; i++) lib[i] = new ENDF(S);

  int nmt = processMF12(&fpin,lib);

  fpin.close();

  /*** generate MF12 */
  for(int i=0 ; i<nmt ; i++){
#ifdef NEUTRON_INELASTIC_ONLY
    if(lib[i]->getENDFmt() >= 91) continue;
#endif
    ENDFWriteMF12(lib[i]);
  }
  ENDFWriteFEND(mat);

  /*** generate MF14, isotropic angular distribution */
  wrk.setENDFmf(14);
  for(int i=0 ; i<nmt ; i++){
    Record cont = lib[i]->rdata[0];
    int nk = cont.n2;
    head.setRecord(za,awr,1,0,nk,0);
    int mt = lib[i]->getENDFmt();
    wrk.setENDFmt(mt);
    wrk.setENDFhead(head);
#ifdef NEUTRON_INELASTIC_ONLY
    if(lib[i]->getENDFmt() >= 91) continue;
#endif
    ENDFWriteMF14(&wrk);
  }
  ENDFWriteFEND(mat);

  /*** clean up */
  for(int i=0 ; i<NSEC ; i++)  delete lib[i];
  return(0);
}


/**********************************************************/
/*      Discrete Level Data Read In Array                 */
/**********************************************************/
int processMF12(ifstream *fp, ENDF **lib)
{
  string dat1, dat2, line;
  double ex[NSUB], xdat[2*NSUB];
  const int lo = 2;
  const int lg = 1;
  const int lp = 0;

  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;
    if(line == "#gammaray") break;
  }

  int np  = 0;
  int nmt = 0;
  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;
    if( line.length() ==0 ) break;
    if(np > 3) break;

    dat1 = line.substr( 5, 4);
    int nk = atoi(&dat1[0]);

    /*** skip capture part */
    if(np == 0){
      for(int k=0 ; k<nk ; k++) getline(*fp,line);
      np ++;
      continue;
    }

    for(int k=0 ; k<NSUB ; k++) ex[k] = 0.0;

    for(int k=0 ; k<nk ; k++){
      getline(*fp,line);
      if(k == 0){
        ex[k] = 0.0;
        continue;
      }

      dat1 = line.substr( 9, 4);  // number of gamma-rays
      dat2 = line.substr(13,13);  // level energy

      int nt = atoi(&dat1[0]);
      ex[k]  = atof(&dat2[0]) * 1e+6;

      int gt=0;
      for(int g=0 ; g<nt ; g++){
        dat1 = line.substr(26+17*g, 4); // final state
        dat2 = line.substr(30+17*g,13); // branching ratio

        /*** eliminate zero transitios */
        double br = atof(&dat2[0]);
        if(br>0.0){
          xdat[2*gt  ] = ex[ atoi(&dat1[0]) ];
          xdat[2*gt+1] = br;
          gt++;
        }
      }
      nt = gt;

      int mt = 0;
      if(     np == 1) mt =  50 + k;
      else if(np == 2) mt = 600 + k;
      else if(np == 3) mt = 800 + k;

      Record head(za,awr,lo,lg,k,0);
      lib[nmt]->setENDFhead(head);
      lib[nmt]->setENDFmat(mat);
      lib[nmt]->setENDFmf(mf);
      lib[nmt]->setENDFmt(mt);

      Record cont(ex[k],0.0,lp,0,2*nt,nt);
      ENDFPackLIST(cont,xdat,lib[nmt]);
      nmt ++;
    }

    np ++;
  }

  return(nmt);
}
