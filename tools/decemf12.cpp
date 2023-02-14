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

//#define NEUTRON_INELASTIC_ONLY

const int NSUB =  50;  // maximum 50 levels for each particle
const int NSEC = 300;  // 50 levels for n, p, a, d, t, and h

int    main       (int, char *[]);
int    processMF12(ifstream *, ENDF **);

static int    mat = 0, mf = 12;
static double za = 0.0, awr = 0.0;

int main(int argc, char *argv[])
{
  if(argc < 2){
    cerr << "usage: decemf12 CoHoutput.dat ENDF_file" << endl;
    cerr << "       CoHoutput generated with -q1 option" << endl;
    exit(-1);
  }

  ifstream fpin;
  string   libname = "", datname = "";
  ENDF     wrk,*lib[NSEC];

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

  /*** read CoH output with -q1 option */
  fpin.open(datname.c_str());
  if(!fpin){
    cerr << "data file cannot open" << endl;  exit(-1);
  }

  for(int i=0 ; i<NSEC ; i++) lib[i] = new ENDF;

  int nmt = processMF12(&fpin,lib);

  fpin.close();

  /*** generate MF12 in the MT order */
  for(int mt = 51 ; mt < 850 ; mt++){

#ifdef NEUTRON_INELASTIC_ONLY
    if(mt >= 91) continue;
#endif

    for(int i=0 ; i<nmt ; i++){
      if(lib[i]->getENDFmt() == mt) ENDFWriteMF12(lib[i]);
    }
  }
  ENDFWriteFEND(mat);

  /*** generate MF14, isotropic angular distribution */
  wrk.setENDFmf(14);

  for(int mt = 51 ; mt < 850 ; mt++){

#ifdef NEUTRON_INELASTIC_ONLY
    if(mt >= 91) continue;
#endif

    for(int i=0 ; i<nmt ; i++){
      if(lib[i]->getENDFmt() == mt){
        Record cont = lib[i]->rdata[0];
        int nk = cont.n2;
        head.setRecord(za,awr,1,0,nk,0);
        int mt = lib[i]->getENDFmt();
        wrk.setENDFmt(mt);
        wrk.setENDFhead(head);
        ENDFWriteMF14(&wrk);
      }
    }
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
  bool   withconv = false;
  string dat1, dat2, dat3, line;
  double *ex = new double [NSUB];
  double *xdat = new double [3 * NSUB];

  const int lo = 2;
  const int lp = 0;

  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;
    if(line == "#gammaray"){
      withconv = false; break;
    }
    else if(line == "#gammaraywithconversion"){
      withconv = true; break;
    }
  }

  int lg = (withconv) ? 2 : 1; // 1: GP is always 1.0, 2: GP (gamma-ray probability) given

  int nmt = 0;
  while(1){
    getline(*fp,line);
    if(fp->eof() != 0) break;
    if(line.length() == 0) break;

    dat1 = line.substr( 5, 4);  int nk = atoi(&dat1[0]); // number of levels
    dat2 = line.substr( 0, 5);  int np = atoi(&dat2[0]); // particle ID

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
      if(k >= 41) continue; // this is because CoH produces 41 levels

      dat1 = line.substr( 9, 4);  // number of gamma-rays
      dat2 = line.substr(13,13);  // level energy

      int nt = atoi(&dat1[0]);
      ex[k]  = atof(&dat2[0]) * 1e+6;

      int gt=0;
      for(int g=0 ; g<nt ; g++){
        if(withconv){
          dat1 = line.substr(26 + 30*g,  4); // final state
          dat2 = line.substr(30 + 30*g, 13); // branching ratio
          dat3 = line.substr(43 + 30*g, 13); // gamma-ray emission probability
        }
        else{
          dat1 = line.substr(26 + 17*g,  4); // final state
          dat2 = line.substr(30 + 17*g, 13); // branching ratio
        }

        /*** eliminate zero transitios */
        double es = ex[ atoi(&dat1[0]) ];
        double br = atof(&dat2[0]);
        double gp = (withconv) ? atof(&dat3[0]) : 1.0;

        if(br > 0.0){
          if(withconv){
            xdat[3*gt  ] = es;
            xdat[3*gt + 1] = br;
            xdat[3*gt + 2] = gp;
          }
          else{
            xdat[2*gt    ] = es;
            xdat[2*gt + 1] = br;
          }
          gt++;
        }
      }
      nt = gt;

      int mt = 0;
      switch(np){
      case 1: mt =  50 + k; break;
      case 2: mt = 600 + k; break;
      case 3: mt = 800 + k; break;
      case 4: mt = 650 + k; break;
      case 5: mt = 700 + k; break;
      case 6: mt = 750 + k; break;
      default: cerr << "particle ID " << np << " out of range" << endl;  exit(-1);
      }

      Record head(za,awr,lo,lg,k,0);
      lib[nmt]->setENDFhead(head);
      lib[nmt]->setENDFmat(mat);
      lib[nmt]->setENDFmf(mf);
      lib[nmt]->setENDFmt(mt);

      Record cont(ex[k],0.0,lp,0,(lg+1)*nt,nt);
      ENDFPackLIST(cont,xdat,lib[nmt]);
      nmt ++;
    }

    np ++;
  }

  delete [] ex;
  delete [] xdat;

  return(nmt);
}
