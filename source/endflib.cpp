/******************************************************************************/
/**                                                                          **/
/**     ENDF LIB : ENDF Formatting Utility Library                           **/
/**                                                                          **/
/******************************************************************************/

#include <cstring>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>

using namespace std;

#include "endflib.h"

static int    ENDFReadArray (ifstream *, int, int, int    *);
static int    ENDFReadArray (ifstream *, int, int, double *);
static int    ENDFReadArray (ifstream *, int, int, int, double *);

static void   ENDFWriteArray (ENDF *, int, int *);
static void   ENDFWriteArray (ENDF *, int, double *);
static void   ENDFWriteArray (ENDF *, int, int, int, double *);

static void   ENDFExceedSubBlock (const string, ENDF *);
static void   ENDFExceedDataSize (const string, ENDF *, const int, const int);

static double ENDFPadExp (string);
static void   ENDFDelExp (double, char *);

static string line;
static string blank = "           ";
static int    seqno = 1;

static bool   LINE_NUMBER = false;


/**********************************************************/
/*      Seek Head                                         */
/**********************************************************/
static int mfs = 0;
static int mts = 0;

int ENDFSeekHead(ifstream *fp, ENDF *lib, const int mfsearch, const int mtsearch)
{
  int    mat=0;
  string s;

  /*** rewind tape, if query section is bigger than previous search */
  if((mfs == 0) && (mts == 0)) fp->seekg(0,ios_base::beg);
  if( (mfsearch < mfs) ||
      ((mfsearch == mfs) && (mtsearch < mts) )) fp->seekg(0,ios_base::beg);

  /*** look for given MF and MT */
  lib->setENDFmat(0);
  lib->setENDFmf(0);
  lib->setENDFmt(0);

  bool found = false;
  while( getline(*fp,line) ){

    /*** if shorter than record witdth,
         assume this is the same section */
    if(line.length() < 75) continue;

    s = line.substr(66,4);  mat = atoi(s.c_str());
    s = line.substr(70,2);  mfs = atoi(s.c_str());
    s = line.substr(72,3);  mts = atoi(s.c_str());

    if( (mfs == mfsearch) && (mts== mtsearch) ){
      found =true;
      mfs = mfsearch;
      mts = mtsearch;

      lib->setENDFmat(mat);
      lib->setENDFmf(mfs);
      lib->setENDFmt(mts);

      Record r = ENDFSplitCONT();
      lib->setENDFhead(r);
      break;
    }

    /*** not found */
    if( (mfs > mfsearch) || ((mfs == mfsearch) && (mts > mtsearch)) ){
      mfs = 0;
      mts = 0;
      break;
    }
  }

  if(found) return( 0);
  else      return(-1);
}


/**********************************************************/
/*      Scan ENDF Tape                                    */
/**********************************************************/
int ENDFScanLibrary(string libname, ENDFDict *dict)
{
  ifstream fp;
  string   s;
  int      mf0=0, mt0=0, mf1=1, mt1=451, overrun = 0;
  ENDF     lib(S);

  /*** open tape */
  fp.open(&libname[0]); if(!fp) return(-1);

  /*** tape ID */
  fp.seekg(0,ios_base::beg);
  getline(fp,line);
  strncpy(dict->tpid,&line[0],66);
  dict->tpid[66] = '\0';

  /*** first HEAD record */
  ENDFSeekHead(&fp,&lib,1,451);
  dict->setDICThead(lib.getENDFhead());
  dict->setMAT(lib.getENDFmat());

  /*** read in 3 CONT fields */
  dict->setDICTcont(0, ENDFNextCONT(&fp));
  dict->setDICTcont(1, ENDFNextCONT(&fp));
  dict->setDICTcont(2, ENDFNextCONT(&fp));

  /*** start at the TEXT field */
  dict->resetSEC();
  int c = 4;
  while( getline(fp,line) ){

    /*** if shorter than record witdth,
         assume this is the same section */
    if(line.length() < 75){
      c++;
      continue;
    }

    s = line.substr(70,2);  mf0 = atoi(s.c_str());
    s = line.substr(72,3);  mt0 = atoi(s.c_str());

    if( (mf0 == mf1) && (mt0 == mt1) ){
      c++;
      continue;
    }else if(mf0 == 0 || mt0 == 0)
      continue;

    dict->addDict(mf1,mt1,c,-1);

    if((mf1 == 1) && (mt1 == 451) && (c >= 9)) dict->setSTDHeader(true);

    mf1 = mf0;
    mt1 = mt0;
    c = 1;

    if(dict->getSEC() >= MAX_SECTION-1){
      overrun = -2;
      break;
    }
  }
  dict->addDict(mf1,mt1,c,-1);
  fp.close();

  /*** when standard header text is given, copy into ENDFDict object */
  if(dict->getSTDHeader()){
    fp.open(&libname[0]);
    getline(fp,line);
    for(int c=0 ; c<4 ; c++) getline(fp,line);
    for(int c=0 ; c<5 ; c++){
      getline(fp,line);
      strncpy(dict->text[c],&line[0],TEXT_WIDTH);
      int nc = line.length();
      if(nc < 66){
        for(int i=nc ; i<66 ; i++) dict->text[c][i] = ' ';
        dict->text[c][66] = '\0';
      }
    }
    fp.close();
  }

  return(overrun);
}


/**********************************************************/
/*      Split CONT Record                                 */
/**********************************************************/
Record ENDFSplitCONT()
{
  Record r;
  int    p=0;
  string s;

  s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;  r.c1 = ENDFPadExp(s);
  s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;  r.c2 = ENDFPadExp(s);
  s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;  r.l1 = atoi(s.c_str());
  s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;  r.l2 = atoi(s.c_str());
  s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;  r.n1 = atoi(s.c_str());
  s = line.substr(p,FIELD_WIDTH);                  r.n2 = atoi(s.c_str());

  return(r);
}


/**********************************************************/
/*      Next CONT Card                                    */
/**********************************************************/
Record ENDFNextCONT(ifstream *fp)
{
  getline(*fp,line);
  return( ENDFSplitCONT() );
}


/**********************************************************/
/*      Read CONT Record                                  */
/**********************************************************/
Record ENDFReadCONT(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadCONT",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);

  /*** keep the same address for arrays */
  lib->iptr[idx+1] = lib->iptr[idx];
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read LIST Record                                  */
/**********************************************************/
Record ENDFReadLIST(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadLIST",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);
  int nc = lib->rdata[idx].n1;

  if( lib->checkMAXDATA(0,nc) ) ENDFExceedDataSize("ReadLIST",lib,0,nc);

  /*** read in double data array */
  ENDFReadArray(fp, 0, nc, lib->xptr[idx]);

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx];
  lib->xptr[idx+1] = lib->xptr[idx] + nc;

  lib->inclPOS();

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read TAB1 Record                                  */
/**********************************************************/
Record ENDFReadTAB1(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadTAB1",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);
  int nr = 2*lib->rdata[idx].n1;
  int np = 2*lib->rdata[idx].n2;

  if( lib->checkMAXDATA(nr,np) ) ENDFExceedDataSize("ReadTAB1",lib,nr,np);

  /*** read in int and double data array */
  ENDFReadArray(fp, 0, nr, lib->iptr[idx]);
  ENDFReadArray(fp, 0, np, lib->xptr[idx]);

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx] + nr;
  lib->xptr[idx+1] = lib->xptr[idx] + np;

  lib->inclPOS();

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read TAB2 Record (LIST Type)                      */
/**********************************************************/
Record ENDFReadTAB2(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadTAB2",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);
  int nr = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  if( lib->checkMAXDATA(2*nr,0) ) ENDFExceedDataSize("ReadTAB2",lib,2*nr,0);

  /*** read in int data array */
  ENDFReadArray(fp, 0, 2*nr, lib->iptr[idx]);

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx] + 2*nr;
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();

  /*** read LIST for each subsection */
  for(int i=0 ; i<np ; i++) ENDFReadLIST(fp, lib);

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read TAB2 Record (TAB1 Type)                      */
/**********************************************************/
Record ENDFReadTAB21(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadTAB21",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);
  int nr = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  if( lib->checkMAXDATA(2*nr,0) ) ENDFExceedDataSize("ReadTAB21",lib,2*nr,0);

  /*** read in int data array */
  ENDFReadArray(fp, 0, 2*nr, lib->iptr[idx]);

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx] + 2*nr;
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();

  /*** read TAB1 for each subsection */
  for(int i=0 ; i<np ; i++) ENDFReadTAB1(fp, lib);

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read TAB2 Record (TAB2 Type)                      */
/**********************************************************/
Record ENDFReadTAB22(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadTAB22",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);
  int nr = lib->rdata[idx].n1;
  int ne = lib->rdata[idx].n2;

  if( lib->checkMAXDATA(2*nr,0) ) ENDFExceedDataSize("ReadTAB22",lib,2*nr,0);

  /*** read in int data array */
  ENDFReadArray(fp, 0, 2*nr, lib->iptr[idx]);

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx] + 2*nr;
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();

  /*** read TAB2 for each subsection */
  for(int i=0 ; i<ne ; i++) ENDFReadTAB21(fp, lib);

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read CONT+INTG Records                            */
/**********************************************************/
Record ENDFReadINTG(ifstream *fp, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("ReadINTG",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = ENDFNextCONT(fp);
  int nd = lib->rdata[idx].l1;
  int nm = lib->rdata[idx].n1;

  static int row[] = {0, 0, 18, 13, 11, 9, 8};

  int nmax = nm * (row[nd] + 2);
  if( lib->checkMAXDATA(0,nmax) ) ENDFExceedDataSize("ReadINTG",lib,0,nmax);

  /***  read in INTG-formatted int data into double array */
  ENDFReadArray(fp,nm,row[nd],nd,lib->xptr[idx]);

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx];
  lib->xptr[idx+1] = lib->xptr[idx] + nmax;

  lib->inclPOS();

  return(lib->rdata[idx]);
}


/**********************************************************/
/*      Read 1-Dim Array (double)                         */
/**********************************************************/
int ENDFReadArray(ifstream *fp, int m, int n, double *x)
{
  string s;

  if( (m == 0) && (n > 0) )     { m = numline(n);   }
  else if( (n == 0) && (m > 0) ){ n = COLUMN_NUMBER * m; }
  else return(0);

  int i=0;
  for(int j=0 ; j<m ; j++){
    getline(*fp,line);
    int p=0;
    for(int k=0 ; k<COLUMN_NUMBER ; k++){
      s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;
      if(s != blank) x[i++] = ENDFPadExp(s);
      else           x[i++] = 0.0;
      if(i >= n) break;
    }
    if(i >= n) break;
  }

  return(i);
}


/**********************************************************/
/*      Read 1-Dim Array (int)                            */
/**********************************************************/
int ENDFReadArray(ifstream *fp, int m, int n, int *x)
{
  string s;

  if( (m == 0) && (n > 0) )     { m = numline(n);   }
  else if( (n == 0) && (m > 0) ){ n = COLUMN_NUMBER * m; }
  else return(0);

  int i=0;
  for(int j=0 ; j<m ; j++){
    getline(*fp,line);
    int p=0;
    for(int k=0 ; k<COLUMN_NUMBER ; k++){
      s = line.substr(p,FIELD_WIDTH); p+=FIELD_WIDTH;
      if(s != blank) x[i++] = atoi(s.c_str());
      else           x[i++] = 0;
      if(i >= n) break;
    }
    if(i >= n) break;
  }

  return(i);
}


/**********************************************************/
/*      Read 1-Dim Array (INTG-type int)                  */
/**********************************************************/
int ENDFReadArray(ifstream *fp, int m, int n, int d, double *x)
{
  string s;

  int i = 0;
  for(int j=0 ; j<m ; j++){
    getline(*fp,line);
    s = line.substr(0,5); x[i++] = atof(s.c_str()); 
    s = line.substr(5,5); x[i++] = atof(s.c_str());

    int p = FIELD_WIDTH;
    if(d == 6) p --;
    for(int k=0 ; k<n ; k++){
      s = line.substr(p,(d+1)); p += d+1;
      x[i++] = (double)atoi(s.c_str());
    }
  }

  return (i);
}


/**********************************************************/
/*      Error Message                                     */
/**********************************************************/
void ENDFExceedSubBlock(const string loc, ENDF *lib)
{
  cerr << "too many sub-block at " << loc;
  cerr << "  MF = " << lib->getENDFmf();
  cerr << "  MT = " << lib->getENDFmt() << endl;
  cerr << "  current block count " << lib->getPOS() << endl;
  exit(-1);
}


void ENDFExceedDataSize(const string loc, ENDF *lib, const int ni, const int nx)
{
  cerr << "too many data-point at " << loc;
  cerr << "  MF = " << lib->getENDFmf();
  cerr << "  MT = " << lib->getENDFmt();
  cerr << "  requested Ni = " << ni << "  Nx = " << nx << endl;
  exit(-1);
}


/**********************************************************/
/*      Write TPID Record                                 */
/**********************************************************/
void ENDFWriteTPID(ENDFDict *dic)
{
  cout << left << setw(66) << dic->tpid;
  ENDFPrintRight(1,0,0);
}


/**********************************************************/
/*      Write HEAD Record                                 */
/**********************************************************/
void ENDFWriteHEAD(ENDF *lib)
{
  seqno = 1;
  Record head = lib->getENDFhead();
  ENDFWriteRecord(head);
  ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());
  lib->resetCTR();
}


/**********************************************************/
/*      Write TEXT Record                                 */
/**********************************************************/
void ENDFWriteTEXT(ENDF *lib, string txt)
{
  cout << txt.substr(0,66);
  ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());
}


/**********************************************************/
/*      Write DICT Record                                 */
/**********************************************************/
void ENDFWriteDICT(ENDF *lib, int mf, int mt, int nc, int mod)
{
  cout << blank << blank;
  cout << setw(11) << mf;
  cout << setw(11) << mt;
  cout << setw(11) << nc;
  cout << setw(11) << mod;
  ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());
}


/**********************************************************/
/*      Write Record                                      */
/**********************************************************/
void ENDFWriteRecord(Record cont)
{
  char   num1[14],num2[14];
  string s1,s2;

  ENDFDelExp(cont.c1,num1);
  ENDFDelExp(cont.c2,num2);

  cout << setw(FIELD_WIDTH) << num1    << setw(FIELD_WIDTH) << num2;
  cout << right;
  cout << setw(FIELD_WIDTH) << cont.l1 << setw(FIELD_WIDTH) << cont.l2;
  cout << setw(FIELD_WIDTH) << cont.n1 << setw(FIELD_WIDTH) << cont.n2;
}


/**********************************************************/
/*      Write CONT Record                                 */
/**********************************************************/
Record ENDFWriteCONT(ENDF *lib)
{
  int idx = lib->getCTR();
  Record cont = lib->rdata[idx];

  ENDFWriteRecord(cont);
  ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());

  lib->inclCTR();

  return(cont);
}


/**********************************************************/
/*      Write LIST Record                                 */
/**********************************************************/
Record ENDFWriteLIST(ENDF *lib)
{
  int idx = lib->getCTR();
  int np  = lib->rdata[idx].n1;

  Record cont = ENDFWriteCONT(lib);
  ENDFWriteArray(lib,np,lib->xptr[idx]);

  return(cont);
}


/**********************************************************/
/*      Write TAB1 Record                                 */
/**********************************************************/
Record ENDFWriteTAB1(ENDF *lib)
{
  int idx = lib->getCTR();
  int nr  = lib->rdata[idx].n1;
  int np  = lib->rdata[idx].n2;

  Record cont = ENDFWriteCONT(lib);
  ENDFWriteArray(lib,2*nr,lib->iptr[idx]);
  ENDFWriteArray(lib,2*np,lib->xptr[idx]);

  return(cont);
}


/**********************************************************/
/*      Write TAB2 Record (LIST Type)                     */
/**********************************************************/
Record ENDFWriteTAB2(ENDF *lib)
{
  int idx = lib->getCTR();
  int nr  = lib->rdata[idx].n1;
  int np  = lib->rdata[idx].n2;

  Record cont = ENDFWriteCONT(lib);
  ENDFWriteArray(lib,2*nr,lib->iptr[idx]);

  for(int n=0 ; n<np ; n++)  ENDFWriteLIST(lib);

  return(cont);
}


/**********************************************************/
/*      Write TAB2 Record (TAB1 Type)                     */
/**********************************************************/
Record ENDFWriteTAB21(ENDF *lib)
{
  int idx = lib->getCTR();
  int nr  = lib->rdata[idx].n1;
  int np  = lib->rdata[idx].n2;

  Record cont = ENDFWriteCONT(lib);
  ENDFWriteArray(lib,2*nr,lib->iptr[idx]);

  for(int n=0 ; n<np ; n++) ENDFWriteTAB1(lib);

  return(cont);
}


/**********************************************************/
/*      Write TAB2 Record (TAB2 Type)                     */
/**********************************************************/
Record ENDFWriteTAB22(ENDF *lib)
{
  int idx = lib->getCTR();
  int nr  = lib->rdata[idx].n1;
  int np  = lib->rdata[idx].n2;

  Record cont = ENDFWriteCONT(lib);
  ENDFWriteArray(lib,2*nr,lib->iptr[idx]);

  for(int n=0 ; n<np ; n++) ENDFWriteTAB21(lib);
  
  return(cont);
}


/**********************************************************/
/*      Write CONT+INTG Records                           */
/**********************************************************/
Record ENDFWriteINTG(ENDF *lib)
{
  int idx = lib->getCTR();
  Record cont = ENDFWriteCONT(lib);
  int nd = cont.l1;
  int nm = cont.n1;

  static int row[] = {0, 0, 18, 13, 11, 9, 8};

  ENDFWriteArray(lib,nm,row[nd],nd,lib->xptr[idx]);

  return(cont);
}


/**********************************************************/
/*      Write SEND Record                                 */
/**********************************************************/
void ENDFWriteSEND(ENDF *lib)
{
  for(int i=0 ; i<COLUMN_NUMBER ; i++) cout << blank;
  seqno = 99999;
  ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),0);
  seqno = 1;
}


/**********************************************************/
/*      Write FEND Record                                 */
/**********************************************************/
void ENDFWriteFEND(int mat)
{
  for(int i=0 ; i<COLUMN_NUMBER ; i++) cout << blank;
  seqno = 0;
  ENDFPrintRight(mat,0,0);
}


/**********************************************************/
/*      Write 1-Dim Data (double)                         */
/**********************************************************/
void ENDFWriteArray(ENDF *lib, int np, double *x)
{
  char num[14];

  int n = numline(np);
  int k = 0;
  for(int i=0 ; i<n ; i++){
    for(int j=0 ; j<COLUMN_NUMBER ; j++){
      if(k < np){
        ENDFDelExp(x[k++],num);
        cout << setw(FIELD_WIDTH) << num;
      }
      else cout << blank;
    }
    ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());
  }
}


/**********************************************************/
/*      Write 1-Dim Data (int)                            */
/**********************************************************/
void ENDFWriteArray(ENDF *lib, int np, int *x)
{
  int n = numline(np);
  int k = 0;
  for(int i=0 ; i<n ; i++){
    for(int j=0 ; j<COLUMN_NUMBER ; j++){
      if(k < np){
        cout << setw(FIELD_WIDTH) << x[k++];
      }
      else cout << blank;
    }
    ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());
  }
}


/**********************************************************/
/*      Write 1-Dim Data (INTG-type int)                  */
/**********************************************************/
void ENDFWriteArray(ENDF *lib, int m, int n, int d, double *x)
{
  int i = 0;
  for(int j=0 ; j<m ; j++){
    cout << setw(5) << (int)x[i++];
    cout << setw(5) << (int)x[i++];
    if(d != 6) cout << " ";

    for(int k=0 ; k<n ; k++){
      if(x[i] == 0.0){
        for(int c=0 ; c <= d ; c++) cout << " ";
      }
      else{
        cout << setw(d+1) << (int)x[i];
      }
      i++;
    }
    int p = 55 - n * (d+1);
    for(int k=0 ; k<p ; k++) cout << " ";
    ENDFPrintRight(lib->getENDFmat(),lib->getENDFmf(),lib->getENDFmt());
  }
}


/**********************************************************/
/*      MAT MT MT and Line Number                         */
/**********************************************************/
void ENDFPrintRight(int mat, int mf, int mt)
{
  cout << right;
  cout << setw(4) <<mat << setw(2) <<mf << setw(3) << mt;

  if(LINE_NUMBER) cout << setw(5) << seqno;
  cout << endl;

  if(seqno == 99999) seqno = 0;
  else seqno++;
}

void ENDFPrintLineNumber(bool x)
{ LINE_NUMBER = x; }


/**********************************************************/
/*      Pack CONT Record                                  */
/**********************************************************/
void ENDFPackCONT(Record cont, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("PackCONT",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = cont;

  /*** keep the same address for arrays */
  lib->iptr[idx+1] = lib->iptr[idx];
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();
}


/**********************************************************/
/*      Pack LIST Record                                  */
/**********************************************************/
void ENDFPackLIST(Record cont, double *xdat, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("PackLIST",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = cont;
  int nc = lib->rdata[idx].n1;

  if( lib->checkMAXDATA(0,nc) ) ENDFExceedDataSize("PackLIST",lib,0,nc);

  /*** store double data array */
  for(int i=0 ; i<nc ; i++) lib->xptr[idx][i] = xdat[i];

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx];
  lib->xptr[idx+1] = lib->xptr[idx] + nc;

  lib->inclPOS();
}


/**********************************************************/
/*      Pack TAB1 Record                                  */
/**********************************************************/
void ENDFPackTAB1(Record cont, int *idat, double *xdat, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("PackTAB1",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = cont;
  int nr = 2*lib->rdata[idx].n1;
  int np = 2*lib->rdata[idx].n2;

  if( lib->checkMAXDATA(nr,np) ) ENDFExceedDataSize("PackTAB1",lib,nr,np);

  /*** store int and double data array */
  for(int i=0 ; i<nr ; i++) lib->iptr[idx][i] = idat[i];
  for(int i=0 ; i<np ; i++) lib->xptr[idx][i] = xdat[i];

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx] + nr;
  lib->xptr[idx+1] = lib->xptr[idx] + np;

  lib->inclPOS();
}


/**********************************************************/
/*      Pack TAB2 Record (LIST Type)                      */
/**********************************************************/
void ENDFPackTAB2(Record cont, Record *cdat, int *idat, double **xtab, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("PackTAB2",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = cont;
  int nr = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  if( lib->checkMAXDATA(2*nr,0) ) ENDFExceedDataSize("PackTAB2",lib,2*nr,0);

  /*** store int data array */
  for(int i=0 ; i<2*nr ; i++) lib->iptr[idx][i] = idat[i];

  lib->iptr[idx+1] = lib->iptr[idx] + 2*nr;
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();

  /*** read LIST for each subsection */
  for(int i=0 ; i<np ; i++) ENDFPackLIST(cdat[i],xtab[i],lib);
}


/**********************************************************/
/*      Pack TAB2 Record (TAB1 Type)                      */
/**********************************************************/
void ENDFPackTAB21(Record cont, int *idat, Record *cdat, int **itab, double **xtab, ENDF *lib)
{
  if( lib->checkSUBBLOCK() ) ENDFExceedSubBlock("PackTAB21",lib);
  int idx = lib->getPOS();

  /*** store CONT record */
  lib->rdata[idx] = cont;
  int nr = lib->rdata[idx].n1;
  int np = lib->rdata[idx].n2;

  if( lib->checkMAXDATA(2*nr,0) ) ENDFExceedDataSize("PackTAB21",lib,2*nr,0);

  /*** store int data array */
  for(int i=0 ; i<2*nr ; i++) lib->iptr[idx][i] = idat[i];

  /*** calculate next address */
  lib->iptr[idx+1] = lib->iptr[idx] + 2*nr;
  lib->xptr[idx+1] = lib->xptr[idx];

  lib->inclPOS();

  /*** read TAB1 for each subsection */
  for(int i=0 ; i<np ; i++) ENDFPackTAB1(cdat[i],itab[i],xtab[i],lib);
}


/**********************************************************/
/*      Pack TAB2 Record (TAB2 Type)                      */
/**********************************************************/
//      not implemented
//Record ENDFPackTAB22(void)
//{
//}


/**********************************************************/
/*      Copy Object                                       */
/**********************************************************/
void ENDFLibCopy(ENDF *libsrc, ENDF *libdst)
{
  if((int)libsrc->getSIZE() > (int)libdst->getSIZE()){
    cerr << "cannot copy to smaller size object" << endl;
    exit(-1);
  }

  int nb = libsrc->getPOS();
  int ni = libsrc->getNI();
  int nx = libsrc->getNX();

  /*** copy HEAD record, reset index */
  libdst->setENDFhead( libsrc->getENDFhead() );
  libdst->setENDFmat( libsrc->getENDFmat() );
  libdst->setENDFmf( libsrc->getENDFmf() );
  libdst->setENDFmt( libsrc->getENDFmt() );

  /*** copy all integer and double data */
  for(int i=0 ; i<ni ; i++){ libdst->idata[i] = libsrc->idata[i]; }
  for(int i=0 ; i<nx ; i++){ libdst->xdata[i] = libsrc->xdata[i]; }

  /*** copy all CONT */
  for(int i=0 ; i<nb ; i++){ libdst->rdata[i] = libsrc->rdata[i]; }

  /*** re-assign pointers */
  libdst->resetPOS();
  int ofset = 0;
  for(int i=1 ; i<=nb ; i++){
    ofset = libsrc->iptr[i] - libsrc->iptr[0];
    libdst->iptr[i] = &libdst->idata[ofset];

    ofset = libsrc->xptr[i] - libsrc->xptr[0];
    libdst->xptr[i] = &libdst->xdata[ofset];

    libdst->inclPOS();
  }
}


/**********************************************************/
/*      Browse Object Contents                            */
/**********************************************************/
void ENDFLibPeek(ENDF *lib)
{
  Record head = lib->getENDFhead();

  /*** HEAD record */
  cout << lib->getENDFmat() << " "
       << lib->getENDFmf()  << " "
       << lib->getENDFmt()  << endl;

  cout <<" HEAD  : ";
  ENDFWriteRecord(head);
  cout << endl;
  cout <<" POS   : " << setw(11) << lib->getPOS() << endl;
  cout <<" CTR   : " << setw(11) << lib->getCTR() << endl;
  /*** each block */
  if(lib->getPOS() == 1){
    cout <<" CONT  : ";
    ENDFWriteRecord(lib->rdata[0]);
    cout << endl;
    cout <<" NINT  : " << setw(11) << lib->getNI() << endl;
    cout <<" NDBL  : " << setw(11) << lib->getNX() << endl;
  }
  else{
    for(int i=0 ; i<lib->getPOS() ; i++){
      cout <<" CONT  : ";
      ENDFWriteRecord(lib->rdata[i]);
      cout << endl;
      cout <<" NINT  : " << setw(11) << lib->iptr[i+1] - lib->iptr[i] << endl;
      cout <<" NDBL  : " << setw(11) << lib->xptr[i+1] - lib->xptr[i] << endl;
    }
  }
}


/**********************************************************/
/*      Extract One Section from ENDF Data and Print      */
/**********************************************************/
void ENDFExtract(ifstream *fp, int mf, int mt)
{
  ENDF   lib(L);
  string s1,s2;

  ENDFSeekHead(fp,&lib,mf,mt);
  ENDFWriteHEAD(&lib);

  while( getline(*fp,line) ){

    int n = line.length();
    if(n < 75){
      cout << left << line;
      n = 65-n;
      while(n >= 0){
        cout <<' ';
        n--;
      }
    }
    else{
      s1 = line.substr(0,66);
      s2 = line.substr(72,3);

      int mt0 = atoi(s2.c_str());
      if(mt0 == 0) break;
      cout << s1;
    }
    ENDFPrintRight(lib.getENDFmat(),mf,mt);
  }

  ENDFWriteSEND(&lib);
}


/**********************************************************/
/*      Interpolate Data                                  */
/**********************************************************/
double ENDFInterpolation(ENDF *lib, double x, bool dupflag, const int idx)
{
  int p = 0;
  if(x < lib->xptr[idx][0]) return(0.0);

  Record cont = lib->rdata[idx];

  /*** check if the max values are different */
  double xmax = lib->xptr[idx][2*(cont.n2-1)];
  if(x > xmax) return(0.0);


  /*** check if no interpolation is needed */
  bool found = false;
  if(dupflag){
    /* duplicated point, take the fisrt one */
    for(int i=0 ; i<cont.n2 ; i++){
      if(x == lib->xptr[idx][2*i]){
        p = 2*i;
        found = true;
        break;
      }
    }
  }
  else{
    /*** take the last one */
    for(int i=0 ; i<cont.n2 ; i++){
      if(x == lib->xptr[idx][2*i]){
        p = 2*i;
        found = true;
      }
    }
  }
  if(found) return(lib->xptr[idx][p+1]);


  p = 0;
  int i = 0, m = 0;
  double y = 0.0;
  for(int ir=0 ; ir<cont.n1 ; ir++){
    m = lib->iptr[idx][2*ir+1];
    for(int ip=i ; ip<lib->iptr[idx][2*ir]-1 ; ip++){
      if(x >= lib->xptr[idx][2*ip] && x<lib->xptr[idx][2*(ip+1)]){
        p = 2*ip;
        break;
      }
    }
    if(p > 0) break;
    i = lib->iptr[idx][2*ir]-1;
  }

  if(m > 0){
    double x1 = lib->xptr[idx][p  ];
    double y1 = lib->xptr[idx][p+1];
    double x2 = lib->xptr[idx][p+2];
    double y2 = lib->xptr[idx][p+3]; 

    if( (x1==x2) || (m==1) ) y = y1;
    else{
      switch(m){
      case 2:
        y =  (y2-y1)/(x2-x1) * (x-x1) + y1;
        break;
      case 3:
        y =  (y2-y1)/(log(x2)-log(x1)) * (log(x)-log(x1)) + y1;
        break;
      case 4:
        if((y1 == 0.0) && (y2 == 0.0)) y = 0.0;
        else{
          y =  (log(y2)-log(y1))/(x2-x1) * (x-x1) + log(y1);
          y = exp(y);
        }
        break;
      case 5:
        if((y1 == 0.0) && (y2 == 0.0)) y = 0.0;
        else{
          y =  (log(y2)-log(y1))/(log(x2)-log(x1)) * (log(x)-log(x1)) + log(y1);
          y = exp(y);
        }
        break;
      default:
        break;
      }
    }
  }
  return(y);
}


/**********************************************************/
/*      Merge X-Data Points Appeared in Two Sections      */
/**********************************************************/
int ENDFMergeXdata(ENDF *lib1, ENDF *lib2, double *z)
{
  int idx1 = 0, idx2 = 0;
  int n1  = lib1->rdata[idx1].n2;
  int n2  = lib2->rdata[idx2].n2;

  int i1 = 0, i2 = 0, i = 0;

  /*** check if the max values are different */
  double x1 = lib1->xptr[idx1][2*(n1-1)];
  double x2 = lib2->xptr[idx2][2*(n2-1)];

  if(x1 > x2){
    lib2->xptr[idx2][2*n2  ] = x1;
    lib2->xptr[idx2][2*n2+1] = 0.0;
    n2 ++;
  }
  else if(x1 < x2){
    lib1->xptr[idx1][2*n1  ] = x2;
    lib1->xptr[idx1][2*n1+1] = 0.0;
    n1 ++;
  }

  do{
    if( lib1->xptr[idx1][i1] <= lib2->xptr[idx2][i2] )
      z[i] = lib1->xptr[idx1][i1];
    else
      z[i] = lib2->xptr[idx2][i2];
    i+=2;

    if( ((i1/2+1) == n1) && ((i2/2+1) == n2) ) break;

    if( lib1->xptr[idx1][i1] < lib2->xptr[idx2][i2] ){
      i1+=2;
    }
    else if( lib1->xptr[idx1][i1] > lib2->xptr[idx2][i2] ){
      i2+=2;
    }
    else{
      i1+=2;
      i2+=2;
    }
  }while(i < MAX_DBLDATA);

  return(i/2);
}


/**********************************************************/
/*      Scan High-Side Boundary of Resonance Region       */
/**********************************************************/
void ENDFMF2boundary(ENDFDict *dict, ENDF *lib)
{
  double emaxRR = 0.0, emaxUR = 0.0;
  int    idx = 0;
  int    ner = lib->rdata[idx++].n1;
  bool   ssf = false;

  for(int i=0 ; i<ner ; i++){
    double eh   = lib->rdata[idx  ].c2;
    int    lru  = lib->rdata[idx  ].l1;
    int    lrf  = lib->rdata[idx  ].l2;
    int    nro  = lib->rdata[idx++].n1;

    int    lssf = lib->rdata[idx  ].l1;
    int    nls  = lib->rdata[idx++].n1;

    if((lru == 2) && (lssf == 1)) ssf = true;

    if( lru == 1 ){
      if(eh > emaxRR) emaxRR = eh;
    }
    else if(lru == 2){
      if(eh > emaxUR) emaxUR = eh;
    }

    if(lru == 1){
      if(lrf <= 3){
        if(nro != 0) idx++;
        idx += nls;
      }
      else if(lrf == 7){
        idx += 2*nls;
      }
    }
    else if(lru == 2 && lrf == 2){
      for(int inls=0 ; inls<nls ; inls++){
        int njs  = lib->rdata[idx++].n1;
        idx += njs;
      }
    }
  }

  /*** actual highest resonance range */
  double emaxRe = 0.0;

  /*** when RRR exists */
  if(emaxRR > 0.0){
    /*** when URR exists too */
    if(emaxUR > 0.0){
      if(ssf) emaxRe = emaxRR;  // if LSSF flag is on, take the RRR boundary
      else    emaxRe = emaxUR;  // otherwise URR high-side
    }
    /*** when only RRR */
    else emaxRe = emaxRR;
  }
  else{
    /*** when only URR exists (no RRR) */
    if(emaxUR > 0.0){
      if(ssf) emaxRe = 0.0;
      else    emaxRe = emaxUR;
    }
    /*** no resonance */
    else emaxRe = 0.0;
  }

  dict->setEboundary(emaxRR,emaxUR,emaxRe);
}


/**********************************************************/
/*      Remove E for FORTRAN Hack                         */
/**********************************************************/
inline void ENDFDelExp(double x, char *num)
{
  double min_val = numeric_limits<double>::min();
  ostringstream os;
  os.setf(ios::scientific, ios::floatfield);

  double z = fabs(x);
  /*** too small numbers will be rounded to zero */
  if(z < min_val){
    strcpy(num," 0.000000+0");
  }
  /*** when the negative exponent has 3 digits */
  else if(z <= 1.0e-99){
    os << setprecision(4) << setw(12) << x;
    string s = os.str();
    strcpy(num,s.c_str());
    for(int i=7 ; i<=10 ; i++) num[i] = num[i+1];
    //  from |+1.1234e-123|
    //  to   |+1.1234-123|
  }
  /*** when the exponent has 2 digits */
  else if( (z < 1.0e-09) || (1.0e+10 <= z && z < 1.0e+100) ){
    os << setprecision(5) << setw(12) << x;
    string s = os.str();
    strcpy(num,s.c_str());
    for(int i=8 ; i<=10 ; i++) num[i] = num[i+1];
    //  from |+1.12345e+12|
    //  to   |+1.12345+12|
  }
  /*** when the exponent has 3 digits */
  else if(z >= 1.0e+100){
    os << setprecision(4) << setw(12) << x;
    string s = os.str();
    strcpy(num,s.c_str());
    for(int i=7 ; i<=10 ; i++) num[i] = num[i+1];
    //  from |+1.1234e+123|
    //  to   |+1.1234+123|
  }
  else{
    os << setprecision(6) << setw(13) << x;
    string s = os.str();
    strcpy(num,s.c_str());
    num[ 9] = num[10];
    num[10] = num[12];
    //  from |+1.123456e+01|
    //  to   |+1.123456+1|
  }
  num[11] = '\0';
}


/**********************************************************/
/*      E Padding for FORTRAN Hack                        */
/**********************************************************/
inline double ENDFPadExp(string str)
{
  /*** if blank, return zero */
  if(str == "           ") return(0.0);

  /*** search for E, if exists return the floating point number */
  bool found = false;
  for(unsigned int i=0 ; i<str.size() ; i++){
    if(toupper(str[i]) == 'E') found = true;
  }
  if(found)  return(atof(str.c_str()));

  int len = str.size();
  int p1 = 0, p2 = 0, p3 = 0;

  /**** sign: check if this number is positive or negative */
  int sig = 1;
  for(p1 = 0; p1<len ; p1++){
    if(str[p1] == ' ') continue;
    else if(str[p1] == '+') sig =  1;
    else if(str[p1] == '-') sig = -1;
    else if( (str[p1] >= '0' && str[p1] <= '9') || str[p1]=='.' ) break;
  }

  /**** find + or - */
  found = false;
  char q = '+';
  for(p2 = p1 ; p2<len ; p2++){
    if(str[p2] == '+' || str[p2] == '-'){
      found = true;
      q = str[p2];
      break;
    }
    else if(str[p2] == ' ') continue;
  }

  /*** reconstruct the real number */
  char num[20];
  int k = 0;
  for(int i=0 ; i<(p2-p1) ; i++){
    if(str[i+p1] == ' ') continue;
    num[k++] = str[i+p1];
  }

  if(found){
    num[k++] = 'e';
    num[k++] = q;
    for(p3 = p2+1 ; p3<len ; p3++){
      if(str[p3] == ' ') continue;
      num[k++] = str[p3];
    }
  }
  num[k] = '\0';

  return(sig*atof(num));
}
