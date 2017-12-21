/*
   endflib.h :
            prototype of ENDF in/out subroutines
            class definition
 */

#define MAX_DBLDATA            50000
#define MAX_INTDATA             1000
#define MAX_SUBBLOCK             500

#define MAX_DBLDATA_LARGE    2000000
#define MAX_INTDATA_LARGE     500000
#define MAX_SUBBLOCK_LARGE    100000

#define MAX_DBLDATA_SMALL        100
#define MAX_INTDATA_SMALL         10
#define MAX_SUBBLOCK_SMALL        10

#define MAX_SECTION          1000

#define FIELD_WIDTH            11
#define COLUMN_NUM              6

enum dataSize{S =0, M =1, L  =2};
typedef enum dataSize DataSize;


/**************************************/
/*      Class Record : Control        */
/**************************************/
class Record{
 public:
    double    c1;
    double    c2;
    int       l1;
    int       l2;
    int       n1;
    int       n2;
    Record(){
      c1 = 0.0; c2 = 0.0;  l1 = 0;  l2 = 0;  n1 = 0;  n2 = 0;
    }
    Record(double a, double b, int c, int d, int e, int f){
      c1 = a  ; c2 = b  ;  l1 = c;  l2 = d;  n1 = e;  n2 = f;
    }
    void setRecord(double a, double b, int c, int d, int e, int f){
      c1 = a  ; c2 = b  ;  l1 = c;  l2 = d;  n1 = e;  n2 = f;
    }
};


/**************************************/
/*      Class ENDF : Tape Data        */
/**************************************/
class ENDF{
 private:
  int       mat    ;     // ENDF MAT number
  int       mf     ;     // ENDF MF number
  int       mt     ;     // ENDF MT number
  Record    head   ;     // Head Record
  DataSize  size   ;     // ENDF data size
  int       pos    ;     // pointer to current block
  int       ni     ;     // current number of integer data
  int       nx     ;     // current number of double data
 public:
  int       *idata ;     // integer data
  double    *xdata ;     // floating point data
  int      **iptr  ;     // int pointer for the 2-dim array
  double   **xptr  ;     // double pointer for the 2-dim array
  Record    *rdata ;     // pointer for 2-dim data CONT

  ENDF(DataSize datasize){
    mat   = 0;
    mf    = 0;
    mt    = 0;
    pos   = 0;
    ni    = 0;
    nx    = 0;

    if(datasize == S){
      idata = new int      [MAX_INTDATA_SMALL];
      xdata = new double   [MAX_DBLDATA_SMALL];
      iptr  = new int    * [MAX_SUBBLOCK_SMALL];
      xptr  = new double * [MAX_SUBBLOCK_SMALL];
      rdata = new Record   [MAX_SUBBLOCK_SMALL];
      size  = S;
      for(int i=0 ; i<MAX_SUBBLOCK_SMALL ; i++){
        iptr[i] = NULL;
        xptr[i] = NULL;
      }
    }
    else if(datasize == L){
      idata = new int      [MAX_INTDATA_LARGE];
      xdata = new double   [MAX_DBLDATA_LARGE];
      iptr  = new int    * [MAX_SUBBLOCK_LARGE];
      xptr  = new double * [MAX_SUBBLOCK_LARGE];
      rdata = new Record   [MAX_SUBBLOCK_LARGE];
      size  = L;
      for(int i=0 ; i<MAX_SUBBLOCK_LARGE ; i++){
        iptr[i] = NULL;
        xptr[i] = NULL;
      }
    }
    else{
      idata = new int      [MAX_INTDATA];
      xdata = new double   [MAX_DBLDATA];
      iptr  = new int    * [MAX_SUBBLOCK];
      xptr  = new double * [MAX_SUBBLOCK];
      rdata = new Record   [MAX_SUBBLOCK];
      size  = M;
      for(int i=0 ; i<MAX_SUBBLOCK ; i++){
        iptr[i] = NULL;
        xptr[i] = NULL;
      }
    }

    iptr[0] = &idata[0];
    xptr[0] = &xdata[0];
  }

  ~ENDF(){
    delete [] idata;
    delete [] xdata;
    delete [] iptr ;
    delete [] xptr ;
    delete [] rdata;
  }
  void   setENDFmat(int n){    mat = n;    }
  void   setENDFmf (int n){    mf  = n;    }
  void   setENDFmt (int n){    mt  = n;    }

  int    getENDFmat()     {    return mat; }
  int    getENDFmf ()     {    return mf;  }
  int    getENDFmt ()     {    return mt;  }

  void   setENDFhead(double c1, double c2, int l1, int l2, int n1, int n2){
    head.setRecord(c1,c2,l1,l2,n1,n2);
    pos = ni = nx = 0;
  }
  void   setENDFhead(Record r){
    head = r;
  }
  Record getENDFhead(){
    return head;
  }
  int getPOS(void){
    return pos;
  }
  int getNI(void){
    return ni;
  }
  int getNX(void){
    return nx;
  }
  void setPOS(int i){
    pos = i;
  }
  void resetPOS(void){
    pos = 0;
  }
  void inclPOS(void){
    pos++;
  }
  void addPOS(int n){
    pos += n;
  }
  Record getENDFcont(){
    return rdata[pos];
  }
  bool checkSUBBLOCK(void){
    if(size == S){
      if(pos >= MAX_SUBBLOCK_SMALL) return true;
    }else if(size == L){
      if(pos >= MAX_SUBBLOCK_LARGE) return true;
    }else{
      if(pos >= MAX_SUBBLOCK      ) return true;
    }
    return false;
  }
  bool checkMAXDATA(int mi, int mx){
    ni += mi;
    nx += mx;
    if(size == S){
      if((ni >= MAX_INTDATA_SMALL) || (nx >= MAX_DBLDATA_SMALL)) return true;
    }
    else if(size == L){
      if((ni >= MAX_INTDATA_LARGE) || (nx >= MAX_DBLDATA_LARGE)) return true;
    }
    else{
      if((ni >= MAX_INTDATA)       || (nx >= MAX_DBLDATA      )) return true;
    }
    return false;
  }

  void copyENDF(ENDF *src){
    int n1,n2,n3;
    if(size == S){
      n1 = MAX_INTDATA_SMALL;
      n2 = MAX_DBLDATA_SMALL;
      n3 = MAX_SUBBLOCK_SMALL;
    }
    else if(size == L){
      n1 = MAX_INTDATA_LARGE;
      n2 = MAX_DBLDATA_LARGE;
      n3 = MAX_SUBBLOCK_LARGE;
    }
    else{
      n1 = MAX_INTDATA;
      n2 = MAX_DBLDATA;
      n3 = MAX_SUBBLOCK;
    }

    for(int i=0 ; i<n1 ; i++) idata[i] = src->idata[i];
    for(int i=0 ; i<n2 ; i++) xdata[i] = src->xdata[i];
    for(int i=0 ; i<n3 ; i++){
      if(src->iptr[i] != NULL){
        int ofset = src->iptr[i] - src->iptr[0];
        iptr[i] = &idata[ofset];
      }
      if(src->xptr[i] != NULL){
        int ofset = src->xptr[i] - src->xptr[0];
        xptr[i] = &xdata[ofset];
      }
      rdata[i] = src->rdata[i];
    }
    ni = src->getNI();
    nx = src->getNX();
  }
};


/**************************************/
/*      Class ENDF : Dictionary       */
/**************************************/
class ENDFDict{
 public:
  int       mat     ;     // ENDF MAT number
  int       sec     ;     // number of sections
  int       *mf     ;     // ENDF MF number
  int       *mt     ;     // ENDF MT number
  int       *nc     ;     // Line count
  int       *mod    ;     // MOD number
  int       *id     ;     // ID for ENDF data on memory
  char      tpid[67];     // Tape ID 
  Record    head    ;     // Tape HEAD Record
  Record    cont[3] ;     // CONT Records
  double    emax    ;     // highest energy in the file
  double    emaxRR  ;     // energy boundary of resolved resonance region
  double    emaxUR  ;     // energy boundary of unresolved resonance region
  double    emaxRe  ;     // either emaxRR or emaxUR depending on LSSF flag
  ENDFDict(){
    mat    = 0;
    sec    = 0;
    emax   = 0.0;
    emaxRR = 0.0;
    emaxUR = 0.0;
    emaxRe = 0.0;
    mf     = new int [MAX_SECTION];
    mt     = new int [MAX_SECTION];
    nc     = new int [MAX_SECTION];
    mod    = new int [MAX_SECTION];
    id     = new int [MAX_SECTION];
  }
  ~ENDFDict(){
    delete [] mf;
    delete [] mt;
    delete [] nc;
    delete [] mod;
    delete [] id;
  }
  int getNWD(){
    return cont[2].n1;
  }
  int getNXC(){
    return cont[2].n2;
  }
  int getProj(){
    int nsub = cont[1].n1;
    return(nsub/10);
  }
  void setNWD(int nwd){
    cont[2].n1 = nwd;
  }
  void setNXC(int nxc){
    cont[2].n2 = nxc;
  }
  /*** if MF and MT are in Dict, return its index in Dict.
       if >=0, unknown if Lib is allocated. */
  int scanDict(int n1, int n2){
    int k = -1;
    for(int i=0 ; i<sec ; i++){
      if( (mf[i]==n1) && (mt[i]==n2) ){ k = i; break; }
    }
    return(k);
  }
  /*** if MF and MT are in Dict, return its Lib number.
       return -1 if not allocated. */
  int getID(int n1, int n2){
    int k = -1;
    for(int i=0 ; i<sec ; i++){
      if( (mf[i]==n1) && (mt[i]==n2) ){ k = id[i]; break; }
    }
    return(k);
  }
  void addDict(int n1, int n2, int c, int k){
    mf[sec]  = n1;
    mt[sec]  = n2;
    nc[sec]  = c ;
    mod[sec] = 0 ;
    id[sec]  = k ;
    sec ++;
  }
  void setID(int i, int k){
    if(0 <= i && i < sec) id[i] = k;
  }
  /*** change MT number */
  void modMT(int i, int t){
    if(0 <= i && i < sec) mt[i] = t;
  }
  /*** save boundary energies of resolved and unresolved resonance regions */
  void setEboundary(double e1, double e2, double e3){
    emaxRR = e1;
    emaxUR = e2;
    emaxRe = e3;
    emax   = cont[1].c2;
  }
  /*** check if fissile */
  bool isFission(){
    return( (head.l2 == 1) ? true : false );
  }
};


static int numline(int);
inline int numline(int n)
{ int m = (int)(n/COLUMN_NUM)+1;
  if(n%COLUMN_NUM == 0) m--;
  return(m); }


/**************************************/
/*      endflib.cpp                   */
/**************************************/
int    ENDFSeekHead       (ifstream *, ENDF *, const int, const int);
int    ENDFScanLibrary    (string, ENDFDict *);
Record ENDFSplitCONT      (void);
Record ENDFNextCONT       (ifstream *);
Record ENDFReadCONT       (ifstream *, ENDF *);
Record ENDFReadLIST       (ifstream *, ENDF *);
Record ENDFReadTAB1       (ifstream *, ENDF *);
Record ENDFReadTAB2       (ifstream *, ENDF *);
Record ENDFReadTAB21      (ifstream *, ENDF *);
Record ENDFReadTAB22      (ifstream *, ENDF *);

int    ENDFReadArray      (ifstream *, int, int, int    *);
int    ENDFReadArray      (ifstream *, int, int, double *);

void   ENDFExceedSubBlock (const string, ENDF *);
void   ENDFExceedDataSize (const string, ENDF *, int, int);

void   ENDFWriteHEAD      (ENDF *);
void   ENDFWriteRecord    (Record *);
void   ENDFWriteTEXT      (ENDF *, string);
void   ENDFWriteDICT      (ENDF *, int, int, int, int);
Record ENDFWriteCONT      (ENDF *);
Record ENDFWriteLIST      (ENDF *);
Record ENDFWriteTAB1      (ENDF *);
Record ENDFWriteTAB2      (ENDF *);
Record ENDFWriteTAB21     (ENDF *);
Record ENDFWriteTAB22     (ENDF *);
void   ENDFWriteSEND      (ENDF *);
void   ENDFWriteFEND      (int);
void   ENDFWriteArray     (ENDF *, int, double *);
void   ENDFWriteArray     (ENDF *, int, int *);
void   ENDFPrintRight     (const int, const int, const int);
void   ENDFPackCONT       (Record, ENDF *);
void   ENDFPackLIST       (Record, double *, ENDF *);
void   ENDFPackTAB1       (Record, int *, double *, ENDF *);
void   ENDFPackTAB2       (Record, Record *, int *, double **, ENDF *);
void   ENDFPackTAB21      (Record, int *, Record *, int **, double **, ENDF *);

void   ENDFLibCopy        (ENDF *, ENDF *);
void   ENDFLibPeek        (ENDF *);
void   ENDFExtract        (ifstream *, const int, const int);

int    ENDFMergeXdata     (ENDF *, ENDF *, double *);
double ENDFInterpolation  (ENDF *, double, bool, const int);
void   ENDFMF2boundary    (ENDFDict *, ENDF *);


/**************************************/
/*      endfio.cpp                    */
/**************************************/
int    ENDFRead           (ifstream *, ENDF *, const int, const int);
void   ENDFWrite          (ENDF *);

int    ENDFReadMF1        (ifstream *, ENDF *, const int);
int    ENDFReadMF2        (ifstream *, ENDF *);
int    ENDFReadMF3        (ifstream *, ENDF *, const int);
int    ENDFReadMF4        (ifstream *, ENDF *, const int);
int    ENDFReadMF5        (ifstream *, ENDF *, const int);
int    ENDFReadMF6        (ifstream *, ENDF *, const int);
int    ENDFReadMF8        (ifstream *, ENDF *, const int);
int    ENDFReadMF9        (ifstream *, ENDF *, const int);
int    ENDFReadMF10       (ifstream *, ENDF *, const int);
int    ENDFReadMF11       (ifstream *, ENDF *, const int);
int    ENDFReadMF12       (ifstream *, ENDF *, const int);
int    ENDFReadMF13       (ifstream *, ENDF *, const int);
int    ENDFReadMF14       (ifstream *, ENDF *, const int);
int    ENDFReadMF15       (ifstream *, ENDF *, const int);
int    ENDFReadMF31       (ifstream *, ENDF *, const int);
int    ENDFReadMF32       (ifstream *, ENDF *);
int    ENDFReadMF33       (ifstream *, ENDF *, const int);
int    ENDFReadMF34       (ifstream *, ENDF *, const int);

void   ENDFWriteMF1       (ENDF *);
void   ENDFWriteMF2       (ENDF *);
void   ENDFWriteMF3       (ENDF *);
void   ENDFWriteMF4       (ENDF *);
void   ENDFWriteMF5       (ENDF *);
void   ENDFWriteMF6       (ENDF *);
void   ENDFWriteMF8       (ENDF *);
void   ENDFWriteMF9       (ENDF *);
void   ENDFWriteMF10      (ENDF *);
void   ENDFWriteMF11      (ENDF *);
void   ENDFWriteMF12      (ENDF *);
void   ENDFWriteMF13      (ENDF *);
void   ENDFWriteMF14      (ENDF *);
void   ENDFWriteMF15      (ENDF *);
void   ENDFWriteMF31      (ENDF *);
void   ENDFWriteMF32      (ENDF *);
void   ENDFWriteMF33      (ENDF *);
void   ENDFWriteMF34      (ENDF *);

void   ENDFPrintLIST      (ENDF *, const int);
void   ENDFPrint1Dim      (ENDF *, const int);
