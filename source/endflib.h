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

#define MAX_SECTION             1000

#define FIELD_WIDTH               11   // width of data field
#define COLUMN_NUMBER              6   // number of data files in one line
#define TEXT_WIDTH                66   // numerical or text data width = FIELD_WIDTH x COLUMN_NUMBER


enum dataSize{S = 0, M = 1, L = 2};
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
  bool      allocated;  // memory allocation flag
  int       mat;        // ENDF MAT number
  int       mf;         // ENDF MF number
  int       mt;         // ENDF MT number
  Record    head;       // Head Record
  DataSize  size;       // ENDF data size
  int       ctr;        // pointer to currently running block
  int       nb;         // total number of blocks
 public:
  int       *idata;     // integer data
  double    *xdata;     // floating point data
  int      **iptr;      // int pointer for the 2-dim array
  double   **xptr;      // double pointer for the 2-dim array
  Record    *rdata;     // pointer for 2-dim data CONT

  ENDF(DataSize datasize){
    mat   = 0;
    mf    = 0;
    mt    = 0;
    allocated = false;
    memalloc(datasize);
    resetPOS();
  }

  ~ENDF(){
    memfree();
  }

  void memalloc(DataSize datasize){
    if(!allocated){
      if(datasize == S){
        idata = new int      [MAX_INTDATA_SMALL];
        xdata = new double   [MAX_DBLDATA_SMALL];
        iptr  = new int    * [MAX_SUBBLOCK_SMALL];
        xptr  = new double * [MAX_SUBBLOCK_SMALL];
        rdata = new Record   [MAX_SUBBLOCK_SMALL];
        size  = S;
      }
      else if(datasize == L){
        idata = new int      [MAX_INTDATA_LARGE];
        xdata = new double   [MAX_DBLDATA_LARGE];
        iptr  = new int    * [MAX_SUBBLOCK_LARGE];
        xptr  = new double * [MAX_SUBBLOCK_LARGE];
        rdata = new Record   [MAX_SUBBLOCK_LARGE];
        size  = L;
      }
      else{
        idata = new int      [MAX_INTDATA];
        xdata = new double   [MAX_DBLDATA];
        iptr  = new int    * [MAX_SUBBLOCK];
        xptr  = new double * [MAX_SUBBLOCK];
        rdata = new Record   [MAX_SUBBLOCK];
        size  = M;
      }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] idata;
      delete [] xdata;
      delete [] iptr ;
      delete [] xptr ;
      delete [] rdata;
      allocated = false;
    }
  }

  /*** reset pointers */
  void resetPOS(){
    if(allocated){
      if(size == S){
        for(int i=0 ; i<MAX_SUBBLOCK_SMALL ; i++){ iptr[i] = NULL; xptr[i] = NULL; }
      }
      else if(size == L){
        for(int i=0 ; i<MAX_SUBBLOCK_LARGE ; i++){ iptr[i] = NULL; xptr[i] = NULL; }
      }
      else{
        for(int i=0 ; i<MAX_SUBBLOCK       ; i++){ iptr[i] = NULL; xptr[i] = NULL; }
      }
      iptr[0] = &idata[0];
      xptr[0] = &xdata[0];
    }
    ctr = nb = 0;
  }

  /*** increment data block pointer */
  void inclPOS (void) { nb ++; }

  /*** number of currently stored blocks */
  int getPOS (void) { return nb; }

  /*** reset / increment block counter */
  void resetCTR (void) { ctr = 0; }
  void inclCTR (void) { ctr ++; }
  int  getCTR (void)   { return ctr; }

  /*** total number of data */
  int getNI (void) {
    int ni = 0;
    if(nb > 0) ni = iptr[nb] - iptr[0];
    return ni;
  }
  int getNX (void) {
    int nx = 0;
    if(nb > 0) nx = xptr[nb] - xptr[0];
    return nx;
  }

  DataSize getSIZE(void) { return size; }
  bool isalloc(void) { return allocated; }

  /** when HEAD is set, block counter reset */
  void setENDFhead(double c1, double c2, int l1, int l2, int n1, int n2){
    head.setRecord(c1,c2,l1,l2,n1,n2);
    resetCTR();
  }
  void setENDFhead(Record r){
    head = r;
    resetCTR();
  }

  /*** setting ENDF data */
  void setENDFmat(int n){ mat = n; }
  void setENDFmf (int n){ mf  = n; }
  void setENDFmt (int n){ mt  = n; }

  /*** inquire ENDF data */
  int getENDFmat(){ return mat; }
  int getENDFmf (){ return mf; }
  int getENDFmt (){ return mt; }

  Record getENDFhead(){ return head; }
  Record getENDFcont(){ return rdata[ctr]; }

  /*** check if more data can be stored */
  bool checkSUBBLOCK(void){
    if(size == S){
      if(nb >= MAX_SUBBLOCK_SMALL) return true;
    }else if(size == L){
      if(nb >= MAX_SUBBLOCK_LARGE) return true;
    }else{
      if(nb >= MAX_SUBBLOCK      ) return true;
    }
    return false;
  }

  bool checkMAXDATA(int mi, int mx){
    int ni0 = getNI() + mi;
    int nx0 = getNX() + mx;
    if(size == S){
      if((ni0 >= MAX_INTDATA_SMALL) || (nx0 >= MAX_DBLDATA_SMALL)) return true;
    }
    else if(size == L){
      if((ni0 >= MAX_INTDATA_LARGE) || (nx0 >= MAX_DBLDATA_LARGE)) return true;
    }
    else{
      if((ni0 >= MAX_INTDATA)       || (nx0 >= MAX_DBLDATA      )) return true;
    }
    return false;
  }
};


/**************************************/
/*      Class ENDF : Dictionary       */
/**************************************/
class ENDFDict{
 private:
  int       mat;          // ENDF MAT number
  int       sec;          // number of sections
  Record    head;         // Tape HEAD Record
  Record    cont[3];      // CONT Records
  char      *cbuf;        // buffer for TEXT lines and TPID
  bool      stdheader;    // flag for standard header
 public:
  int       *mf;          // ENDF MF number
  int       *mt;          // ENDF MT number
  int       *nc;          // Line count
  int       *mod;         // MOD number
  int       *id;          // ID for ENDF data on memory
  char      *tpid;        // Tape ID 
  char      *text[5];     // text data field
  double    emaxRR;       // energy boundary of resolved resonance region
  double    emaxUR;       // energy boundary of unresolved resonance region
  double    emaxRe;       // either emaxRR or emaxUR depending on LSSF flag

  ENDFDict(){
    mat    = 0;
    sec    = 0;
    emaxRR = 0.0;
    emaxUR = 0.0;
    emaxRe = 0.0;
    mf     = new int [MAX_SECTION];
    mt     = new int [MAX_SECTION];
    nc     = new int [MAX_SECTION];
    mod    = new int [MAX_SECTION];
    id     = new int [MAX_SECTION];

    int l = TEXT_WIDTH + 1;
    int p = 0;
    cbuf   = new char [l * 6];
    tpid   = &cbuf[p];  p += l;
    text[0]= &cbuf[p];  p += l;
    text[1]= &cbuf[p];  p += l;
    text[2]= &cbuf[p];  p += l;
    text[3]= &cbuf[p];  p += l;
    text[4]= &cbuf[p];

    stdheader = false;
  }

  ~ENDFDict(){
    delete [] mf;
    delete [] mt;
    delete [] nc;
    delete [] mod;
    delete [] id;
    delete [] cbuf;
  }

  void   setMAT(int n){ mat = n; }
  int    getMAT(){ return mat; }

  void   resetSEC(){ sec = 0; }
  int    getSEC(){ return sec; }

  void   setDICThead(Record r){ head = r; }
  Record getDICThead(){ return head; }

  void   setDICTcont(int i, Record r){ if(0 <= i && i < 3) cont[i] = r; }
  Record getDICTcont(int i){ return cont[i]; }

  void   setZA   (double za  ){ head.c1 = za;   }
  void   setAWR  (double awr ){ head.c2 = awr;  }
  void   setLRP  (int    lrp ){ head.l1 = lrp;  }
  void   setLFI  (int    lfi ){ head.l2 = lfi;  }
  void   setNLIB (int    nlib){ head.n1 = nlib; }
  void   setNMOD (int    nmod){ head.n2 = nmod; }

  double getZA   (){ return head.c1; }
  double getAWR  (){ return head.c2; }
  int    getLRP  (){ return head.l1; }
  int    getLFI  (){ return head.l2; }
  int    getNLIB (){ return head.n1; }
  int    getNMOD (){ return head.n2; }

  void   setELIS (double elis){ cont[0].c1 = elis; }
  void   setSTA  (double sta ){ cont[0].c2 = sta ; }
  void   setLIS  (int    lis ){ cont[0].l1 = lis ; }
  void   setLISO (int    liso){ cont[0].l2 = liso; }
  void   setNFOR (int    nfor){ cont[0].n2 = nfor; }

  double getELIS (){ return cont[0].c1; }
  double getSTA  (){ return cont[0].c2; }
  int    getLIS  (){ return cont[0].l1; }
  int    getLISO (){ return cont[0].l2; }
  int    getNFOR (){ return cont[0].n2; }

  void   setAWI  (double awi ){ cont[1].c1 = awi ; }
  void   setEMAX (double emax){ cont[1].c2 = emax; }
  void   setLREL (int    lrel){ cont[1].l1 = lrel; }
  void   setNSUB (int    nsub){ cont[1].n1 = nsub; }
  void   setNVER (int    nver){ cont[1].n2 = nver; }

  double getAWI  (){ return cont[1].c1; }
  double getEMAX (){ return cont[1].c2; }
  int    getLREL (){ return cont[1].l1; }
  int    getNSUB (){ return cont[1].n1; }
  int    getNVER (){ return cont[1].n2; }

  void   setTEMP (double temp){ cont[2].c1 = temp; }
  void   setLDRV (int    ldrv){ cont[2].l1 = ldrv; }
  void   setNWD  (int    nwd ){ cont[2].n1 = nwd;  }
  void   setNXC  (int    nxc ){ cont[2].n2 = nxc;  }

  double getTEMP (){ return cont[2].c1; }
  int    getLDRV (){ return cont[2].l1; }
  int    getNWD  (){ return cont[2].n1; }
  int    getNXC  (){ return cont[2].n2; }

  int    getProj (){ return(getNSUB()/10); } // no way to get projectile

  /*** if MF and MT are in Dict, return its index in Dict.
       if >=0, unknown if Lib is allocated. */
  int scanDict(int n1, int n2){
    int k = -1;
    for(int i=0 ; i<sec ; i++){
      if( (mf[i] == n1) && (mt[i] == n2) ){ k = i; break; }
    }
    return(k);
  }
  /*** if MF and MT are in Dict, return its Lib number.
       return -1 if not allocated. */
  int getID(int n1, int n2){
    int k = -1;
    for(int i=0 ; i<sec ; i++){
      if( (mf[i] == n1) && (mt[i] == n2) ){ k = id[i]; break; }
    }
    return(k);
  }
  bool addDict(int n1, int n2, int c, int k){
    mf[sec]  = n1;
    mt[sec]  = n2;
    nc[sec]  = c ;
    mod[sec] = 0 ;
    id[sec]  = k ;
    sec ++;
    if(sec >= MAX_SECTION) return true;
    else return false;
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
  }
  /*** check if fissile */
  bool isFission(){
    return( (head.l2 == 1) ? true : false );
  }
  /*** set / get stdheader flag */
  void setSTDHeader(bool c){ stdheader = c; }
  bool getSTDHeader(){ return stdheader; }
};


static int numline(int);
inline int numline(int n)
{ int m = (int)(n/COLUMN_NUMBER)+1;
  if(n%COLUMN_NUMBER == 0) m--;
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


void   ENDFWriteTPID      (ENDFDict *);
void   ENDFWriteHEAD      (ENDF *);
void   ENDFWriteRecord    (Record);
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
void   ENDFPrintLineNumber(const bool);
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
int    ENDFReadMF35       (ifstream *, ENDF *, const int);

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
void   ENDFWriteMF35      (ENDF *);

void   ENDFPrintLIST      (ENDF *, const int);
void   ENDFPrint1Dim      (ENDF *, const int);
