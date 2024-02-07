/******************************************************************************/
/**                                                                          **/
/**     ENDF LIB : ENDF and ENDFDict Class Definition                        **/
/**                                                                          **/
/******************************************************************************/

#ifndef __CSTDLIB__
#define __CSTDLIB__
#include <cstdlib>
#endif


/**************************************/
/*      Initial Memory Size           */
/**************************************/

static const int INIT_SUBBLOCK =      10;  // initially allocated CONT records
static const int INIT_INTDATA  =      10;  // initially allocated INT data buffer
static const int INIT_DBLDATA  =     100;  // initially allocated DBL data buffer
static const int MULT_MEMSIZE  =       2;  // multiplication of memory size


/**************************************/
/*      Limit Memory Size             */
/**************************************/

static const int MAX_SECTION   =    1000;  // max number of sections defined by MF/MT
static const int MAX_SUBBLOCK  =  100000;  // limit max number of sub-blocks (3.2MB)
static const int MAX_INTDATA   =  500000;  // integer data max size of 2 MB
static const int MAX_DBLDATA   = 5000000;  // set double data memory buffer of 400 MB


/**************************************/
/*      Text Field, Don't Change This */
/**************************************/

#define FIELD_WIDTH     11 // width of data field
#define COLUMN_NUMBER    6 // number of data files in one line
#define TEXT_WIDTH      66 // numerical or text data width = FIELD_WIDTH x COLUMN_NUMBER


/**************************************/
/*      Class Record : Control (CONT) */
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


/*  Note on ENDF Class

We assign an ENDF object to data speficied by a unique set of (MF,MT).
In each set, the structure will be:

   HEAD
   CONT + integer data
   CONT + double data
   CONT + (CONT + int) + (CONT + double)
   ... etc

Each ENDF object may consist of several blocks divided by the first
CONT line, and they include numerical data whose structure is
(MF,MT)-dependnet. These data block might be divided into several
sub-blocks (as the third block shown above). By including all the
sub-blocks (3 + 2 = 5), (the total number of CONT records NB = 5 in
the case above), this is the number of blocks currently stored in the
ENDF object. We access NB by getPOS() and inclPOS().

In each block, the numerical data are stored in the buffers, idata,
xdata, and rdata. 2-dimensional pointers, iptr and xptr, allow to
access the numerical data in the i-th sub-block.  xptr[0][0] is the
same as xdata[0]. For example, xptr[1][0] is the first double data in
the second sub-block, and iptr[3][0] is the first integer data in the
third sub-block.

The block counter CTR is used to print an ENDF object.
When the HEAD part is printed, this counter is set to zero,
then when printing each data block, CTR is used as a runing counter
that varies from zero to NB-1.

When the memalloc() method is called, a minimum size of ENDF object
is allocated. The memory size will be increased automatically when 
checkDataSize(mi,mx) says a requested memory area is larger than
currently allocated.
*/

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
  int       isize;      // integer data buffer size
  int       xsize;      // double data buffer size
  int       rsize;      // CONT data buffer size
  int       ctr;        // pointer to currently running block
  int       nb;         // total number of blocks
 public:
  int       *idata;     // integer data buffer
  double    *xdata;     // floating-point (double) data buffer
  int      **iptr;      // int pointer for the 2-dim array
  double   **xptr;      // double pointer for the 2-dim array
  Record    *rdata;     // CONT data buffer

  ENDF(){
    mat   = 0;
    mf    = 0;
    mt    = 0;
    allocated = false;
    memalloc();
    resetPOS();
  }

  ~ENDF(){
    memfree();
  }

  /*** allocate minimum size data buffers for Integer, Double, and CONT data,
       these data buffer will be expanded as needed */
  void memalloc(){
    if(!allocated){
      isize = INIT_INTDATA;  // initial integer data buffer size
      xsize = INIT_DBLDATA;  // initial double data buffer size
      rsize = INIT_SUBBLOCK; // initial CONT data buffer size

      idata = (int    * ) malloc(sizeof(int     ) * isize);
      xdata = (double * ) malloc(sizeof(double  ) * xsize);
      iptr  = (int    **) malloc(sizeof(int    *) * rsize);
      xptr  = (double **) malloc(sizeof(double *) * rsize);
      rdata = (Record * ) malloc(sizeof(Record  ) * rsize);

      allocated = true;
    }
  }

  /*** modify data buffer size by specified sizes */
  void memresize(const int nr, const int ni, const int nx){
    if(allocated){
      rsize = nr;
      isize = ni;
      xsize = nx;

      idata = (int    * ) realloc(idata,sizeof(int     ) * isize);
      xdata = (double * ) realloc(xdata,sizeof(double  ) * xsize);
      iptr  = (int    **) realloc(iptr ,sizeof(int    *) * rsize);
      xptr  = (double **) realloc(xptr ,sizeof(double *) * rsize);
      rdata = (Record * ) realloc(rdata,sizeof(Record  ) * rsize);
    }
  }

  /*** release allocated memories */
  void memfree(){
    if(allocated){
      free(idata);
      free(xdata);
      free(iptr);
      free(xptr);
      free(rdata);
      allocated = false;
    }
  }

  /*** reset all the allocated pointers */
  void resetPOS(){
    if(allocated){
      for(int i=0 ; i<rsize ; i++){ iptr[i] = NULL; xptr[i] = NULL; }
      iptr[0] = &idata[0];
      xptr[0] = &xdata[0];
    }
    ctr = nb = 0;
  }

  /*** increment data block pointer */
  void inclPOS (void) { nb ++; }

  /*** total number of currently stored blocks */
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

  int getISIZE(void) { return isize; }
  int getXSIZE(void) { return xsize; }
  int getRSIZE(void) { return rsize; }

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

  /*** check if more data can be stored,
       return true if memory allocation failed */
  bool checkSUBBLOCK(void){
    if(nb >= rsize-1){
      if(increaseSUBBLOCK() < 0) return true;
    }
    return false;
  }

  /*** check if INT or DBL data need more space,
       return true if memory allocation failed */
  bool checkDataSize(int mi, int mx){
    int ni0 = getNI() + mi;
    int nx0 = getNX() + mx;

    if(ni0 >= isize){
      if(increaseIDATA(mi) < 0) return true;
    }
    if(nx0 >= xsize){
      if(increaseXDATA(mx) < 0) return true;
    }
    return false;
  }

  /*** increase RDATA buffer size */
  int increaseSUBBLOCK(){
    if(allocated){
      /*** reallocate memory */
      rsize *= MULT_MEMSIZE;

      if(rsize >= MAX_SUBBLOCK) return -1;

      iptr  = (int    **) realloc(iptr ,sizeof(int    *) * rsize);
      xptr  = (double **) realloc(xptr ,sizeof(double *) * rsize);
      rdata = (Record  *) realloc(rdata,sizeof(Record  ) * rsize);

      /*** clear new pointers, exept i = NB */
      for(int i=nb+1 ; i<rsize ; i++){ iptr[i] = NULL; xptr[i] = NULL; }
    }
    return rsize;
  }

  /*** increase IDATA buffer size */
  int increaseIDATA(int mi){
    if(allocated){
      int ni0 = getNI() + mi;
      int ni1 = isize;
      if(ni1 > ni0) return isize;

      /*** increase INT buffer size until requested data can be stored */
      while(ni1 < ni0){
        ni1 *= MULT_MEMSIZE;
      }
      isize = ni1;
      if(isize >= MAX_INTDATA) return -1;

      /*** reallocate memory */
      idata = (int *) realloc(idata, sizeof(int) * isize);

      /*** re-assign pointers */
      int *base = iptr[0];
      iptr[0] = &idata[0];
      int ofset = 0;
      for(int i=1 ; i<=nb ; i++){
        ofset = iptr[i] - base;   // calculate number of elements of idata[i][n]
        iptr[i] = &idata[ofset];  // point idata[i][n]+1
      }
    }
    return isize;
  }

  /*** increase XDATA buffer size */
  int increaseXDATA(int mx){
    if(allocated){
      int nx0 = getNX() + mx;
      int nx1 = xsize;
      if(nx1 > nx0) return xsize;

      /*** increase DBL buffer size until requested data can be stored */
      while(nx1 < nx0){
        nx1 *= MULT_MEMSIZE;
      }
      xsize = nx1;
      if(xsize >= MAX_DBLDATA) return -1;

      /*** reallocate memory */
      xdata = (double *) realloc(xdata, sizeof(double) * xsize);

      /*** re-assign pointers */
      double *base = xptr[0];
      xptr[0] = &xdata[0];
      int ofset = 0;
      for(int i=1 ; i<=nb ; i++){
        ofset = xptr[i] - base;   // calculate number of elements of xdata[i][n]
        xptr[i] = &xdata[ofset];  // point xdata[i][n]+1
      }
    }
    return xsize;
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

  bool exist(int n1, int n2){
    if(getID(n1,n2) >= 0) return true;
    else return false;
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
Record ENDFReadTAB2L      (ifstream *, ENDF *);
Record ENDFReadINTG       (ifstream *, ENDF *);

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
Record ENDFWriteTAB2L     (ENDF *);
Record ENDFWriteINTG      (ENDF *);
void   ENDFWriteSEND      (ENDF *);
void   ENDFWriteFEND      (int);
void   ENDFPrintRight     (const int, const int, const int);
void   ENDFPrintLineNumber(const bool);
void   ENDFPackCONT       (Record, ENDF *);
void   ENDFPackLIST       (Record, double *, ENDF *);
void   ENDFPackTAB1       (Record, int *, double *, ENDF *);
void   ENDFPackTAB2       (Record, Record *, int *, double **, ENDF *);
void   ENDFPackTAB21      (Record, Record *, int *, int **, double **, ENDF *);
void   ENDFPackCopyCONT   (ENDF *, ENDF *, int);
void   ENDFPackCopyLIST   (ENDF *, ENDF *, int);
void   ENDFPackCopyTAB1   (ENDF *, ENDF *, int);
void   ENDFPackCopyTAB2   (ENDF *, ENDF *, int);
void   ENDFPackCopyTAB21  (ENDF *, ENDF *, int);

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
int    ENDFReadMF2        (ifstream *, ENDF *, const int);
int    ENDFReadMF3        (ifstream *, ENDF *, const int);
int    ENDFReadMF4        (ifstream *, ENDF *, const int);
int    ENDFReadMF5        (ifstream *, ENDF *, const int);
int    ENDFReadMF6        (ifstream *, ENDF *, const int);
int    ENDFReadMF7        (ifstream *, ENDF *, const int);
int    ENDFReadMF8        (ifstream *, ENDF *, const int);
int    ENDFReadMF9        (ifstream *, ENDF *, const int);
int    ENDFReadMF10       (ifstream *, ENDF *, const int);
int    ENDFReadMF11       (ifstream *, ENDF *, const int);
int    ENDFReadMF12       (ifstream *, ENDF *, const int);
int    ENDFReadMF13       (ifstream *, ENDF *, const int);
int    ENDFReadMF14       (ifstream *, ENDF *, const int);
int    ENDFReadMF15       (ifstream *, ENDF *, const int);
int    ENDFReadMF31       (ifstream *, ENDF *, const int);
int    ENDFReadMF32       (ifstream *, ENDF *, const int);
int    ENDFReadMF33       (ifstream *, ENDF *, const int);
int    ENDFReadMF34       (ifstream *, ENDF *, const int);
int    ENDFReadMF35       (ifstream *, ENDF *, const int);

void   ENDFWriteMF1       (ENDF *);
void   ENDFWriteMF2       (ENDF *);
void   ENDFWriteMF3       (ENDF *);
void   ENDFWriteMF4       (ENDF *);
void   ENDFWriteMF5       (ENDF *);
void   ENDFWriteMF6       (ENDF *);
void   ENDFWriteMF7       (ENDF *);
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
void   ENDFPrintLIST      (ENDF *, const int, string, string);
void   ENDFPrint1Dim      (ENDF *, const int);
void   ENDFPrint1Dim      (ENDF *, const int, string, string);
