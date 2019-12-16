// Length of text field, don't change
const int MAX_FIELD = 66;

// Length of keyword in the patch data file
const int WORD_LENGTH = 9;

/*
   Define text data length and location
   There are five lines starting at the sixth row in ENDF.
   The location is defined as xpos = row - 6, and ypos = column.
   For example, when the text part is like

===========================================================================
 30-Zn- 65 INST       EVAL-Dec09 A.U.Thor                         3028 1451
                      DIST-MAY10                       20191211   3028 1451
----ENDF/B-X          MATERIAL 3028                               3028 1451
-----INCIDENT NEUTRON DATA                                        3028 1451
------ENDF-6 FORMAT                                               3028 1451
===========================================================================

The text for the institution "INST" is (11,0), and the length is 11.
This object will be created by TextField(11,11,0).
*/

class TextField{
  private:
    bool alloc;     // flag for memory allocation
    int  length;    // text data length
    char *text;     // data content
  public:
    int xpos, ypos; // location of the text data

  TextField(const int m, const int x, const int y){
    length = m;
    text = new char [m+1];
    for(int i=0 ; i<length ; i++) text[i] = ' ';
    text[length] = '\0';
    xpos = x;
    ypos = y;
    alloc = true;
  }

  ~TextField(){
    if(alloc){
      delete [] text;
    }
    alloc = false;
  }

  int getlen(){ return length; }

  void read(int c, char *src){
    if(c > length) c = length;
    for(int i=0 ; i<c ; i++) text[i] = src[i];
    if(c < length){
      for(int i=c ; i<length ; i++) text[i] = ' ';
    }
  }

  void copy(int c, char *src){
    if(c >= length + xpos){
      for(int i=0 ; i<length ; i++) text[i] = src[i+xpos];
    }
  }

  void paste(int c, char *dst){
    int ix = length;
    if(c < length + xpos) ix = c - xpos;
    for(int i=0 ; i<ix ; i++){
      char c = text[i];
      if(isprint(c) == 0) c = ' ';
      dst[i+xpos] = c;
    }
  }

  void print(){
    int ix = length-1;
    for(int i=ix ; i>=0 ; i--){
      if(text[i] != ' '){ ix = i; break; }
    }
    for(int i=0 ; i<=ix ; i++) cout << text[i];
  }

  void getText(){
    for(int i=0 ; i<length ; i++) cout << text[i];
  }
};



static int  DeceHeaderScanMF1       (ifstream *, char **, const bool);
static bool DeceHeaderCheckStandard (const int, char *[]);
static void DeceHeaderCopyData      (const int, char *[]);
static void DeceHeaderReadData      (ifstream *);
static void DeceHeaderReplaceData   (ENDFDict *);
static void DeceHeaderCreateLib     (ifstream *, ENDFDict *, const int, char *[]);
static void DeceHeaderPrint         ();


static TextField zsymam(11,0,0), alab(11,11,0), edate(11,22,0), auth(33,33,0);
static TextField refer(22,0,1), ddate(11,22,1), rdate(11,33,1), endate(11,55,1);
static TextField libname(18,4,2);
static TextField sublib(MAX_FIELD-5,5,3);
static TextField format(MAX_FIELD-6,6,4);

