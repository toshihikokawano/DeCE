/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Print / Replace ENDF Header Data                        **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <unistd.h>

using namespace std;

#include "../source/endflib.h"

class TextField{
  private:
    int length;
    char *text;
  public:
    int xpos, ypos;

  TextField(const int m, const int x, const int y){
    length = m;
    text = new char [m+1];
    for(int i=0 ; i<length ; i++) text[i] = ' ';
    text[length] = '\0';
    xpos = x;
    ypos = y;
  }

  ~TextField(){
    delete [] text;
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
    if(c > length + xpos) c = length + xpos;
    for(int i=xpos ; i<c ; i++) dst[i] = text[i-xpos];
  }

  void print(){
    int ix = length-1;
    for(int i=ix ; i>=0 ; i--){
      if(text[i] != ' '){ ix = i; break; }
    }
    for(int i=0 ; i<=ix ; i++) cout << text[i];
  }
};


const int MAX_FIELD = 66;

int main(int, char *[]);
static void readData(ifstream *);
static int checkMF1(ifstream *, char **, const bool);
static bool checkStandardHeader(const int, char *[]);

static void replaceHeader(const int, char *[]);
static void copyHeader(const int, char *[]);
static void printHeader();

static bool stdheader = false;

static TextField zsymam(11,0,0), alab(11,11,0), edate(11,22,0), auth(33,33,0);
static TextField refer(21,1,1), ddate(11,22,1), rdate(11,33,1), endate(11,55,1);
static TextField libname(18,4,2), matnum(4,31,2), sublib(MAX_FIELD-5,5,3);
static TextField format(MAX_FIELD-6,6,4);


int main(int argc, char *argv[])
{
  ENDFDict dict;
  ifstream libin, datin;
  string   libname = "", datname = "";
  char     *cbuf = NULL, **line = NULL;

  /*** command line options */
  int p = 0;
  while((p=getopt(argc,argv,"d:f:"))!=-1){
    switch(p){
    case 'd':
      datname = optarg;
      break;
    case 'f':
      libname = optarg;
      break;
    default:
      break;
    }
  }

  if(ENDFScanLibrary(libname,&dict) < 0){
    cerr << "ENDF file cannot open " << libname << endl;
    exit(-1);
  }


  /*** check MF1 */
  libin.open(libname.c_str());
  int nline = checkMF1(&libin,line,false);

//  cout << nline << " " << dict.getNWD() << " " << dict.getNXC() << endl;

  cbuf = new char [MAX_FIELD * nline];
  line = new char * [nline];
  for(int i=0 ; i<nline ; i++) line[i] = &cbuf[MAX_FIELD * i];

  checkMF1(&libin,line,true);
  stdheader = checkStandardHeader(nline,line);

  if(stdheader) copyHeader(nline,line);

  if(datname != ""){
    datin.open(datname.c_str());
    if(!datin){
      cerr << "Data file cannot open " << datname << endl;
      exit(-1);
    }
    readData(&datin);
    replaceHeader(nline,line);
  }
  else printHeader();


/*
  ENDFWriteTPID(&dict);
  
  for(int mf=1 ; mf <= 40 ; mf++){

    for(int i=0 ; i<dict.getSEC() ; i++){
      if(dict.mf[i] == mf){
        ENDFExtract(&libin,mf,dict.mt[i]);
      }
    }
  }

  ENDFWriteFEND(0);
  ENDFWriteFEND(-1);
*/
  libin.close();
  datin.close();

  delete [] cbuf;
  delete [] line;

  return(0);
}


void readData(ifstream *fpin)
{
  const int WORD_LENGTH = 9;
  string data;
  char s1[WORD_LENGTH+1], s2[MAX_FIELD+1];

  while( getline(*fpin,data) ){
    int n = data.length();
    if(n < WORD_LENGTH) continue;
    else if(n > 66) n = 66;

    for(int i=0 ; i<WORD_LENGTH ; i++) s1[i] = data[i];
    s1[WORD_LENGTH] = '\0';

    int i0, i1;
    for(i0=WORD_LENGTH ; i0<n ; i0++) if(data[i0] != ' ') break;
    for(i1=n-1 ; i1>=-i0 ; i1--) if(data[i1] != ' ') break;
    for(int i=i0 ; i<=i1 ; i++) s2[i-i0] = data[i];
    for(int i=i1-i0+1 ; i<MAX_FIELD ; i++) s2[i] = ' ';
    s2[MAX_FIELD] = '\0';
//  cout << "DEBUB:" << s2 << "| " << i0 <<" " << i1 <<endl;

    if(     !strncmp(s1,"  ZSYMAM:",WORD_LENGTH)){ zsymam.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"    ALAB:",WORD_LENGTH)){ alab.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"    AUTH:",WORD_LENGTH)){ auth.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   REFER:",WORD_LENGTH)){ refer.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   EDATE:",WORD_LENGTH)){ edate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   DDATE:",WORD_LENGTH)){ ddate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   RDATE:",WORD_LENGTH)){ rdate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"  ENDATE:",WORD_LENGTH)){ endate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1," LIBNAME:",WORD_LENGTH)){ libname.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"     MAT:",WORD_LENGTH)){ matnum.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"  SUBLIB:",WORD_LENGTH)){ sublib.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"  FORMAT:",WORD_LENGTH)){ format.read(MAX_FIELD,s2); }
  }
}



int checkMF1(ifstream *fpin, char *line[], bool copy)
{
  ENDF lib(S);
  string data, s;

  fpin->seekg(0,ios_base::beg);
  ENDFSeekHead(fpin,&lib,1,451);

  for(int i=0 ; i < 3 ; i++) getline(*fpin,data);

  int m = 0;
  while( getline(*fpin,data) ){

    int n = data.length();
    if(n >= 75){
      s = data.substr(72,3);
      int mt = atoi(s.c_str());
      if(mt == 0) break;
    }

    if(copy){
      for(int i=0 ; i<66 ; i++) line[m][i] = (i < n) ? data[i] : ' ';
      line[m][66] = '\0';
      //cout << m << " : " << line[m] << endl;
    }

    m ++;
  }

  return m;
}


bool checkStandardHeader(const int nline, char *line[])
{
  if(nline < 5) return false;

  char dash1[5], dash2[6], dash3[7];

  strncpy(dash1, &line[2][0],4); dash1[4] = '\0';
  strncpy(dash2, &line[3][0],5); dash2[5] = '\0';
  strncpy(dash3, &line[4][0],6); dash3[6] = '\0';

  /*** maybe this is in the standard format */
  bool test1 = strncmp(dash1,"----",4);
  bool test2 = strncmp(dash2,"-----",5);
  bool test3 = strncmp(dash3,"------",6);
  bool test4 = (line[2][4] == '-');
  bool test5 = (line[3][5] == '-');
  bool test6 = (line[4][6] == '-');

  if( !test1 && !test2 && !test3 && !test4 && !test5 && !test6) return true;
  else return false;
}


void copyHeader(const int nline, char *line[])
{
  if(nline < 2) return;

  zsymam.copy( MAX_FIELD,line[ zsymam.ypos]);
  alab.copy(   MAX_FIELD,line[   alab.ypos]);
  edate.copy(  MAX_FIELD,line[  edate.ypos]);
  auth.copy(   MAX_FIELD,line[   auth.ypos]);
  refer.copy(  MAX_FIELD,line[  refer.ypos]);
  ddate.copy(  MAX_FIELD,line[  ddate.ypos]);
  rdate.copy(  MAX_FIELD,line[  rdate.ypos]);
  endate.copy( MAX_FIELD,line[ endate.ypos]);
  libname.copy(MAX_FIELD,line[libname.ypos]);
  matnum.copy( MAX_FIELD,line[ matnum.ypos]);
  sublib.copy( MAX_FIELD,line[ sublib.ypos]);
  format.copy( MAX_FIELD,line[ format.ypos]);
}


void replaceHeader(const int nline, char *line[])
{
  if(nline < 2) return;

  zsymam.paste(MAX_FIELD,line[zsymam.ypos]);
  alab.paste(  MAX_FIELD,line[  alab.ypos]);
  edate.paste( MAX_FIELD,line[ edate.ypos]);
  auth.paste(  MAX_FIELD,line[  auth.ypos]);
  refer.paste( MAX_FIELD,line[ refer.ypos]);
  ddate.paste( MAX_FIELD,line[ ddate.ypos]);
  rdate.paste( MAX_FIELD,line[ rdate.ypos]);
  endate.paste(MAX_FIELD,line[endate.ypos]);
}


void printHeader()
{
  cout << "  ZSYMAM: "; zsymam.print(); cout << endl;
  cout << "    ALAB: "; alab.print();   cout << endl;
  cout << "    AUTH: "; auth.print();   cout << endl;
  cout << "   REFER: "; refer.print();  cout << endl;
  cout << "   EDATE: "; edate.print();  cout << endl;
  cout << "   DDATE: "; ddate.print();  cout << endl;
  cout << "   RDATE: "; rdate.print();  cout << endl;
  cout << "  ENDATE: "; endate.print(); cout << endl;

  if(stdheader){
    cout << " LIBNAME: "; libname.print(); cout << endl;
    cout << "     MAT: "; matnum.print();  cout << endl;
    cout << "  SUBLIB: "; sublib.print();  cout << endl;
    cout << "  FORMAT: "; format.print();  cout << endl;
  }
}
